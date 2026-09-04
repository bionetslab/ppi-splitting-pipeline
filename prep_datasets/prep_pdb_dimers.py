#!/usr/bin/env python3
"""Download PDB dimers, prepare their interactions, and export benchmark files.

The script searches for and downloads all matching RCSB biological assemblies,
extracts one sequence and resolved-residue mask per PDB entity, and filters out
interactions with unsuitable sequences. Exact sequences are reclustered before
duplicate cluster-pair interactions are removed. Finally, the script writes the
interaction, sequence, GO-annotation, and species files used by the PPI
splitting benchmark. All artifacts are written below the same ``--outdir``.

"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import gemmi
import pandas as pd


# =============================================================================
# STEP 1: SEARCH RCSB, WRITE CANDIDATE METADATA, AND DOWNLOAD ASSEMBLIES
# =============================================================================

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
FILES_BASE = "https://files.rcsb.org/download"


@dataclass(frozen=True)
class Candidate:
    """Metadata retained for one protein-dimer biological assembly."""

    assembly_id: str
    pdb_id: str
    assembly_number: str
    entity_pair: tuple[str, str]
    uniprot_pair: tuple[str, str]
    species_pair: tuple[str, str]
    taxonomy_pair: tuple[str, str]
    dimer_type: str
    cluster_pair: tuple[str, str]
    resolution: float
    modeled_residue_count: int
    method: str
    oligomeric_count: str
    download_url: str
    local_filename: str


def http_json_post(
    url: str,
    payload: dict[str, Any],
    *,
    timeout: int = 120,
    retries: int = 4,
) -> dict[str, Any]:
    """POST JSON with short exponential retries for transient API failures."""

    body = json.dumps(payload).encode("utf-8")
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    last_error: BaseException | None = None

    for attempt in range(retries):
        request = urllib.request.Request(url, data=body, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                response_body = response.read()
            return json.loads(response_body.decode("utf-8")) if response_body else {}
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            if isinstance(exc, urllib.error.HTTPError) and 400 <= exc.code < 500 and exc.code not in (408, 429):
                detail = exc.read().decode("utf-8", errors="replace")
                raise RuntimeError(f"HTTP {exc.code} from {url}: {detail}") from exc
            time.sleep(2**attempt)

    raise RuntimeError(f"Failed POST {url}: {last_error!r}") from last_error


def download_binary(url: str, path: Path, *, timeout: int = 120, retries: int = 4) -> bool:
    """Download one assembly atomically, returning False for absent RCSB files."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_name(path.name + ".tmp")
    last_error: BaseException | None = None

    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response, temporary_path.open("wb") as output:
                while chunk := response.read(1024 * 1024):
                    output.write(chunk)
            temporary_path.replace(path)
            return True
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            temporary_path.unlink(missing_ok=True)
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (404, 410):
                return False
            if isinstance(exc, urllib.error.HTTPError) and 400 <= exc.code < 500 and exc.code not in (408, 429):
                raise RuntimeError(f"HTTP {exc.code} while downloading {url}") from exc
            time.sleep(2**attempt)

    raise RuntimeError(f"Failed download {url}: {last_error!r}") from last_error


def search_candidate_assemblies(no_ligands: bool) -> list[str]:
    """Return protein-only, globally symmetric assemblies with two polymers."""

    def terminal(attribute: str, operator: str, value: Any) -> dict[str, Any]:
        return {
            "type": "terminal",
            "service": "text",
            "parameters": {"attribute": attribute, "operator": operator, "value": value},
        }

    nodes = [
        terminal("rcsb_entry_info.selected_polymer_entity_types", "exact_match", "Protein (only)"),
        terminal("rcsb_assembly_info.polymer_entity_instance_count", "equals", 2),
        terminal("rcsb_struct_symmetry.kind", "exact_match", "Global Symmetry"),
    ]
    if no_ligands:
        nodes.append(terminal("rcsb_assembly_info.nonpolymer_entity_instance_count", "equals", 0))

    payload = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "assembly",
        "request_options": {
            "return_all_hits": True,
            "results_content_type": ["experimental"],
            "results_verbosity": "compact",
        },
    }
    result_set = http_json_post(SEARCH_URL, payload).get("result_set") or []
    identifiers: set[str] = set()
    for hit in result_set:
        if isinstance(hit, str):
            identifier = hit
        elif isinstance(hit, dict):
            identifier = hit.get("identifier")
        else:
            continue
        if identifier:
            identifiers.add(str(identifier).upper())
    return sorted(identifiers)


ASSEMBLY_QUERY = """
query AssemblyBatch($assembly_ids: [String!]!) {
  assemblies(assembly_ids: $assembly_ids) {
    rcsb_id
    rcsb_assembly_info {
      modeled_polymer_monomer_count
      polymer_entity_instance_count
      polymer_entity_instance_count_protein
    }
    pdbx_struct_assembly {
      oligomeric_count
    }
    pdbx_struct_assembly_gen {
      asym_id_list
    }
    entry {
      rcsb_entry_info {
        experimental_method
        resolution_combined
      }
      polymer_entities {
        rcsb_id
        entity_poly {
          rcsb_entity_polymer_type
        }
        rcsb_entity_source_organism {
          ncbi_taxonomy_id
          ncbi_scientific_name
        }
        rcsb_polymer_entity_container_identifiers {
          asym_ids
          reference_sequence_identifiers {
            database_name
            database_accession
          }
        }
        uniprots {
          rcsb_id
          rcsb_uniprot_container_identifiers {
            uniprot_id
          }
        }
        rcsb_polymer_entity_group_membership {
          aggregation_method
          similarity_cutoff
          group_id
        }
      }
    }
  }
}
"""


def fetch_assembly_metadata(assembly_ids: Sequence[str], batch_size: int) -> list[dict[str, Any]]:
    """Fetch the metadata needed by this and later pipeline steps in batches."""

    metadata: list[dict[str, Any]] = []
    for start in range(0, len(assembly_ids), batch_size):
        batch = assembly_ids[start : start + batch_size]
        payload = {"query": ASSEMBLY_QUERY, "variables": {"assembly_ids": list(batch)}}
        response = http_json_post(GRAPHQL_URL, payload)
        if response.get("errors"):
            raise RuntimeError("GraphQL errors:\n" + json.dumps(response["errors"], indent=2))
        metadata.extend(((response.get("data") or {}).get("assemblies") or []))
        print(f"Fetched metadata batch {start // batch_size + 1}: {len(batch)} assemblies", file=sys.stderr)
    return metadata


def split_csv_ids(value: Any) -> list[str]:
    """Normalize GraphQL ID fields, which may be lists or comma-separated strings."""

    if value is None:
        return []
    if isinstance(value, list):
        return [identifier for item in value for identifier in split_csv_ids(item)]
    return [identifier.strip() for identifier in str(value).split(",") if identifier.strip()]


def cluster_id_for_entity(entity: dict[str, Any]) -> str | None:
    """Read the 100% RCSB sequence cluster, using a stable singleton fallback."""

    for membership in entity.get("rcsb_polymer_entity_group_membership") or []:
        if not isinstance(membership, dict) or membership.get("aggregation_method") != "sequence_identity":
            continue
        try:
            is_exact_cluster = int(membership.get("similarity_cutoff")) == 100
        except (TypeError, ValueError):
            continue
        if is_exact_cluster and membership.get("group_id"):
            return str(membership["group_id"])

    rcsb_id = entity.get("rcsb_id")
    return f"singleton:{rcsb_id}" if rcsb_id else None


def uniprot_ids_for_entity(entity: dict[str, Any]) -> str:
    """Return sorted UniProt accessions, with generic cross-references as fallback."""

    accessions: set[str] = set()
    for uniprot in entity.get("uniprots") or []:
        if not isinstance(uniprot, dict):
            continue
        container = uniprot.get("rcsb_uniprot_container_identifiers") or {}
        values = container.get("uniprot_id") or uniprot.get("rcsb_id")
        if isinstance(values, list):
            accessions.update(str(value) for value in values if value)
        elif values:
            accessions.add(str(values))

    if not accessions:
        container = entity.get("rcsb_polymer_entity_container_identifiers") or {}
        for reference in container.get("reference_sequence_identifiers") or []:
            if not isinstance(reference, dict) or str(reference.get("database_name") or "").lower() != "uniprot":
                continue
            values = reference.get("database_accession")
            if isinstance(values, list):
                accessions.update(str(value) for value in values if value)
            elif values:
                accessions.add(str(values))

    return ";".join(sorted(accessions))


def species_for_entity(entity: dict[str, Any]) -> tuple[str, str]:
    """Return sorted source-species names and NCBI taxonomy IDs."""

    names: set[str] = set()
    taxonomy_ids: set[str] = set()
    for organism in entity.get("rcsb_entity_source_organism") or []:
        if not isinstance(organism, dict):
            continue
        name = organism.get("ncbi_scientific_name")
        taxonomy_id = organism.get("ncbi_taxonomy_id")
        names.update(str(value) for value in (name if isinstance(name, list) else [name]) if value)
        taxonomy_ids.update(
            str(value) for value in (taxonomy_id if isinstance(taxonomy_id, list) else [taxonomy_id]) if value
        )
    return ";".join(sorted(names)), ";".join(sorted(taxonomy_ids))


def participating_protein_entities(assembly: dict[str, Any]) -> list[dict[str, Any]]:
    """Resolve generated assembly asym IDs to unique protein entities."""

    asym_to_entity: dict[str, dict[str, Any]] = {}
    for entity in ((assembly.get("entry") or {}).get("polymer_entities") or []):
        polymer_type = (entity.get("entity_poly") or {}).get("rcsb_entity_polymer_type") or ""
        if "Protein" not in str(polymer_type):
            continue
        container = entity.get("rcsb_polymer_entity_container_identifiers") or {}
        for asym_id in split_csv_ids(container.get("asym_ids")):
            asym_to_entity[asym_id] = entity

    entities_by_id: dict[str, dict[str, Any]] = {}
    for generator in assembly.get("pdbx_struct_assembly_gen") or []:
        for asym_id in split_csv_ids(generator.get("asym_id_list")):
            entity = asym_to_entity.get(asym_id)
            entity_id = str((entity or {}).get("rcsb_id") or "")
            if entity_id:
                entities_by_id[entity_id] = entity
    return [entities_by_id[entity_id] for entity_id in sorted(entities_by_id)]


def candidate_from_assembly(assembly: dict[str, Any]) -> Candidate | None:
    """Validate one API record and convert it to a candidate TSV row."""

    assembly_id = str(assembly.get("rcsb_id") or "").upper()
    if "-" not in assembly_id:
        return None
    pdb_id, assembly_number = assembly_id.rsplit("-", 1)

    info = assembly.get("rcsb_assembly_info") or {}
    if info.get("polymer_entity_instance_count") != 2 or info.get("polymer_entity_instance_count_protein") != 2:
        return None

    entities = participating_protein_entities(assembly)
    entity_records: list[tuple[str, str, str, str, str]] = []
    for entity in entities:
        entity_id = str(entity.get("rcsb_id") or "")
        cluster_id = cluster_id_for_entity(entity)
        if not entity_id or cluster_id is None:
            print(f"Skipping {assembly_id}: missing protein-entity or cluster information", file=sys.stderr)
            return None
        uniprot_id = uniprot_ids_for_entity(entity)
        species, taxonomy = species_for_entity(entity)
        entity_records.append((entity_id, cluster_id, uniprot_id, species, taxonomy))

    if len(entity_records) == 1:
        entity_records *= 2
        dimer_type = "homo"
    elif len(entity_records) == 2:
        entity_records.sort(key=lambda record: record[0])
        dimer_type = "hetero"
    else:
        print(
            f"Skipping {assembly_id}: expected 1 or 2 unique protein entities, got {len(entity_records)}",
            file=sys.stderr,
        )
        return None

    resolution_values = ((assembly.get("entry") or {}).get("rcsb_entry_info") or {}).get(
        "resolution_combined"
    ) or []
    if isinstance(resolution_values, (int, float)):
        resolution_values = [resolution_values]
    numeric_resolutions: list[float] = []
    for value in resolution_values:
        try:
            if value is not None:
                numeric_resolutions.append(float(value))
        except (TypeError, ValueError):
            continue
    resolution = min(numeric_resolutions) if numeric_resolutions else math.inf

    methods = ((assembly.get("entry") or {}).get("rcsb_entry_info") or {}).get("experimental_method") or []
    method = methods if isinstance(methods, str) else ";".join(str(value) for value in methods if value is not None)
    oligomeric_count = str((assembly.get("pdbx_struct_assembly") or {}).get("oligomeric_count") or "")
    local_filename = f"{pdb_id.lower()}-assembly{assembly_number}.cif.gz"

    return Candidate(
        assembly_id=assembly_id,
        pdb_id=pdb_id,
        assembly_number=assembly_number,
        entity_pair=(entity_records[0][0], entity_records[1][0]),
        uniprot_pair=(entity_records[0][2], entity_records[1][2]),
        species_pair=(entity_records[0][3], entity_records[1][3]),
        taxonomy_pair=(entity_records[0][4], entity_records[1][4]),
        dimer_type=dimer_type,
        cluster_pair=(entity_records[0][1], entity_records[1][1]),
        resolution=resolution,
        modeled_residue_count=int(info.get("modeled_polymer_monomer_count") or 0),
        method=method,
        oligomeric_count=oligomeric_count,
        download_url=f"{FILES_BASE}/{local_filename}",
        local_filename=local_filename,
    )


def write_candidate_tsv(path: Path, candidates: Sequence[Candidate]) -> None:
    """Write the complete candidate schema consumed by downstream steps."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(
            [
                "assembly_id",
                "pdb_id",
                "assembly_number",
                "entity_pair",
                "uniprot_pair",
                "uniprot_1",
                "uniprot_2",
                "species_pair",
                "species_1",
                "species_2",
                "taxonomy_pair",
                "taxonomy_1",
                "taxonomy_2",
                "dimer_type",
                "cluster_pair_100pct",
                "resolution_best_angstrom",
                "modeled_polymer_monomer_count",
                "experimental_method",
                "oligomeric_count",
                "download_url",
                "local_filename",
            ]
        )
        for candidate in candidates:
            writer.writerow(
                [
                    candidate.assembly_id,
                    candidate.pdb_id,
                    candidate.assembly_number,
                    ",".join(candidate.entity_pair),
                    "|".join(candidate.uniprot_pair),
                    *candidate.uniprot_pair,
                    "|".join(candidate.species_pair),
                    *candidate.species_pair,
                    "|".join(candidate.taxonomy_pair),
                    *candidate.taxonomy_pair,
                    candidate.dimer_type,
                    ",".join(candidate.cluster_pair),
                    "" if math.isinf(candidate.resolution) else candidate.resolution,
                    candidate.modeled_residue_count,
                    candidate.method,
                    candidate.oligomeric_count,
                    candidate.download_url,
                    candidate.local_filename,
                ]
            )


def write_unavailable_download_tsv(path: Path, unavailable: Sequence[tuple[Candidate, str]]) -> None:
    """Record API candidates whose biological-assembly coordinate file is absent."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(["assembly_id", "pdb_id", "assembly_number", "download_url", "local_filename", "reason"])
        for candidate, reason in unavailable:
            writer.writerow(
                [
                    candidate.assembly_id,
                    candidate.pdb_id,
                    candidate.assembly_number,
                    candidate.download_url,
                    candidate.local_filename,
                    reason,
                ]
            )


def download_dimers(outdir: Path, batch_size: int, no_ligands: bool, no_download: bool, limit: int) -> None:
    """Run step 1 and download parsed candidates, excluding absent coordinate files."""

    print("Step 1/6: searching RCSB for protein-only dimer assemblies", file=sys.stderr)
    assembly_ids = search_candidate_assemblies(no_ligands)
    if limit:
        assembly_ids = assembly_ids[:limit]
    print(f"Search returned {len(assembly_ids)} candidate assemblies", file=sys.stderr)

    metadata = fetch_assembly_metadata(assembly_ids, batch_size)
    candidates = [candidate for assembly in metadata if (candidate := candidate_from_assembly(assembly)) is not None]
    candidates.sort(key=lambda candidate: (candidate.pdb_id, candidate.assembly_number))

    candidate_path = outdir / "dataset_iterations" / "01_all_candidate_assemblies.tsv"
    if no_download:
        write_candidate_tsv(candidate_path, candidates)
        print(f"Wrote {len(candidates)} candidates to {candidate_path}", file=sys.stderr)
        return

    assembly_dir = outdir / "assemblies"
    downloaded_candidates: list[Candidate] = []
    unavailable_downloads: list[tuple[Candidate, str]] = []
    for index, candidate in enumerate(candidates, start=1):
        destination = assembly_dir / candidate.local_filename
        if destination.exists() and destination.stat().st_size > 0:
            print(f"[{index}/{len(candidates)}] exists {destination}", file=sys.stderr)
            downloaded_candidates.append(candidate)
            continue
        print(f"[{index}/{len(candidates)}] downloading {candidate.assembly_id}", file=sys.stderr)
        if download_binary(candidate.download_url, destination):
            downloaded_candidates.append(candidate)
        else:
            reason = "assembly coordinate file returned HTTP 404/410"
            unavailable_downloads.append((candidate, reason))
            print(f"[{index}/{len(candidates)}] skipping {candidate.assembly_id}: {reason}", file=sys.stderr)

    write_candidate_tsv(candidate_path, downloaded_candidates)
    print(f"Wrote {len(downloaded_candidates)} downloaded candidates to {candidate_path}", file=sys.stderr)
    unavailable_path = outdir / "dataset_iterations" / "01_unavailable_candidate_assemblies.tsv"
    write_unavailable_download_tsv(unavailable_path, unavailable_downloads)
    if unavailable_downloads:
        print(f"Wrote {len(unavailable_downloads)} unavailable candidates to {unavailable_path}", file=sys.stderr)


# =============================================================================
# STEP 2: EXTRACT ENTITY SEQUENCES AND RESOLVED-RESIDUE MASKS
# =============================================================================


def extract_sequences(outdir: Path) -> None:
    """Extract one sequence and per-asym residue mask for every unique entity."""

    print("Step 2/6: extracting entity sequences and residue masks", file=sys.stderr)
    candidate_path = outdir / "dataset_iterations" / "01_all_candidate_assemblies.tsv"
    candidates = pd.read_csv(candidate_path, sep="\t", keep_default_na=False)

    records: list[dict[str, str]] = []
    seen_entities: set[str] = set()

    for row in candidates.itertuples(index=False):
        entities = [entity.strip() for entity in str(row.entity_pair).split(",")]
        clusters = [cluster.strip() for cluster in str(row.cluster_pair_100pct).split(",")]
        new_entities: list[tuple[str, str]] = []
        for entity, cluster in zip(entities, clusters):
            if entity in seen_entities:
                continue
            seen_entities.add(entity)
            new_entities.append((entity, cluster))
        if not new_entities:
            continue

        cif_path = outdir / "assemblies" / str(row.local_filename)
        if not cif_path.is_file():
            print(f"Warning: missing {cif_path}", file=sys.stderr)
            records.extend(
                {"entity_name": entity, "cluster_id": cluster, "sequence": "-", "binary_mask": ""}
                for entity, cluster in new_entities
            )
            continue

        # Read the compressed mmCIF once, then derive both sequence and atom-site
        # information from the same block. Entity remapping is needed because an
        # assembly file can renumber the entity IDs stored in ``entity_pair``.
        document = gemmi.cif.read_file(str(cif_path))
        block = document.sole_block()
        structure = gemmi.make_structure_from_block(block)
        structure.setup_entities()

        remapped_entity_ids = {
            str(original_id): str(entity_id)
            for entity_id, original_id in block.find(
                "_pdbx_entity_remapping.",
                ["entity_id", "orig_entity_id"],
            )
        }
        asym_ids_by_entity: dict[str, list[str]] = {}
        for asym_id, entity_id in block.find("_struct_asym.", ["id", "entity_id"]):
            asym_ids_by_entity.setdefault(str(entity_id), []).append(str(asym_id))

        sequence_by_entity = {
            str(entity.name): gemmi.one_letter_code(entity.full_sequence)
            for entity in structure.entities
            if entity.full_sequence
        }
        extraction_info: dict[str, tuple[str, str, list[str]]] = {}
        relevant_asym_ids: set[str] = set()
        for entity_name, cluster_id in new_entities:
            original_id = entity_name.rsplit("_", 1)[-1]
            entity_id = remapped_entity_ids.get(original_id, original_id)
            asym_ids = asym_ids_by_entity.get(entity_id, [])
            extraction_info[entity_name] = (cluster_id, sequence_by_entity.get(entity_id, ""), asym_ids)
            relevant_asym_ids.update(asym_ids)

        # label_seq_id is one-based. A set records each resolved residue only
        # once even though the atom-site table contains many atoms per residue.
        resolved_indices = {asym_id: set() for asym_id in relevant_asym_ids}
        for asym_id, sequence_id in block.find("_atom_site.", ["label_asym_id", "label_seq_id"]):
            asym_id = str(asym_id)
            if asym_id not in resolved_indices:
                continue
            try:
                residue_index = int(sequence_id) - 1
            except (TypeError, ValueError):
                continue
            if residue_index >= 0:
                resolved_indices[asym_id].add(residue_index)

        for entity_name, (cluster_id, sequence, asym_ids) in extraction_info.items():
            if not sequence:
                print(f"Warning: no sequence found for {entity_name} in {cif_path}", file=sys.stderr)
                records.append(
                    {"entity_name": entity_name, "cluster_id": cluster_id, "sequence": "-", "binary_mask": ""}
                )
                continue

            masks: list[str] = []
            for asym_id in asym_ids:
                mask = ["0"] * len(sequence)
                for residue_index in resolved_indices[asym_id]:
                    if residue_index < len(mask):
                        mask[residue_index] = "1"
                    else:
                        print(
                            f"Warning: residue index {residue_index} exceeds the {len(sequence)}-residue sequence "
                            f"for {entity_name} ({asym_id})",
                            file=sys.stderr,
                        )
                masks.append(f"{asym_id}: {''.join(mask)};")

            records.append(
                {
                    "entity_name": entity_name,
                    "cluster_id": cluster_id,
                    "sequence": sequence,
                    "binary_mask": " ".join(masks),
                }
            )

    sequence_path = outdir / "sequences" / "entity_sequences.tsv"
    sequence_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records, columns=["entity_name", "cluster_id", "sequence", "binary_mask"]).to_csv(
        sequence_path,
        sep="\t",
        index=False,
    )
    missing_count = sum(record["sequence"] == "-" for record in records)
    print(f"Wrote {len(records)} entity sequences to {sequence_path} ({missing_count} missing)", file=sys.stderr)


# =============================================================================
# STEP 3: FILTER INTERACTIONS BY SEQUENCE ATTRIBUTES
# =============================================================================


def filter_by_sequence_attributes(outdir: Path) -> None:
    """Keep interactions whose entities are 50-1250 aa with at most 15% X."""

    print("Step 3/6: filtering interactions by sequence attributes", file=sys.stderr)
    dataset_dir = outdir / "dataset_iterations"
    interactions = pd.read_csv(dataset_dir / "01_all_candidate_assemblies.tsv", sep="\t", keep_default_na=False)
    sequences = pd.read_csv(outdir / "sequences" / "entity_sequences.tsv", sep="\t", keep_default_na=False)
    sequence_lengths = sequences["sequence"].str.len()

    def remove_entities(frame: pd.DataFrame, invalid_entities: set[str]) -> pd.DataFrame:
        contains_invalid_entity = frame["entity_pair"].map(
            lambda pair: any(entity.strip() in invalid_entities for entity in str(pair).split(","))
        )
        return frame.loc[~contains_invalid_entity]

    print(f"Interactions before filtering: {len(interactions)}", file=sys.stderr)

    long_entities = set(sequences.loc[sequence_lengths > 1250, "entity_name"])
    filtered = remove_entities(interactions, long_entities)
    print(f"After removing sequences longer than 1250: {len(filtered)}", file=sys.stderr)

    short_entities = set(sequences.loc[sequence_lengths < 50, "entity_name"])
    filtered = remove_entities(filtered, short_entities)
    print(f"After removing sequences shorter than 50: {len(filtered)}", file=sys.stderr)

    unknown_fraction = sequences["sequence"].map(lambda sequence: sequence.count("X") / len(sequence))
    unknown_entities = set(sequences.loc[unknown_fraction > 0.15, "entity_name"])
    filtered = remove_entities(filtered, unknown_entities)
    print(f"After removing sequences with more than 15% X: {len(filtered)}", file=sys.stderr)

    output_path = dataset_dir / "02_filtered_by_seq_attributes.tsv"
    filtered.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote filtered interactions to {output_path}", file=sys.stderr)


# =============================================================================
# STEP 4: RECLUSTER EXACT SEQUENCES AT 100% IDENTITY
# =============================================================================


def cluster_by_100_sequence_identity(outdir: Path) -> None:
    """Assign one deterministic internal cluster ID to each exact sequence."""

    print("Step 4/6: clustering exact sequences", file=sys.stderr)
    dataset_dir = outdir / "dataset_iterations"
    sequence_dir = outdir / "sequences"
    sequences = pd.read_csv(sequence_dir / "entity_sequences.tsv", sep="\t")
    interactions = pd.read_csv(dataset_dir / "02_filtered_by_seq_attributes.tsv", sep="\t")

    # Restrict the reusable entity table to proteins that survived Step 3.
    valid_entity_names = set()
    for _, row in interactions.iterrows():
        entity_pair = row["entity_pair"].split(",")
        valid_entity_names.update(entity.strip() for entity in entity_pair)
    sequences = sequences[sequences["entity_name"].isin(valid_entity_names)]

    # Pandas sorts group keys by default. That ordering intentionally defines
    # the stable fix_<n>_100 IDs used by every later pipeline stage.
    sequence_to_clusters = (
        sequences.groupby("sequence")
        .agg(
            {
                "entity_name": lambda values: list(values.unique()),
                "cluster_id": lambda values: list(values.unique()),
            }
        )
        .reset_index()
    )

    # One exact sequence should have one original RCSB cluster. As in the
    # existing pipeline, remove every interaction from a PDB containing an
    # inconsistent mapping, while retaining all groups in unique_sequences.tsv.
    inconsistent_pdb_ids = []
    for index, row in sequence_to_clusters.iterrows():
        if len(row["cluster_id"]) > 1:
            inconsistent_pdb_ids.extend(entity.split("_")[0] for entity in row["entity_name"])
            print(
                f"Warning: Row {index} has more than one unique cluster ID. "
                f"Entity names: {row['entity_name']}, Cluster IDs: {row['cluster_id']}, "
                f"Sequence: {row['sequence']}",
                file=sys.stderr,
            )

    print(f"Number of interactions before inconsistent-cluster filtering: {len(interactions)}", file=sys.stderr)
    interactions = interactions[~interactions["pdb_id"].isin(inconsistent_pdb_ids)]
    print(f"Number of interactions after inconsistent-cluster filtering: {len(interactions)}", file=sys.stderr)

    sequence_to_clusters["new_cluster_id"] = sequence_to_clusters.index.map(
        lambda index: f"fix_{index + 1}_100"
    )
    sequences = sequences.merge(
        sequence_to_clusters[["sequence", "new_cluster_id"]],
        on="sequence",
        how="left",
    )

    # Preserve entity-pair orientation. Step 5 canonicalizes a temporary copy
    # only when it checks whether two interactions represent the same pair.
    new_cluster_pairs = []
    for _, row in interactions.iterrows():
        cluster_ids = []
        for entity in row["entity_pair"].split(","):
            cluster_id = sequences[sequences["entity_name"] == entity]["new_cluster_id"].values
            if len(cluster_id) > 0:
                cluster_ids.append(cluster_id[0])
            else:
                print(f"Warning: Entity {entity} not found in sequence table", file=sys.stderr)
                cluster_ids.append("unknown_cluster")
        new_cluster_pairs.append(",".join(cluster_ids))

    interactions.insert(
        interactions.columns.get_loc("cluster_pair_100pct") + 1,
        "new_cluster_pair",
        new_cluster_pairs,
    )

    interaction_output = dataset_dir / "03_new_clusters.tsv"
    sequence_output = sequence_dir / "unique_sequences.tsv"
    interactions.to_csv(interaction_output, sep="\t", index=False)
    sequence_to_clusters.to_csv(sequence_output, sep="\t", index=False)
    print(f"Wrote clustered interactions to {interaction_output}", file=sys.stderr)
    print(f"Wrote unique sequences to {sequence_output}", file=sys.stderr)


# =============================================================================
# STEP 5: REMOVE DUPLICATE SEQUENCE-PAIR INTERACTIONS
# =============================================================================


def filter_duplicate_interactions(outdir: Path) -> None:
    """Keep the first interaction for each unordered exact-sequence pair."""

    print("Step 5/6: removing duplicate sequence-pair interactions", file=sys.stderr)
    dataset_dir = outdir / "dataset_iterations"
    interactions = pd.read_csv(dataset_dir / "03_new_clusters.tsv", sep="\t")

    interactions["tmp_cluster_pair"] = interactions["new_cluster_pair"].apply(
        lambda pair: ",".join(sorted(pair.split(",")))
    )
    filtered = interactions.drop_duplicates(subset=["tmp_cluster_pair"], keep="first").copy()

    print(f"Number of interactions before filtering duplicates: {len(interactions)}", file=sys.stderr)
    print(f"Number of interactions after filtering duplicates: {len(filtered)}", file=sys.stderr)

    filtered.drop(columns=["tmp_cluster_pair"], inplace=True)
    output_path = dataset_dir / "04_removed_duplicates.tsv"
    filtered.to_csv(output_path, sep="\t", index=False)
    print(f"Wrote deduplicated interactions to {output_path}", file=sys.stderr)


# =============================================================================
# STEP 6: EXPORT FILES FOR THE PPI SPLITTING BENCHMARK
# =============================================================================


def fetch_uniprot_go_annotations(
    uniprot_ids: Sequence[str],
    batch_size: int = 100,
    sleep_seconds: float = 0.2,
) -> tuple[dict[str, tuple[str, str, str]], list[str]]:
    """Fetch biological-process, molecular-function, and component GO terms."""

    fields = "accession,go_p,go_f,go_c"

    def fetch_tsv(accessions: Sequence[str]) -> str:
        query = urllib.parse.urlencode(
            {
                "accessions": ",".join(accessions),
                "fields": fields,
                "format": "tsv",
            }
        )
        with urllib.request.urlopen(
            f"https://rest.uniprot.org/uniprotkb/accessions?{query}",
            timeout=120,
        ) as response:
            return response.read().decode("utf-8")

    go_by_uniprot: dict[str, tuple[str, str, str]] = {}
    missing_uniprot_ids: list[str] = []
    total_batches = (len(uniprot_ids) + batch_size - 1) // batch_size

    for start in range(0, len(uniprot_ids), batch_size):
        batch = uniprot_ids[start : start + batch_size]
        print(
            f"UniProt GO batch {start // batch_size + 1}/{total_batches}",
            file=sys.stderr,
        )

        try:
            response_text = fetch_tsv(batch)
        except urllib.error.HTTPError:
            # A single invalid accession can make a whole batch fail. Retry the
            # accessions separately so valid proteins still receive GO terms.
            response_lines = [
                "Entry\tGene Ontology (biological process)\t"
                "Gene Ontology (molecular function)\t"
                "Gene Ontology (cellular component)"
            ]
            for uniprot_id in batch:
                try:
                    response_lines.extend(fetch_tsv([uniprot_id]).splitlines()[1:])
                except urllib.error.HTTPError:
                    missing_uniprot_ids.append(uniprot_id)
            response_text = "\n".join(response_lines)

        for row in csv.DictReader(response_text.splitlines(), delimiter="\t"):
            go_by_uniprot[row["Entry"]] = (
                row["Gene Ontology (biological process)"],
                row["Gene Ontology (molecular function)"],
                row["Gene Ontology (cellular component)"],
            )

        if sleep_seconds and start + batch_size < len(uniprot_ids):
            time.sleep(sleep_seconds)

    return go_by_uniprot, missing_uniprot_ids


def create_benchmark_files(outdir: Path) -> None:
    """Write interaction, FASTA, GO, species, and missing-accession files."""

    print("Step 6/6: exporting PPI splitting benchmark files", file=sys.stderr)
    interactions_path = outdir / "dataset_iterations" / "04_removed_duplicates.tsv"
    sequences_path = outdir / "sequences" / "entity_sequences.tsv"
    benchmark_dir = outdir / "ppi_splitting_pipeline"

    splitter = re.compile(r"[;|,\s]+")
    missing_values = {"", "-", "na", "n/a", "nan", "none", "null"}
    interactions: list[tuple[str, str, str, str]] = []
    entity_to_uniprots: defaultdict[str, Counter[str]] = defaultdict(Counter)
    entity_to_taxa: defaultdict[str, Counter[str]] = defaultdict(Counter)

    with interactions_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"entity_pair", "uniprot_1", "uniprot_2", "taxonomy_1", "taxonomy_2"}
        missing_columns = required_columns - set(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(f"{interactions_path} is missing columns: {', '.join(sorted(missing_columns))}")

        for row_number, row in enumerate(reader, start=2):
            entity_ids = [part.strip() for part in row["entity_pair"].split(",")]
            if len(entity_ids) != 2:
                raise ValueError(f"{interactions_path}:{row_number} has an invalid entity_pair")

            uniprot_ids = []
            for column in ("uniprot_1", "uniprot_2"):
                uniprot_ids.append(
                    next(
                        (
                            value.strip().upper()
                            for value in splitter.split(row[column])
                            if value.strip().lower() not in missing_values
                        ),
                        "",
                    )
                )

            entity_1, entity_2 = entity_ids
            uniprot_1, uniprot_2 = uniprot_ids
            interactions.append((entity_1, entity_2, uniprot_1, uniprot_2))

            if uniprot_1:
                entity_to_uniprots[entity_1][uniprot_1] += 1
            if uniprot_2:
                entity_to_uniprots[entity_2][uniprot_2] += 1
            if row["taxonomy_1"].strip():
                entity_to_taxa[entity_1][row["taxonomy_1"].strip()] += 1
            if row["taxonomy_2"].strip():
                entity_to_taxa[entity_2][row["taxonomy_2"].strip()] += 1

    protein_ids = sorted({entity_id for interaction in interactions for entity_id in interaction[:2]})
    entity_to_uniprot = {
        entity_id: counts.most_common(1)[0][0]
        for entity_id, counts in entity_to_uniprots.items()
    }
    entity_to_taxon = {
        entity_id: counts.most_common(1)[0][0]
        for entity_id, counts in entity_to_taxa.items()
    }
    uniprot_ids = sorted(set(entity_to_uniprot.values()))
    protein_id_set = set(protein_ids)

    print(
        f"Benchmark input: {len(interactions)} interactions, {len(protein_ids)} proteins, "
        f"{len(uniprot_ids)} UniProt IDs",
        file=sys.stderr,
    )

    benchmark_dir.mkdir(parents=True, exist_ok=True)
    temporary_interactions = benchmark_dir / "pdb_interactions.csv.tmp"
    temporary_sequences = benchmark_dir / "sequences.fasta.tmp"
    temporary_go = benchmark_dir / "go_annotations.tsv.tmp"
    temporary_species = benchmark_dir / "species.tsv.tmp"
    temporary_missing = benchmark_dir / "missing_uniprot_go_ids.txt.tmp"

    with temporary_interactions.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["pdbid_1", "pdbid_2", "protein1", "protein2"])
        writer.writerows(interactions)

    sequences: dict[str, str] = {}
    with sequences_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["entity_name"] in protein_id_set:
                sequences[row["entity_name"]] = row["sequence"]

    with temporary_sequences.open("w") as handle:
        for protein_id in protein_ids:
            if protein_id not in sequences:
                print(f"Warning: missing sequence for {protein_id}", file=sys.stderr)
                continue
            handle.write(f">{protein_id}\n{sequences[protein_id]}\n")

    go_by_uniprot, missing_uniprot_ids = fetch_uniprot_go_annotations(uniprot_ids)

    with temporary_go.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id", "go_bp", "go_mf", "go_cc"])
        for protein_id in protein_ids:
            go_terms = go_by_uniprot.get(entity_to_uniprot.get(protein_id, ""), ("", "", ""))
            writer.writerow([protein_id, *go_terms])

    with temporary_species.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id", "taxon_id"])
        for protein_id in protein_ids:
            writer.writerow([protein_id, entity_to_taxon.get(protein_id, "")])

    with temporary_missing.open("w") as handle:
        for uniprot_id in sorted(missing_uniprot_ids):
            handle.write(f"{uniprot_id}\n")

    temporary_interactions.replace(benchmark_dir / "pdb_interactions.csv")
    temporary_sequences.replace(benchmark_dir / "sequences.fasta")
    temporary_go.replace(benchmark_dir / "go_annotations.tsv")
    temporary_species.replace(benchmark_dir / "species.tsv")
    temporary_missing.replace(benchmark_dir / "missing_uniprot_go_ids.txt")

    print(f"Wrote benchmark files to {benchmark_dir}", file=sys.stderr)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run six PDB-dimer preparation and benchmark-export steps.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("rcsb_protein_dimers"),
        help="Output root [default: %(default)s]",
    )
    parser.add_argument("--batch-size", type=int, default=200, help="RCSB GraphQL batch size [default: %(default)s]")
    parser.add_argument("--no-ligands", action="store_true", help="Require assemblies with no non-polymer instances")
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Skip coordinate downloads; steps 2-6 require assemblies already present",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Limit assemblies after search for debugging [default: all]",
    )
    args = parser.parse_args()
    if args.batch_size < 1:
        parser.error("--batch-size must be positive")
    if args.limit < 0:
        parser.error("--limit cannot be negative")
    return args


def main():
    args = parse_args()
    download_dimers(args.outdir, args.batch_size, args.no_ligands, args.no_download, args.limit)
    extract_sequences(args.outdir)
    filter_by_sequence_attributes(args.outdir)
    cluster_by_100_sequence_identity(args.outdir)
    filter_duplicate_interactions(args.outdir)
    create_benchmark_files(args.outdir)


if __name__ == "__main__":
    main()
