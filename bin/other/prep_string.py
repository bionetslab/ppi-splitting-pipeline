# download from https://string-db.org/cgi/download?sessionId=bI0YqS2EiTEj
# -> organism: Homo sapiens. Advanced options: Includes AB pairs only.
# -> download 9606.protein.physical.links.full.v12.0.onlyAB.txt.gz
# -> download 9606.protein.links.full.v12.0.onlyAB.txt.gz
# -> download 9606.protein.aliases.v12.0.txt
# -> download https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/sec_ac.txt

import json
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

from bin.fetch_data import _next_page_url, _ssl_context

IDMAPPING_RUN_URL     = "https://rest.uniprot.org/idmapping/run"
IDMAPPING_STATUS_URL  = "https://rest.uniprot.org/idmapping/status/{job_id}"
IDMAPPING_RESULTS_URL = "https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}"
IDMAPPING_BATCH_SIZE  = 5000  # UniProt's id mapping jobs accept far more ids per job than /uniprotkb/search
path_to_string = "~/Downloads"

def submit_id_mapping_job(ids, from_db="STRING", to_db="UniProtKB"):
    data = urllib.parse.urlencode({"from": from_db, "to": to_db, "ids": ",".join(ids)}).encode()
    req = urllib.request.Request(IDMAPPING_RUN_URL, data=data)
    with urllib.request.urlopen(req, timeout=30, context=_ssl_context()) as resp:
        return json.loads(resp.read().decode())["jobId"]


def wait_for_id_mapping_job(job_id, poll_interval=1.0, timeout=300):
    url = IDMAPPING_STATUS_URL.format(job_id=job_id)
    deadline = time.time() + timeout
    while True:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=30, context=_ssl_context()) as resp:
            body = json.loads(resp.read().decode())
        if "results" in body or "failedIds" in body or body.get("jobStatus") == "FINISHED":
            return
        if time.time() > deadline:
            raise RuntimeError(f"ID mapping job {job_id} did not finish within {timeout}s")
        time.sleep(poll_interval)


def fetch_id_mapping_results(job_id):
    """Return {from_id: [primary_accession, ...]} for a finished STRING->UniProtKB job."""
    resolved = {}
    url = f"{IDMAPPING_RESULTS_URL.format(job_id=job_id)}?fields=accession&size=500"
    while url:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=60, context=_ssl_context()) as resp:
            body = json.loads(resp.read().decode())
            next_url = _next_page_url(resp.headers.get("Link"))
        for r in body.get("results", []):
            resolved.setdefault(r["from"], []).append(r["to"]["primaryAccession"])
        url = next_url
    return resolved


def resolve_via_id_mapping(string_ids, batch_size=IDMAPPING_BATCH_SIZE):
    """Map STRING protein IDs directly to their current UniProt accession via
    UniProt's ID Mapping service (from=STRING, to=UniProtKB). STRING IDs are
    already unambiguous per-organism identifiers, so this resolves demerges
    correctly with no candidate/organism disambiguation needed -- unlike
    going through sec_ac.txt alone.

    Returns ({string_id: accession}, [(string_id, [0 or >1 hits]), ...]).
    """
    string_ids = sorted(string_ids)
    resolved, failed = {}, []
    for i in range(0, len(string_ids), batch_size):
        batch = string_ids[i:i + batch_size]
        job_id = submit_id_mapping_job(batch)
        wait_for_id_mapping_job(job_id)
        hits = fetch_id_mapping_results(job_id)
        for string_id in batch:
            accs = hits.get(string_id, [])
            if len(accs) == 1:
                resolved[string_id] = accs[0]
            else:
                failed.append((string_id, accs))
    return resolved, failed

string = pd.read_csv(f"{path_to_string}/9606.protein.physical.links.full.v12.0.onlyAB.txt", sep=" ")
string_all = pd.read_csv(f"{path_to_string}/9606.protein.links.full.v12.0.onlyAB.txt", sep=" ")
string_info = pd.read_csv(f"{path_to_string}/9606.protein.aliases.v12.0.txt", sep="\t")

string_info = string_info[string_info["source"] == "UniProt_AC"]
secondary_acc_mapping = pd.read_csv(f"{path_to_string}/sec_ac.txt",
                                    sep="\\s+",
                                    names=["secondary_ac", "primary_ac"],
                                    skiprows=31
                                    )
string_info_merged = string_info.merge(secondary_acc_mapping, how="left", left_on="alias", right_on="secondary_ac")
string_info_merged["primary_ac"] = string_info_merged["primary_ac"].fillna(string_info_merged["alias"])
string_info_merged = string_info_merged[["#string_protein_id", "primary_ac"]]
string_info_merged = string_info_merged.drop_duplicates()
# show all remaining duplicates in #string_protein_id
dups = string_info_merged[string_info_merged.duplicated(subset="#string_protein_id", keep=False)]
lookup_with_api = set(dups["#string_protein_id"])

print(
    f"Resolving {len(lookup_with_api)} ambiguous STRING protein(s) via "
    f"UniProt's STRING->UniProtKB ID mapping...",
    file=sys.stderr,
)
resolved_map, failed_ids = resolve_via_id_mapping(lookup_with_api)

if failed_ids:
    print(
        f"Warning: ID mapping did not return exactly one UniProt accession for "
        f"{len(failed_ids)} STRING protein(s):",
        file=sys.stderr,
    )
    for string_id, accs in failed_ids:
        print(f"  {string_id}: {accs or '(no hit)'}", file=sys.stderr)

resolved_df = pd.DataFrame(
    [{"#string_protein_id": k, "primary_ac": v} for k, v in resolved_map.items()],
    columns=["#string_protein_id", "primary_ac"],
)

# Replace the ambiguous rows with their single ID-mapping-resolved match;
# STRING IDs left unresolved above are dropped here and only surfaced via
# the warning.
string_info_merged = pd.concat([
    string_info_merged[~string_info_merged["#string_protein_id"].isin(lookup_with_api)],
    resolved_df,
], ignore_index=True)

# construct final csv protein1,protein2
string = string.rename(columns={"protein1": "string1", "protein2": "string2"})
string_all = string_all.rename(columns={"protein1": "string1", "protein2": "string2"})

for df, filename in [(string_all, "string_all"), (string, "string")]:
    string_out = df.merge(string_info_merged, left_on="string1", right_on="#string_protein_id")
    string_out = string_out.rename(columns={"primary_ac": "protein1"})
    string_out = string_out.merge(string_info_merged, left_on="string2", right_on="#string_protein_id")
    string_out = string_out.rename(columns={"primary_ac": "protein2"})
    string_out = string_out[["protein1", "protein2", "experiments", "experiments_transferred", "database", "database_transferred", "textmining", "textmining_transferred", "combined_score"]]
    string_out = string_out.drop_duplicates()
    string_out.to_csv(f"../data/{filename}.csv", index=False)

string900 = string_out[string_out["combined_score"] >= 900]
string900.to_csv("../data/string_900.csv", index=False)

string_exp = string_out[string_out["experiments"] > 0]
string_exp.to_csv("../data/string_experimental.csv", index=False)
string_db = string_out[string_out["database"] > 0]
string_db.to_csv("../data/string_database.csv", index=False)
string_text = string_out[string_out["textmining"] > 0]
string_text.to_csv("../data/string_textmining.csv", index=False)
