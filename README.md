# PPI Splitting Pipeline

Automated leakage-aware splitting of a protein–protein interaction (PPI) dataset into train, validation, and test sets, with redundancy removal, negative sampling, embedding-based classification, and bias analysis.

## Quick Start

### 1. Install dependencies

```bash
conda env create -f environment.yml
conda activate ppi-splitting-pipeline
```

### 2. Prepare your input

Create a CSV file with at least two columns (`protein1`, `protein2`) containing UniProt accession IDs. Additional columns (e.g. STRING evidence scores) are preserved throughout the pipeline.

```
protein1,protein2
P45985,Q14315
Q86TC9,P35609
O14836-2,P12345
...
```

### 3. Run the pipeline

```bash
nextflow run main.nf --ppis ppis.csv --outdir results
```

If you have a GPU, additionally specify `-profile gpu`, which will submit only the embedding steps to the GPU:

```bash
nextflow run main.nf --ppis ppis.csv --outdir results -profile gpu -c my_config.config
```

Config example for an HPC with slurm and a dedicated GPU queue: https://nf-co.re/configs/daisybio/. Important part:

```bash
 profiles {
     ...
     gpu {
            docker.runOptions       = '-u $(id -u):$(id -g) --gpus all'
            apptainer.runOptions    = '--nv'
            singularity.runOptions  = '--nv'
        process{
                withLabel:process_gpu {
                    queue = 'shared-gpu'
                    clusterOptions = '--qos=limitgpus --gpus=a40:1 --exclude compms-gpu-1.exbio.wzw.tum.de'
                }
            }
        }
        }
} 
```

Key parameters (all optional):

| Parameter | Default | Description |
|---|---|---|
| `--ppis` | `ppis.csv` | Input PPI CSV file |
| `--outdir` | `results` | Output directory |
| `--embedding_model` | `esm2` | Embedding model: `none` (one-hot), `esm2`, `prot_t5`, or path to a pre-computed `.npz` file |
| `--edge_weight` | `normalized_bitscore` | BLAST edge weight for the similarity graph: `bitscore` or `normalized_bitscore` |
| `--kahip_k` | `3` | Number of partitions (train / val / test) |
| `--kahip_seed` | `1234` | KaHIP random seed |
| `--kahip_preconfiguration` | `strong` | KaHIP mode: `strong`, `eco`, `fast`, `ultrafast` |
| `--cdhit_identity` | `0.4` | CD-HIT sequence identity threshold for redundancy removal |
| `--seed` | `42` | Random seed for negative sampling, classification, and bias analysis |

Example with custom parameters:

```bash
nextflow run main.nf \
    --ppis         string900_ppis.csv \
    --outdir       results/string900 \
    --embedding_model prot_t5 \
    --kahip_k      3
```

### 4. View the report

Open `results/multiqc/multiqc_report.html` in a browser.

---

## Workflow

![Pipeline overview](pipeline_overview.png)

```
ppis.csv
   │
   ├─ FETCH_DATA ──────────────────────────────────── sequences.fasta
   │       │                                           go_annotations.tsv
   │       │                                           species.tsv
   │       │
   ├─ GET_LENGTHS ──────────────────────────────────── lengths.tsv
   │
   ├─ RUN_BLAST ────────────────────────────────────── all_vs_all.tsv
   │
   ├─ MAKE_METIS ───────────────────────────────────── similarity.graph
   │                                                   node_mapping.tsv
   ├─ RUN_KAHIP ────────────────────────────────────── partitioned_proteome.txt
   │
   ├─ SORT_PPIS ────────────────────────────────────── train/val/test .csv + .fasta
   │
   ├─ CDHIT (train↔val, train↔test)
   │
   ├─ REMOVE_REDUNDANT ─────────────────────────────── train_nr/val_nr/test_nr .csv + .fasta
   │
   ├─ SAMPLE_NEGATIVES ─────────────────────────────── train/val/test_balanced/test_realistic .csv
   │
   ├─ EMBED_SEQUENCES ──────────────────────────────── embeddings.npz
   │
   ├─ TRAIN_CLASSIFIER ─────────────────────────────── classifier_metrics_mqc.tsv
   │
   ├─ BIAS_ANALYSIS (×6–7 attributes, parallel) ────── *_bias_mqc.tsv
   │
   ├─ COLLECT_BIAS ─────────────────────────────────── bias_scatter_mqc.html
   │
   ├─ SIMILARITY_HEATMAP ───────────────────────────── similarity_heatmap_mqc.html
   │
   └─ MULTIQC ──────────────────────────────────────── multiqc_report.html
```

### Step descriptions

**FETCH_DATA** — Queries UniProt for all proteins in the input CSV. Retrieves sequences (canonical + isoform-specific via the FASTA endpoint), GO annotations (biological process, molecular function, cellular component), and NCBI taxon IDs. Outputs `sequences.fasta`, `go_annotations.tsv`, and `species.tsv`.

**GET_LENGTHS** — Computes per-protein sequence lengths for length-normalised BLAST scores.

**RUN_BLAST** — Runs all-against-all BLASTp with `makeblastdb` + `blastp` to quantify pairwise sequence similarity.

**MAKE_METIS** — Converts the BLAST results into a weighted similarity graph in METIS format. Edge weights are either raw bitscore or bitscore normalised by the geometric mean of protein lengths.

**RUN_KAHIP** — Partitions the similarity graph into `k` parts using KaHIP's `kaffpa`. Proteins within the same partition are kept together; cross-partition PPIs are discarded. The largest partition becomes train, the second largest val, the smallest test.

**SORT_PPIS** — Assigns each PPI to a split based on the KaHIP partition. Writes per-split CSV and FASTA files.

**CDHIT** — Runs CD-HIT-2D between train↔val and train↔test to identify cross-split similar sequences.

**REMOVE_REDUNDANT** — Removes proteins from val and test that are too similar to any training protein (above the CD-HIT identity threshold).

**SAMPLE_NEGATIVES** — Samples random negative pairs for each split. Negatives are drawn such that each protein's degree distribution is approximately preserved. Produces a balanced test set (1:1 positive:negative) and a realistic test set (1:10 ratio).

**EMBED_SEQUENCES** — Computes per-protein embeddings using the selected model:
- `none` — 21-dimensional mean-pooled one-hot amino acid composition
- `esm2` — ESM-2 650M (dimension 1280), mean-pooled over residues
- `prot_t5` — ProtT5-XL (dimension 1024), mean-pooled over residues
- A path to a pre-computed `.npz` file skips this step entirely.

**TRAIN_CLASSIFIER** — Trains a Random Forest classifier on concatenated pair embeddings. Hyperparameters are tuned on the validation AUROC over 3 configurations (max_depth 5/10/30, max_samples 0.2), then the best model is retrained on train+val and evaluated on the balanced and realistic test sets.

**BIAS_ANALYSIS** — Runs in parallel for each attribute, computing:
- *Utility* — NMI(A; Y) = MI / √(H(A)·H(Y)): how much the attribute is correlated with the PPI label
- *Detectability* — Spearman ρ of a Ridge regressor predicting the attribute from pair embeddings

Attributes analysed:
| Attribute | Description |
|---|---|
| `sequence_similarity` | BLASTp pident between the two proteins, normalised to [0, 1] |
| `embedding_similarity` | Cosine similarity of the two individual protein embeddings |
| `functional_relatedness_BP/MF/CC` | Jaccard similarity of GO term sets (biological process / molecular function / cellular component) |
| `self_interactions` | 1 if both proteins are identical, 0 otherwise |
| `same_species` | 1 if both proteins share the same NCBI taxon ID, 0 otherwise (only included if the dataset contains proteins from more than one species) |

**COLLECT_BIAS** — Aggregates all per-attribute TSVs into a single interactive Plotly scatter plot (NMI vs detectability, coloured by attribute, shaped by split).

**SIMILARITY_HEATMAP** — Plots a heatmap of pairwise BLASTp similarity between proteins in different splits, to visualise the degree of leakage.

**MULTIQC** — Collects all `*_mqc.tsv` and `*_mqc.html` files into a single MultiQC report.

---

## Outputs

```
results/
├── multiqc/
│   └── multiqc_report.html       # Main report
├── data/
│   └── embeddings.npz            # Pre-computed embeddings (reusable)
├── train.csv                     # Final labelled splits (positives + negatives)
├── val.csv
├── test_balanced.csv
└── test_realistic.csv
```

---

## Standalone STRING channel analysis

To investigate which STRING evidence channels explain classifier performance differences between datasets, use the standalone script (not part of the Nextflow pipeline):

```bash
python bin/analyse_string_channels.py \
    --train      results/train.csv \
    --test       results/test_balanced.csv \
    --embeddings results/data/embeddings.npz \
    --out        string_channel_analysis.tsv
```

This fits a Ridge regressor (on positive pairs only) to predict each STRING evidence channel score from pair embeddings, and reports train and test Spearman ρ per channel. `combined_score` is excluded since it is derived from the individual channels.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.10
- Conda (for the environment) — or install the packages in `environment.yml` manually
- Internet access for the initial UniProt fetch (subsequent runs use cached Nextflow work directories)
- A GPU is recommended but not required for `esm2` and `prot_t5` embedding models
