![Logo](logo.png)

Automated leakage-aware splitting of a protein–protein interaction (PPI) dataset into train, validation, and test sets, with redundancy removal, negative sampling, embedding-based classification, and bias analysis.

Have a look at the [Wiki](https://github.com/bionetslab/ppi-splitting-pipeline/wiki) for more information.


![Pipeline overview](metro_map.svg)

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 26
- Conda (for the environment) — or install the packages in `environment.yml` manually
- Internet access for the initial UniProt fetch (subsequent runs use cached Nextflow work directories)
- A GPU is recommended but not required for `esm2` and `prot_t5` embedding models

---

## Quick Start

### Input PPI File

To custom-split your PPI dataset, you need to provide it to the pipeline as a CSV file with at least two columns (`protein1`, `protein2`) containing UniProt accession IDs. Additional columns (e.g., STRING evidence scores) are preserved throughout the pipeline.

```
protein1,protein2
P45985,Q14315
Q86TC9,P35609
O14836-2,P12345
...
```

### Samplesheet preparation

You provide all parameters for the pipeline via a samplesheet CSV where one row corresponds to one run. E.g., :

| id        | ppis             | split_method | negative_sampling_method | neg_ilp_solver | gurobi_license     |
|-----------|------------------|--------------|--------------------------|----------------|--------------------|
| fast-run  | data/my_ppis.csv | kahip        | default                  |                |                    |
| split-ilp | data/my_ppis.csv | ilp          | default                  |                |                    |
| all-ilp   | data/my_ppis.csv | ilp          | ilp                      | gurobi         | path/to/gurobi.lic |

### Run the pipeline

If you have a GPU, `-profile gpu` will submit the embedding step to a GPU, as specified by your nextflow config.

```
nextflow run main.nf --samplesheet samplesheet.csv --outdir results -profile gpu -c my_config.config
```

In `my_config.config`, you have to specify where the GPU is located. SLURM example:

```
profiles {
    ...
    gpu {
        process {
            withLabel:process_gpu {
                queue = 'shared-gpu'
                clusterOptions = '--qos=limitgpus --gpus=a40:1 --exclude small-gpu'
            }
        }
    }
}
```

### Seeds and robustness test

Each source of randomness has its own seed:

| Parameter         | Default | Controls                                                                  |
|-------------------|---------|---------------------------------------------------------------------------|
| `--split_seed`    | 1234    | KaHIP partitioning, the split ILP solver, and `split_method=random`       |
| `--sampling_seed` | 42      | Negative sampling (degree-preserving, uniform, and the ILP sampler)       |
| `--training_seed` | 42      | Baseline RandomForest classifier                                          |
| `--bias_seed`     | 42      | Bias diagnostics and heatmap subsampling (analysis only)                  |

`--robustness_test` runs every combination of three hard-coded seeds per source
(split: 1234, 1, 2; sampling: 42, 1, 2; training: 42, 1, 2), i.e. 3 partitions, 9 negative
samplings and 27 classifier trainings per dataset. BLAST, METIS and the embeddings are computed
only once. The first seed of each list is the default, so one combination reproduces a normal run.
Bias diagnostics and QC are skipped; the MultiQC report instead shows a per-dataset table with the
mean and standard deviation of every classifier metric, split into the contribution of each seed
source. Outputs land in `results/<id>/robustness/split<S>/sampling<N>/training<T>/`, with the
summary in `results/<id>/robustness/robustness_runs.tsv` and `robustness_summary_mqc.tsv`.
`--robustness_test` can't be combined with `--split_only` or `--bias_only`.

```
nextflow run main.nf --samplesheet samplesheet.csv --outdir results_robustness --robustness_test
```

### View the report

The MultiQC report can be found at `results/multiqc/multiqc_report.html`, which you can view in a browser.

