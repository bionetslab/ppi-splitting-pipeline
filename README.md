![Logo](logo.png)

Automated leakage-aware splitting of a protein–protein interaction (PPI) or domain-domain interaction (DDI) dataset into train, validation, and test sets, with redundancy removal, negative sampling, embedding-based classification, and bias analysis.

Have a look at the [Wiki](https://github.com/bionetslab/ppi-splitting-pipeline/wiki) for more information.


![Pipeline overview](metro_map.svg)

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 26
- Conda (for the environment), install the packages in `environment.yml` manually, or use `-profile docker` and the image built from `docker/Dockerfile`
- Internet access for the initial UniProt/Pfam fetch (subsequent runs use cached Nextflow work directories)
- A GPU is recommended but not required for `esm2` and `prot_t5` embedding models. 

---

## Quick Start

### Input PPI File

To custom-split your PPI dataset, you need to provide it to the pipeline as a CSV file with at least two columns (`protein1`, `protein2`) containing UniProt accession IDs. Additional columns (e.g., STRING evidence scores) are preserved throughout the pipeline. In DDI mode the file has exactly the same shape, but the two columns hold Pfam family accessions; the column names stay the same. 

```
protein1,protein2
P45985,Q14315
Q86TC9,P35609
O14836-2,P12345
...
```

### Samplesheet preparation

You provide all parameters for the pipeline via a samplesheet CSV where one row corresponds to one run. E.g., :

| id        | ppis             | split_method | negative_sampling_method |
|-----------|------------------|--------------|--------------------------|
| fast-run  | data/my_ppis.csv | kahip        | default                  |
| split-ilp | data/my_ppis.csv | ilp          | default                  |
| all-ilp   | data/my_ppis.csv | ilp          | ilp                      |


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

### View the report

The MultiQC report can be found at `results/multiqc/multiqc_report.html`, which you can view in a browser.
