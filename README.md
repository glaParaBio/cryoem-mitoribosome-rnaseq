<!-- vim-markdown-toc GFM -->

* [Description](#description)
* [Set up](#set-up)
* [Run](#run)

<!-- vim-markdown-toc -->

# Description

This repository contains the code for the analysis of the sequencing data in Shikha
*et al.* (manuscript in preparation). 

# Set up

Set up the software dependencies using the package manager
[mamba](https://github.com/mamba-org/mamba) in a dedicated environment. See
[requirements.txt](requirements.txt) for the list of dependencies.

```
mamba create -n cryoem-mitoribosome-rnaseq
mamba activate cryoem-mitoribosome-rnaseq
mamba install --yes --file requirements.txt
```

Download fastq file from ENA (study PRJEB72258, currently private) and edit the
tab-separated file [sample_sheet.tsv](sample_sheet.tsv) to have the columns
`fastq_r1/2` listing the full path to the respective local fastq file.

# Run

The analysis steps are executed by snakemake. The option `--directory` sets the
output directory (edit as desired). Options `--cluster` and `--cluster-cancel`
are specific for the  execution on HPC with the job manager slurm and the
should be omitted for local execution.

Option `--dry-run` shows the commands that would be executed without doing
anything. Omit this option for actual execution. 

For details of the analysis see the [Snakefile](Snakefile).

```
snakemake -j 10 -p \
    --dry-run \
    --rerun-trigger mtime \
    --directory output \
    --config ss=$PWD/sample_sheet.tsv \
        tgmito=$PWD/rRNA_sequences_231221.fa \
        chains=$PWD/chains.tsv \
    --latency-wait 60 \
    --default-resources "mem='1G'" "cpus_per_task='1'" \
    --cluster 'sbatch --cpus-per-task={resources.cpus_per_task} --mem={resources.mem} --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err -A none' \
    --cluster-cancel scancel
```

The workflow on test data should complete in a a few minutes.

To report bugs and questions, please submit an [issue](https://github.com/glaParaBio/cryoem-mitoribosome-rnaseq/issues).
