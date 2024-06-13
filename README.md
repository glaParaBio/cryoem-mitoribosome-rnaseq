<!-- vim-markdown-toc GFM -->

* [Set up](#set-up)
* [Run](#run)

<!-- vim-markdown-toc -->

# Set up

```
mamba create -n cryoem-mitoribosome-rnaseq
mamba activate cryoem-mitoribosome-rnaseq
mamba install --yes --file requirements.txt
```

# Run

```
snakemake -j 10 -p -n \
    --rerun-trigger mtime \
    --directory /users/db291g/sharedscratch/projects/cryoem-mitoribosome-rnaseq \
    --config ss=$PWD/sample_sheet.tsv \
        tgmito=$PWD/rRNA_sequences_231221.fa \
        chains=$PWD/chains.tsv \
    --latency-wait 60 \
    --default-resources "mem='1G'" "cpus_per_task='1'" \
    --cluster 'sbatch --cpus-per-task={resources.cpus_per_task} --mem={resources.mem} --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err -A none' \
    --cluster-cancel scancel
```

---

Dev note: code for this repository is originally from
https://github.com/glaParaBio/dariober/tree/2bc723dda9a154351ea1b53027066deb0fb4b1a0/cryoem-mitoribosome-rnaseq
