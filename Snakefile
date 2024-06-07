import pandas

os.makedirs('slurm', exist_ok=True)

MONKEY = 'Chlorocebus_sabaeus.ChlSab1.1'
TGONDII = 'ToxoDB-66_TgondiiGT1'
TGMITO = os.path.basename(config['tgmito'])
#MITO_GENES = os.path.basename(config['mito_genes'])

REF = ['tgmito']

PID = ['0.90', '0.95', '0.98']

ss = pandas.read_csv(config['ss'], sep='\t', comment='#')
assert len(list(set(ss.library_id)))

mtrna = []
with open(config['tgmito']) as fin:
    for line in fin:
        if line.strip().startswith('>'):
            ctg = line.strip().lstrip('>').strip().split()[0]
            mtrna.append(ctg)

wildcard_constraints:
    TGONDII = re.escape(TGONDII),
    MONKEY = re.escape(MONKEY),
    library_id = '|'.join([re.escape(str(x)) for x in ss.library_id]),
    mtrna = '|'.join([re.escape(x) for x in mtrna]),
    r = '|'.join([re.escape(str(x)) for x in ['1', '2', 'mrg']]),
    pid = '|'.join([str(x) for x in PID])

rule all:
    input:
        # expand('vsearch/reads.{pid}.uc', pid=['0.90', '0.95', '0.98']),
        # expand('bowtie2/genomes/unmapped/{library_id}.fa.gz', library_id=ss.library_id),
        # expand('bowtie2/{library_id}.align_intersect.tsv', library_id=ss.library_id),
        # expand('consensus/{mtrna}/{library_id}.cons.fa', mtrna=mtrna, library_id=ss.library_id),
        # expand('vsearch/pid_{pid}/clusters_tgmito.tsv', pid=PID),
        expand('bowtie2/{ref}.align_summary.tsv', ref=REF),
        expand('bowtie2/{ref}.align_summary.pdf', ref=REF),
        os.path.join(workflow.basedir, 'results/polyA.pdf'),
        # expand('bwa/tgmito/{library_id}.bam', library_id=ss.library_id),
        # expand('bwa/genomes/{library_id}.bam', library_id=ss.library_id),
        #expand('bwa/{library_id}.bam', library_id=ss.library_id),
        #expand('bwa/{library_id}.bam.bai', library_id=ss.library_id),
        #expand('hisat2/{library_id}.bam', library_id=ss.library_id),
        #expand('hisat2/{library_id}.bam.bai', library_id=ss.library_id),
        # expand('mergepairs/{library_id}.fastq', library_id=ss.library_id, r=['1', '2']),


rule download_monkey:
    output:
        gff=f'ref/{MONKEY}.gff3',
        fa=f'ref/{MONKEY}.dna_sm.toplevel.fa',
        fai=f'ref/{MONKEY}.dna_sm.toplevel.fa.fai',
    shell:
        r"""
        curl -L https://ftp.ensembl.org/pub/release-110/gff3/chlorocebus_sabaeus/{MONKEY}.110.gff3.gz \
        | gunzip > {output.gff}

        curl -L https://ftp.ensembl.org/pub/release-110/fasta/chlorocebus_sabaeus/dna/{MONKEY}.dna_sm.toplevel.fa.gz \
        | gunzip > {output.fa}
        samtools faidx {output.fa}
        """


rule download_tgondii:
    output:
        gff=f'ref/{TGONDII}.gff',
        fa=f'ref/{TGONDII}.fasta',
        fai=f'ref/{TGONDII}.fasta.fai',
    shell:
        r"""
        curl -L https://toxodb.org/common/downloads/release-66/TgondiiGT1/gff/data/{TGONDII}.gff > {output.gff}
        curl -L https://toxodb.org/common/downloads/release-66/TgondiiGT1/fasta/data/{TGONDII}_Genome.fasta > {output.fa}
        samtools faidx {output.fa}
        """


rule concat_genomes:
    input:
        fa=['ref/{TGONDII}.fasta', 'ref/{MONKEY}.dna_sm.toplevel.fa'],
    output:
        fa='ref/{TGONDII}_{MONKEY}.fa',
    shell:
        r"""
        cat {input.fa} > {output.fa}
        """


rule trim_reads:
    input:
        fastq_r1=lambda wc: ss[ss.library_id == wc.library_id].fastq_r1.iloc[0],
        fastq_r2=lambda wc: ss[ss.library_id == wc.library_id].fastq_r2.iloc[0],
    output:
        fastq_r1='cutadapt/{library_id}_1.fastq.gz',
        fastq_r2='cutadapt/{library_id}_2.fastq.gz',
        xlog='cutadapt/{library_id}.log',
    resources:
        cpus_per_task='4',
    shell:
        r"""
        cutadapt --quality-cutoff 0 --minimum-length 5 --cores {resources.cpus_per_task} \
            -a AGATCGGAAGAGC -A GATCGTCGGACT -o {output.fastq_r1} -p {output.fastq_r2} {input.fastq_r1} {input.fastq_r2} > {output.xlog}
        """


rule count_polya:
    input:
        bam='bowtie2/tgmito/{library_id}_2.bam',
    output:
        cnt='bowtie2/tgmito/polyA/{library_id}.cnt',
    run:
        import pysam

        samfile = pysam.AlignmentFile(input.bam, "rb")

        counts = {}
        for rec in samfile:
            if rec.is_unmapped or rec.is_reverse or rec.is_read2:
                continue
            if rec.reference_name not in counts:
                counts[rec.reference_name] = {}
            polylen = len(rec.seq) - len(rec.seq.rstrip('A'))
            if polylen not in counts[rec.reference_name]:
                counts[rec.reference_name][polylen] = 0
            counts[rec.reference_name][polylen] += 1
        
        with open(output.cnt, 'w') as fout:
            fout.write('\t'.join(['library_id', 'reference_name', 'polya_length', 'count']) + '\n')
            for ref in counts.keys():
                polylen = list(counts[ref].keys())
                polylen.sort()
                for l in polylen:
                    out = [wildcards.library_id, ref, l, counts[ref][l]]
                    out = [str(x) for x in out]
                    fout.write('\t'.join(out) + '\n')


rule plot_polya:
    input:
        cnt=expand('bowtie2/tgmito/polyA/{library_id}.cnt', library_id=ss.library_id),
        ss=config['ss'],
        rrna=config['tgmito'],
    output:
        ridge=os.path.join(workflow.basedir, 'results/polyA.pdf'),
        tsv='bowtie2/tgmito//polyA/polyA.tsv.gz',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)
library(Biostrings)
library(ggridges)

ss <- fread('{input.ss}')
fasta <- readDNAStringSet('{input.rrna}')
fasta <- as.character(fasta)
fasta <- data.table(reference_name= sub(' .*', '', names(fasta)), seq=fasta)
fasta[, polyA := nchar(seq) - nchar(sub('A*$', '', seq))]

ff <- strsplit('{input.cnt}', split=' ')[[1]]

cnt <- list()
for (x in ff) {{
    fin <- fread(x)
    fin[, cpm := 1000000 * (count/sum(count))]
    cnt[[length(cnt) + 1]] <- fin
}}
cnt <- rbindlist(cnt)
cnt <- merge(cnt, ss[, list(library_id, type, line)], by='library_id')
cnt <- merge(cnt, fasta[, list(reference_name, polyA)], by='reference_name')
cnt[, polya_pct := 100 * (count / sum(count)), list(reference_name, library_id)]

fwrite(cnt, '{output.tsv}', sep='\t')

keep <- cnt[, list(size=sum(cpm)/length(unique(library_id))), reference_name][size > 2000]
dat <- merge(cnt, keep, by='reference_name')
dat <- dat[order(library_id, reference_name, polya_length)][, list(polya_length, cumpct=cumsum(polya_pct), polya_pct, size, polyA, type, line), by=list(library_id, reference_name)]

dat[, title := sprintf('%s | size=%.1f', gsub('toxo_|rRNA_', '', reference_name), size)]
xord <- unique(dat[order(-polyA, title), list(polyA, title)])$title
dat[, title := factor(title, xord)]
dat[, type := factor(type, c('control', 'PolyA-KO', 'Untreated', 'pull_down'))]
dat[, line := factor(line, c("SDHB", "57-WC", "PAP1", 'L7'))]
xord <- unique(dat[order(line, type, library_id)]$library_id)
dat[, library_id := factor(library_id, rev(xord))]

gg <- ggplot(data=dat[library_id %in% c('57-WC-noATc', '57-WC-48ATc', 'L7-12-01', 'L7-12-02') & cumpct < 99], aes(x=polya_length, y=type, colour=type)) +
    geom_density_ridges(aes(height=polya_pct, group=library_id), stat = "identity", scale = 1, fill='grey80') +
    geom_vline(data=unique(dat[, list(polyA, title)]), aes(xintercept=polyA), colour='grey30', linetype='dashed') +
    # scale_size(range=c(0, 6)) +
    facet_wrap(~title, ncol=6) +
    scale_y_discrete(expand = c(0, 0)) +
    xlab('PolyA length') +
    ylab('% reads') +
    theme_light() +
    theme(strip.text=element_text(colour='grey30')) +
    theme(legend.position="none")
ggsave('{output.ridge}', width=40, height=30, units='cm')
EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule fastqc_trimmed:
    input:
        fastq='cutadapt/{library_id}_{r}.fastq.gz',
    output:
        qc='cutadapt/fastqc/{library_id}_{r}_fastqc.html',
    shell:
        r"""
        fastqc --outdir cutadapt/fastqc {input.fastq}
        """


rule merge_pairs:
    input:
        fastq_r1='cutadapt/{library_id}_1.fastq.gz',
        fastq_r2='cutadapt/{library_id}_2.fastq.gz',
    output:
        fastq_mrg='mergepairs/{library_id}_mrg.fastq',
        fastq_r1='mergepairs/{library_id}_1.fastq',
        fastq_r2='mergepairs/{library_id}_2.fastq',
        xlog='mergepairs/{library_id}.log',
    shell:
        r"""
        vsearch --fastq_mergepairs {input.fastq_r1} --reverse {input.fastq_r2} \
            --fastqout {output.fastq_mrg} \
            --fastqout_notmerged_fwd {output.fastq_r1} \
            --fastqout_notmerged_rev {output.fastq_r2} 2> {output.xlog}
        """


# rule reads_to_fasta:
#     input:
#         fastq='mergepairs/{library_id}_{r}.fastq',
#     output:
#         fa='mergepairs/fasta/{library_id}_{r}.fa',
#     shell:
#         r"""
#         cat {input.fastq} \
#         | paste - - - - \
#         | awk '{{print ">" $1 "\n" $2}}' > {output.fa}
#         """


rule bowtie_index_tgmito:
    input:
        fa=config['tgmito'],
    output:
        fa=f'ref/{TGMITO}',
        idx=f'ref/{TGMITO}.rev.2.bt2',
    shell:
        r"""
        cp {input.fa} {output.fa}
        bowtie2-build {output.fa} {output.fa} 
        """


rule bowtie_index_tgmitogenome_marialuisa:
    input:
        fa='ref/TgMitogenome_MariaLuisa.fa',
    output:
        idx='ref/TgMitogenome_MariaLuisa.fa.rev.2.bt2',
    shell:
        r"""
        bowtie2-build --seed 1234 {input.fa} {input.fa} 
        """


rule bowtie_index_genomes:
    input:
        fa='ref/{TGONDII}_{MONKEY}.fa',
    output:
        idx='ref/{TGONDII}_{MONKEY}.fa.rev.2.bt2',
    resources:
        mem='10G',
        cpus_per_task='8'
    shell:
        r"""
        bowtie2-build --threads {resources.cpus_per_task} --seed 1234 {input.fa} {input.fa} 
        """


# rule bowtie_index_mito_genes:
#     input:
#         fa=config['mito_genes'],
#     output:
#         fa=f'ref/{MITO_GENES}',
#         idx=f'ref/{MITO_GENES}.rev.2.bt2',
#     shell:
#         r"""
#         cp {input.fa} {output.fa}
#         bowtie2-build {output.fa} {output.fa} 
#         """

rule bowtie_read2:
    input:
        fastq_r2='cutadapt/{library_id}_2.fastq.gz',
        idx=f'ref/{TGMITO}.rev.2.bt2',
    output:
        bam= 'bowtie2/{ref}/{library_id}_2.bam',
        smry= 'bowtie2/{ref}/{library_id}_2.log',
    params:
        idx_prefix=lambda wc, input: re.sub('.rev.2.bt2$', '', input.idx),
    resources:
        mem='20G',
        cpus_per_task='32'
    shell:
        # We don;t really to need to revcomp, but it makes it easier to work with
        r"""
        vsearch --fastx_revcomp {input.fastq_r2} --fastqout - \
        | bowtie2 -x {params.idx_prefix} -U - -N 1 -L 8 --local -p {resources.cpus_per_task} --score-min L,5,1 2> {output.smry} \
        | samtools sort -@ 16 > {output.bam}
        """


rule bowtie_align:
    input:
        fastq_mrg='mergepairs/{library_id}_mrg.fastq',
        fastq_r1='mergepairs/{library_id}_1.fastq',
        idx=f'ref/{TGMITO}.rev.2.bt2',
        # idx=lambda wc: f'ref/{TGONDII}_{MONKEY}.fa.rev.2.bt2' if wc.ref == 'genomes' else 'ref/TgMitogenome_MariaLuisa.fa.rev.2.bt2' if wc.ref == 'TgMitogenome_MariaLuisa' else f'ref/{TGMITO}.rev.2.bt2' if wc.ref == 'tgmito' else f'ref/{MITO_GENES}.rev.2.bt2' if wc.ref == 'mito_genes' else 'NA',
    output:
        bam= 'bowtie2/{ref}/{library_id}.bam',
        smry= 'bowtie2/{ref}/{library_id}.log',
    params:
        idx_prefix=lambda wc, input: re.sub('.rev.2.bt2$', '', input.idx),
    resources:
        mem='10G',
        cpus_per_task='16'
    shell:
        r"""
        bowtie2 -x {params.idx_prefix} -U {input.fastq_mrg} {input.fastq_r1} -N 1 -L 8 --local -p {resources.cpus_per_task} --score-min L,5,1 2> {output.smry} \
        | samtools sort -@ 4 > {output.bam}
        """


rule bowtie_align_index:
    input:
        bam= 'bowtie2/{ref}/{library_id}.bam',
    output:
        bai= 'bowtie2/{ref}/{library_id}.bam.bai',
    shell:
        r"""
        samtools index {input.bam}
        """


rule bowtie_summary:
    input:
        bam= expand('bowtie2/{{ref}}/{library_id}.bam', library_id=ss.library_id),
        bai= expand('bowtie2/{{ref}}/{library_id}.bam.bai', library_id=ss.library_id),
    output:
        smry='bowtie2/{ref}.align_summary.tsv',
        plot='bowtie2/{ref}.align_summary.pdf',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)

bams <- strsplit('{input.bam}', ' ')[[1]]

smry <- list()
for(bam in bams) {{
    print(bam)
    r <- fread(cmd=sprintf('samtools idxstat %s', bam), header=FALSE, col.names=c('ref', 'length', 'count', 'other'))
    r[, library_id := sub('\\.bam$', '', basename(bam))]
    r[ref == '*']$count <- r[ref == '*']$other
    r[ref == '*']$length <- NA
    r[ref == '*']$ref <- 'unmapped'
    r[, other := NULL]
    r[, cpm := 1000000 * count / sum(count)]
    smry[[length(smry) + 1]] <- r
}}
smry <- rbindlist(smry)

dat <- dcast(data=smry, ref + length ~ library_id, value.var='cpm')
dat[, ref := gsub('.*_rRNA_', '', ref)]
dat <- melt(data=dat[, list(ref, `L7-12-01`, `L7-12-02`, `SDHB-NC`)], id.vars=c('ref'), variable.name='library_id', value.name='cpm')
xord <- dat[, list(cpm=sum(cpm)), ref][order(cpm)]$ref
dat[, ref := factor(ref, xord)]

gg <- ggplot(data=dat[ref != 'unmapped'], aes(x=ref, y=cpm, colour=library_id, shape=library_id)) +
    geom_point() +
    coord_flip() +
    xlab('')  +
    ylab('Abundance as counts per million') +
    scale_y_log10() +
    theme_light() 
ggsave('{output.plot}', width=16, height=18, units='cm')

smry[, cpm := sprintf('%.1f', cpm)]
out <- dcast(data=smry, ref + length ~ library_id, value.var='cpm')[order(length)]
write.table(out, '{output.smry}', sep='\t', row.names=FALSE, quote=FALSE)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule unmapped:
    input:
        bam='bowtie2/{ref}/{library_id}.bam',
        bai='bowtie2/{ref}/{library_id}.bam.bai',
    output:
        fa='bowtie2/{ref}/unmapped/{library_id}.fa.gz',
    shell:
        # NB: Don't use this for mapped reads. Use samtools fastq 
        r"""
        samtools view -f 4 {input.bam} \
        | cut -f 10 \
        | LC_ALL=C sort \
        | uniq -c \
        | awk '{{printf(">id%08d:count=%d:length=%d\t%s\t%s\n", NR, $1, length($2), $1, $2)}}' \
        | sort -k2,2nr \
        | awk '{{print $1 "\n" $3}}' \
        | gzip > {output.fa}
        """


rule align_intersect:
    input:
        bams=expand('bowtie2/{ref}/{{library_id}}.bam', ref=REF),
        bai=expand('bowtie2/{ref}/{{library_id}}.bam.bai', ref=REF),
    output:
        tsv='bowtie2/{library_id}.align_intersect.tsv',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)

bams <- strsplit('{input.bams}', ' ')[[1]]

reads <- list()
for(bam in bams) {{
    r <- fread(cmd=sprintf('samtools view -F4 %s | cut -f 1', bam), header=FALSE, col.names='rname')
    r[, ref := basename(dirname(bam))]
    reads[[length(reads) + 1]] <- r
}}
reads <- rbindlist(reads)
nreads <- reads[, list(.N, ref=paste(ref, collapse='|')), by='rname']
nreads <- nreads[, .N, ref]
nreads[, pct := 100 * (N/sum(N))]
nreads <- nreads[order(-pct)]
nreads[, pct := sprintf('%.1f', pct)]
fwrite(nreads, '{output.tsv}', sep='\t')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule cluster_tgmito:
    input:
        bam='bowtie2/tgmito/{library_id}.bam',
        bai='bowtie2/tgmito/{library_id}.bam.bai',
    output:
        cons='consensus/{mtrna}/{library_id}.cons.fa',
        msa='consensus/{mtrna}/{library_id}.msa.fa',
        uc='consensus/{mtrna}/{library_id}.uc.tsv',
        centr='consensus/{mtrna}/{library_id}.centroids.fa',
    shell:
        r"""
        samtools view -F 2304 -u {input.bam} {wildcards.mtrna} \
        | samtools fastq \
        | vsearch \
            --msaout {output.msa} \
            --fasta_width 0 \
            --centroids {output.centr} \
            --consout {output.cons} \
            --uc {output.uc} \
            --threads 8 \
            --clusterout_sort \
            --minseqlength 1 \
            --iddef 1 \
            --minwordmatches 10 \
            --wordlength 6 \
            --id 0.95 \
            --cluster_size -
        """



rule add_library_to_reads:
    input:
        fq='mergepairs/{library_id}_{r}.fastq',
    output:
        fa=temp('mergepairs/{library_id}_{r}.fa'),
    shell:
        r"""
        cat {input.fq} \
        | sed 's/ .*//' \
        | paste - - - - \
        | sed 's/^@/>/' \
        | awk -v id={wildcards.library_id} '{{print $1 "::" id "\n" $2}}' > {output.fa}
        """


rule cat_reads:
    input:
        fa=expand('mergepairs/{library_id}_{r}.fa', library_id=ss.library_id, r=['mrg', '1']),
    output:
        fa=temp('mergepairs/reads.fa'),
    shell:
        r"""
        cat {input.fa} > {output.fa}
        """


# rule vsearch_cluster_tgmito:
#     input:
#         fq='mergepairs/{library_id}_mrg.fastq',
#         tg=f'ref/{TGMITO}'
#     output:
#         msa='tgmito/{library_id}.{pid}.fa', 
#         centroids='tgmito/{library_id}.centroids.{pid}.fa', 
#         uc='tgmito/{library_id}.{pid}.uc', 
#         xlog='tgmito/{library_id}.{pid}.log',
#     resources:
#         cpus_per_task='4',
#         mem='4G',
#     shell:
#         r"""
#         cat {input.fq} {input.tg} \
#         | vsearch --msaout {output.msa} \
#             --fasta_width 0 \
#             --centroids {output.centroids} \
#             --uc {output.uc} \
#             --threads {resources.cpus_per_task} \
#             --clusterout_sort \
#             --minseqlength 1 \
#             --gapopen 0E/5I \
#             --gapext 0E/2I \
#             --mismatch -8 \
#             --match 6 \
#             --iddef 1 \
#             --minwordmatches 10 \
#             --wordlength 6 \
#             --id {wildcards.pid} \
#             --cluster_fast {input.fa} 2> {output.xlog}
#         """


rule vsearch_cluster_all:
    input:
        fa='mergepairs/reads.fa',
    output:
        msa='vsearch/reads.{pid}.fa', 
        centroids='vsearch/reads.centroids.{pid}.fa', 
        uc='vsearch/reads.{pid}.uc', 
        xlog='vsearch/reads.{pid}.log',
    resources:
        cpus_per_task='4',
        mem='10G',
    shell:
        r"""

        vsearch --msaout {output.msa} \
            --fasta_width 0 \
            --centroids {output.centroids} \
            --uc {output.uc} \
            --threads {resources.cpus_per_task} \
            --clusterout_sort \
            --minseqlength 1 \
            --gapopen 0E/5I \
            --gapext 0E/2I \
            --mismatch -8 \
            --match 6 \
            --iddef 1 \
            --minwordmatches 10 \
            --wordlength 6 \
            --id {wildcards.pid} \
            --cluster_fast {input.fa} 2> {output.xlog}
        """


rule cluster_summary:
    input:
        centroids='vsearch/reads.centroids.{pid}.fa', 
        uc='vsearch/reads.{pid}.uc', 
    output:
        clst=temp('vsearch/pid_{pid}/clusters.tsv'),
    resources:
        mem='8G',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R

library(data.table)

centroids <- fread(cmd="awk -v RS='>' '$1 != \"\" {{print $1, $2}}' {input.centroids}", col.names=c('centroid', 'sequence'), header=FALSE)
uc <- fread('{input.uc}', sep='\t', header=FALSE, col.names=c('type', 'seqlen', 'pct_id', 'strand', 'query', 'centroid'), select=c(1, 3, 4, 5, 9, 10))

uc[, library_id := sub('.*:', '', query)]
uc[, centroid := ifelse(type == 'S', query, centroid)]
uc[, pct_id := ifelse(query == centroid, 100, pct_id)]
uc <- uc[type != 'C']
uc[, type := NULL]
uc[, pct_id := as.numeric(pct_id)]
size <- uc[, list(cluster_size=.N), centroid]
uc <- merge(uc, size, by='centroid')

smry <- uc[, list(.N), list(library_id, centroid, cluster_size)][order(-cluster_size, centroid, library_id)]
libsize <- uc[, list(libsize=.N), library_id]
smry <- merge(smry, libsize, by='library_id', sort=FALSE)
smry[, rpm := (N/libsize)*1000000]
smry[, tot_rpm := sum(rpm), by=centroid]
smry[, pct_count := 100 * (rpm / tot_rpm), by=centroid]
smry[, libsize := NULL]
smry <- merge(smry, centroids)

smry[, pct_count := sprintf('%.2f', pct_count)]
ctrl_depl <- smry[tot_rpm > 10]
xdat <- dcast(data=smry[centroid %in% ctrl_depl$centroid], centroid + sequence + tot_rpm ~ library_id, value.var='pct_count')
xdat <- xdat[order(-tot_rpm)]
xdat[, tot_rpm := sprintf('%.1f', tot_rpm)]
xdat[, sequence_length := nchar(sequence)]
fwrite(xdat, '{output.clst}', sep='\t')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


# rule align_denovo_clusters:
#     input:
#         clst='vsearch/pid_{pid}/clusters.tsv',
#         fa=f'ref/{TGMITO}',
#     output:
#         tsv=temp('vsearch/pid_{pid}/clusters.tmp'),
#     shell:
#         r"""
#         awk 'NR > 1 {{print ">"$1 "\n" $2}}' {input.clst} \
#         | bowtie2 -f -x {input.fa} -U - -N 1 -L 8 --local --score-min L,5,0.2 \
#         | samtools view \
#         | python -c "
# import sys 
# print('centroid best_tgmito cigar n_mismatch'.replace(' ', '\t'))
# for line in sys.stdin:
#     line = line.split('\t')  
#     nm = [x for x in line if x.startswith('NM:i:')]
#     if len(nm) == 1:
#         nm = nm[0].replace('NM:i:', '')
#     elif len(nm) == 0:
#         nm = '-1'
#     else:
#         raise Exception()
#     line[2] = 'NA' if line[2] == '*' else line[2]
#     line[5] = 'NA' if line[5] == '*' else line[5]
#     print('\t'.join([line[0], line[2], line[5], nm])) 
#         " > {output.tsv}
#         """


# rule add_ref_to_clusters:
#     input:
#         tgmito='vsearch/pid_{pid}/clusters.tmp',
#         fa=f'ref/{TGMITO}',
#         clst='vsearch/pid_{pid}/clusters.tsv',
#     output:
#         clst='vsearch/pid_{pid}/clusters_tgmito.tsv',
#     shell:
#         r"""
# cat <<'EOF' > {rule}.$$.tmp.R
# library(data.table)
# 
# clst <- fread('{input.clst}')
# aln <- fread('{input.tgmito}')
# fa <- fread('{input.fa}', header=FALSE)
# 
# fx <- cbind(fa[grepl('>', V1), list(best_tgmito=sub('>', '', V1))], fa[!grepl('>', V1), list(tgmito_seq=V1)])
# stopifnot(nrow(fx) == nrow(fa)/2)
# 
# dat <- merge(merge(clst, aln, by='centroid', sort=FALSE), fx, by='best_tgmito', sort=FALSE, all=TRUE)
# dat <- dat[order(best_tgmito, -tot_rpm)]
# write.table(dat, '{output.clst}', sep='\t', quote=FALSE, row.names=FALSE)
# 
# EOF
# Rscript {rule}.$$.tmp.R
# rm {rule}.$$.tmp.R
#         """
