Timepoints = ['TYPD', 'T6h', 'T8h25']

S_minimap2 = '/nas/longleaf/home/wuhao/anaconda3/envs/nanopore/bin/minimap2' # minimap2=2.1.1-r341
S_chopper = '/nas/longleaf/home/wuhao/anaconda3/envs/nanopore/bin/chopper' # chopper=0.8.0
S_pychopper = '/nas/longleaf/home/wuhao/anaconda3/envs/nanopore/bin/pychopper' # pychopper==2.7.0 hmmer-3.4
S_samtools = '/nas/longleaf/rhel8/apps/samtools/1.21/bin/samtools' # samtools=1.21
S_spliced_bam2gff = '/nas/longleaf/home/wuhao/anaconda3/envs/nanopore/bin/spliced_bam2gff' # pinfish-200316-2bc39d4_1
S_cluster_gff = '/nas/longleaf/home/wuhao/anaconda3/envs/nanopore/bin/cluster_gff'
S_polish_clusters = '/nas/longleaf/home/wuhao/anaconda3/envs/nanopore/bin/polish_clusters'
S_collapse_partials = '/nas/longleaf/home/wuhao/anaconda3/envs/nanopore/bin/collapse_partials'
S_stringtie = '/nas/longleaf/home/wuhao/anaconda3/envs/nanopore/bin/stringtie' # S_stringtie=2.2.3
S_gffcompare = '/nas/longleaf/home/wuhao/anaconda3/envs/nanopore/bin/gffcompare' # gffcompare=0.12.6
S_featureCounts = '/nas/longleaf/home/wuhao/Packages/subread-2.0.6-Linux-x86_64/bin/featureCounts' # featureCounts=2.0.6
# S_gffcompare = '/nas/longleaf/home/wuhao/anaconda3/envs/nanopore/gffcompare/gffcompare/gffcompare' # gffcompare v0.12.9
R_fa = '/nas/longleaf/home/wuhao/Reference/sk1_revised_100331_filtered.fasta'
R_gtf = '/nas/longleaf/home/wuhao/Reference/sk1_gene.gtf'
R_gff = '/work/users/w/u/wuhao/Nanopore_antisense/sk1_gene.gff'
wd = '/work/users/w/u/wuhao/Nanopore_antisense/'

import pandas as pd
import csv

def filterGtfTPM2(file_in, file_out):
    col=['seqname', 'source', 'feature', 'start', 'end', 
    'score', 'strand', 'frame', 'attribute']
    gtf=pd.read_table(file_in, comment='#', header=None,names=col)
    gtf['TPM'] = gtf['attribute'].str.extract(r'TPM "([^"]+)"')
    gtf['TPM'] = gtf['TPM'].astype(float)

    gtf['cov'] = gtf['attribute'].str.extract(r'cov "([^"]+)"')
    gtf['cov'] = gtf['cov'].astype(float)

    gtf['transcript_id'] = gtf['attribute'].str.extract(r'transcript_id "([^"]+)"')
    gtf=gtf[gtf.transcript_id.isin(gtf.loc[gtf['TPM']>=2,'transcript_id'])]
    gtf[col].to_csv(file_out, header=None,index=False,sep='\t',
                    quotechar='"', quoting=csv.QUOTE_NONE, escapechar='\\')

rule all:
    input:
        expand(wd + '1_chopper/{timepoint}_chopped.fq.gz', timepoint = Timepoints),
        expand(wd + '2_pychopper/{timepoint}_full_length.fq', timepoint = Timepoints),
        expand(wd + '2_pychopper/{timepoint}_chopped_full_length.fq.gz', timepoint = Timepoints),
        expand(wd + '3_bam/{timepoint}.sam', timepoint = Timepoints),
        expand(wd + '3_bam/{timepoint}_samstat.txt', timepoint = Timepoints),
        expand(wd + '3_bam/{timepoint}_sorted.bam', timepoint = Timepoints),
        # expand(wd + '4_pinfish/{timepoint}_raw_transcripts.gff', timepoint = Timepoints),
        # expand(wd + '4_pinfish/{timepoint}_clustered_transcripts.gff', timepoint = Timepoints),  
        # expand(wd + '4_pinfish/{timepoint}_consensus.fas', timepoint = Timepoints),
        # expand(wd + '4_pinfish/{timepoint}_collapsed.gff', timepoint = Timepoints),
        # expand(wd + '5_stringtie/{timepoint}_stringtie_TPM2.gtf',  timepoint = Timepoints),
        # wd + '5_stringtie/A14201_ONT_sense_count.table',
        # wd + '5_stringtie/A14201_ONT_anti_count.table',
        # wd + '5_stringtie/Merge_stringtie.gtf'
        # wd + '6_stringtieMerge/Merge_stringtie.gtf',
        expand(wd + '3_bam/{timepoint}_split.fwd.bam',timepoint = Timepoints),
        expand(wd + '3_bam/{timepoint}_split.rev.bam',timepoint = Timepoints),
        expand(wd + '3_bam/{timepoint}_start_sites.bed',timepoint = Timepoints),

rule chopper:
    input:
        wd + 'fastq/{timepoint}_pass.fq.gz'
    output:
        fq_gz=wd + '1_chopper/{timepoint}_chopped.fq.gz'
    log:
        wd + '1_chopper/logs/{timepoint}_chopper.log'
    shell:
        """
        gunzip -c {input} | {S_chopper} -q 7 -l 50 | gzip > {output.fq_gz} 2> {log}
        """

rule pychopper:
    input:
        fq_gz=wd + '1_chopper/{timepoint}_chopped.fq.gz'
    output:
        report=wd + '2_pychopper/{timepoint}_report.pdf',
        unclassified=wd + '2_pychopper/{timepoint}_unclassified.fq',
        rescued=wd + '2_pychopper/{timepoint}_rescued.fq',
        full_length=wd + '2_pychopper/{timepoint}_full_length.fq',
        log=wd + '2_pychopper/{timepoint}_pychopper.log'
        
    shell:
        """
        {S_pychopper} -t 16 -r {output.report} -k PCS109  \
        -Y 10000 -u {output.unclassified} -w {output.rescued} \
        -B 100000  {input[0]} {output.full_length} \
        > {output.log} 
        """

rule chopper2:
    input:
        wd + '2_pychopper/{timepoint}_full_length.fq'
    output:
        fq_gz=wd + '2_pychopper/{timepoint}_chopped_full_length.fq.gz'
    log:
        wd + '2_pychopper/logs/{timepoint}_chopper.log'
    shell:
        """
        {S_chopper} -q 7 -l 50 -i {input} | gzip > {output.fq_gz} 2> {log}
        """

rule minimap2:
    input:
        wd + '2_pychopper/{timepoint}_chopped_full_length.fq.gz'
    output:
        wd + '3_bam/{timepoint}.sam' 
    shell:
        '{S_minimap2} \
        -ax splice \
        -uf -k14 {R_fa} \
        {input[0]} > {output[0]}'
        
rule samstat:
    input:
        wd + '3_bam/{timepoint}.sam'
    output:
        wd + '3_bam/{timepoint}_samstat.txt' 
    shell:
        '{S_samtools} flagstat {input[0]} > {output[0]}'

rule bamsort:
    input:
        wd + '3_bam/{timepoint}.sam'
    output:
        wd + '3_bam/{timepoint}_sorted.bam' 
    shell:
        """
        {S_samtools} view -buS {input[0]} | {S_samtools} sort -o {output[0]} 
        {S_samtools} index {output[0]}
        """

# rule bam2gff:
#     input:
#         wd + '3_bam/{timepoint}_sorted.bam'
#     output:
#         wd + '4_pinfish/{timepoint}_raw_transcripts.gff' 
#     shell:
#         """
#         {S_spliced_bam2gff} -M {input[0]} > {output[0]}
        # """

# rule cluster_gff:
#     input:
#         wd + '4_pinfish/{timepoint}_raw_transcripts.gff'
#     output:
#         gff=wd + '4_pinfish/{timepoint}_clustered_transcripts.gff',
#         tsv=wd + '4_pinfish/{timepoint}_clustered_transcripts.tsv'
#     shell:
#         """
#         {S_cluster_gff} -t 16 -a {output.tsv} {input[0]} > {output.gff}
#         """

# rule polish_clusters:
#     input:
#         bam=wd + '3_bam/{timepoint}_sorted.bam',
#         tsv=wd + '4_pinfish/{timepoint}_clustered_transcripts.tsv'
#     output:
#         fas=wd + '4_pinfish/{timepoint}_consensus.fas'
#     shell:
#         """
#         {S_polish_clusters} -t 16 -m \
#         -a {input.tsv} \
#         -o {output.fas} {input.bam}
        # """

# rule collapse_partials:
#     input:
#         gff=wd + '4_pinfish/{timepoint}_clustered_transcripts.gff'
#     output:
#         gff=wd + '4_pinfish/{timepoint}_collapsed.gff'
#     shell:
#         """
#         {S_collapse_partials} -t 16 {input.gff} > {output.gff}
#         """

# rule map_consensus:
#     input:
#         fas=wd + '4_pinfish/{timepoint}_consensus.fas'
#     output:
#         wd + '5_consensus/{timepoint}_consensus.bam', 
#     shell:
#         """
#         {S_minimap2} -ax splice -uf -k14 {R_fa} {input.fas} > {output[0]}
#         {S_samtools} view -buS {input[0]} | {S_samtools} sort -o {output[0]} 
#         {S_samtools} index {output[0]}
#         """

# rule stringtie:
#     input:
#         wd + '3_bam/{timepoint}_sorted.bam',
#         gtf={R_gff}
#     output:
#         wd + '5_stringtie/{timepoint}_stringtie.gtf',
#     params:
#         '{timepoint}'
#     shell:
#         """
#         {S_stringtie} --conservative -j 10 -l {params} -L -R -o {output[0]} {input[0]}
#         """     

# rule filterGtf:
#     input:
#         wd + '5_stringtie/{timepoint}_stringtie.gtf'
#     output:
#         wd + '5_stringtie/{timepoint}_stringtie_TPM2.gtf'
#     run:
#         filterGtfTPM2(input[0], output[0])

# rule merge_sort_gtf:
#     input:
#         expand(wd + '5_stringtie/{timepoint}_stringtie.gtf', timepoint = Timepoints)
#     params:
#         ' '.join(expand(wd + '5_stringtie/{timepoint}_stringtie.gtf', timepoint = Timepoints))
#     output:
#         wd + '5_stringtie/Merge_stringtie.gtf',
#         wd + '5_stringtie/Merge_stringtie.sorted.gtf'
#     shell:
#         """
#         cat {params} > {output[0]}
#         """
        # sort -k1,1 -k4,4n -k5,5n {output[0]} > {output[1]}
        


# rule featureCounts:
#     input:
#         expand(wd + '3_bam/{timepoint}_sorted.bam', timepoint = Timepoints)
#     params:
#         ' '.join(expand(wd + '3_bam/{timepoint}_sorted.bam', timepoint = Timepoints))
#     output:
#         sense=wd + '5_stringtie/A14201_ONT_sense_count.table',
#         anti=wd + '5_stringtie/A14201_ONT_anti_count.table'
#     shell:
#         """
#         featureCounts --minOverlap 40 -T 12 -s 1 -O -L -a {R_gtf}  -o {output.sense} {params[0]}
#         featureCounts --minOverlap 40 -T 12 -s 2 -O -L -a {R_gtf}  -o {output.anti} {params[0]}
#         """

# rule collapse_partials2:
#     input:
#         gff=wd + '4_pinfish/{timepoint}_clustered_transcripts.gff'
#     output:
#         gff=wd + '4_pinfish/{timepoint}_collapsed.gff'
#     shell:
#         """
#         {S_collapse_partials} -t 16 {input.gff} > {output.gff}
#         """


# rule stringtieMerge:
#     input:
#         expand(wd + '5_stringtie/{timepoint}_stringtie.gtf', timepoint = Timepoints)
#     params:
#         ' '.join(expand(wd + '5_stringtie/{timepoint}_stringtie.gtf', timepoint = Timepoints))
#     output:
#         wd + '5_stringtie/Merge_stringtie.gtf'
#     shell:
#         """
#         {S_stringtie} --merge -o {output} -m 50 -T 10 {params}
#         """

# rule gffcompare:
#     input:
#         expand(wd + '5_stringtie/{timepoint}_stringtie_TPM2.gtf', timepoint = Timepoints),
#         gtf = {R_gff}
#     params:
#         ' '.join(expand(wd + '5_stringtie/{timepoint}_stringtie_TPM2.gtf', timepoint = Timepoints))
#     output:
#         wd + '6_stringtieMerge/Merge_stringtie.gtf'
#     shell:
#         """
#         {S_gffcompare} -r {input.gtf} -p merge -R -D -S -d 50  -o {output} {params}
#         """


rule splitbam:
    input:
        wd + '3_bam/{timepoint}_sorted.bam'
    output:
        fwd = wd + '3_bam/{timepoint}_split.fwd.bam',
        rev = wd + '3_bam/{timepoint}_split.rev.bam'
    shell:
        """
        {S_samtools} view -h -F 0x10 {input[0]} | {S_samtools} sort -o {output.fwd}
        {S_samtools} view -h -f 0x10 {input[0]} | {S_samtools} sort -o {output.rev}
        {S_samtools} index {output.fwd}
        {S_samtools} index {output.rev}
        """
        
rule fetch_start_sites:	
    input:
        bam=wd + '3_bam/{timepoint}_sorted.bam'
    output:
        bed=wd + '3_bam/{timepoint}_start_sites.bed'
    shell:
        """
        /work/users/w/u/wuhao/Nanopore_antisense/fetch_start_sites.sh {S_samtools} {input.bam} {output.bed}
        """