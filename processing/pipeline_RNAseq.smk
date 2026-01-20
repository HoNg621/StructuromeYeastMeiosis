Timepoints = ['YPD', '0h', '6h', '7h', '7h25', '7h5', '8h25', '8h75', '9h5']
Repeats = ['A', 'B']
Reads = ['R1', 'R2']
Strands  = ['fwd', 'rev']
Modifications = ['plus', 'minus']

rule all:
    input:
        expand('1_fastq/unabDNA/unabDNA_pb_{timepoint}{modification}_{repeat}_{read}.fq.gz', timepoint = Timepoints, modification = Modifications, repeat = Repeats, read = Reads),
        '2_processing/RNAexpression/pb_A14201_count.table'
        

rule cutadapt:
    input:
        '1_fastq/{timepoint}{modification}_{repeat}_R1.fastq.gz',
        '1_fastq/{timepoint}{modification}_{repeat}_R2.fastq.gz'
    output:
        '1_fastq/pb_{timepoint}{modification}_{repeat}_R1.fq.gz',
        '1_fastq/pb_{timepoint}{modification}_{repeat}_R2.fq.gz'

    shell:
        'cutadapt -j 0 -q 20 -m 15 --pair-filter=any -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} && rm {input[0]} && rm {input[1]}'
    

rule remove_abDNA:
    input:
        '1_fastq/pb_{timepoint}{modification}_{repeat}_R1.fq.gz',
        '1_fastq/pb_{timepoint}{modification}_{repeat}_R2.fq.gz'
    output:
        temp('1_fastq/abDNA/pb_{timepoint}{modification}_{repeat}.bam'),
        '1_fastq/unabDNA/unabDNA_pb_{timepoint}{modification}_{repeat}_R1.fq.gz',
        '1_fastq/unabDNA/unabDNA_pb_{timepoint}{modification}_{repeat}_R2.fq.gz'
        
    shell:
        'bowtie2 -x indexfile/bt2_sk1_revised_100331_abDNA/bt2_sk1_revised_100331_abDNA -1 {input[0]}  -2 {input[1]} -t --un-conc-gz 1_fastq/unabDNA/unabDNA_pb_{wildcards.timepoint}{wildcards.modification}_{wildcards.repeat}_R%.fq.gz -p 10 | samtools view -bSF4 - > {output[0]}'


    
rule starAlign:
    input:
        expand('1_fastq/unabDNA/unabDNA_pb_{timepoint}{modification}_{repeat}_{read}.fq.gz', timepoint = Timepoints, modification = Modifications, repeat = Repeats, read = Reads)

    output:
        expand('2_processing/RNAexpression/pb_{timepoint}{modification}_{repeat}Aligned.out.bam', timepoint = Timepoints, modification = Modifications, repeat = Repeats)
    run:
        for timepoint in Timepoints:
            for modification in Modifications:
                for repeat in Repeats:
                    fastq1 = '1_fastq/unabDNA/unabDNA_pb_{timepoint}{modification}_{repeat}_R1.fq.gz'
                    fastq2 = '1_fastq/unabDNA/unabDNA_pb_{timepoint}{modification}_{repeat}_R2.fq.gz'
                    shell('STAR --runThreadN 20 --alignEndsType Local  --runMode alignReads --outMultimapperOrder Random  --genomeDir indexfile/star_sk1_revised_100331_filtered/ --readFilesIn <(gunzip -c ' + fastq1 + ') <(gunzip -c ' + fastq2 +') --limitOutSJcollapsed 10000000 --limitIObufferSize 1400000000 --outFileNamePrefix 2_processing/RNAexpression/pb_'+ timepoint + modification + '_'+ repeat +' --outSAMtype BAM Unsorted --outSAMattrIHstart 0 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical') 

rule sort_index:
    input:
        '2_processing/RNAexpression/pb_{timepoint}{modification}_{repeat}Aligned.out.bam'
    output:
        '2_processing/RNAexpression/pb_{timepoint}{modification}_{repeat}Aligned.sorted.bam',
        '2_processing/RNAexpression/pb_{timepoint}{modification}_{repeat}Aligned.sorted.bam.bai'
    shell:
        'samtools sort {input[0]} -o {output[0]}  && samtools index {output[0]} && rm {input[0]}' 


rule RNAcount:
    input:
        expand('2_processing/RNAexpression/pb_{timepoint}{modification}_{repeat}Aligned.sorted.bam', timepoint = Timepoints, modification = Modifications, repeat = Repeats),
        GTF = 'reference/sk1_gene.gtf'

    params:
        ' '.join(expand('2_processing/RNAexpression/pb_{timepoint}{modification}_{repeat}Aligned.sorted.bam', timepoint = Timepoints, modification = Modifications, repeat = Repeats))
    output:
        '2_processing/RNAexpression/pb_A14201_count.table'
    shell:
        'featureCounts -a {input.GTF}  -o 2_processing/RNAexpression/pb_A14201_count.table {params[0]}'