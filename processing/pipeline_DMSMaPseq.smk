Timepoints = ['YPD', '0h', '6h', '7h', '7h5', '7h25', '8h25', '8h75', '9h5'] 
Repeats = ['A', 'B']
Reads = ['R1', 'R2']
Strands  = ['fwd', 'rev']
Modifications = ['plus', 'minus']

rule all:
    input:
        expand('stranded/RNAstructure/{timepoint}{repeat}/Pipeline_sk1_transcriptome_profile.txt', timepoint = Timepoints, repeat = Repeats)

rule bowtie2_Align:
    input:
        '1_fastq/unabDNA/unabDNA_pb_{timepoint}{modification}_{repeat}_R1.fq.gz',
        '1_fastq/unabDNA/unabDNA_pb_{timepoint}{modification}_{repeat}_R2.fq.gz',
        'indexfile/sk1_geneUTR150_transcriptome_bowtie2Index.1.bt2'
    output:
        temp('stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.bam'),
        'stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.sorted.bam'
    
    shell:
        'bowtie2 --wrapper basic-0 --fr --no-mixed --no-discordant --sensitive-local --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 --maxins 800 --ignore-quals -x indexfile/sk1_geneUTR150_transcriptome_bowtie2Index -1 {input[0]}  -2 {input[1]} -t -p 40 | samtools view -bSF4 - > {output[0]} && samtools sort {output[0]} -o {output[1]}'

rule split_bam:
    input:
        'stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.sorted.bam'
    output:
        temp('stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.fwd1.bam'),
        temp('stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.fwd2.bam'),
        temp('stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.rev1.bam'),
        temp('stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.rev2.bam')
    run:
        shell('samtools view -b -f 128 -F 16 {input} > {output[0]}')
        shell('samtools view -b -f 64 -F 32 {input} > {output[1]}')
        shell('samtools view -b -f 144 {input} > {output[2]}')
        shell('samtools view -b -f 96 {input} > {output[3]}')
        # shell('samtools view -b -f 128 -F 16 {input} > {output[0]}') # fwd2
        # shell('samtools view -b -f 64 -F 32 {input} > {output[1]}') # fwd1
        # shell('samtools view -b -f 80 {input} > {output[2]}')        # rev1
        # shell('samtools view -b -f 96 {input} > {output[3]}')        # rev2     
rule merge_bam:
    input:
        'stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.fwd1.bam',
        'stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.fwd2.bam',
        'stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.rev1.bam',
        'stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.rev2.bam'
    output:
        fwd = temp('stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.fwd.bam'),
        rev = temp('stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.rev.bam')
    run:
        shell('samtools merge -f {output.fwd} {input[0]} {input[1]}')
        shell('samtools merge -f {output.rev} {input[2]} {input[3]}')

rule sort_byQuery:
    input:
        'stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.{strand}.bam'
    output:
        temp('stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.{strand}.sorted.bam')
    shell:
        'samtools sort -n {input} -o {output}'

rule split_fastq:
    input:
        'stranded/RNAstructure/bt2/pb_{timepoint}{modification}_{repeat}.{strand}.sorted.bam'
    output:
        'stranded/RNAstructure/fastq/pb_{timepoint}{modification}_{repeat}_{strand}_R1.fq',
        'stranded/RNAstructure/fastq/pb_{timepoint}{modification}_{repeat}_{strand}_R2.fq'
    run:
        shell('bedtools bamtofastq -i {input} -fq {output[0]} -fq2 {output[1]}')

rule pigz:
    input:
        'stranded/RNAstructure/fastq/pb_{timepoint}{modification}_{repeat}_{strand}_R1.fq',
        'stranded/RNAstructure/fastq/pb_{timepoint}{modification}_{repeat}_{strand}_R2.fq'
    output:
        'stranded/RNAstructure/fastq/pb_{timepoint}{modification}_{repeat}_{strand}_R1.fq.gz',
        'stranded/RNAstructure/fastq/pb_{timepoint}{modification}_{repeat}_{strand}_R2.fq.gz'
    run:
        shell('pigz {input[0]}  {input[1]}')
        
rule shapemapper_fwd:
    input:
        'stranded/RNAstructure/fastq/pb_{timepoint}plus_{repeat}_fwd_R1.fq.gz',
        'stranded/RNAstructure/fastq/pb_{timepoint}plus_{repeat}_fwd_R2.fq.gz',
        'stranded/RNAstructure/fastq/pb_{timepoint}minus_{repeat}_fwd_R1.fq.gz',
        'stranded/RNAstructure/fastq/pb_{timepoint}minus_{repeat}_fwd_R2.fq.gz'
    output:
        'stranded/RNAstructure/{timepoint}{repeat}/Pipeline_sk1_transcriptome_profile.txt'
    shell:
        'shapemapper --nproc 20 --min-depth 10 --overwrite --name {wildcards.timepoint}fwd{wildcards.repeat} --target reference/sk1_geneUTR150_transcriptome.fasta  --modified --R1 {input[0]} --R2 {input[1]} --untreated --R1 {input[2]} --R2 {input[3]} --random-primer-len 6 --out stranded/RNAstructure/{wildcards.timepoint}{wildcards.repeat}/  --log stranded/RNAstructure/{wildcards.timepoint}{wildcards.repeat}/{wildcards.timepoint}{wildcards.repeat}.log --temp shapemapper_tem{wildcards.timepoint}fwd{wildcards.repeat}/' 
