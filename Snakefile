#!/usr/bin/env python3

import multiprocessing


#############
# FUNCTIONS #
#############

def input_reads_resolver(wildcards):
    my_sample_number = sample_id_key[wildcards.sample]
    return({
        'r1': ('data/raw_reads/CDWA0ANXX-3848-'
               f'{my_sample_number:02}-49-01_'
               f'S{my_sample_number}_L001_R1_001.fastq.gz'),
        'r2': ('data/raw_reads/CDWA0ANXX-3848-'
               f'{my_sample_number:02}-49-01_'
               f'S{my_sample_number}_L001_R2_001.fastq.gz')})


###########
# GLOBALS #
###########

sample_id_key = {
    's0_A': 1,
    's1_A': 2,
    's2_A': 3,
    's3_A': 4,
    'sQMP_A': 5,
    's0_B': 6,
    's1_B': 7,
    's2_B': 8,
    's3_B': 9,
    'sQMP_B': 10,
    's0_C': 11,
    's1_C': 12,
    's2_C': 13,
    's3_C': 14,
    'sQMP_C': 15}

# singularity containers
bbmap = 'shub://TomHarrop/singularity-containers:bbmap_38.50b'
salmon = 'local_containers/salmon_0.14.1.sif'
salmontools = 'local_containers/salmontools_23eac84.sif'
bioconductor = 'local_containers/bioconductor_3.9.sif'

#########
# RULES #
#########

rule target:
    input:
        'output/030_deseq/lrt_results.csv'


rule deseq2_lrt:
    input:
        dds = 'output/030_deseq/dds.Rds'
    output:
        pca_plot = 'output/030_deseq/pca.pdf',
        lrt_results = 'output/030_deseq/lrt_results.csv'
    params:
        alpha = 0.01
    log:
        'output/logs/030_deseq/deseq2_lrt.log'
    singularity:
        bioconductor
    script:
        'src/deseq2_lrt.R'


rule generate_deseq2_object:
    input:
        quant_files = expand('output/020_salmon/{sample}/quant.sf',
                             sample=list(sample_id_key.keys())),
        gff = 'data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.gff'
    output:
        'output/030_deseq/dds.Rds'
    log:
        'output/logs/030_deseq/generate_deseq2_object.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        bioconductor
    script:
        'src/generate_deseq2_object.R'

rule salmon_quant:
    input:
        r1 = 'output/010_bbduk-trim/{sample}_r1.fq.gz',
        r2 = 'output/010_bbduk-trim/{sample}_r2.fq.gz',
        index = 'output/005_index'
    output:
        'output/020_salmon/{sample}/quant.sf'
    log:
        'output/logs/020_salmon/{sample}.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        salmon
    shell:
        'salmon quant '
        '--libType ISR '
        '--index {input.index} '
        '--mates1 {input.r1} '
        '--mates2 {input.r2} '
        '--output {output} '
        '--threads {threads} '
        '--validateMappings '
        '--gcBias '
        '&> {log}'


rule bbduk_trim:
    input:
        unpack(input_reads_resolver)
    output:
        r1 = 'output/010_bbduk-trim/{sample}_r1.fq.gz',
        r2 = 'output/010_bbduk-trim/{sample}_r2.fq.gz'
    params:
        adapters = '/adapters.fa'
    log:
        'output/logs/010_bbduk-trim/{sample}.log'
    threads:
        1
    singularity:
        bbmap
    shell:
        'bbduk.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'

rule generate_index:
    input:
        transcriptome = 'output/000_ref/gentrome.fa',
        decoys = 'output/000_ref/decoys.txt'
    output:
        directory('output/005_index')
    log:
        'output/logs/005_index/generate_index.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        salmon
    shell:
        'salmon index '
        '--transcripts {input.transcriptome} '
        '--index {output} '
        '--threads {threads} '
        '--decoys {input.decoys} '
        '&> {log}'

rule generate_decoy_trancriptome:
    input:
        fasta = 'data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.fna',
        transcriptome = 'data/ref/GCF_003254395.2_Amel_HAv3.1_rna.fna',
        gff = 'data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.gff'
    output:
        'output/000_ref/gentrome.fa',
        'output/000_ref/decoys.txt'
    params:
        outdir = 'output/000_ref'
    log:
        'output/logs/000_ref/generate_decoy_trancriptome.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        salmontools
    shell:
        'generateDecoyTranscriptome.sh '
        '-j {threads} '
        '-b /usr/bin/bedtools '
        '-m /usr/local/bin/mashmap '
        '-a {input.gff} '
        '-g {input.fasta} '
        '-t {input.transcriptome} '
        '-o {params.outdir} '
        '&> {log}'
