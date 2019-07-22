#!/usr/bin/env python3


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


#########
# RULES #
#########

rule target:
    input:
        expand('output/010_bbduk-trim/{sample}_r1.fq.gz',
               sample=list(sample_id_key.keys()))

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


