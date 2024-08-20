## PacBio MAS-Seq/Kinnex scRNA-seq data processing pipeline
# Author: Moe
# PART 1

# Workflow description:
# skera: perform segmentation on MAS-Seq reads using MAS-Seq primers list
# lima: remove 10x primers from segmented reads (barcoded samples require barcodes.fa and barcode orientation added to parameters [lima_FAQ](https://lima.how/barcode-design.html))
# tag: clip UMI and barcodes to use for deduplication a later stage (barcode design: https://isoseq.how/umi/umi-barcode-design.html)
# refine: clip polyA tail if persent in sample

# Sample wildcards
IDS, = glob_wildcards("../data/{id}.bam")

# Define local variable to adjust rule parameters globally
# Number of worker cores for general jobs excluding mapping and memory intensive tasks
general_job_j = 64
lib_design = "T-12U-16B"
primers = "../misc/10x_3kit_primers.fasta"

# Final outputs for the pipeline
rule all:
    input:
        expand("../refine/{id}.fltnc.bam", id=IDS),

# Perform segmentation on MAS-Seq reads
rule skera:
    input:
        bam = "../data/{id}.bam"
    output:
        segmented = "../skera/{id}.segmented.bam"
    conda: "pbskera_conda.yaml"
    benchmark: "../benchmarks/{id}_skera.benchmark"
    params:
        log = "../logs/{id}_skera.log",
        log_level = "TRACE",
        mas_primers = "../misc/mas16_primers.fasta",
        j = general_job_j,
    shell:
        '''
        skera split -j {params.j} \
        --log-level {params.log_level} \
        --log-file {params.log} \
        {input.bam} {params.mas_primers} {output.segmented}
        '''

# Remove 10x primers and re-orient reads to 5'->3'
rule lima:
    input:
        segmented = "../data/{id}.bam",
        # segmented = "../skera/{id}.segmented.bam"
    output:
        lima = "../lima/{id}.fl.5p--3p.bam"
    conda: "lima_conda.yaml"
    benchmark: "../benchmarks/{id}_lima.benchmark"
    params: 
        log = "../logs/{id}_lima.log",
        log_level = "TRACE",
        primers = primers,
        outfile = "../lima/{id}.fl.bam",
        j = general_job_j,
    shell:
        '''
        lima --isoseq {input.segmented} \
        {params.primers} \
        {params.outfile} \
        --per-read \
        --peek-guess \
        --log-level {params.log_level} \
        --log-file {params.log} \
        -j {params.j}
        '''

# Clip UMI/barcodes and store for downstream deduplication (barcode design should be specified: https://isoseq.how/umi/umi-barcode-design.html e.g. 10X 3' Chromium kit design = T-12U-16B)
rule tag:
    input:
        lima = "../lima/{id}.fl.5p--3p.bam"
    output:
        tag = "../tag/{id}.flt.bam"
    conda: "isoseq_conda.yaml"
    benchmark: "../benchmarks/{id}_tag.benchmark"
    params:
        log = "../logs/{id}_tag.log",
        log_level = "TRACE",
        # Adjust design based on experiment
        design = lib_design,
        j = general_job_j,
    shell:
        '''
        isoseq tag --design {params.design} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        -j {params.j} \
        {input.lima} \
        {output.tag}
        '''

# Trim polyA and remove concatenated reads
# ***Note that output of this step is an optional input for `collapse` (bulk ISO-Seq data)
rule refine:
    input:
        tag = "../tag/{id}.flt.bam"
    output:
        refine = "../refine/{id}.fltnc.bam"
    conda: "isoseq_conda.yaml"
    benchmark: "../benchmarks/{id}_refine.benchmark"
    params:
        primer = primers,
        log = "../logs/{id}_refine.log",
        log_level = "TRACE",
        j = general_job_j,
    shell:
        '''
        isoseq refine \
        --require-polya \
        --log-file {params.log} \
        --log-level {params.log_level} \
        -j {params.j} \
        {input.tag} \
        {params.primer} \
        {output.refine}
        '''