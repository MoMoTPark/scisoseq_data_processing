## PacBio MAS-Seq/Kinnex scRNA-seq data processing pipeline
# Author: Moe
# PART 2

# Workflow description:
# correct: perform barcode correction to recover more reads and make deduplication more effective
# bcstats: get cell calling statistics to adjust correction methods if required
# plot_knees: generate cell calling QC knee plots

# Sample wildcards
IDS, = glob_wildcards("../data/{id}.bam")

# Define local variable to adjust rule parameters globally
# Number of worker cores for general jobs excluding mapping and memory intensive tasks
general_job_j = 64
# Barcode whitelist file (note 10x V4 would require the latest barcode whitelist file supplied at ../misc/3M-3pgex-may-2023_rc.txt.gz)
whitelist = "../misc/10x_barcodes/3M-february-2018-REVERSE-COMPLEMENTED.txt.gz"

# Final outputs for the pipeline
rule all:
    input:
        expand("../bcstats/{id}.bcstats.tsv", id=IDS),
        expand("../bcstats/{id}.knee.png", id=IDS),

# Barcode correction to recover about 5% of reads in most cases and make deduplication process more efficient (option to choose between `percentile` and `knee` methods [knee])
rule correct:
    input:
        refine = "../refine/{id}.fltnc.bam"
    output:
        corrected = "../barcode_corrected/{id}.fltncc.bam"
    conda: "isoseq_conda.yaml"
    benchmark: "../benchmarks/{id}_correct.benchmark"
    params:
        log = "../logs/{id}_correct.log",
        log_level = "TRACE",
        barcodes = whitelist,
        j = general_job_j,
        # method = "knee",
    shell:
        '''
        isoseq correct -j {params.j} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        --barcodes {params.barcodes} \
        {input.refine} \
        {output.corrected}
        ''' 
# --method percentile \
# --percentile 99 \

# Calculate cell calling statistics (choose method according to `rule correct` either `knee` by default or change to `--method percentile` and define value with `--percentile 99`)
rule bcstats:
    input:
        corrected = "../barcode_corrected/{id}.fltncc.bam"
    output:
        bcstats = "../bcstats/{id}.bcstats.tsv"
    conda: "isoseq_conda.yaml"
    benchmark: "../benchmarks/{id}_bcstats.benchmark"
    params:
        log = "../logs/{id}_bcstats.log",
        log_level = "TRACE",
        json = "../bcstats/{id}.bcstats.json",
        j = general_job_j,
        # method = "knee",
    shell:
        '''
        isoseq bcstats -j {params.j} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        -o {output.bcstats} \
        --json {params.json} \
        {input.corrected}
        ''' 
# --method percentile \
# --percentile 99 \

# Generate cell calling QC knee plots
rule plot_knees:
    input:
        bcstats = "../bcstats/{id}.bcstats.tsv"
    output:
        knee_plot = "../bcstats/{id}.knee.png"
    conda: "plotknees_conda.yaml"
    benchmark: "../benchmarks/{id}_plot_knees.benchmark"
    params:
        log = "../logs/{id}_plot_knees.log",
        percentile = 95,
        output = "../bcstats/{id}",
    shell:
        '''
        python plot_knees.py \
        --tsv {input.bcstats} \
        --output {params.output} \
        --estimate_percentile {params.percentile} 2> {params.log}
        '''