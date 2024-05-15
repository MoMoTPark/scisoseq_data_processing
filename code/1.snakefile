## PacBio MAS-Seq/Kinnex scRNA-seq data processing pipeline

# Overall workflow:
# skera: perform segmentation on MAS-Seq reads using MAS-Seq primers list
# lima: remove 10x primers from segmented reads (barcoded samples require barcodes.fa and barcode orientation added to parameters [lima_FAQ](https://lima.how/barcode-design.html))
# tag: clip UMI and barcodes to use for deduplication a later stage (barcode design: https://isoseq.how/umi/umi-barcode-design.html)
# refine: clip polyA tail if persent in sample
# correct: perform barcode correction to recover more reads and make deduplication more effective
# bcstats: get cell calling statistics to adjust correction methods if required
# sort: sort bam files with corrected barcodes by `CB` tag to prepare input as groupdedup
# groupdedup: deduplicate reads based on UMIs grouped by corrected cell barcodes
# pbmm2: map reads to indexed reference genome
# collapse: collapse redundant transcripts
# pigeon_prepare: prepare input files for pigeon to classify and filter transcripts
# pigeon_classify: classification of transcripts into SQANTI3 defined categories
# pigeon_filter: filter transcripts to reduce false positives
# pigeon_report: prepare a report for classification and filtering steps
# make_seurat: create a barcode-feature matrix compatible with Seurat as input

# Sample wildcards
IDS, = glob_wildcards("../data/{id}.bam")

# Define local variable to adjust rule parameters globally
# Number of worker cores for general jobs excluding mapping and memory intensive tasks
general_job_j = 64

# Final outputs for the pipeline
rule all:
    input:
        expand("../dedup/{id}.fltncc.sorted.dedup.bam", id=IDS),
        expand("../pbmm2/{id}.dedup.mapped.bam", id=IDS),
        expand("../bcstats/{id}.bcstats.tsv", id=IDS),
        expand("../collapse/{id}.collapsed.sorted.gff", id=IDS),
        expand("../pigeon_report/{id}_pigeon_report.txt", id=IDS),
        expand("../seurat/{id}.info.csv", id=IDS),

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
        segmented = "../skera/{id}.segmented.bam"
    output:
        lima = "../lima/{id}.fl.5p--3p.bam"
    conda: "lima_conda.yaml"
    benchmark: "../benchmarks/{id}_lima.benchmark"
    params: 
        log = "../logs/{id}_lima.log",
        log_level = "TRACE",
        primers = "../misc/10x_3kit_primers.fasta",
        outfile = "../lima/{id}.fl.bam",
        j = general_job_j,
    shell:
        '''
        lima --isoseq {input.segmented} \
        {params.primers} \
        {params.outfile} \
        --per-read \
        --peak-guess \
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
        design = "T-12U-16B",
        j = general_job_j,
    shell:
        '''
        isoseq tag --design {params.design} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        -j {params.j} \
        {input.lima} \
        {output.tag};
        rm {input.lima}
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
        primer = "../misc/10x_3kit_primers.fasta",
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
        {output.refine};
        rm {input.tag}
        '''

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
        barcodes = "../misc/10x_barcodes/3M-february-2018-REVERSE-COMPLEMENTED.txt.gz",
        j = general_job_j,
    shell:
        '''
        isoseq correct -j {params.j} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        --barcodes {params.barcodes} \
        {input.refine} \
        {output.corrected};
        rm {input.refine}
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

# Sort corrected bam file by CB tag 
rule sort:
    input:
        corrected = "../barcode_corrected/{id}.fltncc.bam",
        # `bcstats` input added to force snakemake execute `bcstats` rule prior to this rule to prevent pre-mature deletion of intermediate files
        bcstats = "../bcstats/{id}.bcstats.tsv",
    output:
        sort = "../sorted/{id}.fltncc.sorted.bam"
    conda: "samtools_conda.yaml"
    benchmark: "../benchmarks/{id}_sort.benchmark"
    params:
        log = "../logs/{id}_sort.log",
        j = general_job_j,
    shell:
        '''
        samtools sort -@ {params.j} \
        -t CB \
        {input.corrected} \
        -O BAM \
        -o {output.sort} 2> {params.log};
        rm {input.corrected}
        '''

# Deduplicate reads based on UMI grouped by cell barcodes
rule dedup:
    input:
        sort = "../sorted/{id}.fltncc.sorted.bam"
    output:
        dedup = "../dedup/{id}.fltncc.sorted.dedup.bam"
    conda: "isoseq_conda.yaml"
    benchmark: "../benchmarks/{id}_dedup.benchmark"
    params:
        log = "../logs/{id}_dedup.log",
        log_level = "TRACE",
        j = general_job_j,
    shell:
        '''
        isoseq groupdedup -j {params.j} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        {input.sort} \
        {output.dedup}
        '''

# Map reads to reference with pbmm2
rule pbmm2:
    input:
        dedup = "../dedup/{id}.fltncc.sorted.dedup.bam"
    output:
        mapped = "../pbmm2/{id}.dedup.mapped.bam"
    conda: "pbmm2_conda.yaml"
    benchmark: "../benchmarks/{id}_pbmm2.benchmark"
    params:
        log = "../logs/{id}_pbmm2.log",
        log_level = "TRACE",
        ref = "../ref/GRCh38.p14.genome.mmi",
        j = 32,
    shell:
        '''
        pbmm2 align --preset ISOSEQ --sort \
        --log-file {params.log} \
        --log-level {params.log_level} \
        -j {params.j} \
        {params.ref} \
        {input.dedup} \
        {output.mapped}
        '''

# STARlong step, benchmark step (pending)
# rule starlong:

# Collapse reduntant transcripts into unique isoforms
# Alternative tool for this step would be [TAMA](https://github.com/GenomeRIK/tama/wiki/Tama-Collapse) or [Cupcake](https://github.com/Magdoll/cDNA_Cupcake) collapse script 
rule collapse:
    input:
        mapped = "../pbmm2/{id}.dedup.mapped.bam"
    output:
        collapsed = "../collapse/{id}.collapsed.gff"
    conda: "isoseq_conda.yaml"
    benchmark: "../benchmarks/{id}_collapse.benchmark"
    params:
        log = "../logs/{id}_collapse.log",
        log_level = "TRACE",
        j = general_job_j,
    shell:
        '''
        isoseq collapse -j {params.j} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        --do-not-collapse-extra-5exons \
        {input.mapped} \
        {output.collapsed}
        '''
# Alternative workflow: Sqanti3 QC and Filter 
# Sort collapse output gff
rule pigeon_prepare:
    input:
        collapsed = "../collapse/{id}.collapsed.gff"
    output:
        sorted = "../collapse/{id}.collapsed.sorted.gff"
    conda: "pbpigeon_conda.yaml"
    benchmark: "../benchmarks/{id}_pigeon_prepare.benchmark"
    params:
        log = "../logs/{id}_pigeon_prepare.log",
        log_level = "TRACE",
    shell:
        '''
        pigeon prepare --log-file {params.log} \
        --log-level {params.log_level} \
        {input.collapsed}
        '''

# Classify transcripts into categories: https://isoseq.how/classification/categories
rule pigeon_classify:
    input:
        sorted = "../collapse/{id}.collapsed.sorted.gff"
    output:
        classification = "../classify/{id}_classification.txt"
    conda: "pbpigeon_conda.yaml"
    benchmark: "../benchmarks/{id}_pigeon_classify.benchmark"
    params:
        log = "../logs/{id}_pigeon_classify.log",
        log_level = "TRACE",
        j = general_job_j,
        out_dir = "../classify",
        prefix = "{id}",
        anno = "../ref/gencode.v45.primary_assembly.annotation.sorted.gtf",
        ref = "../ref/GRCh38.p14.genome.fa",
        polyA = "../ref/polyA.list.txt",
        cage_peaks = "../ref/hg38_liftover_CAGE_peaks_phase1and2.sorted.bed",
        coverage = "../ref/intropolis.v1.hg19_with_liftover_to_hg38.min_count_10.sorted.tsv",
    shell:
        '''
        pigeon classify -j {params.j} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        -d {params.out_dir} \
        -o {params.prefix} \
        --cage-peak {params.cage_peaks} \
        --poly-a {params.polyA} \
        --coverage {params.coverage} \
        --gene-id \
        {input.sorted} \
        {params.anno} \
        {params.ref}
        '''

# Filter classified transcripts
rule pigeon_filter:
    input:
        classification = "../classify/{id}_classification.txt"
    output:
        filtered = "../classify/{id}_classification.filtered_lite_classification.txt"
    conda: "pbpigeon_conda.yaml"
    benchmark: "../benchmarks/{id}_pigeon_filter.benchmark"
    params:
        log = "../logs/{id}_pigeon_filter.log",
        log_level = "TRACE",
        j = general_job_j,
        i = "../collapse/{id}.collapsed.sorted.gff"
    shell:
        '''
        pigeon filter -j {params.j} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        -i {params.i} \
        {input.classification}
        '''

# Generate report from filtered transcripts output
rule pigeon_report:
    input:
        filtered = "../classify/{id}_classification.filtered_lite_classification.txt"
    output:
        report = "../pigeon_report/{id}_pigeon_report.txt"
    conda: "pbpigeon_conda.yaml"
    benchmark: "../benchmarks/{id}_pigeon_report.benchmark"
    params:
        log = "../logs/{id}_pigeon_report.log",
        log_level = "TRACE",
        j = general_job_j,
    shell:
        '''
        pigeon report -j {params.j} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        {input.filtered} \
        {output.report}
        '''

# Generate a downstream Seurat input matrix
rule make_seurat:
    input:
        filtered = "../classify/{id}_classification.filtered_lite_classification.txt"
    output:
        seurat = "../seurat/{id}.info.csv"
    conda: "pbpigeon_conda.yaml"
    benchmark: "../benchmarks/{id}_make_seurat.benchmark"
    params:
        log = "../logs/{id}_make_seurat.log",
        log_level = "TRACE",
        j = general_job_j,
        dedup = "../dedup/{id}.fltncc.sorted.dedup.fasta",
        group = "../collapse/{id}.collapsed.group.txt",
        out_prefix = "{id}",
        out_dir = "../seurat",
    shell:
        '''
        pigeon make-seurat -j {params.j} \
        --log-file {params.log} \
        --log-level {params.log_level} \
        --dedup {params.dedup} \
        --group {params.group} \
        -o {params.out_prefix} \
        -d {params.out_dir} \
        {input.filtered}
        '''
