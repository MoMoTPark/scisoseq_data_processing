## PacBio MAS-Seq/Kinnex scRNA-seq data processing pipeline
# Author: Moe
# PART 3

# Workflow description:
# sort: sort bam files with corrected barcodes by `CB` tag to prepare input as groupdedup
# dedup: deduplicate reads based on UMIs grouped by corrected cell barcodes
# pbmm2: map reads to indexed reference genome
# filter_bam: filter mapped bam file to generate two separate outputs with 'cell' and 'non-cell' reads
# collapse: collapse redundant transcripts (only 'cell' reads)

# Sample wildcards
IDS, = glob_wildcards("../data/{id}.bam")

# Define local variable to adjust rule parameters globally
# Number of worker cores for general jobs excluding mapping and memory intensive tasks
general_job_j = 64

# ** TODO: add as config ** 
# Human pbmm2 index
ref_mmi = "../ref/GRCh38.primary_assembly.genome.mmi"
ref_fa = "../ref/GRCh38.primary_assembly.genome.fa"
anno_bed = "../ref/gencode.v46.primary_assembly.annotation.bed"

# Final outputs for the pipeline
rule all:
    input:
        expand("../dedup/{id}.fltncc.sorted.dedup.bam", id=IDS),
        expand("../pbmm2/{id}.dedup.mapped.bam", id=IDS),
        expand("../pbmm2/{id}.dedup.mapped.cell.bam", id=IDS),
        expand("../pbmm2/{id}.dedup.mapped.non_cell.bam", id=IDS),
        expand("../collapse/{id}.collapsed.gff", id=IDS),
        expand("../dedup/{id}.fltncc.sorted.dedup.cell.bam", id=IDS),

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
        -o {output.sort} 2> {params.log}
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
        --keep-non-real-cells \
        {input.sort} \
        {output.dedup}
        '''

# Filter reads identified as "rc:i:1" aka "cell" as input for `make_seurat`
rule filter_dedup:
    input:
        dedup = "../dedup/{id}.fltncc.sorted.dedup.bam"
    output:
        dedup_cell = "../dedup/{id}.fltncc.sorted.dedup.cell.bam"
    conda: "samtools_conda.yaml"
    benchmark: "../benchmarks/{id}_filter_dedup.benchmark"
    params:
        log = "../logs/{id}_filter_dedup.log",
        j = general_job_j,
    shell:
        '''
        samtools view -@ {params.j} -h --tag rc:1 {input.dedup} -o {output.dedup_cell}
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
        ref = ref_mmi,
        j = 64,
    shell:
        '''
        pbmm2 align --preset ISOSEQ \
        --sort \
        --log-file {params.log} \
        --log-level {params.log_level} \
        -j {params.j} \
        {params.ref} \
        {input.dedup} \
        {output.mapped}
        '''

# Filter reads with 'cell' tag: `rc:i:1` and reads with 'non-cell' tag: `rc:i:0` into two separate bam files
rule filter_bam:
    input:
        mapped = "../pbmm2/{id}.dedup.mapped.bam"
    output:
        cell = "../pbmm2/{id}.dedup.mapped.cell.bam",
        non_cell = "../pbmm2/{id}.dedup.mapped.non_cell.bam"
    conda: "samtools_conda.yaml"
    benchmark: "../benchmarks/{id}_filter_bam.benchmark"
    params:
        j = general_job_j,
    shell:
        '''
        samtools view -@ {params.j} -h --tag rc:1 {input.mapped} -o {output.cell};
        samtools index -@ {params.j} {output.cell};
        samtools view -@ {params.j} -h --tag rc:0 {input.mapped} -o {output.non_cell};
        samtools index -@ {params.j} {output.non_cell}
        '''

# Correct microexon misallignment with MisER
rule miser:
    input:
        cell = "../pbmm2/{id}.dedup.mapped.cell.bam"
    output:
        corrected = "../pbmm2/{id}.dedup.mapped.cell.miser.corrected.bam"
    conda: "miser_conda.yaml"
    benchmark: "../benchmarks/{id}_miser.benchmark"
    log: "../logs/{id}_miser.log"
    params:
        out_region = "../pbmm2/{id}.dedup.mapped.cell.miser.regions.out.txt",
        j = general_job_j,
        ref = ref_fa,
        bed = anno_bed
    shell:
        '''
        MisER -c {params.j} --strandSpecific --setTag --outBam {output.corrected} {input.cell} {params.ref} {params.bed} {params.out_region} 2> {log}
        '''

# Replace "M" with "=" in CIGAR strings to generate PacBio compatible bam
rule modify_bam:
    input:
        corrected = "../pbmm2/{id}.dedup.mapped.cell.miser.corrected.bam"
    output:
        modified = "../pbmm2/{id}.dedup.mapped.cell.miser.corrected.modified.sorted.bam"
    conda: "samtools_conda.yaml"
    log: "../logs/{id}_modify_bam.log"
    benchmark: "../benchmarks/{id}_modify_bam.benchmark"
    params:
        sam = "../pbmm2/{id}.dedup.mapped.cell.miser.corrected.modified.sam",
        j = general_job_j
    shell:
        '''
        samtools index {input.corrected};
        samtools view -H {input.corrected} > {params.sam};
        samtools view {input.corrected} | awk 'BEGIN {{OFS="\\t"}} {{$6 = gensub("M", "=", "g", $6); print}}' >> {params.sam};
        samtools sort -@ {params.j} -O BAM -o {output.modified} {params.sam};
        samtools index {output.modified}
        '''

# Collapse reduntant transcripts into unique isoforms
# Alternative tool for this step would be [TAMA](https://github.com/GenomeRIK/tama/wiki/Tama-Collapse) or [Cupcake](https://github.com/Magdoll/cDNA_Cupcake) collapse script 
rule collapse:
    input:
        modified = "../pbmm2/{id}.dedup.mapped.cell.miser.corrected.modified.sorted.bam"
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
        {input.modified} \
        {output.collapsed}
        '''
# REMOVE UNDER DEV
# --do-not-collapse-extra-5exons \