## Single-cell PacBio long-read data processing pipeline (compatible with Kinnex)

This pipeline can be utilised to process single-cell PacBio long-read RNA-seq data generated with Kinnex (MAS-Seq) libraries. It can be used as a modular pipeline to process data end-to-end in order to generate a gene/transcript count matrix or partially to generate intermediate outputs to be used with other pipelines. Snakefiles are organised with `1. -> 5.` prefix to indicate sequential execution of scripts. Each snakefile contains informative comments about included rules and methods used, and users are required to interact with snakefiles at various stages of pipeline execution to provide appropriate orthogonal inputs for correct processing of their data. Additionally each rule has an associated conda environment to ensure reproducibility across platforms. Please note that this is not a fully automated pipeline at the moment, and user input is required for when specified in instructions. Nonetheless, a fully automated pipeline in compliance with [Snakemake best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html) is under development.  

### Pipeline features

- Use of isolated conda environments per Snakemake rule
- Snakemake utilisation for efficient and scalable processing of large datasets
- Extensive customization options
- Bespoke pipeline with included supplementary Python and R scripts
- Multiple novel transcript discovery and annotation options (Pigeon and SQANTI3 with custom filtering)
- Multiple quantitation options with various methods (Isosceles and Pigeon)

### Requirements

* [Snakemake](https://snakemake.github.io)
* [Miniconda](https://docs.anaconda.com/miniconda/)
* Python (latest)
* R (latest)
* Isosceles (Bioconductor)
* Any other requirement is downloaded, installed, and managed by defined conda environments specified per rule.

### Instructions

Root working directory structure post pipeline completion:

```
.
├── analysis
│   ├── m84014_pbmc_0.01
│   └── scripts
├── barcode_corrected
├── bcstats
├── benchmarks
├── classify
├── code 
│   ├── SQANTI3-5.2.2 
│   └── ** ALL PATH VARIABLES DEFINED RELATIVE TO THIS DIRECTORY **
├── collapse
├── data
├── dedup
├── gffcompare
├── lima
├── logs
├── misc
│   └── 10x_barcodes
├── pbmm2
├── pigeon_report
├── ref ** INCLUDE YOUR REFERENCE **
│   └── ** INCLUDE INDEX AND REFERENCE FILES HERE ** 
├── refine
├── seurat
│   ├── m84014_pbmc_0.01_pigeon
│   │   ├── genes_seurat
│   │   └── isoforms_seurat
│   ├── m84014_pbmc_0.01_sq3
│   │   ├── genes_seurat
│   │   └── isoforms_seurat
├── sorted
├── sqanti3_filter
│   └── rules
│       └── m84014_pbmc_0.01
├── sqanti3_qc
│   └── m84014_pbmc_0.01
└── tag
```

Execute the following commands from `./code/` directory in the specified order, and provide required input files for each processing step by modifying the corresponding snakefile directly (note that variables are specified under 'Input requirements'):

1. `snakemake -j <n> --rerun-incomplete --use-conda -s 1.snakefile`  
    Workflow description:  
        - skera: perform segmentation on MAS-Seq reads using MAS-Seq primers list  
        - lima: remove 10x primers from segmented reads (barcoded samples require barcodes.fa and barcode orientation added to parameters [lima_FAQ](https://lima.how/barcode-design.html))  
        - tag: clip UMI and barcodes to use for deduplication a later stage [barcode design](https://isoseq.how/umi/umi-barcode-design.html)  
        - refine: clip polyA tail if persent in sample  
    Input requirements:  
        - Barcode design: `lib_design` consult [umi-barcode-design](https://isoseq.how/umi/umi-barcode-design.html) to choose an appropriate design  
        - Primer sequences: `primers` refers to 10X primer sequences used during library prep (default provided file is suitable for 10X 3' library prep kits)
2. `snakemake -j <n> --rerun-incomplete --use-conda -s 2.snakefile`  
    Input requirements:  
        - Barcodes whitelist: `whitelist` refers to 10X barcode whitelist file  
    Notes:  
            Cell calling method can be adjusted in `rule correct` and `rule bcstats` steps to use percentile instead of default knee method.
3. `snakemake -j <n> --rerun-incomplete --use-conda -s 3.snakefile`  
    Input requirements:  
        - Minimap2 index: `ref_mmi`  
        - Reference genome fasta: `ref_fa`  
        - Reference annotation bed: `anno_bed`
4. `snakemake -j <n> --rerun-incomplete --use-conda -s 4.snakefile_pigeon`  
    Input requirements:  
        - Reference genome fasta: `ref_fa`  
        - Reference annotation gtf: `anno_gtf`  
        - Transcript start sites: `tss`  
        - PolyA motif list: `polya_motif`  
5. `snakemake -j <n> --rerun-incomplete --use-conda -s 4.snakefile_sqanti3`  
    Input requirements:  
        - Reference genome fasta: `general_job_j`  
        - Reference annotation gtf: `anno_gtf`  
        - CAGE peaks: `cage_peak`  
        - PolyA motif list: `polyA_motif`  
        - PolyA sites: `polyA_sites`  
6. `./isosceles_gene_level.R --help`  
    ```
    -r ROOT_DIR, --root_dir=ROOT_DIR
                Root directory of the analysis

        -f REF_FA, --ref_fa=REF_FA
                Reference genome fasta sequence file

        -a REF_ANNO, --ref_anno=REF_ANNO
                Reference annotation gtf file

        -b BAM, --bam=BAM
                Reference mapped bam file

        -w WLD, --wld=WLD
                Sample wildcard/prefix

        -o OUT_DIR, --out_dir=OUT_DIR
                Path for output directory

        -n N_CORES, --n_cores=N_CORES
                Number of cores for multi-processing

        -h, --help
                Show this help message and exit
    ```

    Note:  
        Isosceles package must be installed locally and available through `R_LIBS_USER` environmental variable.

In addition all scripts have `general_job_j` specifying number of processing cores for resource management per rule. `n` in `-j` specifies the number of parallel jobs and in most cases samples to be processed concurrently.  

The outputs then can be used for further downstream analysis with `Seurat` or `Scanpy` for instance.

### Notes

Current directory structure should be preserved on your local environment to ensure correct execution of this pipeline. In summary, `./code/*snakefile*` contain the workflow scripts. All required conda environment config files are also stored in `./code/*.yaml` directory. Raw sequencing data (i.e., initial input to the workflow) should be stored in `./data/` directory. Snakemake will take `./data/*.bam` as wildcards therefore, rename samples appropriately in `./data/` directory if original sample names are uninformative or too long (e.g., `565ytur74.bam` -> `pbmc_kinnex.bam`). Additional supplementary input files such as sequences for primers and barcodes are stored in `./misc/` directory.