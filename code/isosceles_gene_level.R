#!/usr/bin/env Rscript

# Author: Moe
# Script performs initial `isosceles` single-cell processing that is compute intensive and may take a long time to run in an interactive session; it generates a number of `rds` output files that are ready for further downstream analysis via on-demand R sessions.
# Input: takes mapped reads, annotation and reference input file.
# Output: multiple `rds` objects described by in-line comments.

library(optparse)
library(Isosceles)

option_list <- list(
    optparse::make_option(c("-r", "--root_dir"), type="character", default=NULL, help="Root directory of the analysis"),
    optparse::make_option(c("-f", "--ref_fa"), type="character", default=NULL, help="Reference genome fasta sequence file"),
    optparse::make_option(c("-a", "--ref_anno"), type="character", default=NULL, help="Reference annotation gtf file"),
    optparse::make_option(c("-b", "--bam"), type="character", default=NULL, help="Reference mapped bam file"),
    optparse::make_option(c("-w", "--wld"), type="character", default=NULL, help="Sample wildcard/prefix"),
    optparse::make_option(c("-o", "--out_dir"), type="character", default=NULL, help="Path for output directory"),
    optparse::make_option(c("-n", "--n_cores"), type="integer", default=NULL, help="Number of cores for multi-processing")
)


opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser)

message(opt$root_dir)
message(getwd())

setwd(opt$root_dir)
ref_fa <- opt$ref_fa
message(ref_fa)
ref_anno <- opt$ref_anno
message(ref_anno)
bam <- opt$bam
message(bam)
n <- opt$n_cores
message(n)

bam_file <- c(Sample = bam)
# Prepare parsed bam input for isosceles downstream functions
bam_parsed <- bam_to_read_structures(bam_file)
wld <- opt$wld
message(wld)
out_dir <- opt$out_dir
message(out_dir)
out_file <- paste0(wld, "_bam_filtered_parsed.rds")
message(file.path(out_dir, out_file))
# Store parsed bam output as rds
saveRDS(bam_parsed, file = file.path(out_dir, out_file))

# Prepare transcript data for further analysis
message("\n\t*** Preparing transcripts ***\n")
transcript_data <- prepare_transcripts(
  gtf_file = ref_anno,
  genome_fasta_file = ref_fa,
  bam_parsed = bam_parsed,
  known_intron_motifs = c("GT-AG", "GC-AG", "AT-AC"),
  min_bam_splice_read_count = 1,
  min_bam_splice_fraction = 0.01,
  rescue_annotated_introns = TRUE,
  )
out_file <- paste0(wld, "_transcript_data.rds")
saveRDS(transcript_data, file = file.path(out_dir, out_file))
message("\n\t*** Done ***\n")

# run_mode: strict - only annotated reference transcripts will be quantified
# de_novo_strict - de novo transcripts enabled, but all splice sites must be known (i.e. found in reference annotations or provided by the user).
# de_novo_loose - de novo transcripts enabled, but all splice sites must be known or reproducibly passing filters in the aligned reads
message("\n\t*** Preparing bam file ***\n")
se_object <- bam_to_tcc(
      bam_files = bam_file,
      transcript_data = transcript_data,
      run_mode = "strict",
      min_read_count = 1,
      min_relative_expression = 0,
      extend_spliced_transcripts = 100,
      is_single_cell = TRUE,
      barcode_tag = "CB",
      # adjust chunk size in case of memory error
      chunk_size = 1e+06,
      ncpu = n
      )
out_file <- paste0(wld, "_se_object.rds")
saveRDS(se_object, file = file.path(out_dir, out_file))
message("\n\t*** Done ***\n")

# Generate transcript level quantitation
message("\n\t*** Convering bam to transcript-level quantitation ***\n")
se_transcript <- tcc_to_transcript(
      se_tcc = se_object,
      em.maxiter = 100,
      em.conv = 0.01,
      use_length_normalization = FALSE,
      ncpu = n
      )
out_file <- paste0(wld, "_se_transcript.rds")
saveRDS(se_transcript, file = file.path(out_dir, out_file))
message("\n\t*** Done ***\n")

# Collapse transcripts to gene-level quantitation
message("\n\t*** Collapsing transcript-level to gene-level quantitation ***\n")
se_gene <- tcc_to_gene(se_tcc = se_object)
out_file <- paste0(wld, "_se_gene.rds")
saveRDS(se_gene, file = file.path(out_dir, out_file))
message("\n\t*** Done ***\n")