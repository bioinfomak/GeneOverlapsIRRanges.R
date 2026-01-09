library(VariantAnnotation)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

setwd("PATH")

# List all filtered VCFs
vcf_files <- list.files(pattern = "_filtered_variants\\.vcf$")

# Loop over each VCF
for (vcf_file in vcf_files) {
  
  message("Processing: ", vcf_file)
  
  # Read VCF
  vcf <- readVcf(vcf_file, genome = "hg38")
  
  # Keep only PASS variants
  vcf_pass <- vcf[fixed(vcf)$FILTER == "PASS", ]
  
  # Define output name
  out_file <- sub("_filtered_variants\\.vcf$", "_output_filtered.vcf", vcf_file)
  
  # Write VCF
  writeVcf(vcf_pass, out_file)
}


############################################################
## Map variants to genes for multiple VCFs (with variant_id)
############################################################
# Input and output folders
input_dir  <- "PATH"
output_dir <- "PATH"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

setwd(input_dir)

# Load gene annotation once
gr_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys       = as.character(gr_genes$gene_id),
  column     = "SYMBOL",
  keytype    = "ENTREZID",
  multiVals  = "first"
)
gr_genes$gene_name <- gene_symbols

# Get all filtered VCFs
vcf_files <- list.files(pattern = "_output_filtered\\.vcf$", full.names = TRUE)
stopifnot(length(vcf_files) > 0)

for (vcf_file in vcf_files) {
  message("Processing: ", basename(vcf_file))
  
  vcf <- readVcf(vcf_file, genome = "hg38")
  gr_variants <- rowRanges(vcf)
  hits <- findOverlaps(gr_variants, gr_genes, ignore.strand = TRUE)
  
  if (length(hits) == 0) {
    message("  No gene overlaps found. Skipping.")
    next
  }
  
  alt_alleles <- sapply(
    mcols(gr_variants)$ALT[queryHits(hits)],
    function(x) paste(as.character(x), collapse = ",")
  )
  
  variant_id_vec <- paste0(
    as.character(seqnames(gr_variants)[queryHits(hits)]), ":",
    start(gr_variants)[queryHits(hits)], "_",
    as.character(mcols(gr_variants)$REF[queryHits(hits)]), "/",
    alt_alleles
  )
  
  overlap_table <- data.frame(
    variant_id = variant_id_vec,
    chr        = as.character(seqnames(gr_variants)[queryHits(hits)]),
    pos        = start(gr_variants)[queryHits(hits)],
    ref        = as.character(mcols(gr_variants)$REF[queryHits(hits)]),
    alt        = alt_alleles,
    gene_id    = gr_genes$gene_id[subjectHits(hits)],
    gene_name  = gr_genes$gene_name[subjectHits(hits)],
    stringsAsFactors = FALSE
  )
  
  out_csv <- file.path(
    output_dir,
    sub("_output_filtered\\.vcf$", "_VCF_Gene.csv", basename(vcf_file))
  )
  
  write.csv(overlap_table, out_csv, row.names = FALSE)
  message("  Saved: ", out_csv)
}
