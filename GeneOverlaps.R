library(VariantAnnotation)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#----1) convert filtered VCF File to output.filtered VCF -----
setwd ("D:\ALS_26122025\BAM")
list.files()
# Load the VCF file
vcf_file <- "SS91_filtered_variants.vcf"
vcf <- readVcf(vcf_file, genome = "hg38")
vcf_filter<- vcf[(vcf@fixed$FILTER=="PASS"),]
# Save the VCF file
writeVcf(vcf_filter, "SS91_output_filterd.vcf")

# --- 1) Load your VCF ---
list.files()
vcf_file <- "SS91_output_filterd.vcf"
vcf <- readVcf(vcf_file, genome="hg38")
gr_variants <- rowRanges(vcf)  # variant GRanges
info(vcf)
# --- 2) Load gene annotation from TxDb ---
gr_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
#cds(TxDb.Hsapiens.UCSC.hg38.knownGene,columns="cds_id", filter=NULL, use.names=FALSE)

# Map ENTREZ IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys=as.character(gr_genes$gene_id),
                       column="SYMBOL",
                       keytype="ENTREZID",
                       multiVals="first")
gr_genes$gene_name <- gene_symbols

# Quick check
head(gr_genes)

# --- 3) Find overlaps ---
hits <- findOverlaps(gr_variants, gr_genes, ignore.strand=TRUE)

# --- 4) Create overlap table ---
# Collapse ALT alleles so lengths match
alt_alleles <- sapply(mcols(gr_variants)$ALT[queryHits(hits)], function(x) {
  paste(as.character(x), collapse=",")
})

overlap_table <- data.frame(
  variant_id = names(gr_variants)[queryHits(hits)],
  chr        = as.character(seqnames(gr_variants)[queryHits(hits)]),
  pos        = start(gr_variants)[queryHits(hits)],
  ref        = as.character(mcols(gr_variants)$REF[queryHits(hits)]),
  alt        = alt_alleles,
  gene_id    = gr_genes$gene_id[subjectHits(hits)],
  gene_name  = gr_genes$gene_name[subjectHits(hits)],
  stringsAsFactors = FALSE
)

head(overlap_table)

write.csv(overlap_table,file= "SS91_VCF_Gene.csv")



