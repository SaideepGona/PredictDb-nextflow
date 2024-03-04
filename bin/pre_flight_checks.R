## Prechecks

# Prechecks should be:

# 1.) Formatting of each file should match documentation
# 1.1) TSV files
# 1.2) Expression with NAME column

# 2.) Gene overlap between expression data and gene annotations should be high relative to the number of genes in the expression data.

# 3.) Sample overlap between expression data and genotype data should be high relative to the number of samples in the expression data.

# 4.) SNP overlap between genotype data and SNP annotations should be high relative to the number of SNPs in the genotype data.

suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages((library(reshape2)))
suppressMessages(library(methods))
suppressMessages(library(tools))
suppressMessages(library(glue))
"%&%" <- function(a,b) paste(a,b, sep = "")

#! /usr/bin/env Rscript
# suppressMessages(library(optparse))

# # create options
# option_list <- list(
#   make_option(c("--snp_annotation"), type="character", default=NULL,
#               help="The snp annotation file for the chromosome you are processing",
#               metavar="character"),
#   make_option(c("--gene_annotation"), type="character", default=NULL,
#               help="The gene annotation file for the chromosome you are processing",
#               metavar="character"),
#   make_option(c("--genotype_file"), type="character", default=NULL,
#               help="The genotype file for the chromosome you are processing",
#               metavar="character"),
#   make_option(c("--gene_expression"), type="character", default=NULL,
#               help="The gene expression file",
#               metavar="character")
# )

# opt_parser <- OptionParser(option_list=option_list)
# args <- parse_args(opt_parser)


# snp_annot_file <- args$snp_annotation
# gene_annot_file <- args$gene_annotation
# genotype_file <- args$genotype_file
# expression_file <- args$gene_expression

snp_annot_file <- "/beagle3/haky/users/saideep/projects/aracena_modeling/Inputs/ref/rsIDs.for.predixcan_b37.txt"
gene_annot_file <- "/beagle3/haky/users/saideep/projects/aracena_modeling/Inputs/ref/gencode.v44lift37.annotation.gtf.parsed.txt"
genotype_file <- "/beagle3/haky/users/saideep/projects/aracena_modeling/Inputs/QTL_mapping/SNP_genotypes_b37.txt"
expression_file <- "/beagle3/haky/users/saideep/projects/aracena_modeling/Inputs/counts_matrices/final_testing/fully_preprocessed_flu_expression_for_predictDB_3_7_lnormTRUE_beforeTRUE_scale1.txt"

maf=0.01
n_times=3
n_k_folds=10
cis_window=1e6
alpha=0.5


snp_annot <- read.table(snp_annot_file, header = T, stringsAsFactors = F)
print("Dimensions of SNP annotation file:")
print(dim(snp_annot))

gene_types=c('protein_coding', 'pseudogene', 'lincRNA')
gene_df <- read.table(gene_annot_file, header = TRUE, stringsAsFactors = FALSE)
print("Dimensions of gene annotation file:")
print(dim(gene_df))

gt_df <- read.table(genotype_file, header = T, stringsAsFactors = F, row.names = 1)
print("Dimensions of genotype file:")
print(dim(gt_df))

expr_df <- as.data.frame(read.table(expression_file, header = T, stringsAsFactors = F, row.names = 1))
colnames(expr_df) <- file_path_sans_ext(colnames(expr_df)) # remove gene version number
print("Dimensions of gene expression file:")
print(dim(expr_df))


# 1.) Formatting of each file should match documentation

if (length(colnames(snp_annot)) == 1) {
  stop("SNP annotation file only reading a single column. Check file format.")
} 
if (length(colnames(gene_df)) == 1) {
  stop("Gene annotation file only reading a single column. Check file format.")
}
if (length(colnames(gt_df)) == 1) {
  stop("Genotype file only reading a single column. Check file format.")
}
if (length(colnames(expr_df)) == 1) {
  stop("Gene expression file only reading a single column. Check file format.")
}


# 2.) Gene overlap between expression data and gene annotations should be high relative to the number of genes in the expression data.

gene_overlap <- length(intersect(rownames(expr_df), gene_df$gene_id)) / length(rownames(expr_df))
print("Proportio of expression genes in gene annotation file:")
print(gene_overlap)

# 3.) Sample overlap between expression data and genotype data should be high relative to the number of samples in the expression data.

sample_overlap <- length(intersect(colnames(expr_df), colnames(gt_df))) / length(colnames(expr_df))
print("Proportion of expression samples in genotype file:")
print(sample_overlap)

# 4.) SNP overlap between genotype data and SNP annotations should be high relative to the number of SNPs in the genotype data.

# print(head(gt_df))
# print(head(snp_annot))

snp_overlap <- length(intersect(rownames(gt_df), snp_annot$varID)) / length(rownames(gt_df))
print("Proportion of SNPs in genotype file in SNP annotation file:")
print(snp_overlap)


