# load libraries  ---------------------------------------------------------
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(here)


# load data ---------------------------------------------------------------
# gene x sample tpm expression df
gene_tpm_all <- read.csv('work/data/nextflow_nut_healthy_lung_RNA/star_salmon/salmon.merged.gene_tpm.tsv', sep = "\t", check.names = F) %>% select(-gene_id)

# curated set of CTAs
cta_df <- read.csv('work/data/cell_lists/cancer-testis-antigens.csv')
CTAs <- cta_df$Symbol

# sum together genes with same names
gene_tpm_all <- gene_tpm_all %>%
  group_by(gene_name) %>%
  reframe(across(everything(), ~ if (is.numeric(.)) sum(., na.rm = TRUE) else unique(.))) %>% 
  distinct(gene_name, .keep_all = TRUE)

# move gene_name column to rownames
gene_tpm_all <- gene_tpm_all %>% tibble::column_to_rownames(var='gene_name')

# remove the Faegerberg healthy lung samples
gene_tpm_nc <- gene_tpm_all[ , grepl('scram|PDX', colnames(gene_tpm_all))]

colnames(gene_tpm_nc) <- gsub('X', "", colnames(gene_tpm_nc))
rename(gene_tpm_nc, "PDX" = "PD")


gene_tpm_nc_cta <- gene_tpm_nc[rownames(gene_tpm_nc) %in% CTAs, ]
gene_tpm_nc_cta_filtered <- gene_tpm_nc_cta[rowMeans(gene_tpm_nc_cta) > 1,]

cell_lines <- factor(
  c(rep('10-15', 4), rep('14169', 4), rep('TC-797', 4), rep('JCM1', 3), 'PDX', rep('PER403', 4)),
  levels = c('PER403', '14169', '10-15', 'JCM1', 'TC-797', 'PDX'))

col_fun <- colorRamp2(c(0, max(gene_tpm_nc_cta_filtered)), c("blue", "red"))

pdf('work/results/rna_seq/nut_cell_lines/Figure2c_nut_carcinoma_cell_lines_top_cta_expression.pdf', width = 6.5, height = 4)
heatmap_rna_nc_cell_lines <- Heatmap(as.matrix(gene_tpm_nc_cta_filtered),
        name = 'TPM',
        col = col_fun,
        show_column_names = FALSE,
        column_split = cell_lines,
        column_title_side = "bottom",
        cluster_columns = TRUE,
        row_title = paste0('Top ' , nrow(gene_tpm_nc_cta_filtered) , ' Expressed CTAs'),
        row_title_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 8, fontface = "italic"),
        column_title_gp = gpar(fontsize = 11))
dev.off()



