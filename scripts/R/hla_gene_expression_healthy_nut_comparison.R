library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(NatParksPalettes)

# read in normalized gene x sample matrix from salmon
hlung_nut <- read.csv('data/nextflow_nut_heathy_lung_RNA/star_salmon/salmon.merged.gene_tpm.tsv', sep = '\t')
hlung_nut <- hlung_nut %>% select(-gene_id)
hlung_nut <- hlung_nut %>% select(-'PDX')

nutsample_list <- colnames(hlung_nut)[grepl('scram', colnames(hlung_nut))]
hsample_list <- colnames(hlung_nut)[grepl('ERS', colnames(hlung_nut))]

nutsamples <- length(grep('scram', colnames(hlung_nut))) %>% as.numeric()
hsamples <- length(grep('ERS', colnames(hlung_nut))) %>% as.numeric()

# sum together genes with the same name 
hlung_nut <- hlung_nut %>%
  group_by(gene_name) %>%
  reframe(across(everything(), ~ if (is.numeric(.)) sum(., na.rm = TRUE) else unique(.))) %>% 
  distinct(gene_name, .keep_all = TRUE)


# gene lists for genes of interest for comparison
HLACLI_genes <- c('B2M', 'HLA-A', 'HLA-B', 'HLA-C')

antigenprocessing_genes <- c('ERAP1', 'ERAP2', 'TAP1', 'TAPBP')

HLACLII_genes <- c('HLA-DMA','HLA-DMB','HLA-DOA','HLA-DOB','HLA-DPA1','HLA-DPA2',
                   'HLA-DPA3','HLA-DPB1','HLA-DPB2','HLA-DQA1', 'HLA-DQA2','HLA-DQB1',
                   'HLA-DQB1-AS1','HLA-DQB2','HLA-DQB3','HLA-DRA','HLA-DRB1','HLA-DRB5',
                   'HLA-DRB6','HLA-DRB9')

genes_of_interest <- c(HLACLI_genes, antigenprocessing_genes, HLACLII_genes)


# to do things by expression ----------------------------------------------

# pull out genes of interest, and convert data to long format for plotting
# do this one list as a time


# Process HLA Class I genes
long_CLI <- data.frame()
for (gene in HLACLI_genes) {
  temp <- hlung_nut %>% filter(gene_name == gene)
  long_temp <- pivot_longer(temp, cols = -c(gene_name), values_to = 'TPM', names_to = 'sample')
  # Assign group based on whether the sample is in nutsamples_list
  long_temp$group <- ifelse(long_temp$sample %in% nutsample_list, 'NC', 'Normal')
  long_temp$gene_group <- "HLA-Class I"
  long_CLI <- rbind(long_CLI, long_temp)
}

# Process APM genes
long_APM <- data.frame()
for (gene in antigenprocessing_genes) {
  temp <- hlung_nut %>% filter(gene_name == gene)
  long_temp <- pivot_longer(temp, cols = -c(gene_name), values_to = 'TPM', names_to = 'sample')
  long_temp$group <- ifelse(long_temp$sample %in% nutsample_list, 'NC', 'Normal')
  long_temp$gene_group <- "Antigen Processing" 
  long_APM <- rbind(long_APM, long_temp)
}

# Process HLA Class II genes
long_CLII <- data.frame()
for (gene in HLACLII_genes) {
  temp <- hlung_nut %>% filter(gene_name == gene)
  long_temp <- pivot_longer(temp, cols = -c(gene_name), values_to = 'TPM', names_to = 'sample')
  long_temp$group <- ifelse(long_temp$sample %in% nutsample_list, 'NC', 'Normal')
  long_temp$gene_group <- "HLA-Class II" 
  long_CLII <- rbind(long_CLII, long_temp)
}

# lump together HLA Class II genes for now
long_CLII_grouped <- long_CLII %>% mutate(gene_name = 'HLA-Class II') 

# Combine dataframe for each cell list into single df for plotting
long_all_genes_grouped <- rbind(long_APM, long_CLI, long_CLII_grouped)

# set order for the gene group for the plot 
long_all_genes_grouped$gene_group <- factor(long_all_genes_grouped$gene_group, levels = c("Antigen Processing", "HLA-Class I", "HLA-Class II"))

# set order for the gene_names for the plot
HLACLII_genes <- 'HLA-Class II'
long_all_genes_grouped$gene_name <- factor(long_all_genes_grouped$gene_name, levels = c(HLACLI_genes, antigenprocessing_genes, HLACLII_genes) )

long_all_genes_grouped <- long_all_genes_grouped %>% mutate(sample = gsub("^X|scram_rep\\d+|NUT_rep\\d+", "", sample))

long_all_genes_grouped_agg <- long_all_genes_grouped %>% 
  group_by(sample, gene_name, group) %>% 
  summarize(
    mean = mean(TPM)
  ) %>%
  ungroup()

long_all_genes_grouped_agg$group <- factor(long_all_genes_grouped_agg$group, levels = c('Normal', 'NC'))

colors <- NatParksPalettes::natparks.pals('Banff', n = 2)

# Visualize the combined data
healthy_HLA_APM_comp <- ggplot(long_all_genes_grouped_agg, aes(x = gene_name, y = mean, color = group)) +
  geom_point(position = position_dodge(width = 0.7), size = 3, alpha = 0.6) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.7, position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_y_log10() +  # Log transformation while keeping actual values on the axis
  scale_color_manual(values = colors,
                     labels = c("Normal" = "Healthy Lung", "NC" = "NUT Carcinoma"),
                     breaks = c("Normal", "NC")) +
  labs(title = "", x = "", y = "Mean TPM") +  
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

#ggsave('results/rna_seq_figures/Fagerberg/Fagerberg_RNA_classI_APM_Healthy_Lung_NC_TPM.pdf', plot = healthy_HLA_APM_comp, dpi = 500)

wilcox_res <- long_all_genes_grouped %>% group_by(gene_name) %>% summarise(p.val = wilcox.test(TPM ~ group)$p.value)
wilcox_res <- wilcox_res %>% mutate(adj.p.val = p.adjust(p.val, method = "BH"))  

# view HLA-Class II genes separated out by gene_name
hla_cl2_separated <- ggplot(long_CLII, aes(x = gene_name, y = log(TPM + 1), color = group)) +
  geom_point(position = position_dodge(width = 0.7), size = 3, alpha = 0.6) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.7, position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_color_manual(values = colors,
                     labels = c("Normal" = "Healthy Lung", "NC" = "NUT Carcinoma" ),
                     breaks = c("Normal", "NC")) +
  labs(title = "", x = "", y = "log(TPM + 1)") + 
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
