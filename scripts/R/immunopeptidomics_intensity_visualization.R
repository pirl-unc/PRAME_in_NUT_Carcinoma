library(dplyr)
library(ggplot2)
library(here)
library(UpSetR)
library(tidyr)


# western quantification dataframe

densitometry <- data.frame(cell_line = c('10-15','10-15','14169','14169','PER403','PER403','JCM1','JCM1'),
                           condition = rep(c('scram','NUT'), 4),
                           PRAME_GAPDH_ratio = c(1.36, 0.078, 3.22, 2.23, 0.63, 0.16, 2.32, 0.015))

densitometry$condition <- factor(densitometry$condition, levels = c("scram", "NUT"))


ggplot(densitometry, aes(x = cell_line, y = PRAME_GAPDH_ratio, fill = condition)) +
  geom_col(position = "dodge") + 
  scale_fill_manual(values = c("#35b779", "#31688e"),
                    labels = c('siRNA NUT', 'siRNA scram')) + 
  theme_classic() +
  labs(x = "Cell Line", y = "PRAME/GAPDH Ratio", fill = "Condition") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "lines"),
        legend.position = c(0.85, 0.85),
        plot.margin = margin(7, 7, 7, 10))  


# MS intensity of peptides

ctas <- read.csv('data/source_data/cell_lists/cancer-testis-antigens.csv')

ms_set1 <- read.csv('data/source_data/proteomics/cell_lines_1015_797_14169_biognosys_intensity_results.csv', na.strings = c("", "NA", "#N/A"), check.names = F)
ms_set2 <- read.csv('data/source_data/proteomics/cell_lines_PER403_JCM1_PDX_biognosys_intensity_results.csv', na.strings = c("", "NA", "#N/A"), check.names = F)


# update NUTM1 gene name to be NUTM1 since it is currently differentiated by breakpoint but we don't need to know that for this purpose
ms_set1$Gene_ID <- sub("NUTM1.*", "NUTM1", ms_set1$Gene_ID)
ms_set2$Gene_ID <- sub("NUTM1.*", "NUTM1", ms_set2$Gene_ID)

# update NUTM1 protein name to be NUT since it is currently differentiated by breakpoint but we don't need to know that for this purpose
ms_set1$Protein_Name <- ifelse(grepl("BRD4_HUMAN|BRD3_HUMAN|BRD2_HUMAN", ms_set1$Protein_Name), "NUT", ms_set1$Protein_Name)
ms_set2$Protein_Name <- ifelse(grepl("BRD4_HUMAN|BRD3_HUMAN|BRD2_HUMAN", ms_set2$Protein_Name), "NUT", ms_set2$Protein_Name)

# Identify metadata columns
metadata_cols <- intersect(names(ms_set1), names(ms_set2))

# Perform full join, keeping only one copy of metadata columns but all intensity values
ms_all <- full_join(ms_set1, ms_set2, by = metadata_cols)


# Update Gene_ID in nut_peptides based on Peptide_Sequence using mutate and case_when
ms_all <- ms_all %>%
  mutate(Gene_ID = case_when(
    Peptide_Sequence == 'VFDPIGHF' ~ 'BRD4',
    Peptide_Sequence == 'SPPALHNAL' ~ 'BRD4',
    Peptide_Sequence == 'RVVLKTLWK' ~ 'BRD4',
    Peptide_Sequence == 'RLAELQEQL' ~ 'BRD4',
    Peptide_Sequence == 'KMPDEPVEA' ~ 'BRD3',
    Peptide_Sequence == 'FAADVRLMF' ~ 'BRD3',
    Peptide_Sequence == 'DVYENFRQW' ~ 'NUTM1',
    Peptide_Sequence == 'AVVSPPALHNA' ~ 'BRD4',
    Peptide_Sequence == 'AAYAWPFYK' ~ 'BRD4',
    TRUE ~ Gene_ID  # Retain existing Gene_ID if no match is found
  ))

# Check the result
head(nut_peptides)


# convert NAs to 0's as these are peptide that were not detected
# there are true NA in the dataset and also #N/As so need to replace both
ms_all[is.na(ms_all)] <- 0

length(unique(ms_all$Peptide_ID))


# make list of peptides per cell line
cells797 <- ms_all %>% filter(`TC-797` > 0.001) %>% select(Peptide_ID, `TC-797`)
cells14169 <- ms_all %>% filter(`14169` > 0.001) %>% select(Peptide_ID, `14169`)
cells1015 <- ms_all %>% filter(`1015` > 0.001) %>% select(Peptide_ID, `1015`)
cellsPER403 <- ms_all %>% filter(`PER-403` > 0.001) %>% select(Peptide_ID, `PER-403`)
cellsJCM1 <- ms_all %>% filter(JCM1 > 0.001) %>% select(Peptide_ID, JCM1)
cellsPDX <- ms_all %>% filter(`PDX` > 0.001) %>% select(Peptide_ID, `PDX`)


proteins <- list(
  "clTC797" = cells797$Peptide_ID, 
  "cl14169" = cells14169$Peptide_ID, 
  "cl1015" = cells1015$Peptide_ID, 
  "clPER403" = cellsPER403$Peptide_ID, 
  "clJCM1" = cellsJCM1$Peptide_ID, 
  "PDX" = cellsPDX$Peptide_ID
)

upset_plt <- upset(fromList(proteins),
                   order.by = "freq",
                   nsets = 6,
                   intersections = list(
                     list("clTC797", "cl14169", "cl1015", "clPER403","clJCM1", "PDX"),
                     list("clTC797"),
                     list("cl14169"),
                     list("cl1015"),
                     list("clPER403"),
                     list("clJCM1"),
                     list("PDX")),
                   mainbar.y.label = "Number of Unique Peptides", 
                   sets.x.label = "Peptides Per Sample",
                   decreasing = TRUE,
                   main.bar.color = 'black',
                   sets.bar.color = 'grey',
                   matrix.color = 'black',
                   shade.color =  'grey',
                   text.scale = c(1.5, 1.3, 1, 1.5, 1.3, 2))


# if they are BRD3 or BRD4 peptide, only keep if the peptide matches the fusion (canonical BRD expression it not targetable)

# PDX is BRD3, remove BRD4 peptide from it
ms_all[ms_all$Peptide_Sequence == 'AVVSPPALHNA','PDX'] <- 0

# 797 and PER403 are BRD4, remove BRD3 peptide
ms_all[ms_all$Peptide_Sequence == 'KMPDEPVEA', c('TC-797', 'PER-403')] <- 0

# JCM1 is BRD4, remove BRD3 peptide
ms_all[ms_all$Peptide_Sequence == 'FAADVRLMF', 'JCM1'] <- 0

ctas <- ctas$Symbol %>% as.vector() %>% c('BRD4', 'BRD3')

cta_peptides <- ms_all[ms_all$Gene_ID %in% ctas, ]

cta_peptides <- cta_peptides[rowSums(cta_peptides[8:13] > 0) > 0, ]


# visualize all CTAs found

cta_peptides_long <- cta_peptides %>%
  select(Gene_ID, Peptide_Sequence, `TC-797`, `1015`, `14169`, `PER-403`, `JCM1`, `PDX`) %>%
  pivot_longer(
    cols = c(`TC-797`, `1015`, `14169`, `PER-403`, `JCM1`, `PDX`), 
    names_to = "sample",  
    values_to = "intensity"  
  )

cta_peptides_long$log2Intensity <- log2(cta_peptides_long$intensity)
cta_peptides_long$alpha <- ifelse(is.finite(cta_peptides_long$log2Intensity), 1, 0)
cta_peptides_long$log2Intensity <- ifelse(is.finite(cta_peptides_long$log2Intensity), cta_peptides_long$log2Intensity, 0)

gene_id_order <- c('PRAME','NUTM1','BRD4','BRD3','CTAGE1','MAGEA1','CT83','ACTL8')
cta_peptides_long$Gene_ID <- factor(cta_peptides_long$Gene_ID, levels = gene_id_order)
cta_peptides_long <- cta_peptides_long %>%
  mutate(Peptide_Sequence = factor(Peptide_Sequence, 
                                   levels = unique(Peptide_Sequence[order(Gene_ID, Peptide_Sequence, decreasing=TRUE)])))
samples_order <- c('PER-403', '1015','14169', 'JCM1', 'TC-797', 'PDX')
cta_peptides_long$sample <- factor(cta_peptides_long$sample, levels = samples_order)

ggplot(cta_peptides_long, aes(x = factor(sample), y = Peptide_Sequence)) +
  geom_point(aes(fill = log2Intensity, alpha = as.factor(alpha)), 
             shape = 21, size = 4, stroke = 1,
             show.legend = c(alpha = FALSE)) +  
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "blue", high = "red", 
                      name = "Log2(Intensity)",
                      limits = c(min(cta_peptides_long$log2Intensity, na.rm = TRUE), 
                                 max(cta_peptides_long$log2Intensity, na.rm = TRUE))) + 
  theme_minimal() +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, color='black'),
        axis.text.x = element_text(size = 12, color='black'),
        axis.text.y = element_text(size = 10, color='black')) +
  labs(x = "", y = "Peptide", fill = "Log2(Intensity)")
