library(tidyverse)
library(readxl)
library(microbiomics)
library(vegan)
library(RColorBrewer)
library(cowplot)
library(ggbeeswarm)
library(Maaslin2)
library(rstatix)
library(openxlsx)

metadata_dna <- 
  read_xlsx("nipper_sample_metadata.xlsx", sheet = 2) %>%
  rename(sampleID = `Sequencing name`,
         stool_sample_id = `Sample ID`,
         country = `Stool Site`,
         mother_baby = `Mother/Baby`,
         preservation_method = `Sample preservation method`,
         extraction_protocol = `Extraction protocol`,
         extraction_site = `Extraction site`,
         family_id = `Participant ID`) %>%
  select(sampleID, stool_sample_id, country, mother_baby, preservation_method, extraction_protocol, extraction_site, family_id) %>%
  mutate(sample_origin = paste(family_id, mother_baby)) %>%
  filter(preservation_method == "RNA later",
         mother_baby != "NA")

metacyc_dna <- 
  read_tsv("nipper_metacyc_pathway_abundance_dna.txt") %>%
  rename(pathway = `# Pathway`) %>%
  pivot_longer(-pathway, names_to = "sampleID", values_to = "cpm") %>%
  mutate(sampleID = str_extract(sampleID, "neslig_[0-9]+")) %>%
  filter(sampleID %in% metadata_dna$sampleID) %>%
  filter(!(grepl("\\|", pathway))) 


metacyc_for_permanova <- 
  metacyc_dna %>%
  pivot_wider(names_from = pathway, values_from = cpm) %>%
  column_to_rownames("sampleID") %>%
  as.data.frame()

adonis2(metacyc_for_permanova ~ mother_baby + family_id, 
        metadata_dna,
        by = "margin")
# variance explained:
# mother baby = 57.9% (more than taxonomic profiles)

metacyc_for_maaslin <- 
  metacyc_dna %>%
  filter(!(grepl("^UN", pathway))) %>%
  pivot_wider(names_from = pathway, values_from = cpm) %>%
  column_to_rownames("sampleID") %>%
  as.data.frame()

metadata_dna_maaslin <-
  as.data.frame(metadata_dna)
rownames(metadata_dna_maaslin) <- metadata_dna_maaslin$sampleID

# Maaslin2 v0.99.18
maaslin_dna <- 
  Maaslin2(input_data = metacyc_for_maaslin, 
           input_metadata = metadata_dna_maaslin,
           output = "maaslin_dna",
           fixed_effects = c("mother_baby", "country"),
           random_effects = "family_id",
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4,
           transform = "none",
           normalization = "none")

dna_results <- 
  maaslin_dna$results %>%
  tibble() %>%
  filter(qval < 0.05,
         metadata == "mother_baby")
# 152 pathways differ between moms and babies in DNA

# total features tested:
dim(metacyc_for_maaslin)
# 466 pathways
# 112 filtered out (not present in 10% of samples)
466-112
# 354 pathways tested 

# 152 differentially abundant (q < 0.05)
152 / 354

dna_results %>%
  mutate(direction = case_when(coef > 0 ~ "Mother up",
                               coef < 0 ~ "Infant up")) %>%
  group_by(direction) %>%
  summarise(n = n())
# 48 pathways more abundant in moms
# 104 pathways more abundance in infants
dna_pathways_tested <- unique(maaslin_dna$results$feature)

tibble(metacyc_diversity = diversity(metacyc_for_permanova),
       sampleID = rownames(metacyc_for_permanova)) %>%
  left_join(metadata_dna) %>%
  ggplot(aes(x=mother_baby, y=metacyc_diversity)) +
  geom_boxplot() +
  geom_quasirandom() +
  theme_bw()

metacyc_dna %>%
  filter(grepl("UNMAPPED", pathway)) %>%
  left_join(metadata_dna) %>%
  wilcox_test(cpm ~ mother_baby)

metacyc_dna %>%
  filter(grepl("UNMAPPED", pathway)) %>%
  left_join(metadata_dna) %>%
  ggplot(aes(x=mother_baby, y= (1e6-cpm) / 1e6 )) +
  geom_boxplot(outlier.colour = NA, size = 0.2) +
  geom_quasirandom(dodge.width = 1, size = 1) +
  theme_bw() +
  ylab("Relative abundance") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("mapped_dna.pdf", w=1.5, h=2.5, useDingbats = F)
# p = 0.000325

metacyc_dna %>%
  filter(grepl("^UN", pathway)) %>%
  pivot_wider(names_from = pathway, values_from = cpm) %>%
  left_join(metadata_dna) %>%
  mutate(INTEGRATED = 1e6 - UNMAPPED - UNINTEGRATED,
         prop_integrated = INTEGRATED / (INTEGRATED + UNINTEGRATED)) %>%
  wilcox_test(prop_integrated ~ mother_baby)
# p = 0.165

metacyc_dna %>%
  filter(grepl("^UN", pathway)) %>%
  pivot_wider(names_from = pathway, values_from = cpm) %>%
  left_join(metadata_dna) %>%
  mutate(INTEGRATED = 1e6 - UNMAPPED - UNINTEGRATED) %>%
  ggplot(aes(x=mother_baby, y= INTEGRATED / ( INTEGRATED + UNINTEGRATED))) +
  geom_boxplot(outlier.colour = NA, size = 0.2) +
  geom_quasirandom(dodge.width = 1, size = 1) +
  theme_bw() +
  ylab("Proportion of the mapped reads") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("integrated_dna.pdf", w=1.5, h=2.5, useDingbats = F)

metacyc_dna %>%
  filter(grepl("^UN", pathway)) %>%
  pivot_wider(names_from = pathway, values_from = cpm) %>%
  left_join(metadata_dna) %>%
  mutate(INTEGRATED = 1e6 - UNMAPPED - UNINTEGRATED) %>%
  select(sampleID, mother_baby, INTEGRATED, UNINTEGRATED) %>%
  pivot_longer(cols = ends_with("ED"), names_to = "type", values_to = "value") %>%
  mutate(value = value / 1e6) %>%
  ggplot(aes(x=sampleID, y=value, fill = type)) +
  geom_bar(stat = "identity") + 
  facet_grid(~mother_baby, scales = "free_x") +
  theme_bw() +
  ylab("Proportion of mapped reads")

# RNA
metadata_rna <- 
  read_xlsx("nipper_sample_metadata.xlsx", sheet = 1) %>%
  rename(sampleID = `Sequencing name`,
         stool_sample_id = `Sample ID`,
         country = `Stool Site`,
         mother_baby = `Mother/Baby`,
         preservation_method = `Sample preservation method`,
         extraction_protocol = `Extraction protocol`,
         extraction_site = `Extraction site`,
         family_id = `Participant ID`) %>%
  select(sampleID, stool_sample_id, country, mother_baby, preservation_method, extraction_protocol, extraction_site, family_id) %>%
  mutate(sample_origin = paste(family_id, mother_baby)) %>%
  filter(mother_baby != "NA")

metacyc_rna <- 
  read_tsv("nipper_metacyc_pathway_abundance_rna.txt") %>%
  rename(pathway = `# Pathway`) %>%
  pivot_longer(-pathway, names_to = "sampleID", values_to = "cpm") %>%
  mutate(sampleID = str_extract(sampleID, "neslig_[0-9]+")) %>%
  filter(sampleID %in% metadata_rna$sampleID) %>%
  filter(!(grepl("^UN", pathway)),
         !(grepl("\\|", pathway)))

metacyc_rna_for_permanova <- 
  metacyc_rna %>%
  pivot_wider(names_from = pathway, values_from = cpm) %>%
  column_to_rownames("sampleID") %>%
  as.data.frame()

metadata_rna <- 
  metadata_rna %>%
  filter(sampleID %in% rownames(metacyc_rna_for_permanova))

metacyc_rna_for_permanova <- metacyc_rna_for_permanova[ metadata_rna$sampleID , ]
metadata_rna$sampleID == rownames(metacyc_rna_for_permanova)

adonis2(metacyc_rna_for_permanova ~ mother_baby + family_id, 
        metadata_rna,
        by = "margin")
# variance explained:
# mother baby = 27.0% (less than DNA profiles; functional or taxonomic)

metadata_rna_maaslin <-
  metadata_rna %>%
  column_to_rownames("sampleID") %>%
  as.data.frame()

# Maaslin2 v0.99.18
maaslin_rna <- 
  Maaslin2(input_data = metacyc_rna_for_permanova, 
           input_metadata = metadata_rna_maaslin,
           output = "maaslin_rna",
           fixed_effects = c("mother_baby", "country"),
           random_effects = "family_id",
           plot_scatter = T,
           min_prevalence = 0.1,
           cores = 4,
           transform = "none",
           normalization = "none")

rna_results <- 
  maaslin_rna$results %>%
  tibble() %>%
  filter(qval < 0.05,
         metadata == "mother_baby") 
# 40 pathways differ between moms and babies in DNA

sum(rna_results$feature %in% dna_results$feature)
# 38 of 40 pathways are also different in DNA

# 152 pathways differ between moms and babies in DNA

# total features tested:
dim(metacyc_rna_for_permanova)
# 494 pathways
# 109 filtered out (not present in 10% of samples)
494-109
# 385 pathways tested 

rna_pathways_tested <- unique(maaslin_rna$results$feature)

maaslin_results_dna <-
  maaslin_dna$results %>%
  tibble() %>%
  filter(qval < 0.05,
         metadata == "mother_baby") %>%
  select(feature, coef, stderr, pval, qval) %>%
  mutate(feature = str_replace(feature, "\\.\\.", "-"),
         feature = str_extract(feature, "[[:alnum:]\\.]+")) 


metacyc_dna_with_metadata %>%
  group_by(pathway, mother_baby) %>%
  summarise(mean_cpm = mean(cpm)) %>%
  pivot_wider(names_from = mother_baby, values_from = mean_cpm) %>%
  rename(mean_in_infants = Baby,
         mean_in_mothers = Mother) %>%
  mutate(feature = str_extract(pathway, "[[:alnum:]_-]+:"),
         feature = str_replace_all(feature, "-", "."),
         feature = str_remove(feature, ":")) %>%
  right_join(maaslin_results_dna) %>%
  select(-feature) %>%
  filter(!(is.na(pathway))) %>%
  arrange(qval) %>%
  write.xlsx(file = "Table_S1_1.xlsx")
  
maaslin_results_rna <-
  maaslin_rna$results %>%
  tibble() %>%
  filter(qval < 0.05,
         metadata == "mother_baby") %>%
  select(feature, coef, stderr, pval, qval) %>%
  mutate(feature = str_replace(feature, "\\.\\.", "-"),
         feature = str_extract(feature, "[[:alnum:]\\.]+")) 


metacyc_rna_with_metadata %>%
  group_by(pathway, mother_baby) %>%
  summarise(mean_cpm = mean(cpm)) %>%
  pivot_wider(names_from = mother_baby, values_from = mean_cpm) %>%
  rename(mean_in_infants = Baby,
         mean_in_mothers = Mother) %>%
  mutate(feature = str_extract(pathway, "[[:alnum:]_-]+:"),
         feature = str_replace_all(feature, "-", "."),
         feature = str_remove(feature, ":")) %>%
  right_join(maaslin_results_rna) %>%
  select(-feature) %>%
  filter(!(is.na(pathway))) %>%
  arrange(qval) %>%
  write.xlsx(file = "Table_S1_2.xlsx")



