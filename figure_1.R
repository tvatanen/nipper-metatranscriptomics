library(tidyverse)
library(readxl)
library(microbiomics)
library(vegan)
library(RColorBrewer)
library(cowplot)
library(TSP)
library(ggbeeswarm)

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

metaphlan_species <- read_metaphlan_table("nipper_metaphlan_profiles.tsv")

metaphlan_species_long <- 
  metaphlan_species %>%
  rownames_to_column("sampleID") %>%
  pivot_longer(-sampleID, names_to = "taxon", values_to = "relative_abundance") %>%
  mutate(sampleID = str_extract(sampleID, "neslig_[0-9]+"),
         species = str_extract(taxon, "s__.+"),
         species = str_remove(species, "s__"),
         species = str_replace(species, "_", " "))

metaphlan_species_rna_later <- 
  metaphlan_species_long %>%
  filter(sampleID %in% metadata_dna$sampleID)

metaphlan_for_mds <- 
  metaphlan_species_rna_later %>%
  select(-taxon) %>%
  pivot_wider(names_from = species, values_from = relative_abundance) %>%
  as.data.frame()

rownames(metaphlan_for_mds) <- metaphlan_for_mds$sampleID
metaphlan_for_mds <- metaphlan_for_mds[ , -1 ]
set.seed(1523345)
mds_res <- metaMDS(metaphlan_for_mds)

tibble(sampleID = rownames(mds_res$points),
       MDS1 = mds_res$points[ , 1],
       MDS2 = mds_res$points[ , 2]) %>%
  inner_join(metadata_dna) %>%
  ggplot(aes(x=MDS1, y=MDS2, fill = mother_baby)) + 
  geom_point(shape = 21) + 
  coord_equal() +
  theme_bw() +
  scale_fill_manual(values = brewer.pal(3, "Set2")[c(2,3)]) +
  xlab("PCo1") +
  ylab("PCo2") + 
  scale_x_continuous(breaks = c()) +
  scale_y_continuous(breaks = c()) +
  theme(legend.position = c(0.1, 0.9),
        legend.background = element_rect(color = "black", size=0.2, linetype="solid"),
        legend.title = element_blank())
ggsave("pcoa.pdf", w = 2, h = 2, useDingbats = F)

metadata_dna$sampleID == rownames(mds_res$points)

adonis2(metaphlan_for_mds ~ mother_baby + family_id, 
        metadata_dna,
        by = "margin")
# variance explained:
# mother baby = 31.4%

# infant bar plots
infant_species <-
  metaphlan_species_rna_later %>%
  inner_join(metadata_dna) %>%
  filter(mother_baby == "Baby") %>%
  group_by(species) %>%
  summarise(mean_rel_ab = mean(relative_abundance)) %>%
  top_n(14, wt = mean_rel_ab) %>%
  pull(species)

metaphlan_profiles_infants <- 
  metaphlan_species_rna_later %>%
  inner_join(metadata_dna) %>%
  filter(mother_baby == "Baby",
         species %in% infant_species) %>%
  group_by(sampleID) %>%
  summarise(relative_abundance = 1 - sum(relative_abundance)) %>%
  mutate(species = "Other") %>%
  bind_rows(metaphlan_species_rna_later %>%
              filter(species %in% infant_species)) %>%
  inner_join(metadata_dna) %>%
  filter(mother_baby == "Baby") 

metaphlan_wide <-
  metaphlan_profiles_infants %>%
  select(sampleID, species, relative_abundance) %>%
  pivot_wider(names_from = species, values_from = relative_abundance) %>%
  column_to_rownames("sampleID")

sample_dist <- vegdist(metaphlan_wide)
tsp <- TSP(sample_dist)
order <- solve_TSP(tsp)


metaphlan_profiles_infants %>%
  mutate(sampleID = factor(sampleID, levels = labels(order)),
         species = factor(species, levels = c(infant_species, "Other"))) %>%
  ggplot(aes(x=sampleID, y = relative_abundance, fill = species)) +
  geom_bar(stat = "identity") +
  theme_cowplot() + 
  scale_fill_manual(values = c(brewer.pal(7, "Set2"), brewer.pal(7, "Set3"), "gray"),
                    name = "") +
  ylab("Relative abundance") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("infant_profiles.pdf", w = 6, h = 5)
  
# maternal bar plots
maternal_species <-
  metaphlan_species_rna_later %>%
  inner_join(metadata_dna) %>%
  filter(mother_baby == "Mother") %>%
  group_by(species) %>%
  summarise(mean_rel_ab = mean(relative_abundance)) %>%
  top_n(14, wt = mean_rel_ab) %>%
  pull(species)

metaphlan_profiles_moms <- 
  metaphlan_species_rna_later %>%
  inner_join(metadata_dna) %>%
  filter(mother_baby == "Mother",
         species %in% maternal_species) %>%
  group_by(sampleID) %>%
  summarise(relative_abundance = 1 - sum(relative_abundance)) %>%
  mutate(species = "Other") %>%
  bind_rows(metaphlan_species_rna_later %>%
              filter(species %in% maternal_species)) %>%
  inner_join(metadata_dna) %>%
  filter(mother_baby == "Mother") 

metaphlan_wide <-
  metaphlan_profiles_moms %>%
  select(sampleID, species, relative_abundance) %>%
  pivot_wider(names_from = species, values_from = relative_abundance) %>%
  column_to_rownames("sampleID")

sample_dist <- vegdist(metaphlan_wide)
tsp <- TSP(sample_dist)
order <- solve_TSP(tsp)

metaphlan_profiles_moms %>%
  mutate(sampleID = factor(sampleID, levels = labels(order)),
         species = factor(species, levels = c(maternal_species, "Other"))) %>%
  ggplot(aes(x=sampleID, y = relative_abundance, fill = species)) +
  geom_bar(stat = "identity") +
  theme_cowplot() + 
  scale_fill_manual(values = c(brewer.pal(7, "Set2"), brewer.pal(7, "Set3"), "gray"),
                    name = "") +
  ylab("Relative abundance") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("maternal_profiles.pdf", w = 6.5, h = 5)

alpha_div <- 
  metaphlan_species_rna_later %>%
  select(-taxon) %>%
  pivot_wider(names_from = species, values_from = relative_abundance) %>%
  column_to_rownames("sampleID") %>%
  vegan::diversity()

tibble(alpha_div = alpha_div,
       sampleID = names(alpha_div)) %>%
  inner_join(metadata_dna) %>%
  ggplot(aes(y=alpha_div, x =mother_baby)) +
  geom_boxplot(size = 0.2) +
  geom_quasirandom(size = 1) + 
  theme_bw() + 
  ylab("Shannon's diversity") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("alpha_div.pdf", w=1.5, h=2.5, useDingbats = F)