library(tidyverse)
library(readxl)
library(microbiomics)
library(vegan)
library(RColorBrewer)
library(cowplot)
library(ggbeeswarm)
library(phangorn)

metadata_dna <- 
  read_xlsx("../metadata//nipper_sample_metadata.xlsx", sheet = 2) %>%
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

load("strain_haplotype_dist.RData")

comparisons_metadata <-
  comparisons %>%
  filter(sample1 %in% metadata_dna$sampleID,
         sample2 %in% metadata_dna$sampleID) %>%
  inner_join(metadata_dna %>% 
               rename(sample1 = sampleID,
                      sample1_origin = sample_origin,
                      sample1_family = family_id) %>%
               select(sample1, sample1_origin, sample1_family)) %>%
  inner_join(metadata_dna %>% 
               rename(sample2 = sampleID,
                      sample2_origin = sample_origin,
                      sample2_family = family_id) %>%
               select(sample2, sample2_origin, sample2_family)) %>%
  mutate(comparison_category = case_when(sample1_family == sample2_family ~ "mother-infant comparison",
                                         TRUE ~ "other comparisons")) %>%
  mutate(species = str_remove(species, ".fasta"))
  

comparisons_norm <- 
  comparisons_metadata %>%
  group_by(species) %>%
  summarise(median_dist = median(dist)) %>%
  left_join(comparisons_metadata) %>%
  mutate(norm_dist = dist / median_dist)

comparisons_norm %>%
  filter(norm_dist < 3) %>%
  ggplot(aes(x=norm_dist, fill = comparison_category)) + 
  geom_density(alpha = 0.5) +
  theme_cowplot() +
  xlab("Normalized strain haplotype difference") +
  geom_vline(xintercept = 0.2, linetype = "dashed")
ggsave("strain_haplotype_densities.pdf", h = 4, w= 6)

# conservative threshold 0.2
transmission_events <- 
  comparisons_norm %>%
  filter(comparison_category == "mother-infant comparison",
         norm_dist <= 0.2) 
View(transmission_events)
# 52 transmission events
# 51 when excluding 1 phage

transmission_events %>%
  filter(!(grepl("phage", species))) %>%
  distinct(sample1_family, species) %>%
  mutate(species = str_remove(species, "s__"),
         species = str_replace_all(species, "_", " ")) %>%
  group_by(species) %>%
  summarise(`Number of transmission events` = n()) %>%
  arrange(-`Number of transmission events`) %>%
  write.xlsx(file = "Table_1.xlsx")
# 30 bacterial species

transmission_events %>%
  filter(!(grepl("phage", species))) %>%
  mutate(genus = str_extract(species, "s__[:alnum:]+")) %>%
  group_by(genus) %>%
  summarise(n = n()) %>%
  arrange(-n) 
# from 7 genera. Most transmission events (36 of 52) in Bacteroides
35/51
# 6 bifidobacteria
6/51

transmission_events %>%
  filter(grepl("phage", species)) %>%
  select(species)
# one tentative phage transmission; Enterobacteria phage HK022

transmission_events %>%
  select(species, sample1_family) %>%
  write_delim(delim = "\t", path = "maternally_transmitted_strains.tsv")

# Phylogenetic tree of B. ovatus
seq <- read.dna("s__Bacteroides_ovatus.fasta", format="fasta")

ovatus_phyDat <- phyDat(seq, type = "DNA", levels = NULL)
ovatus_phyDat <- subset(ovatus_phyDat, subset = str_extract(names(ovatus_phyDat), "neslig_[0-9]+") %in% metadata_dna$sampleID)

dna_dist <- dist.ml(ovatus_phyDat, model="JC69")
dim(dna_dist)

dna_dist <- filter_dist_outliers(as.matrix(dna_dist), quant=0.8, n_neighbour = 2)
dim(dna_dist)

# NJ phylogeny
ovatus_nj <- nj(dna_dist)
fit_nj <- pml(ovatus_nj, data=ovatus_phyDat, model = "K80")
fit_nj_optim <- optim.pml(fit_nj, model = "K80", optNni = T)

fit_nj_optim$logLik

library(ggtree)
p <- ggtree(fit_nj_optim$tree, ladderize = T) +
  geom_treescale(x=0, y=0, width = 0.001)

p <- p + geom_point(data = p$data %>%
                      filter(!(is.na(label))) %>%
                      mutate(sampleID = str_extract(label, "neslig_[0-9]+")) %>%
                      inner_join(metadata_dna),
                    aes(x=x,y=y, fill=family_id, shape=mother_baby), size = 3)+
  scale_fill_manual(values = c(brewer.pal(7,"Set2"),brewer.pal(5,"Set1"))) + 
  scale_shape_manual(values = c(21,23)) +
  coord_flip() +
  scale_x_reverse()
p
ggsave("ovatus_phylogeny.pdf", h=4, w=2.5, useDingbats = F)

