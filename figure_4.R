library(tidyverse)
library(RColorBrewer)

# load gene mgx and mtx profiles of maternally transmitted strains
load("combined_gene_quantifications.RData")

combined_genes %>%
  filter(!(grepl("unknown", gene_family))) %>% # this filters out UniRef90_unknown rows
  group_by(species, stool_sample_id) %>%
  summarise(n_genes = n()) %>%
  arrange(-n_genes)
# Total of 48 sample-species pairs with mtx and mgx data

combined_genes %>%
  filter(cpm_mtx > 0,
         cpm_mgx > 0,
         !(grepl("unknown", gene_family))) %>%
  group_by(species, family_id) %>%
  summarise(n_genes = n()) %>%
  arrange(-n_genes) %>%
  View()

combined_genes %>%
  filter(cpm_mtx > 0,
         cpm_mgx > 0,
         !(grepl("unknown", gene_family))) %>%
  group_by(species, stool_sample_id) %>%
  summarise(n_genes = n()) %>%
  arrange(-n_genes) %>%
  View()
# 45 sample-species pairs with non-zero mgx and mtx data
# 959 - 6,393 transcribed genes detected per population

combined_genes %>%
  filter(cpm_mtx > 0,
         cpm_mgx > 0,
         !(grepl("unknown", gene_family))) %>%
  group_by(species, family_id, mother_baby) %>%
  summarise(n_genes = n()) %>%
  arrange(-n_genes) %>%
  pivot_wider(names_from = mother_baby, values_from = n_genes) %>%
  mutate(difference = Mother - Baby) %>%
  View()
# Three strains where the infant data is missing
# Filter these out for the final strain list for transcription change experiments

final_strain_list <-
  combined_genes %>%
  filter(cpm_mtx > 0,
         cpm_mgx > 0,
         !(grepl("unknown", gene_family))) %>%
  group_by(species, family_id, mother_baby) %>%
  summarise(n_genes = n()) %>%
  arrange(-n_genes) %>%
  pivot_wider(names_from = mother_baby, values_from = n_genes) %>%
  mutate(difference = Mother - Baby) %>%
  filter(!(is.na(difference))) %>%
  select(species, family_id)

combined_genes_final <- 
  inner_join(combined_genes, final_strain_list) %>%
  filter(!(grepl("unknown", gene_family)))

combined_genes_final %>%
  group_by(species, family_id) %>%
  summarise(n_genes = n()) %>%
  arrange(-n_genes) %>%
  View()
# 21 strains in the final analysis

combined_genes_final %>%
  filter(cpm_mtx > 0,
         cpm_mgx > 0) %>%
  group_by(species, stool_sample_id) %>%
  summarise(n_genes = n()) %>%
  arrange(-n_genes) %>%
  View()
# 42 sample-species pairs with non-zero mgx and mtx data
# 959 - 6,393 transcribed genes detected per population

# median expression per genome per sample
median_experession_ratio <- 
  combined_genes_final %>%
  filter(cpm_mtx > 0,
         cpm_mgx > 0) %>%
  mutate(cpm_ratio = cpm_mtx / cpm_mgx) %>%
  group_by(species, stool_sample_id) %>%
  summarise(median_expression = median(cpm_ratio))

normalized_expression <- 
  combined_genes_final %>%
  filter(cpm_mgx > 0) %>% # allow mtx data to contain zeroes
  left_join(median_experession_ratio) %>%
  mutate(cpm_ratio = cpm_mtx / cpm_mgx,
         norm_expression = cpm_ratio / median_expression)

# check that everything looks good for b longum in sample NP46 
panel_a_data <- 
  normalized_expression %>%
  filter(stool_sample_id == "NP46",
         grepl("longum", species)) 

normalized_expression %>%
  select(species, stool_sample_id, median_expression)
# median expression ranges between 0.34 (B dorei) and 1.08

ggplot(panel_a_data, aes(x=cpm_mgx, y=cpm_mtx)) +
  geom_point(alpha = 0.3, shape = 16) +
  theme_bw() +
  geom_abline(slope = 1, 
              intercept =  log10(panel_a_data$median_expression[1]),
              size = 1,
              color = "darkblue") +
  coord_equal() +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Gene DNA abundance (CPM)") +
  ylab("Gene RNA abundance (CPM)")
ggsave("b_longum_np46_expression.pdf", w=4, h=4, useDingbats = F)

normalized_expression2 <- 
  normalized_expression %>%
  select(gene_family, species, mother_baby, family_id, norm_expression) %>%
  pivot_wider(names_from = mother_baby, values_from = norm_expression) %>%
  mutate(difference_baby_mom = Baby / Mother) %>%
  filter(!(is.na(difference_baby_mom))) %>%
  mutate(upregulated =  ifelse(Baby > 2 & difference_baby_mom > 2, T, F),
         downregulated = ifelse(Baby < 0.5 & difference_baby_mom < 0.5, T, F),
         category = case_when(upregulated ~ "Upregulated",
                              downregulated ~ "Downregulated",
                              TRUE ~ "No change"),
         species = str_remove(species, "s__"),
         species = str_replace(species, "_", " "))


dim(normalized_expression2)
# 68850 genes

normalized_expression2 %>%
  group_by(category) %>%
  summarise(n = n())
12574 / 68850
14844 / 68850

ggplot(normalized_expression2, aes(y=Baby, x = difference_baby_mom, color = category)) +
  geom_point(size = 1) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() + 
  coord_equal() +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  scale_color_manual(values = brewer.pal(9, "Set1")[c(2,9,1)]) +
  xlab("Expression; Baby / mother") +
  ylab("Relative expression in Baby")
ggsave("baby_vs_mom_expression.pdf", h = 4, w = 6, useDingbats = F)

n_strains_per_species <- 
  normalized_expression2 %>%
  select(species, family_id) %>%
  distinct() %>%
  group_by(species) %>%
  summarise(n_strains = n())

normalized_expression2 <-
  normalized_expression2 %>%
  left_join(n_strains_per_species) %>%
  mutate(species_label = str_c(species, " (", n_strains,")"))

labeling <- 
  normalized_expression2 %>%
  select(gene_family, species_label, category) %>%
  distinct() %>%
  group_by(species_label, category) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = category, values_from = n) %>%
  mutate(Upregulated = str_c("N = ", Upregulated),
         Downregulated = str_c("N = ", Downregulated))

# alternative labeling with the average number of over/underexpressed genes per family
labeling_alt <- 
  normalized_expression2 %>%
  distinct(species, family_id) %>%
  group_by(species) %>%
  summarise(n_strains = n()) %>%
  full_join(normalized_expression2) %>%
  group_by(species, n_strains) %>%
  summarise(n_over = str_c("N = ", as.integer(sum(upregulated) / max(n_strains))),
            n_under = str_c("N = ", as.integer(sum(downregulated) / max(n_strains))))

ggplot(normalized_expression2, aes(y=Baby, x = difference_baby_mom, color=category)) +
  geom_point(size = 0.5) +
  scale_y_log10() +
  scale_x_log10(labels = c("a","0.01","1","100","e")) +
  theme_bw() + 
  facet_wrap(~species_label, nrow = 3) +
  geom_text(data = labeling, check_overlap = TRUE,
             x = 2, y = 3.2 , color = brewer.pal(9, "Set1")[1], size = 2,
             aes(label = Upregulated)) +
  geom_text(data = labeling, check_overlap = TRUE,
            x = -2, y = -2.55 , color = brewer.pal(9, "Set1")[2], size = 2,
            aes(label = Downregulated)) +
  scale_color_manual(values = brewer.pal(9, "Set1")[c(2,9,1)]) + 
  theme(legend.position = "none",
        strip.text.x = element_text(size = 6)) +
  coord_equal() +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("Relative expression; Baby / mother") +
  ylab("Relative expression in Baby") 
ggsave("baby_vs_mom_expression_all_transmitted_strains.pdf", h = 8, w = 8, useDingbats = F)

upregulated_genes <- 
  normalized_expression2 %>%
  filter(upregulated) %>%
  left_join(gene_names)
# 12,574 genes upregulated (relative expression > 2, baby mom-ratio > 2)

upregulated_genes %>%
  select(gene_name,
         gene_family,
         species,
         family_id,
         Mother,
         Baby,
         difference_baby_mom) %>%
  arrange(species) %>%
  rename(relative_expression_infant = Baby,
         relative_expression_mother = Mother,
         infant_mother_expression_ratio = difference_baby_mom) %>%
  write.xlsx(file = "Table_S2_1.xlsx")
  
upregulated_genes %>%
  group_by(species, family_id) %>%
  summarise(num_upregulated_genes = n()) %>%
  ungroup() %>%
  group_by(species) %>%
  summarise(avg_genes = sum(num_upregulated_genes) / n()) %>%
  arrange(-avg_genes)
# Bacteroides species have most genes that change expression

normalized_expression2 %>%
  left_join(gene_names) %>%
  filter(gene_name == "NO_NAME")
37156 / 68850

filter(upregulated_genes, gene_name == "NO_NAME")
# 5,764 with no annotation
5764 / 12564
filter(upregulated_genes, grepl("[T|t]ransport", gene_name))
# 358 transport related
filter(upregulated_genes, grepl("[T|t]ransferase", gene_name)) 
# 673 Transferases
filter(upregulated_genes, grepl("[H|h]ydrolase", gene_name)) 
# 138 Hydrolases
filter(upregulated_genes, grepl("[P|p]hage", gene_name)) 
# 17 phage proteins
filter(upregulated_genes, grepl("[S|s]us[A-Z]", gene_name))
# 83 SusC/D related genes
filter(upregulated_genes, grepl("[T|t]on[B|b]", gene_name))
# 71 TonB related genes
filter(upregulated_genes, grepl("Two_component_system", gene_name))
# 17 Two component system (quorum sensing) related genes


downregulated_genes <- 
  normalized_expression2 %>%
  filter(downregulated) %>%
  left_join(gene_names)
# 14,844 genes downregulated (relative expression < 0.5, baby mom-ratio < 0.5)

downregulated_genes %>%
  select(gene_name,
         gene_family,
         species,
         family_id,
         Mother,
         Baby,
         difference_baby_mom) %>%
  arrange(species) %>%
  rename(relative_expression_infant = Baby,
         relative_expression_mother = Mother,
         infant_mother_expression_ratio = difference_baby_mom) %>%
  write.xlsx(file = "Table_S2_2.xlsx")

filter(downregulated_genes, gene_name == "NO_NAME")
# 9,935 with no annotation
9935 / 14844
filter(downregulated_genes, grepl("[T|t]ransport", gene_name))
# 395 transport related
filter(downregulated_genes, grepl("[T|t]ransferase", gene_name)) 
# 334 Transferases
filter(downregulated_genes, grepl("[H|h]ydrolase", gene_name)) 
# 257 Hydrolases
filter(downregulated_genes, grepl("[P|p]hage", gene_name)) 
# 68 phage proteins
filter(downregulated_genes, grepl("[S|s]us[A-Z]", gene_name))
# 317 SusC/D related genes
filter(downregulated_genes, grepl("[T|t]on[B|b]", gene_name))
# 179 TonB related genes
filter(downregulated_genes, grepl("Two_component_system", gene_name))
# 71 Two component system (quorum sensing) related genes

# Similar:
# Transport: 358 vs. 395
# Transferases: 673 vs. 334
# Hydrolases: 138 vs. 257
# Phage proteins: 17 vs. 68
# SusC/D: 83 vs. 317
# TonB: 71 vs. 179
# Quorum sensing: 17 vs. 71

gene_replication <- upregulated_genes %>%
  group_by(gene_family, species) %>%
  summarise(n_replicates = n()) %>%
  ungroup() 

# Check for replication across families
gene_replication %>%
  group_by(species) %>%
  summarise(n_genes = n()) %>%
  right_join(gene_replication) %>%
  group_by(species, n_replicates, n_genes) %>%
  summarise(n_genes_per_replication = n()) %>%
  filter(n_replicates == 1) %>%
  mutate(n_replicated = n_genes - n_genes_per_replication) %>%
  mutate(perc_replicated = n_replicated / n_genes) %>%
  select(-n_genes_per_replication)