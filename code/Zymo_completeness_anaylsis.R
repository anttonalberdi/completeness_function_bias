# Raphael Eisenhofer, 2023

library(tidyverse)
library(janitor)

cdb <- read_delim("data/zymo_MAG_data/Cdb.csv")
genomeinfo <- read_delim("data/zymo_MAG_data/genomeInfo.csv")
ncontigs <- read_delim("data/zymo_MAG_data/n_contigs.tsv", col_names = c("genome", "ncontigs"))

combined <- cdb %>%
  inner_join(., genomeinfo, by = "genome") %>%
  left_join(., ncontigs, by = "genome") %>%
  select(!threshold & !cluster_method & !comparison_algorithm)

write_tsv(combined, "data/zymo_MAG_data/zymo_MAG_stats.tsv")

# I selected these MAGs for subsampling:
# "SRR12324251_bin.3.fa.gz" = "Lactobacillus_fermentum" (65)
# "SRR12324251-5Mss_bin.4.fa.gz" = "Bacillus_subtilis" (98)
# "SRR12324251-5Mss_bin.1.fa.gz" = "Enterococcus_faecalis" (42)
# "SRR12324251-2Mss_bin.4.fa.gz" = "Staphylococcus_aureus" (215)
# "SRR12324251-5Mss_bin.6.fa.gz" = "Listeria_monocytogenes" (59)
# "SRR12324251-1Mss_bin.1.fa.gz" = "Pseudomonas_aeruginosa" (80)
# "SRR12324251-2Mss_bin.3.fa.gz" = "Salmonella_enterica" (72)
# "SRR12324251-2Mss_bin.5.fa.gz" = "Escherichia_coli" (89)


checkm <- read_delim("data/zymo_MAG_data/subsampled_MAGs_checkm.tsv") %>%
  clean_names() %>%
  select(bin_id, completeness, contamination) %>%
  mutate(source = case_when(bin_id = str_detect(bin_id, "SRR12324251_bin.3") ~ "Lactobacillus_fermentum",
                            bin_id = str_detect(bin_id, "SRR12324251-5Mss_bin.4") ~ "Bacillus_subtilis",
                            bin_id = str_detect(bin_id, "SRR12324251-5Mss_bin.1") ~ "Enterococcus_faecalis",
                            bin_id = str_detect(bin_id, "SRR12324251-2Mss_bin.4") ~ "Staphylococcus_aureus",
                            bin_id = str_detect(bin_id, "SRR12324251-5Mss_bin.6") ~ "Listeria_monocytogenes",
                            bin_id = str_detect(bin_id, "SRR12324251-1Mss_bin.1") ~ "Pseudomonas_aeruginosa",
                            bin_id = str_detect(bin_id, "SRR12324251-2Mss_bin.3") ~ "Salmonella_enterica",
                            bin_id = str_detect(bin_id, "SRR12324251-2Mss_bin.5") ~ "Escherichia_coli"
                            ),
         subsample_rate = case_when(bin_id = str_detect(bin_id, "RATE0.6") ~ 60,
                                    bin_id = str_detect(bin_id, "RATE0.7") ~ 70,
                                    bin_id = str_detect(bin_id, "RATE0.8") ~ 80,
                                    bin_id = str_detect(bin_id, "RATE0.9") ~ 90,
                                    )
         )

write_tsv(checkm, "data/zymo_MAG_data/checkM_stats_zymo_subsampled.tsv")

figure1 <- checkm %>%
  ggplot(aes(x = source, y = completeness, colour = source)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.7) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, size = 12, face = "italic"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 16)
  ) +
  labs(x = "Source genome", y = "CheckM completeness")
  
ggsave(x = figure1, filename = "figures/extras/simulated_completeness_ranges.png")

figure2 <- checkm %>%
  ggplot(aes(x = subsample_rate, y = completeness, colour = source)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 16)
  ) +
  labs(x = "Contig sampling rate (%)", y = "CheckM completeness (%)")

ggsave(x = figure1, filename = "figures/extras/contig_sampling_rate.png")

