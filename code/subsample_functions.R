##
##
##

library(tidyverse)

#Import data
test <- read.delim("../data/DRAM_annotate/Apodemus_annotations.tsv.gz")
genome_stats <- read.delim("../data/MAG_info/final_bins_Info_Apodemus.csv", sep = ',')
dram_genome_stats <- read.delim("../data/DRAM/DRAM_genome_stats_Apodemus.tsv")
taxonomy <- read.delim("../data/GTDB_tk/gtdbtk.bac120.summary_Apodemus.tsv")

#Join data, figure out which MAGs to subsample
derep_mags <- inner_join(test, genome_stats, by = c("fasta" = "genome")) %>%
  inner_join(., taxonomy, by = c("fasta" = "user_genome")) %>%
  inner_join(., dram_genome_stats, by = c("fasta" = "genome")) %>%
  select(fasta, completeness, contamination, classification, number.of.scaffolds) %>%
  distinct()


#Subsample genes function:

test1 <- test %>%
  filter(fasta == "bin_m10.mtb70")

test2 <- test %>%
  filter(fasta == "bin_m12.mtb81")

#
set.seed(1337)
proportions <- c("0.95", "0.9", "0.85", "0.8", "0.75", "0.7", "0.65", "0.6",
                 "0.55", "0.5")

sample_test <- slice_sample(test1, prop = 0.95)

test <- replicate(100, slice_sample(test1, prop = 0.1))


#Subsample contigs function
test1 %>%
  group_by(scaffold) %>%
  summarize(n = n())

test2 %>%
  group_by(scaffold) %>%
  summarise(n = n()) %>%
  summarise(range = range(n), mean = mean(n))

test1 %>%
  group_by(scaffold) %>%
  summarise(scaffold_len = max(end_position)) %>%
  summarise(total_len = sum(scaffold_len))

test2 %>%
  group_by(scaffold) %>%
  summarise(scaffold_len = max(end_position)) %>%
  summarise(total_len = sum(scaffold_len), range = range(scaffold_len),
            mean = mean(scaffold_len))
