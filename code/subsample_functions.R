##
##
##

library(tidyverse)

## Import data
test <- read.delim("data/DRAM_annotate/Apodemus_annotations.tsv.gz")
genome_stats <- read.delim("data/MAG_info/final_bins_Info_Apodemus.csv", sep = ',')
dram_genome_stats <- read.delim("data/DRAM/Hacked/genome_stats_Apodemus.tsv")
taxonomy <- read.delim("data/GTDB_tk/gtdbtk.bac120.summary_Apodemus.tsv")

## Join data, figure out which MAGs to subsample
derep_mags_100_0 <- inner_join(test, genome_stats, by = c("fasta" = "genome")) %>%
  inner_join(., taxonomy, by = c("fasta" = "user_genome")) %>%
  inner_join(., dram_genome_stats, by = c("fasta" = "genome")) %>%
  select(fasta, completeness, contamination, classification, number.of.scaffolds) %>%
  distinct() %>%
  filter(completeness == 100 & contamination == 0)

## Write function to automate this import and filtering of tables



# Subsample genes function:

test1 <- test %>%
  filter(fasta == "bin_m10.mtb70")

test2 <- test %>%
  filter(fasta == "bin_m12.mtb81")


### Set the seed for reproducibility
### Create vector containing the desired proportions
set.seed(1337)
props <- c(0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5)


#Sanity test slice_sample function
sample_test <- slice_sample(test1, prop = 0.95)

#Sanity test replicate function
test <- replicate(100, slice_sample(test1, prop = 0.1))

## Custom function to run through n proportion randomly 100 times


listo <- lapply(props, function(prop) {
  replicate(100, slice_sample(test1, prop = prop))
}
)



## CheckM completion estimates
#Load the list of CheckM marker sets, select universal Bacterial set (104 markers)
#Then split into 1 marker per row.
#Actually, I think it makes sense to split each 'set', as these are colocated???
bacterial_marker <- read_delim("data/MAG_info/taxon_marker_sets.tsv", delim = '\t',
                        col_names = c("A", "B", "C", "D", "E", "F", "G")) %>%
  filter(C == "Bacteria") %>%
  select(G) %>%
  str_split(., pattern = ",") %>%
  as.data.frame(., col.names = "marker")
  
#Need to clean the marker genes: "'", ")", etc..
bacterial_marker <- bacterial_marker %>%
  mutate(marker = str_replace_all(
    marker, 
    pattern = c("\\'|\\(|\\)|\\[|\\]|set| |\\..."), 
    replacement = "")
  )





# Subsample contigs function
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
