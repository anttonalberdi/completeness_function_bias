library(tidyverse)
library(readxl)
library(janitor)

df <- read_xlsx("data/mag_data/REDUX/selection.xlsx") %>%
  mutate(Path = word(Path, sep = "/", -1))

checkm <- df %>%
  select(Path, Completeness, Redundancy) %>%
  rename("Bin Id" = Path, "Contamination" = Redundancy)

write.csv(checkm, "data/mag_data/REDUX/checkm.csv")


gtdb <- df %>%
  select(Path, Taxonomy) %>%
  rename("user_genome" = Path, "classification" = Taxonomy)

write.csv(gtdb, "data/mag_data/REDUX/gtdb.csv")
