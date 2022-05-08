#
#
#

#Load packages
library(tidyverse)
library(ggforce)
library(ggdist)
library(ggtext)

#Import data
apodemus <- read_delim("data/MAG_info/final_bins_Info_Apodemus.csv") %>%
  filter(completeness > 70 & contamination < 10) %>%
  mutate(source = "apodemus")
crocidura <- read_delim("data/MAG_info/final_bins_Info_Crocidura.csv") %>%
  filter(completeness > 70 & contamination < 10) %>%
  mutate(source = "crocidura")
felis <- read_delim("data/MAG_info/final_bins_Info_Felis.csv") %>%
  filter(completeness > 70 & contamination < 10)  %>%
  mutate(source = "felis")
gallus <- read_delim("data/MAG_info/final_bins_Info_Gallus.csv") %>%
  filter(completeness > 70 & contamination < 10)  %>%
  mutate(source = "gallus")
mus <- read_delim("data/MAG_info/final_bins_Info_Mus.csv") %>%
  filter(completeness > 70 & contamination < 10)  %>%
  mutate(source = "mus")

df <- rbind(apodemus, crocidura, felis, gallus, mus)

#Plot
ggplot(df, aes(x = source, y = completeness, colour = source)) +
  stat_halfeye(
    adjust = .5,
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA
  ) + 
  geom_point(
    aes(colour = source),
    size = 2,
    alpha = 0.4,
    position = position_jitter(seed = 1, width = .1)
  ) +
  theme_minimal() +
  theme(
    legend.position = "",
    panel.grid.major.x = element_blank(),
  )
  labs(y = "CheckM completeness", x = "")
