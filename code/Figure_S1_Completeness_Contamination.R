#
#
#

##########################################################################################
## Load packages
library(tidyverse)
library(ggforce)
library(ggdist)
library(ggtext)
library(patchwork)

##########################################################################################
## Import data
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

##########################################################################################
## Plot em

#Setup theme
theme_RE <- theme(
  legend.position = "",
  panel.grid.major.x = element_blank(),
  axis.text = element_text(face = "bold", size = 14, family = "Aharoni"),
  axis.text.x = element_text(angle = 45),
  axis.title = element_text(face = "bold", size = 16, family = "Aharoni")
)


#Fig XA
figXA <- ggplot(df, aes(x = source, y = completeness, colour = source)) +
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
  theme_RE +
  scale_y_continuous(breaks = seq(70, 100, by = 5)) +
  labs(y = "CheckM completeness", x = "")


#Fig XB
figXB <- ggplot(df, aes(x = source, y = contamination, colour = source)) +
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
  theme_RE +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  labs(y = "CheckM contamination", x = "")


#Stich together
figXA | figXB 

#Save
ggsave("FigureX.png", width = 10, height = 8)

