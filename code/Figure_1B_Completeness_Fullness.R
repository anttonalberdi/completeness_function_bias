### Completeness by (average of fullness | each fullness) -- colour by species
### Raphael Eisenhofer May 2022

##########################################################################################
## Load packages
library(tidyverse)
library(ggforce)
library(ggdist)
library(ggtext)
library(patchwork)

##########################################################################################
## Import MAG info
apodemus <- read_delim("data/mag_data/MAG_info/final_bins_Info_Apodemus.csv") %>%
  filter(completeness > 70 & contamination < 10) %>%
  mutate(source = "Apodemus")
crocidura <- read_delim("data/mag_data/MAG_info/final_bins_Info_Crocidura.csv") %>%
  filter(completeness > 70 & contamination < 10) %>%
  mutate(source = "Crocidura")
felis <- read_delim("data/mag_data/MAG_info/final_bins_Info_Felis.csv") %>%
  filter(completeness > 70 & contamination < 10)  %>%
  mutate(source = "Felis")
gallus <- read_delim("data/mag_data/MAG_info/final_bins_Info_Gallus.csv") %>%
  filter(completeness > 70 & contamination < 10)  %>%
  mutate(source = "Gallus")
mus <- read_delim("data/mag_data/MAG_info/final_bins_Info_Mus.csv") %>%
  filter(completeness > 70 & contamination < 10)  %>%
  mutate(source = "Mus")

#Combine into one dataframe
df_genomeinfo <- rbind(apodemus, crocidura, felis, gallus, mus)

#Import DRAM product files
apodemus_dram <- read_delim("data/mag_data/DRAM/Hacked/product_Apodemus.tsv.gz")
crocidura_dram <- read_delim("data/mag_data/DRAM/Hacked/product_Crocidura.tsv.gz")
felis_dram <- read_delim("data/mag_data/DRAM/Hacked/product_Felis.tsv.gz")
gallus_dram <- read_delim("data/mag_data/DRAM/Hacked/product_Gallus.tsv.gz")
mus_dram <- read_delim("data/mag_data/DRAM/Hacked/product_Mus.tsv.gz")

#Combine into one DF
df_dram <- rbind(apodemus_dram, crocidura_dram, felis_dram, gallus_dram, mus_dram)

#Merge genome info and DRAM info
df_joined <- left_join(df_genomeinfo, df_dram, by = "genome") %>%
  mutate(mean_module_fullness = rowMeans(select(., -genome, -completeness, -contamination, -source) * 100))

##########################################################################################
#Define colour palette
colours <- c("Apodemus" = "#DACC69", "Crocidura" = "#E5685C", "Felis" = "#F08D49",
             "Gallus" = "#7BCCC2", "Mus" = "#869AB2")

##########################################################################################
##Plot
df_joined %>%
  ggplot(aes(x = completeness, y = mean_module_fullness, colour = source)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  scale_colour_manual(values = colours) +
  theme_classic() + 
  theme(
        legend.position = c(0.35, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(size = 0.5),
        panel.grid.major.x = element_blank(),
#        text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(face = "bold", size = 14),
#        axis.title = element_text(face = "bold", size = 16)
       ) +
  labs(x = "Completeness (%)", y = "Mean module fullness (%)")

ggsave("figures/Figure1B.pdf", width = 10, height = 6)
