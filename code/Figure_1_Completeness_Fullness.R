### Code for reproducing figure 1
### Raphael Eisenhofer Sept 2022

##########################################################################################
## Load packages
library(tidyverse)
library(ggforce)
library(ggdist)
library(ggtext)
library(patchwork)

##########################################################################################
## Import MAG info
genome_stats <- read_delim("data/genome_stats.csv") %>%
  mutate(phylum = word(Taxonomy, 2, sep = ";") %>%
           str_replace(., "p__", "") %>%
           str_replace(., "_.", "")) %>%
  mutate(Accession = str_replace(Accession, ".._", ""))


## Import DRAM product files
dram <- read_delim("data/mag_data/DRAM/product.tsv.gz")

## Merge genome info and DRAM info
df_joined <- left_join(genome_stats, dram, by = c("Accession" = "genome")) %>%
  mutate(mean_module_fullness = rowMeans(select(., -Accession, -Completeness, -phylum,
                                                -Redundancy, -Taxonomy) * 100))

##########################################################################################
## Define colour palette
colours <- c("Firmicutes" = "#DACC69", "Proteobacteria" = "#E5685C", 
             "Actinobacteriota" = "#F08D49", "Bacteroidota" = "#7BCCC2")

##########################################################################################
# Plot 1A

## Setup theme
theme_RE <- theme(
  legend.position = "",
  panel.grid.major.x = element_blank(),
  axis.text = element_text(face = "bold", size = 14),
  axis.text.x = element_blank(),
  axis.title = element_text(size = 14)
)

fig1a <- ggplot(df_joined, aes(x = phylum, y = Completeness, colour = phylum)) +
  stat_halfeye(
    adjust = .6,
    width = .5, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA
  ) + 
  geom_point(
    aes(colour = phylum),
    size = 2,
    alpha = 0.4,
    position = position_jitter(seed = 1, width = .1)
  ) +
  scale_colour_manual(values = colours) +
  theme_minimal() +
  theme_RE +
  scale_y_continuous(breaks = seq(70, 100, by = 5)) +
  labs(y = "CheckM completeness", x = "")

ggsave(plot = fig1a, filename = "figures/Fig_1A.pdf", height = 8, width = 10, unit = "in")



#Plot 1B
fig1b <- ggplot(df_joined, aes(x = phylum, y = Redundancy, colour = phylum)) +
  stat_halfeye(
    adjust = .6,
    width = .5, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA
  ) + 
  geom_point(
    aes(colour = phylum),
    size = 2,
    alpha = 0.4,
    position = position_jitter(seed = 1, width = .1)
  ) +
  scale_colour_manual(values = colours) +
  theme_minimal() +
  theme_RE +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  labs(y = "CheckM contamination", x = "")


ggsave(plot = fig1b, filename = "figures/Fig_1B.pdf", height = 8, width = 10, unit = "in")



# Plot 1C
fig1c <- df_joined %>%
  ggplot(aes(x = Completeness, y = mean_module_fullness, colour = phylum)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess", size = 2) +
  scale_colour_manual(values = colours) +
  theme_classic() + 
  theme(
    #        legend.position = c(0.1, 0.93),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.box.background = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(face = "bold", size = 14),
  ) +
  coord_cartesian(xlim = c(69.5, 101), ylim = c(0, 35), expand = FALSE) +
  labs(x = "Completeness (%)", y = "Mean module fullness (%)")

ggsave(plot = fig1c, filename = "figures/Fig_1C.pdf", height = 8, width = 10, unit = "in")




