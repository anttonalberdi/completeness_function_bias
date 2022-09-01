library(tidyverse)
library(ggplot2)

## Load datasets
## *************

dram_product <- read.csv("data/dram_product_fullness_data.csv",header = TRUE)
rownames(dram_product) <- dram_product$X
dram_product <- dram_product[,-1]
Explanatory_dataset <- read.csv("data/explanatory_dataset.csv",header = TRUE)
models_list <- readRDS(file="data/models_list.rds")
Results_table <- read.csv("data/Results_table.csv")
module_hierarchy <- read.csv("data/module_hierarchy.csv")
n_steps <- read.csv("data/number_of_steps_in_modules.csv",header = TRUE)

## Exploratory analyses
## ********************
domain_palette=c("#DB93A9","#E89D74","#69AD86","#E9CA80","#D7BFAF","#99796D","#2D758C","#A3DBDD","#ABA5D2","#8E8E8E")
phylum_palette=c("#C66154","#B6DCDD","#F1DD7C","#98CC6B")

# Multiply 0-1 bound data by 100 to turn into %
Results_table <- data.frame(Results_table,Domain=module_hierarchy$Domain,
                            Category=module_hierarchy$Category,Steps=n_steps$Steps)
Results_table$Actinobacteriota_pred_70 <- Results_table$Actinobacteriota_pred_70*100
Results_table$Actinobacteriota_pred_80 <- Results_table$Actinobacteriota_pred_80*100
Results_table$Actinobacteriota_pred_90 <- Results_table$Actinobacteriota_pred_90*100
Results_table$Actinobacteriota_pred_100 <- Results_table$Actinobacteriota_pred_100*100
Results_table$Bacteroidota_pred_70 <- Results_table$Bacteroidota_pred_70*100
Results_table$Bacteroidota_pred_80 <- Results_table$Bacteroidota_pred_80*100
Results_table$Bacteroidota_pred_90 <- Results_table$Bacteroidota_pred_90*100
Results_table$Bacteroidota_pred_100 <- Results_table$Bacteroidota_pred_100*100
Results_table$Firmicutes_pred_70 <- Results_table$Firmicutes_pred_70*100
Results_table$Firmicutes_pred_80 <- Results_table$Firmicutes_pred_80*100
Results_table$Firmicutes_pred_90 <- Results_table$Firmicutes_pred_90*100
Results_table$Firmicutes_pred_100 <- Results_table$Firmicutes_pred_100*100
Results_table$Proteobacteria_pred_70 <- Results_table$Proteobacteria_pred_70*100
Results_table$Proteobacteria_pred_80 <- Results_table$Proteobacteria_pred_80*100
Results_table$Proteobacteria_pred_90 <- Results_table$Proteobacteria_pred_90*100
Results_table$Proteobacteria_pred_100 <- Results_table$Proteobacteria_pred_100*100

# Percentage of positive slope parameter estimates for Fullness-Completeness relationship
mean(Results_table$Actinobacteriota_slope>0)  # 92%
mean(Results_table$Bacteroidota_slope>0)      # 95%
mean(Results_table$Firmicutes_slope>0)        # 94%
mean(Results_table$Proteobacteria_slope>0)    # 93%

# Computation of the predicted differences in fullness of each module between completeness
# values of 70%, 80%, 90% and 100%, separated by Phylum.

Diff_80_70_Actino <- Results_table[,7]-Results_table[,6]
Diff_90_80_Actino <- Results_table[,8]-Results_table[,7]
Diff_100_90_Actino <- Results_table[,9]-Results_table[,8]
Diff_80_70_Bactero <- Results_table[,11]-Results_table[,10]
Diff_90_80_Bactero <- Results_table[,12]-Results_table[,11]
Diff_100_90_Bactero <- Results_table[,13]-Results_table[,12]
Diff_80_70_Firmi <- Results_table[,15]-Results_table[,14]
Diff_90_80_Firmi <- Results_table[,16]-Results_table[,15]
Diff_100_90_Firmi <- Results_table[,17]-Results_table[,16]
Diff_80_70_Proteo <- Results_table[,19]-Results_table[,18]
Diff_90_80_Proteo <- Results_table[,20]-Results_table[,19]
Diff_100_90_Proteo <- Results_table[,21]-Results_table[,20]

data_10_ranges_Actinobacteriota <- data.frame(Diff_80_70_Actino,Diff_90_80_Actino,Diff_100_90_Actino) %>%
  pivot_longer(everything())
data_10_ranges_Actinobacteriota$Completeness_range <- rep(c("70-80","80-90","90-100"),nrow(data_10_ranges_Actinobacteriota)/3)
data_10_ranges_Actinobacteriota$Module <- rep(Results_table$Function,each=3)
data_10_ranges_Bacteroidota <- data.frame(Diff_80_70_Bactero,Diff_90_80_Bactero,Diff_100_90_Bactero) %>%
  pivot_longer(everything())
data_10_ranges_Bacteroidota$Completeness_range <- rep(c("70-80","80-90","90-100"),nrow(data_10_ranges_Bacteroidota)/3)
data_10_ranges_Bacteroidota$Module <- rep(Results_table$Function,each=3)
data_10_ranges_Firmicutes <- data.frame(Diff_80_70_Firmi,Diff_90_80_Firmi,Diff_100_90_Firmi) %>%
  pivot_longer(everything())
data_10_ranges_Firmicutes$Completeness_range <- rep(c("70-80","80-90","90-100"),nrow(data_10_ranges_Firmicutes)/3)
data_10_ranges_Firmicutes$Module <- rep(Results_table$Function,each=3)
data_10_ranges_Proteobacteria <- data.frame(Diff_80_70_Proteo,Diff_90_80_Proteo,Diff_100_90_Proteo) %>%
  pivot_longer(everything())
data_10_ranges_Proteobacteria$Completeness_range=rep(c("70-80","80-90","90-100"),nrow(data_10_ranges_Proteobacteria)/3)
data_10_ranges_Proteobacteria$Module=rep(Results_table$Function,each=3)

## Figure 2A in main manuscript

data_10_ranges_Actinobacteriota$Phylum <- "Actinobacteriota"
data_10_ranges_Bacteroidota$Phylum <- "Bacteroidota"
data_10_ranges_Firmicutes$Phylum <- "Firmicutes"
data_10_ranges_Proteobacteria$Phylum <- "Proteobacteria"
data_10_ranges <- data.frame(rbind(data_10_ranges_Actinobacteriota,data_10_ranges_Bacteroidota,
                                data_10_ranges_Firmicutes,data_10_ranges_Proteobacteria))

ggplot() +
  geom_point(aes(x=Completeness_range,y=value, color=Phylum),
             position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.7),
             alpha=0.3,
             data = data_10_ranges)+
  scale_color_manual(values = phylum_palette)+
  geom_boxplot(aes(x=Completeness_range,y=value),
               data = data_10_ranges,
               color="grey50",
               fill=NA,
               size=1)+
  theme_classic()+
  ylab("Change in Fullness (%)")+
  xlab("Completeness range (%)")+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        title = element_text(size=18))

ggsave("code/Fig_2A.pdf",width = 8,height = 8)

# Diff Fullness between Completeness 70%-100% by Phylum
Diff_100_70_Actino <- Results_table[,9]-Results_table[,6]
Diff_100_70_Bactero <- Results_table[,13]-Results_table[,10]
Diff_100_70_Firmi <- Results_table[,17]-Results_table[,14]
Diff_100_70_Proteo <- Results_table[,21]-Results_table[,18]

# 16% average increase in Fullness with 30% increase in Completeness
mean(rowMeans(data.frame(Diff_100_70_Actino,Diff_100_70_Bactero,
                         Diff_100_70_Firmi,Diff_100_70_Proteo)))
sd(rowMeans(data.frame(Diff_100_70_Actino,Diff_100_70_Bactero,
                       Diff_100_70_Firmi,Diff_100_70_Proteo)))
