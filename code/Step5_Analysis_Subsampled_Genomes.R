library(lme4)
library(vegan)
library(ggplot2)
library(lemon)

# Load datasets
# *************
subsampled_table <- read.delim("data/simulated_zymo_all_distillate_product.tsv")
rownames(subsampled_table)<-subsampled_table$genome
subsampled_table <- subsampled_table[,grepl("M0",colnames(subsampled_table))]
head(subsampled_table)
dim(subsampled_table)

stats_data<-read.delim("data/checkM_stats_zymo_subsampled.tsv")
head(stats_data)
dim(stats_data)
stats_data$source[1]<-"Bacillus_subtilis"
stats_data$source[2]<-"Enterococcus_faecalis"
stats_data$source[3]<-"Escherichia_coli"
stats_data$source[4]<-"Lactobacillus_fermentum"
stats_data$source[5]<-"Listeria_monocytogenes"
stats_data$source[6]<-"Pseudomonas_aeruginosa"
stats_data$source[335]<-"Salmonella_enterica"
stats_data$source[336]<-"Staphylococcus_aureus"
stats_data$phylum<-NA
stats_data$phylum[grepl("Bacillus",stats_data$source)]<-"Firmicutes"
stats_data$phylum[grepl("Enterococcus",stats_data$source)]<-"Firmicutes"
stats_data$phylum[grepl("Escherichia",stats_data$source)]<-"Proteobacteria"
stats_data$phylum[grepl("Lactobacillus",stats_data$source)]<-"Firmicutes"
stats_data$phylum[grepl("Listeria",stats_data$source)]<-"Firmicutes"
stats_data$phylum[grepl("Pseudomonas",stats_data$source)]<-"Proteobacteria"
stats_data$phylum[grepl("Salmonella",stats_data$source)]<-"Proteobacteria"
stats_data$phylum[grepl("Staphylococcus",stats_data$source)]<-"Firmicutes"

subsampled_table<-subsampled_table[stats_data$subsample_rate!="60"|is.na(stats_data$subsample_rate),]
stats_data<-stats_data[stats_data$subsample_rate!="60"|is.na(stats_data$subsample_rate),]
stats_data$subsample_rate[is.na(stats_data$subsample_rate)]<-100
stats_data$subsample_rate<-factor(stats_data$subsample_rate)

mean(rownames(subsampled_table)==stats_data$bin_id)
dim(subsampled_table)

dram_product <- read.csv("data/dram_product_fullness_data.csv",header = TRUE)
rownames(dram_product) <- dram_product$X
dram_product <- dram_product[,-1]

subsampled_table <- subsampled_table[,colnames(subsampled_table)%in%colnames(dram_product)]

dim(subsampled_table)

## Ordination of subsampled data
set.seed(1)
subsampled_data_pcoa<-cmdscale(dist(subsampled_table))
subsampled_data_pcoa_sc<-scores(subsampled_data_pcoa,display = "sites")
subsampled_data_pcoa_sc_df<-data.frame(subsampled_data_pcoa_sc,
                                       subsample_rate=stats_data$subsample_rate,
                                       genome=stats_data$source,
                                       phylum=stats_data$phylum)
ggplot(subsampled_data_pcoa_sc_df)+
  geom_point(aes(x=Dim1,Dim2,col=genome,shape=subsample_rate), size=c(2,2,2,4)[stats_data$subsample_rate])+
  scale_shape_manual(values = c(1,2,5,15))+
  theme_bw()

## Apply full model to correct subsampled genomes
models_list <- readRDS(file="data/models_list.rds")

subsampled_table_corrected <- matrix(0,nrow = nrow(subsampled_table),ncol = ncol(subsampled_table))
colnames(subsampled_table_corrected)<-colnames(subsampled_table)
rownames(subsampled_table_corrected)<-rownames(subsampled_table)

n_steps <- read.csv("data/number_of_steps_in_modules.csv",header = TRUE)
n_steps<-n_steps[n_steps$ModuleID%in%colnames(subsampled_table),]


for(i in 1:ncol(subsampled_table)){
  Model <- models_list[[i]]
  for(j in 1:nrow(subsampled_table)){
    # Model prediction of fullness if completeness was 100%
    pred_100 <- round(predict(Model,newdata = data.frame(Phylum=stats_data$phylum[j],Completeness=100),
                                type = "response"),1)
    # Model prediction of fullness for actual completeness of the focal genome
    pred_focal <- round(predict(Model,newdata = data.frame(Phylum=stats_data$phylum[j],Completeness=stats_data$completeness[j]),
                                type = "response"),1)
    # The expected change in function fullness if focal MAG was 100% complete
    pred_diff <- pred_100-pred_focal
    subsampled_table_corrected[j,i]=subsampled_table[j,i]+pred_diff
  }
}
# If corrected fullness >1, convert it to 1. If corrected fullness <0, convert it to 0.
subsampled_table_corrected[subsampled_table_corrected>1] <- 1
subsampled_table_corrected[subsampled_table_corrected<0] <- 0

subsampled_table_corrected[stats_data$subsample_rate=="100",]<-as.matrix(subsampled_table[stats_data$subsample_rate=="100",])


## Ordination of whole data (Full genomes and raw and corrected sub sampled genomes together)

subsampled_table_whole<-data.frame(rbind(subsampled_table,subsampled_table_corrected[!stats_data$subsample_rate=="100",]))
stats_data_corrected<-stats_data
stats_data$status<-NA
stats_data_corrected$status<-NA
stats_data$status[stats_data$subsample_rate=="100"]<-"full"
stats_data$status[!stats_data$subsample_rate=="100"]<-"raw"
stats_data_corrected$status[stats_data_corrected$subsample_rate=="100"]<-"full"
stats_data_corrected$status[!stats_data_corrected$subsample_rate=="100"]<-"corrected"

stats_data_whole<-data.frame(rbind(stats_data,stats_data_corrected[!stats_data_corrected$subsample_rate=="100",]))
set.seed(1)
subsampled_data_whole_pcoa<-cmdscale(dist(subsampled_table_whole))
subsampled_data_whole_pcoa_sc<-scores(subsampled_data_whole_pcoa,display = "sites")
subsampled_data_whole_pcoa_sc_df<-data.frame(subsampled_data_whole_pcoa_sc,
                                             subsample_rate=stats_data_whole$subsample_rate,
                                             genome=stats_data_whole$source,
                                             phylum=stats_data_whole$phylum,
                                             status=stats_data_whole$status)
ggplot(subsampled_data_whole_pcoa_sc_df)+
  geom_point(aes(x=Dim1,Dim2,col=genome,shape=subsample_rate), size=c(2,2,2,4)[stats_data_whole$subsample_rate])+
  scale_shape_manual(values = c(1,2,5,15))+
  theme_bw()+
  scale_color_manual(values = c('#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))

corrected_90<-ggplot()+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="90"&stats_data_whole$status=="corrected",],
             aes(x=Dim1,Dim2,col=genome),
             size=2,shape=1)+
  ggtitle("Corrected 90%")+
  stat_ellipse(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="90"&stats_data_whole$status=="corrected",],
               aes(x=Dim1,Dim2,col=genome))+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="100",],
             aes(x=Dim1,Dim2,col=genome),
             size=4,shape=16)+
  theme_bw()+
  coord_cartesian(xlim=c(-4.5,6),ylim = c(-2.5,3))+
  scale_color_manual(values = c('#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))


raw_90<-ggplot()+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="90"&stats_data_whole$status=="raw",],
             aes(x=Dim1,Dim2,col=genome),
             size=2,shape=1)+
  ggtitle("Raw 90%")+
  stat_ellipse(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="90"&stats_data_whole$status=="raw",],
               aes(x=Dim1,Dim2,col=genome))+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="100",],
             aes(x=Dim1,Dim2,col=genome),
             size=4,shape=16)+
  theme_bw()+
  coord_cartesian(xlim=c(-4.5,6),ylim = c(-2.5,3))+
  scale_color_manual(values = c('#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))


corrected_80<-ggplot()+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="80"&stats_data_whole$status=="corrected",],
             aes(x=Dim1,Dim2,col=genome),
             size=2,shape=1)+
  ggtitle("Corrected 80%")+
  stat_ellipse(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="80"&stats_data_whole$status=="corrected",],
               aes(x=Dim1,Dim2,col=genome))+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="100",],
             aes(x=Dim1,Dim2,col=genome),
             size=4,shape=16)+
  theme_bw()+
  coord_cartesian(xlim=c(-4.5,6),ylim = c(-2.5,3))+
  scale_color_manual(values = c('#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))


raw_80<-ggplot()+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="80"&stats_data_whole$status=="raw",],
             aes(x=Dim1,Dim2,col=genome),
             size=2,shape=1)+
  ggtitle("Raw 80%")+
  stat_ellipse(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="80"&stats_data_whole$status=="raw",],
               aes(x=Dim1,Dim2,col=genome))+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="100",],
             aes(x=Dim1,Dim2,col=genome),
             size=4,shape=16)+
  theme_bw()+
  coord_cartesian(xlim=c(-4.5,6),ylim = c(-2.5,3))+
  scale_color_manual(values = c('#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))


corrected_70<-ggplot()+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="70"&stats_data_whole$status=="corrected",],
             aes(x=Dim1,Dim2,col=genome),
             size=2,shape=1)+
  ggtitle("Corrected 70%")+
  stat_ellipse(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="70"&stats_data_whole$status=="corrected",],
               aes(x=Dim1,Dim2,col=genome))+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="100",],
             aes(x=Dim1,Dim2,col=genome),
             size=4,shape=16)+
  theme_bw()+
  coord_cartesian(xlim=c(-4.5,6),ylim = c(-2.5,3))+
  scale_color_manual(values = c('#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))


raw_70<-ggplot()+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="70"&stats_data_whole$status=="raw",],
             aes(x=Dim1,Dim2,col=genome),
             size=2,shape=1)+
  ggtitle("Raw 70%")+
  stat_ellipse(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="70"&stats_data_whole$status=="raw",],
               aes(x=Dim1,Dim2,col=genome))+
  geom_point(data=subsampled_data_whole_pcoa_sc_df[stats_data_whole$subsample_rate=="100",],
             aes(x=Dim1,Dim2,col=genome),
             size=4,shape=16)+
  theme_bw()+
  coord_cartesian(xlim=c(-4.5,6),ylim = c(-2.5,3))+
  scale_color_manual(values = c('#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))


grid_arrange_shared_legend(raw_90,corrected_90,
                           raw_80,corrected_80,
                           raw_70,corrected_70,
                           ncol = 2,nrow = 3)


