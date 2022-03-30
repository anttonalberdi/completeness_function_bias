########################################################################
### Association between MAG completeness and Functional capabilities ###
########################################################################
library(tidyverse)
library(sjPlot)

# 1. Load and prepare the datasets 
# ********************************

### 1.1. Apodemus dataset ###

# DRAM product

Apodemus_DRAM_product=read.delim("data/DRAM/Hacked/product_Apodemus.tsv",sep="\t")
rownames(Apodemus_DRAM_product)=Apodemus_DRAM_product$genome
Apodemus_DRAM_product=Apodemus_DRAM_product[,-1]
# Convert True/False character to 1/0 binary numeric data
for(i in 1:ncol(Apodemus_DRAM_product)){
  if(is.character(Apodemus_DRAM_product[,i])){
    Apodemus_DRAM_product[,i][Apodemus_DRAM_product[,i]=="True"]=1
    Apodemus_DRAM_product[,i][Apodemus_DRAM_product[,i]=="False"]=0
    Apodemus_DRAM_product[,i]=as.integer(Apodemus_DRAM_product[,i])
  }
}
# Apodemus_DRAM_product$CAZy_richness=rowSums(Apodemus_DRAM_product[,grepl("CAZy",colnames(Apodemus_DRAM_product))])
# Apodemus_DRAM_product$CAZy_richness=as.integer(Apodemus_DRAM_product$CAZy_richness)
str(Apodemus_DRAM_product)

n_pathways=read.csv("data/DRAM/Hacked/pathway_numbers.csv",sep=";",header = TRUE)

# Discard completely absent functions from further exploration
Apodemus_dram=Apodemus_DRAM_product[,colSums(Apodemus_DRAM_product)>9]
n_pathways=n_pathways[colSums(Apodemus_DRAM_product)>9,]

mean(colnames(Apodemus_dram)==n_pathways$Function)

## BIN quality
Apodemus_quality=read.delim("data/MAG_info/final_bins_Info_Apodemus.csv",sep=",")
Apodemus_quality=Apodemus_quality[Apodemus_quality$genome%in%rownames(Apodemus_dram),]
Apodemus_quality=Apodemus_quality[match(rownames(Apodemus_dram),Apodemus_quality$genome),]

## Taxonomy
Apodemus_taxonomy=read.delim("data/GTDB_tk/gtdbtk.bac120.summary_Apodemus.tsv",sep="\t")
mean(Apodemus_taxonomy$user_genome==rownames(Apodemus_DRAM_product))

Apodemus_taxonomyclean=Apodemus_taxonomy %>% 
  mutate(classification = str_replace_all(classification, ".__", "")) %>%
  # Split the classification string into columns for each taxonnomic rank
  separate(col = classification, sep = ";", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
# Remove first column
Apodemus_taxmatrix=Apodemus_taxonomyclean[,c(2:6)]
# Set MAG name as row names of the dataframe
rownames(Apodemus_taxmatrix)=Apodemus_taxonomyclean[,1]

# Pool Firmicutes A,B C into Firmicutes
# taxmatrix_mod=data.frame(taxmatrix)
Apodemus_taxmatrix$Phylum[grepl("Firmi",Apodemus_taxmatrix$Phylum)]="Firmicutes"

mean(rownames(Apodemus_taxmatrix)==rownames(Apodemus_dram))
mean(rownames(Apodemus_taxmatrix)==Apodemus_quality$genome)

Explanatory_dataset=data.frame(Phylum=Apodemus_taxmatrix$Phylum,
                            completeness=Apodemus_quality$completeness-100)

Apodemus_dram_phy=Apodemus_dram[Explanatory_dataset$Phylum=="Bacteroidota"|
                                  Explanatory_dataset$Phylum=="Firmicutes"|
                                  Explanatory_dataset$Phylum=="Proteobacteria",]
Explanatory_dataset_phy=Explanatory_dataset[Explanatory_dataset$Phylum=="Bacteroidota"|
                                    Explanatory_dataset$Phylum=="Firmicutes"|
                                    Explanatory_dataset$Phylum=="Proteobacteria",]
Apodemus_table=data.frame(matrix(ncol=12,nrow = ncol(Apodemus_dram_phy)))
colnames(Apodemus_table)=c("Dataset","Function","Bacteroidota_beta0","Firmicutes_beta0","Proteobacteria_beta0","Bacteroidota_beta1","Firmicutes_beta1","Proteobacteria_beta1","P_int","mean_beta0","mean_beta1","P_compl")
for(i in 1:ncol(Apodemus_dram_phy)){
    response_int=Apodemus_dram_phy[,i]
    response=Apodemus_dram[,i]
    weights_int=rep(n_pathways[i,2],nrow(Apodemus_dram_phy))
    weights=rep(n_pathways[i,2],nrow(Apodemus_dram))
    M_int=glm(response_int~Phylum*completeness,family = binomial,weights = weights_int,data = Explanatory_dataset_phy)
    M=glm(response~completeness,family = binomial,weights = weights,data = Explanatory_dataset)
    Apodemus_table[i,1]="Apodemus"
    Apodemus_table[i,2]=n_pathways[i,1]
    Apodemus_table[i,3]=exp(coef(M_int)[1])/(1+exp(coef(M_int)[1]))
    Apodemus_table[i,4]=exp(coef(M_int)[1]+coef(M_int)[2])/(1+exp(coef(M_int)[1]+coef(M_int)[2]))
    Apodemus_table[i,5]=exp(coef(M_int)[1]+coef(M_int)[3])/(1+exp(coef(M_int)[1]+coef(M_int)[3]))
    Apodemus_table[i,6]=coef(M_int)[4]
    Apodemus_table[i,7]=coef(M_int)[4]+coef(M_int)[5]
    Apodemus_table[i,8]=coef(M_int)[4]+coef(M_int)[6]
    Apodemus_table[i,9]=anova(M_int,test = "Chisq")[4,5]
    Apodemus_table[i,10]=exp(coef(M)[1])/(1+exp(coef(M)[1]))
    Apodemus_table[i,11]=coef(M)[2]
    Apodemus_table[i,12]=anova(M,test = "Chisq")[2,5]
    }

## The completeness variable was transformed so that 100 completeness = 0.
## With this transformation, the intercepts (beta0) can be interpreted as: expected log-odds
## if completeness was 100%. 
## In the table, I have transformed the log-odds into probabilities, ranging between 0-1.
## Then, the intercept is the average fullness of a function when completeness = 100%.
## I left the slopes (beta1) in the log-odds scale, hence the slope can be interpreted
## as change in log-odds with 1% change in completeness.
## log-odds change linearly with explanatory variables what allows interpreting the slope
## following the rules of the linear model. We cannot simply transform the slopes to the
## probability scale, since then the relationship with completeness is expected to be non-linear.
## IDEA: we could predict the expected fullness (in probability scale) of a function when 
## completeness = 70% and when completeness = 100%. The difference between them would result in 
## a variable indicating % increase (decrease) in fullness, when changing completeness from 70% to 100%.

head(Apodemus_table)

### 1.2. Crocidura dataset ###

# DRAM product

Crocidura_DRAM_product=read.delim("data/DRAM/Hacked/product_Crocidura.tsv",sep="\t")
rownames(Crocidura_DRAM_product)=Crocidura_DRAM_product$genome
Crocidura_DRAM_product=Crocidura_DRAM_product[,-1]
# Convert True/False character to 1/0 binary numeric data
for(i in 1:ncol(Crocidura_DRAM_product)){
  if(is.character(Crocidura_DRAM_product[,i])){
    Crocidura_DRAM_product[,i][Crocidura_DRAM_product[,i]=="True"]=1
    Crocidura_DRAM_product[,i][Crocidura_DRAM_product[,i]=="False"]=0
    Crocidura_DRAM_product[,i]=as.integer(Crocidura_DRAM_product[,i])
  }
}
# Crocidura_DRAM_product$CAZy_richness=rowSums(Crocidura_DRAM_product[,grepl("CAZy",colnames(Crocidura_DRAM_product))])
# Crocidura_DRAM_product$CAZy_richness=as.integer(Crocidura_DRAM_product$CAZy_richness)
str(Crocidura_DRAM_product)

n_pathways=read.csv("data/DRAM/Hacked/pathway_numbers.csv",sep=";",header = TRUE)

# Discard completely absent functions from further exploration
Crocidura_dram=Crocidura_DRAM_product[,colSums(Crocidura_DRAM_product)>9]
n_pathways=n_pathways[colSums(Crocidura_DRAM_product)>9,]

mean(colnames(Crocidura_dram)==n_pathways$Function)

## BIN quality
Crocidura_quality=read.delim("data/MAG_info/final_bins_Info_Crocidura.csv",sep=",")
Crocidura_quality=Crocidura_quality[Crocidura_quality$genome%in%rownames(Crocidura_dram),]
Crocidura_quality=Crocidura_quality[match(rownames(Crocidura_dram),Crocidura_quality$genome),]

## Taxonomy
Crocidura_taxonomy=read.delim("data/GTDB_tk/gtdbtk.bac120.summary_Crocidura.tsv",sep="\t")
Crocidura_taxonomy=Crocidura_taxonomy[match(rownames(Crocidura_DRAM_product),Crocidura_taxonomy$user_genome),]
mean(Crocidura_taxonomy$user_genome==rownames(Crocidura_DRAM_product))

Crocidura_taxonomyclean=Crocidura_taxonomy %>% 
  mutate(classification = str_replace_all(classification, ".__", "")) %>%
  # Split the classification string into columns for each taxonnomic rank
  separate(col = classification, sep = ";", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
# Remove first column
Crocidura_taxmatrix=Crocidura_taxonomyclean[,c(2:6)]
# Set MAG name as row names of the dataframe
rownames(Crocidura_taxmatrix)=Crocidura_taxonomyclean[,1]

# Pool Firmicutes A,B C into Firmicutes
# taxmatrix_mod=data.frame(taxmatrix)
Crocidura_taxmatrix$Phylum[grepl("Firmi",Crocidura_taxmatrix$Phylum)]="Firmicutes"

mean(rownames(Crocidura_taxmatrix)==rownames(Crocidura_dram))
mean(rownames(Crocidura_taxmatrix)==Crocidura_quality$genome)

Explanatory_dataset=data.frame(Phylum=Crocidura_taxmatrix$Phylum,
                               completeness=Crocidura_quality$completeness-100)

Crocidura_dram_phy=Crocidura_dram[Explanatory_dataset$Phylum=="Bacteroidota"|
                                  Explanatory_dataset$Phylum=="Firmicutes"|
                                  Explanatory_dataset$Phylum=="Proteobacteria",]
Explanatory_dataset_phy=Explanatory_dataset[Explanatory_dataset$Phylum=="Bacteroidota"|
                                              Explanatory_dataset$Phylum=="Firmicutes"|
                                              Explanatory_dataset$Phylum=="Proteobacteria",]
Crocidura_table=data.frame(matrix(ncol=12,nrow = ncol(Crocidura_dram_phy)))
colnames(Crocidura_table)=c("Dataset","Function","Bacteroidota_beta0","Firmicutes_beta0","Proteobacteria_beta0","Bacteroidota_beta1","Firmicutes_beta1","Proteobacteria_beta1","P_int","mean_beta0","mean_beta1","P_compl")
for(i in 1:ncol(Crocidura_dram_phy)){
  response_int=Crocidura_dram_phy[,i]
  response=Crocidura_dram[,i]
  weights_int=rep(n_pathways[i,2],nrow(Crocidura_dram_phy))
  weights=rep(n_pathways[i,2],nrow(Crocidura_dram))
  M_int=glm(response_int~Phylum*completeness,family = binomial,weights = weights_int,data = Explanatory_dataset_phy)
  M=glm(response~completeness,family = binomial,weights = weights,data = Explanatory_dataset)
  Crocidura_table[i,1]="Crocidura"
  Crocidura_table[i,2]=n_pathways[i,1]
  Crocidura_table[i,3]=exp(coef(M_int)[1])/(1+exp(coef(M_int)[1]))
  Crocidura_table[i,4]=exp(coef(M_int)[1]+coef(M_int)[2])/(1+exp(coef(M_int)[1]+coef(M_int)[2]))
  Crocidura_table[i,5]=exp(coef(M_int)[1]+coef(M_int)[3])/(1+exp(coef(M_int)[1]+coef(M_int)[3]))
  Crocidura_table[i,6]=coef(M_int)[4]
  Crocidura_table[i,7]=coef(M_int)[4]+coef(M_int)[5]
  Crocidura_table[i,8]=coef(M_int)[4]+coef(M_int)[6]
  Crocidura_table[i,9]=anova(M_int,test = "Chisq")[4,5]
  Crocidura_table[i,10]=exp(coef(M)[1])/(1+exp(coef(M)[1]))
  Crocidura_table[i,11]=coef(M)[2]
  Crocidura_table[i,12]=anova(M,test = "Chisq")[2,5]
}

head(Crocidura_table)
