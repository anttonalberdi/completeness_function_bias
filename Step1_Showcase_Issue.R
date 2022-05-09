########################################################################
### Association between MAG completeness and Functional capabilities ###
########################################################################
library(tidyverse)

# 1. Load and prepare the datasets 
# ********************************

### 1.1. Apodemus dataset ###

# DRAM product

Apodemus_DRAM_product=read.delim(gzfile("data/DRAM/Hacked/product_Apodemus.tsv.gz"),sep="\t")
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

n_steps=read.csv("data/DRAM/Hacked/step_numbers.csv",sep=";",header = TRUE)
n_steps=n_steps[match(colnames(Apodemus_DRAM_product),n_steps$ModuleID),]

mean(colnames(Apodemus_DRAM_product)==n_steps$ModuleID)

# Discard completely absent functions from further exploration
Apodemus_dram=Apodemus_DRAM_product[,colSums(Apodemus_DRAM_product)>=0]
n_steps=n_steps[colSums(Apodemus_DRAM_product)>=0,]

mean(colnames(Apodemus_dram)==n_steps$ModuleID)

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

Apodemus_Explanatory_dataset=data.frame(Host="Apodemus",Phylum=Apodemus_taxmatrix$Phylum,
                               Completeness=Apodemus_quality$completeness)

### 1.2. Crocidura dataset ###

# DRAM product

Crocidura_DRAM_product=read.delim(gzfile("data/DRAM/Hacked/product_Crocidura.tsv.gz"),sep="\t")
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

n_steps=read.csv("data/DRAM/Hacked/step_numbers.csv",sep=";",header = TRUE)
n_steps=n_steps[match(colnames(Crocidura_DRAM_product),n_steps$ModuleID),]

mean(colnames(Crocidura_DRAM_product)==n_steps$ModuleID)

# Discard completely absent functions from further exploration
Crocidura_dram=Crocidura_DRAM_product[,colSums(Crocidura_DRAM_product)>=0]
n_steps=n_steps[colSums(Crocidura_DRAM_product)>=0,]

mean(colnames(Crocidura_dram)==n_steps$ModuleID)

## BIN quality
Crocidura_quality=read.delim("data/MAG_info/final_bins_Info_Crocidura.csv",sep=",")
Crocidura_quality=Crocidura_quality[Crocidura_quality$genome%in%rownames(Crocidura_dram),]
Crocidura_quality=Crocidura_quality[match(rownames(Crocidura_dram),Crocidura_quality$genome),]

## Taxonomy
Crocidura_taxonomy=read.delim("data/GTDB_tk/gtdbtk.bac120.summary_Crocidura.tsv",sep="\t")
Crocidura_taxonomy=Crocidura_taxonomy[match(rownames(Crocidura_dram),Crocidura_taxonomy$user_genome),]
mean(Crocidura_taxonomy$user_genome==rownames(Crocidura_dram))

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

Crocidura_Explanatory_dataset=data.frame(Host="Crocidura",Phylum=Crocidura_taxmatrix$Phylum,
                                        Completeness=Crocidura_quality$completeness)

### 1.3. Felis dataset ###

# DRAM product

Felis_DRAM_product=read.delim(gzfile("data/DRAM/Hacked/product_Felis.tsv.gz"),sep="\t")
rownames(Felis_DRAM_product)=Felis_DRAM_product$genome
Felis_DRAM_product=Felis_DRAM_product[,-1]
# Convert True/False character to 1/0 binary numeric data
for(i in 1:ncol(Felis_DRAM_product)){
  if(is.character(Felis_DRAM_product[,i])){
    Felis_DRAM_product[,i][Felis_DRAM_product[,i]=="True"]=1
    Felis_DRAM_product[,i][Felis_DRAM_product[,i]=="False"]=0
    Felis_DRAM_product[,i]=as.integer(Felis_DRAM_product[,i])
  }
}
# Felis_DRAM_product$CAZy_richness=rowSums(Felis_DRAM_product[,grepl("CAZy",colnames(Felis_DRAM_product))])
# Felis_DRAM_product$CAZy_richness=as.integer(Felis_DRAM_product$CAZy_richness)
str(Felis_DRAM_product)

n_steps=read.csv("data/DRAM/Hacked/step_numbers.csv",sep=";",header = TRUE)
n_steps=n_steps[match(colnames(Felis_DRAM_product),n_steps$ModuleID),]

mean(colnames(Felis_DRAM_product)==n_steps$ModuleID)

# Discard completely absent functions from further exploration
Felis_dram=Felis_DRAM_product[,colSums(Felis_DRAM_product)>=0]
n_steps=n_steps[colSums(Felis_DRAM_product)>=0,]

mean(colnames(Felis_dram)==n_steps$ModuleID)

## BIN quality
Felis_quality=read.delim("data/MAG_info/final_bins_Info_Felis.csv",sep=",")
Felis_quality=Felis_quality[Felis_quality$genome%in%rownames(Felis_dram),]
Felis_dram=Felis_dram[rownames(Felis_dram)%in%Felis_quality$genome,]
Felis_quality=Felis_quality[match(rownames(Felis_dram),Felis_quality$genome),]
mean(rownames(Felis_dram)==Felis_quality$genome)

## Taxonomy
Felis_taxonomy=read.delim("data/GTDB_tk/gtdbtk.bac120.summary_Felis.tsv",sep="\t")
Felis_taxonomy=Felis_taxonomy[match(rownames(Felis_dram),Felis_taxonomy$user_genome),]
mean(Felis_taxonomy$user_genome==rownames(Felis_dram))

Felis_taxonomyclean=Felis_taxonomy %>% 
  mutate(classification = str_replace_all(classification, ".__", "")) %>%
  # Split the classification string into columns for each taxonnomic rank
  separate(col = classification, sep = ";", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
# Remove first column
Felis_taxmatrix=Felis_taxonomyclean[,c(2:6)]
# Set MAG name as row names of the dataframe
rownames(Felis_taxmatrix)=Felis_taxonomyclean[,1]

# Pool Firmicutes A,B C into Firmicutes
# taxmatrix_mod=data.frame(taxmatrix)
Felis_taxmatrix$Phylum[grepl("Firmi",Felis_taxmatrix$Phylum)]="Firmicutes"

mean(rownames(Felis_taxmatrix)==rownames(Felis_dram))
mean(rownames(Felis_taxmatrix)==Felis_quality$genome)

Felis_Explanatory_dataset=data.frame(Host="Felis",Phylum=Felis_taxmatrix$Phylum,
                                        Completeness=Felis_quality$completeness)

### 1.4. Gallus dataset ###

# DRAM product

Gallus_DRAM_product=read.delim(gzfile("data/DRAM/Hacked/product_Gallus.tsv.gz"),sep="\t")
rownames(Gallus_DRAM_product)=Gallus_DRAM_product$genome
Gallus_DRAM_product=Gallus_DRAM_product[,-1]
# Convert True/False character to 1/0 binary numeric data
for(i in 1:ncol(Gallus_DRAM_product)){
  if(is.character(Gallus_DRAM_product[,i])){
    Gallus_DRAM_product[,i][Gallus_DRAM_product[,i]=="True"]=1
    Gallus_DRAM_product[,i][Gallus_DRAM_product[,i]=="False"]=0
    Gallus_DRAM_product[,i]=as.integer(Gallus_DRAM_product[,i])
  }
}
# Gallus_DRAM_product$CAZy_richness=rowSums(Gallus_DRAM_product[,grepl("CAZy",colnames(Gallus_DRAM_product))])
# Gallus_DRAM_product$CAZy_richness=as.integer(Gallus_DRAM_product$CAZy_richness)
str(Gallus_DRAM_product)

n_steps=read.csv("data/DRAM/Hacked/step_numbers.csv",sep=";",header = TRUE)
n_steps=n_steps[match(colnames(Gallus_DRAM_product),n_steps$ModuleID),]

mean(colnames(Gallus_DRAM_product)==n_steps$ModuleID)

# Discard completely absent functions from further exploration
Gallus_dram=Gallus_DRAM_product[,colSums(Gallus_DRAM_product)>=0]
n_steps=n_steps[colSums(Gallus_DRAM_product)>=0,]

mean(colnames(Gallus_dram)==n_steps$ModuleID)

## BIN quality
Gallus_quality=read.delim("data/MAG_info/final_bins_Info_Gallus.csv",sep=",")
Gallus_quality=Gallus_quality[Gallus_quality$completeness>=70,]
Gallus_quality=Gallus_quality[Gallus_quality$genome%in%rownames(Gallus_dram),]
Gallus_dram=Gallus_dram[rownames(Gallus_dram)%in%Gallus_quality$genome,]
Gallus_quality=Gallus_quality[match(rownames(Gallus_dram),Gallus_quality$genome),]
mean(rownames(Gallus_dram)==Gallus_quality$genome)

## Taxonomy
Gallus_taxonomy=read.delim("data/GTDB_tk/gtdbtk.bac120.summary_Gallus.tsv",sep="\t")
Gallus_taxonomy=Gallus_taxonomy[match(rownames(Gallus_dram),Gallus_taxonomy$user_genome),]
mean(Gallus_taxonomy$user_genome==rownames(Gallus_dram))

Gallus_taxonomyclean=Gallus_taxonomy %>% 
  mutate(classification = str_replace_all(classification, ".__", "")) %>%
  # Split the classification string into columns for each taxonnomic rank
  separate(col = classification, sep = ";", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
# Remove first column
Gallus_taxmatrix=Gallus_taxonomyclean[,c(2:6)]
# Set MAG name as row names of the dataframe
rownames(Gallus_taxmatrix)=Gallus_taxonomyclean[,1]

# Pool Firmicutes A,B C into Firmicutes
# taxmatrix_mod=data.frame(taxmatrix)
Gallus_taxmatrix$Phylum[grepl("Firmi",Gallus_taxmatrix$Phylum)]="Firmicutes"

mean(rownames(Gallus_taxmatrix)==rownames(Gallus_dram))
mean(rownames(Gallus_taxmatrix)==Gallus_quality$genome)

Gallus_Explanatory_dataset=data.frame(Host="Gallus",Phylum=Gallus_taxmatrix$Phylum,
                                        Completeness=Gallus_quality$completeness)

### 1.5. Mus dataset ###

# DRAM product

Mus_DRAM_product=read.delim(gzfile("data/DRAM/Hacked/product_Mus.tsv.gz"),sep="\t")
rownames(Mus_DRAM_product)=Mus_DRAM_product$genome
Mus_DRAM_product=Mus_DRAM_product[,-1]
# Convert True/False character to 1/0 binary numeric data
for(i in 1:ncol(Mus_DRAM_product)){
  if(is.character(Mus_DRAM_product[,i])){
    Mus_DRAM_product[,i][Mus_DRAM_product[,i]=="True"]=1
    Mus_DRAM_product[,i][Mus_DRAM_product[,i]=="False"]=0
    Mus_DRAM_product[,i]=as.integer(Mus_DRAM_product[,i])
  }
}
# Mus_DRAM_product$CAZy_richness=rowSums(Mus_DRAM_product[,grepl("CAZy",colnames(Mus_DRAM_product))])
# Mus_DRAM_product$CAZy_richness=as.integer(Mus_DRAM_product$CAZy_richness)
str(Mus_DRAM_product)

n_steps=read.csv("data/DRAM/Hacked/step_numbers.csv",sep=";",header = TRUE)
n_steps=n_steps[match(colnames(Mus_DRAM_product),n_steps$ModuleID),]

mean(colnames(Mus_DRAM_product)==n_steps$ModuleID)

# Discard completely absent functions from further exploration
Mus_dram=Mus_DRAM_product[,colSums(Mus_DRAM_product)>=0]
n_steps=n_steps[colSums(Mus_DRAM_product)>=0,]

mean(colnames(Mus_dram)==n_steps$ModuleID)

## BIN quality
Mus_quality=read.delim("data/MAG_info/final_bins_Info_Mus.csv",sep=",")
Mus_quality=Mus_quality[Mus_quality$genome%in%rownames(Mus_dram),]
Mus_dram=Mus_dram[rownames(Mus_dram)%in%Mus_quality$genome,]
Mus_quality=Mus_quality[match(rownames(Mus_dram),Mus_quality$genome),]
mean(rownames(Mus_dram)==Mus_quality$genome)

## Taxonomy
Mus_taxonomy=read.delim("data/GTDB_tk/gtdbtk.bac120.summary_Mus.tsv",sep="\t")
Mus_taxonomy=Mus_taxonomy[match(rownames(Mus_dram),Mus_taxonomy$user_genome),]
mean(Mus_taxonomy$user_genome==rownames(Mus_dram))

Mus_taxonomyclean=Mus_taxonomy %>% 
  mutate(classification = str_replace_all(classification, ".__", "")) %>%
  # Split the classification string into columns for each taxonnomic rank
  separate(col = classification, sep = ";", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
# Remove first column
Mus_taxmatrix=Mus_taxonomyclean[,c(2:6)]
# Set MAG name as row names of the dataframe
rownames(Mus_taxmatrix)=Mus_taxonomyclean[,1]

# Pool Firmicutes A,B C into Firmicutes
# taxmatrix_mod=data.frame(taxmatrix)
Mus_taxmatrix$Phylum[grepl("Firmi",Mus_taxmatrix$Phylum)]="Firmicutes"

mean(rownames(Mus_taxmatrix)==rownames(Mus_dram))
mean(rownames(Mus_taxmatrix)==Mus_quality$genome)

Mus_Explanatory_dataset=data.frame(Host="Mus",Phylum=Mus_taxmatrix$Phylum,
                                        Completeness=Mus_quality$completeness)


### 1.6. Combine datasets ###

Explanatory_dataset=data.frame(rbind(Apodemus_Explanatory_dataset,
                                     Crocidura_Explanatory_dataset,
                                     Felis_Explanatory_dataset,
                                     Gallus_Explanatory_dataset,
                                     Mus_Explanatory_dataset))
dram_product=data.frame(rbind(Apodemus_dram,
                              Crocidura_dram,
                              Felis_dram,
                              Gallus_dram,
                              Mus_dram))

# Keep Phyla represented in all hosts
dram_product_RowRed=dram_product[Explanatory_dataset$Phylum=="Actinobacteriota"|
                                   Explanatory_dataset$Phylum=="Bacteroidota"|
                                   Explanatory_dataset$Phylum=="Firmicutes"|
                                   Explanatory_dataset$Phylum=="Proteobacteria",]
# Keep Phyla represented in all hosts
Explanatory_dataset_RowRed=Explanatory_dataset[Explanatory_dataset$Phylum=="Actinobacteriota"|
                                                 Explanatory_dataset$Phylum=="Bacteroidota"|
                                                 Explanatory_dataset$Phylum=="Firmicutes"|
                                                 Explanatory_dataset$Phylum=="Proteobacteria",]

# Keep KEGG modules only
module_hierarchy=read.delim("data/DRAM/Hacked/module_hierarchy.tsv",sep = "\t")
dram_product_RowRed=dram_product_RowRed[,colnames(dram_product_RowRed)%in%module_hierarchy$Module]
n_steps=n_steps[n_steps$ModuleID%in%module_hierarchy$Module,]
module_hierarchy=module_hierarchy[module_hierarchy$Module%in%colnames(dram_product_RowRed),]
module_hierarchy=module_hierarchy[match(colnames(dram_product_RowRed),module_hierarchy$Module),]
mean(module_hierarchy$Module==colnames(dram_product_RowRed))==1
mean(module_hierarchy$Module==n_steps$ModuleID)==1

# Keep functions present in at least 3 MAGs in all hosts in dram table
dram_product_RowRed_ColRed=dram_product_RowRed[,as.logical((colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Apodemus",]>0)>4)*
                                                            (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Crocidura",]>0)>4)*
                                                            (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Felis",]>0)>4)*
                                                            (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Gallus",]>0)>4)*
                                                            (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Mus",]>0)>4))]
# Keep functions present in all hosts in n_steps table
n_steps_RowRed=n_steps[as.logical((colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Apodemus",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Crocidura",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Felis",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Gallus",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Mus",]>0)>4)),]
# Keep functions present in all hosts in hierarchy table
module_hierarchy_RowRed=module_hierarchy[as.logical((colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Apodemus",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Crocidura",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Felis",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Gallus",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Mus",]>0)>4)),]


# Remove functions with 100% fullness for all species in all hosts in n_steps table
n_steps_RowRed=n_steps_RowRed[!colMeans(dram_product_RowRed_ColRed)==1,]
# Remove functions with 100% fullness for all species in all hosts in n_steps table
module_hierarchy_RowRed=module_hierarchy_RowRed[!colMeans(dram_product_RowRed_ColRed)==1,]
# Remove functions with 100% fullness for all species in all hosts
dram_product_RowRed_ColRed=dram_product_RowRed_ColRed[,!colMeans(dram_product_RowRed_ColRed)==1]


# 2. Run Fullnes vs. Completeness models 
# **************************************

Results_table=data.frame(matrix(ncol=26,nrow = ncol(dram_product_RowRed_ColRed)))
colnames(Results_table)=c("Function","Random_effect","Actinobacteriota_slope","Bacteroidota_slope",
                          "Firmicutes_slope","Proteobacteria_slope",
                          "Actinobacteriota_pred_70","Actinobacteriota_pred_80",
                          "Actinobacteriota_pred_90","Actinobacteriota_pred_100",
                          "Bacteroidota_pred_70","Bacteroidota_pred_80",
                          "Bacteroidota_pred_90","Bacteroidota_pred_100",
                          "Firmicutes_pred_70","Firmicutes_pred_80",
                          "Firmicutes_pred_90","Firmicutes_pred_100",
                          "Proteobacteria_pred_70","Proteobacteria_pred_80",
                          "Proteobacteria_pred_90","Proteobacteria_pred_100",
                          "R2_m","R2_c","pvalue","pvalue_int")
models_list=vector(mode="list",length = ncol(dram_product_RowRed_ColRed))
for(i in 1:ncol(dram_product_RowRed_ColRed)){
  response=dram_product_RowRed_ColRed[,i]
  weights=rep(n_steps_RowRed[i,3],nrow(dram_product_RowRed_ColRed))
  M_slope=glmmTMB::glmmTMB(response~Phylum*Completeness+(Completeness|Host),data = Explanatory_dataset_RowRed,
                           weights = weights,family = binomial)
  
  M_intercept=glmmTMB::glmmTMB(response~Phylum*Completeness+(1|Host),data = Explanatory_dataset_RowRed,
                           weights = weights,family = binomial)
  M=stats::glm(response~Phylum*Completeness,data = Explanatory_dataset_RowRed,
        weights = weights,family = binomial)
  AIC_slope=AIC(M_slope)
  AIC_intercept=AIC(M_intercept)
  AIC_m=AIC(M)
  cond_slope=!(any(is.nan(summary(M_slope)$coefficients$cond))|grepl("singular",M_slope$fit$message)|(AIC_slope>AIC_intercept))
  if(is.na(cond_slope)){
    cond_slope=FALSE
  }
  cond_intercept=!(any(is.nan(summary(M_intercept)$coefficients$cond))|grepl("singular",M_intercept$fit$message))
  if(is.na(cond_intercept)){
    cond_intercept=FALSE
  }
  if(cond_slope){
    models_list[[i]]=M_slope
    Results_table[i,1]=n_steps_RowRed$ModuleID[i]
    Results_table[i,2]="Random_slope"
    summ=summary(M_slope)
    Results_table[i,3]=summ$coefficients$cond[5,1]
    Results_table[i,4]=summ$coefficients$cond[5,1]+summ$coefficients$cond[6,1]
    Results_table[i,5]=summ$coefficients$cond[5,1]+summ$coefficients$cond[7,1]
    Results_table[i,6]=summ$coefficients$cond[5,1]+summ$coefficients$cond[8,1]
    Results_table[i,7]=predict(M_slope,newdata = data.frame(Phylum="Actinobacteriota",
                                                            Completeness=70,Host=NA,
                                                            weights=weights[1]),
                               type = "response")
    Results_table[i,8]=predict(M_slope,newdata = data.frame(Phylum="Actinobacteriota",
                                                            Completeness=80,Host=NA,
                                                            weights=weights[1]),
                               type = "response")
    Results_table[i,9]=predict(M_slope,newdata = data.frame(Phylum="Actinobacteriota",
                                                            Completeness=90,Host=NA,
                                                            weights=weights[1]),
                               type = "response")
    Results_table[i,10]=predict(M_slope,newdata = data.frame(Phylum="Actinobacteriota",
                                                             Completeness=100,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,11]=predict(M_slope,newdata = data.frame(Phylum="Bacteroidota",
                                                             Completeness=70,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,12]=predict(M_slope,newdata = data.frame(Phylum="Bacteroidota",
                                                             Completeness=80,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,13]=predict(M_slope,newdata = data.frame(Phylum="Bacteroidota",
                                                             Completeness=90,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,14]=predict(M_slope,newdata = data.frame(Phylum="Bacteroidota",
                                                             Completeness=100,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,15]=predict(M_slope,newdata = data.frame(Phylum="Firmicutes",
                                                             Completeness=70,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,16]=predict(M_slope,newdata = data.frame(Phylum="Firmicutes",
                                                             Completeness=80,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,17]=predict(M_slope,newdata = data.frame(Phylum="Firmicutes",
                                                             Completeness=90,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,18]=predict(M_slope,newdata = data.frame(Phylum="Firmicutes",
                                                             Completeness=100,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,19]=predict(M_slope,newdata = data.frame(Phylum="Proteobacteria",
                                                             Completeness=70,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,20]=predict(M_slope,newdata = data.frame(Phylum="Proteobacteria",
                                                             Completeness=80,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,21]=predict(M_slope,newdata = data.frame(Phylum="Proteobacteria",
                                                             Completeness=90,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,22]=predict(M_slope,newdata = data.frame(Phylum="Proteobacteria",
                                                             Completeness=100,Host=NA,
                                                             weights=weights[1]),
                                type = "response")
    Results_table[i,23]=MuMIn::r.squaredGLMM(M_slope)[2,1]
    Results_table[i,24]=MuMIn::r.squaredGLMM(M_slope)[2,2]
    Results_table[i,25]=car::Anova(M_slope,test="Chisq")[2,3]
    Results_table[i,26]=car::Anova(M_slope,test="Chisq")[3,3]
  }else if(cond_intercept){
    models_list[[i]]=M_intercept
    Results_table[i,1]=n_steps_RowRed$ModuleID[i]
    Results_table[i,2]="Random_intercept"
    summ=summary(M_intercept)
    Results_table[i,3]=summ$coefficients$cond[5,1]
    Results_table[i,4]=summ$coefficients$cond[5,1]+summ$coefficients$cond[6,1]
    Results_table[i,5]=summ$coefficients$cond[5,1]+summ$coefficients$cond[7,1]
    Results_table[i,6]=summ$coefficients$cond[5,1]+summ$coefficients$cond[8,1]
    Results_table[i,7]=predict(M_intercept,newdata = data.frame(Phylum="Actinobacteriota",
                                                                Completeness=70,Host=NA,
                                                                weights=weights[1]),
                               type = "response")
    Results_table[i,8]=predict(M_intercept,newdata = data.frame(Phylum="Actinobacteriota",
                                                                Completeness=80,Host=NA,
                                                                weights=weights[1]),
                               type = "response")
    Results_table[i,9]=predict(M_intercept,newdata = data.frame(Phylum="Actinobacteriota",
                                                                Completeness=90,Host=NA,
                                                                weights=weights[1]),
                               type = "response")
    Results_table[i,10]=predict(M_intercept,newdata = data.frame(Phylum="Actinobacteriota",
                                                                 Completeness=100,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,11]=predict(M_intercept,newdata = data.frame(Phylum="Bacteroidota",
                                                                 Completeness=70,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,12]=predict(M_intercept,newdata = data.frame(Phylum="Bacteroidota",
                                                                 Completeness=80,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,13]=predict(M_intercept,newdata = data.frame(Phylum="Bacteroidota",
                                                                 Completeness=90,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,14]=predict(M_intercept,newdata = data.frame(Phylum="Bacteroidota",
                                                                 Completeness=100,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,15]=predict(M_intercept,newdata = data.frame(Phylum="Firmicutes",
                                                                 Completeness=70,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,16]=predict(M_intercept,newdata = data.frame(Phylum="Firmicutes",
                                                                 Completeness=80,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,17]=predict(M_intercept,newdata = data.frame(Phylum="Firmicutes",
                                                                 Completeness=90,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,18]=predict(M_intercept,newdata = data.frame(Phylum="Firmicutes",
                                                                 Completeness=100,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,19]=predict(M_intercept,newdata = data.frame(Phylum="Proteobacteria",
                                                                 Completeness=70,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,20]=predict(M_intercept,newdata = data.frame(Phylum="Proteobacteria",
                                                                 Completeness=80,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,21]=predict(M_intercept,newdata = data.frame(Phylum="Proteobacteria",
                                                                 Completeness=90,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,22]=predict(M_intercept,newdata = data.frame(Phylum="Proteobacteria",
                                                                 Completeness=100,Host=NA,
                                                                 weights=weights[1]),
                                type = "response")
    Results_table[i,23]=MuMIn::r.squaredGLMM(M_intercept)[2,1]
    Results_table[i,24]=MuMIn::r.squaredGLMM(M_intercept)[2,2]
    Results_table[i,25]=car::Anova(M_intercept,test="Chisq")[2,3]
    Results_table[i,26]=car::Anova(M_intercept,test="Chisq")[3,3]
    }else{
      models_list[[i]]=M
      Results_table[i,1]=n_steps_RowRed$ModuleID[i]
      Results_table[i,2]="No_random"
      summ=summary(M)
      Results_table[i,3]=summ$coefficients[5,1]
      Results_table[i,4]=summ$coefficients[5,1]+summ$coefficients[6,1]
      Results_table[i,5]=summ$coefficients[5,1]+summ$coefficients[7,1]
      Results_table[i,6]=summ$coefficients[5,1]+summ$coefficients[8,1]
      Results_table[i,7]=predict(M,newdata = data.frame(Phylum="Actinobacteriota",
                                                                  Completeness=70,Host=NA,
                                                                  weights=weights[1]),
                                 type = "response")
      Results_table[i,8]=predict(M,newdata = data.frame(Phylum="Actinobacteriota",
                                                                  Completeness=80,Host=NA,
                                                                  weights=weights[1]),
                                 type = "response")
      Results_table[i,9]=predict(M,newdata = data.frame(Phylum="Actinobacteriota",
                                                                  Completeness=90,Host=NA,
                                                                  weights=weights[1]),
                                 type = "response")
      Results_table[i,10]=predict(M,newdata = data.frame(Phylum="Actinobacteriota",
                                                                   Completeness=100,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,11]=predict(M,newdata = data.frame(Phylum="Bacteroidota",
                                                                   Completeness=70,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,12]=predict(M,newdata = data.frame(Phylum="Bacteroidota",
                                                                   Completeness=80,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,13]=predict(M,newdata = data.frame(Phylum="Bacteroidota",
                                                                   Completeness=90,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,14]=predict(M,newdata = data.frame(Phylum="Bacteroidota",
                                                                   Completeness=100,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,15]=predict(M,newdata = data.frame(Phylum="Firmicutes",
                                                                   Completeness=70,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,16]=predict(M,newdata = data.frame(Phylum="Firmicutes",
                                                                   Completeness=80,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,17]=predict(M,newdata = data.frame(Phylum="Firmicutes",
                                                                   Completeness=90,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,18]=predict(M,newdata = data.frame(Phylum="Firmicutes",
                                                                   Completeness=100,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,19]=predict(M,newdata = data.frame(Phylum="Proteobacteria",
                                                                   Completeness=70,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,20]=predict(M,newdata = data.frame(Phylum="Proteobacteria",
                                                                   Completeness=80,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,21]=predict(M,newdata = data.frame(Phylum="Proteobacteria",
                                                                   Completeness=90,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,22]=predict(M,newdata = data.frame(Phylum="Proteobacteria",
                                                                   Completeness=100,Host=NA,
                                                                   weights=weights[1]),
                                  type = "response")
      Results_table[i,23]=MuMIn::r.squaredGLMM(M)[2,1]
      Results_table[i,24]=MuMIn::r.squaredGLMM(M)[2,2]
      Results_table[i,25]=car::Anova(M,test="Wald")[2,3]
      Results_table[i,26]=car::Anova(M,test="Wald")[3,3]
  }
}
# saveRDS(models_list,file = "models_list.rds")
# write.csv(Results_table,file = "Results_table.csv",row.names = FALSE)
models_list=readRDS(file="models_list.rds")
Results_table=read.csv("Results_table.csv")

# 49% of the functions benefit from a random slope (with respect to random intercept) model  
mean(Results_table$Random_effect=="Random_slope")
mean(Results_table$Random_effect=="Random_intercept")
# 8% of the functions did not get convergence using mixed models
mean(Results_table$Random_effect=="No_random")

models_list_red=models_list[Results_table$Random_effect!="No_random"]
Results_table_red=Results_table[Results_table$Random_effect!="No_random",]
module_hierarchy_RowRed=module_hierarchy_RowRed[Results_table$Random_effect!="No_random",]


toplot_re_slope=which(Results_table_red$Random_effect=="Random_slope"&(Results_table_red$pvalue_int<0.05|Results_table_red$pvalue<0.05))
set.seed(1)
plot_list=list()
for(i in 1:16){
  index=sample(toplot_re_slope,1)
  p=sjPlot::plot_model(models_list_red[[index]],
                     type = "pred",
                     terms = c("Completeness","Host"),
                     pred.type = "re",
                     title = paste(Results_table_red$Function[index]))
  plot_list[[i]]=p
}
for (i in 1:16) {
  file_name = paste("model_exploration_plots/re_slope_plot_", i, ".tiff", sep="")
  tiff(file_name)
  print(plot_list[[i]])
  dev.off()
}

# Significant Phylum*Completeness 25% of the functions
mean(Results_table_red$pvalue_int<0.05)
# Significant Completeness 25% of the functions
mean(Results_table_red$pvalue<0.05)
# Positive parameter value for Completeness
mean(Results_table_red$Actinobacteriota_slope[(Results_table_red$pvalue_int<0.05|Results_table_red$pvalue<0.05)]>0)
mean(Results_table_red$Bacteroidota_slope[(Results_table_red$pvalue_int<0.05|Results_table_red$pvalue<0.05)]>0)
mean(Results_table_red$Firmicutes_slope[(Results_table_red$pvalue_int<0.05|Results_table_red$pvalue<0.05)]>0)
mean(Results_table_red$Proteobacteria_slope[(Results_table_red$pvalue_int<0.05|Results_table_red$pvalue<0.05)]>0)

toplot_int_Phylum=which(Results_table_red$pvalue_int<0.05)
set.seed(1)
plot_list=list()
for(i in 1:16){
  index=sample(toplot_int_Phylum,1)
  p=sjPlot::plot_model(models_list_red[[index]],
                       type = "pred",
                       terms = c("Completeness","Phylum"),
                       title = paste(Results_table_red$Function[index]))
  plot_list[[i]]=p
}
for (i in 1:16) {
  file_name = paste("model_exploration_plots/int_Phylum_plot_", i, ".tiff", sep="")
  tiff(file_name)
  print(plot_list[[i]])
  dev.off()
}


# Distribution of parameter values across Phyla
boxplot(Results_table_red$Actinobacteriota_slope,Results_table_red$Bacteroidota_slope,
        Results_table_red$Firmicutes_slope,Results_table_red$Proteobacteria_slope)
abline(h=0)
boxplot(Results_table_red$Actinobacteriota_slope,Results_table_red$Bacteroidota_slope,
        Results_table_red$Firmicutes_slope,Results_table_red$Proteobacteria_slope,
        ylim=c(-0.2,0.2))
abline(h=0)

# Fullness values when completeness = 70%
par(mfrow=c(1,2))
boxplot(Results_table_red$Actinobacteriota_pred_70,Results_table_red$Bacteroidota_pred_70,
        Results_table_red$Firmicutes_pred_70,Results_table_red$Proteobacteria_pred_70,
        ylim=c(0,1),xaxt="n",main="Fullness at 70% Completeness")
axis(side=1,at=c(1,2,3,4),labels = c("Actino","Bactero","Firmi","Proteo"))
# Fullness values when completeness = 100%
boxplot(Results_table_red$Actinobacteriota_pred_100,Results_table_red$Bacteroidota_pred_100,
        Results_table_red$Firmicutes_pred_100,Results_table_red$Proteobacteria_pred_100,
        ylim=c(0,1),xaxt="n",main="Fullness at 100% Completeness")
axis(side=1,at=c(1,2,3,4),labels = c("Actino","Bactero","Firmi","Proteo"))

# Diff Fullness between Completeness 70%-100% by Phylum
par(mfrow=c(2,2))
Diff_100_70_Actino=Results_table_red[,10]-Results_table_red[,7]
Diff_100_70_Bactero=Results_table_red[,14]-Results_table_red[,11]
Diff_100_70_Firmi=Results_table_red[,18]-Results_table_red[,15]
Diff_100_70_Proteo=Results_table_red[,22]-Results_table_red[,19]
boxplot(Diff_100_70_Actino,Diff_100_70_Bactero,
        Diff_100_70_Firmi,Diff_100_70_Proteo,
        xaxt="n",main="Diff in Fullness 70%-100% Completeness",
        ylim=c(-0.1,0.4))
axis(side=1,at=c(1,2,3,4),labels = c("Actino","Bactero","Firmi","Proteo"))
abline(h=0)
# Diff Fullness between Completeness 70%-80% by Phylum
Diff_80_70_Actino=Results_table_red[,8]-Results_table_red[,7]
Diff_80_70_Bactero=Results_table_red[,12]-Results_table_red[,11]
Diff_80_70_Firmi=Results_table_red[,16]-Results_table_red[,15]
Diff_80_70_Proteo=Results_table_red[,20]-Results_table_red[,19]
boxplot(Diff_80_70_Actino,Diff_80_70_Bactero,
        Diff_80_70_Firmi,Diff_80_70_Proteo,
        xaxt="n",main="Diff in Fullness 70%-80% Completeness",
        ylim=c(-0.1,0.4))
axis(side=1,at=c(1,2,3,4),labels = c("Actino","Bactero","Firmi","Proteo"))
abline(h=0)
# Diff Fullness between Completeness 80%-90% by Phylum
Diff_90_80_Actino=Results_table_red[,9]-Results_table_red[,8]
Diff_90_80_Bactero=Results_table_red[,13]-Results_table_red[,12]
Diff_90_80_Firmi=Results_table_red[,17]-Results_table_red[,16]
Diff_90_80_Proteo=Results_table_red[,21]-Results_table_red[,20]
boxplot(Diff_90_80_Actino,Diff_90_80_Bactero,
        Diff_90_80_Firmi,Diff_90_80_Proteo,
        xaxt="n",main="Diff in Fullness 80%-90% Completeness",
        ylim=c(-0.1,0.4))
axis(side=1,at=c(1,2,3,4),labels = c("Actino","Bactero","Firmi","Proteo"))
abline(h=0)
# Diff Fullness between Completeness 80%-90% by Phylum
Diff_100_90_Actino=Results_table_red[,10]-Results_table_red[,9]
Diff_100_90_Bactero=Results_table_red[,14]-Results_table_red[,13]
Diff_100_90_Firmi=Results_table_red[,18]-Results_table_red[,17]
Diff_100_90_Proteo=Results_table_red[,22]-Results_table_red[,21]
boxplot(Diff_100_90_Actino,Diff_100_90_Bactero,
        Diff_100_90_Firmi,Diff_100_90_Proteo,
        xaxt="n",main="Diff in Fullness 90%-100% Completeness",
        ylim=c(-0.1,0.4))
axis(side=1,at=c(1,2,3,4),labels = c("Actino","Bactero","Firmi","Proteo"))
abline(h=0)

# Fullness vs. Completeness by Phyla
par(mfrow=c(2,2))
Actinobacteriota_preds=Results_table_red[,c(7:10)]
Actinobacteriota_preds_long=Actinobacteriota_preds %>%
  pivot_longer(cols = everything())
Actinobacteriota_preds_long$completeness=rep(c(70,80,90,100),nrow(Actinobacteriota_preds_long)/4)
Actinobacteriota_preds_long$module=rep(Results_table_red$Function,each=4)
plot(value~completeness,data = Actinobacteriota_preds_long[Actinobacteriota_preds_long$module=="M00001",],
     type="l",ylim=c(0,1),ylab = "Fullness",main="Actinobacteriota Fullness vs. Completeness",col="red")
for(i in unique(Actinobacteriota_preds_long$module)){
  data=Actinobacteriota_preds_long[Actinobacteriota_preds_long$module==i,]
  if(data$value[4]-data$value[1]>0){
    lines(value~completeness,data = data,col="red")
  }else{
    lines(value~completeness,data = data,col="blue")
  }
}
Bacteroidota_preds=Results_table_red[,c(11:14)]
Bacteroidota_preds_long=Bacteroidota_preds %>%
  pivot_longer(cols = everything())
Bacteroidota_preds_long$completeness=rep(c(70,80,90,100),nrow(Bacteroidota_preds_long)/4)
Bacteroidota_preds_long$module=rep(Results_table_red$Function,each=4)
plot(value~completeness,data = Bacteroidota_preds_long[Bacteroidota_preds_long$module=="M00001",],
     type="l",ylim=c(0,1),ylab = "Fullness",main="Bacteroidota Fullness vs. Completeness",col="red")
for(i in unique(Bacteroidota_preds_long$module)){
  data=Bacteroidota_preds_long[Bacteroidota_preds_long$module==i,]
  if(data$value[4]-data$value[1]>0){
    lines(value~completeness,data = data,col="red")
  }else{
    lines(value~completeness,data = data,col="blue")
  }
}
Firmicutes_preds=Results_table_red[,c(15:18)]
Firmicutes_preds_long=Firmicutes_preds %>%
  pivot_longer(cols = everything())
Firmicutes_preds_long$completeness=rep(c(70,80,90,100),nrow(Firmicutes_preds_long)/4)
Firmicutes_preds_long$module=rep(Results_table_red$Function,each=4)
plot(value~completeness,data = Firmicutes_preds_long[Firmicutes_preds_long$module=="M00001",],
     type="l",ylim=c(0,1),ylab = "Fullness",main="Firmicutes Fullness vs. Completeness",col="red")
for(i in unique(Firmicutes_preds_long$module)){
  data=Firmicutes_preds_long[Firmicutes_preds_long$module==i,]
  if(data$value[4]-data$value[1]>0){
    lines(value~completeness,data = data,col="red")
  }else{
    lines(value~completeness,data = data,col="blue")
  }
}
Proteobacteria_preds=Results_table_red[,c(19:22)]
Proteobacteria_preds_long=Proteobacteria_preds %>%
  pivot_longer(cols = everything())
Proteobacteria_preds_long$completeness=rep(c(70,80,90,100),nrow(Proteobacteria_preds_long)/4)
Proteobacteria_preds_long$module=rep(Results_table_red$Function,each=4)
plot(value~completeness,data = Proteobacteria_preds_long[Proteobacteria_preds_long$module=="M00001",],
     type="l",ylim=c(0,1),ylab = "Fullness",main="Proteobacteria Fullness vs. Completeness",col="red")
for(i in unique(Proteobacteria_preds_long$module)){
  data=Proteobacteria_preds_long[Proteobacteria_preds_long$module==i,]
  if(data$value[4]-data$value[1]>0){
    lines(value~completeness,data = data,col="red")
  }else{
    lines(value~completeness,data = data,col="blue")
  }
}

boxplot(Results_table_red$R2_m,Results_table_red$R2_c-Results_table_red$R2_m)
par(mfrow=c(1,2))
boxplot(Results_table_red$R2_m[Results_table_red$Random_effect=="Random_slope"],Results_table_red$R2_c[Results_table_red$Random_effect=="Random_slope"]-Results_table_red$R2_m[Results_table_red$Random_effect=="Random_slope"])
boxplot(Results_table_red$R2_m[Results_table_red$Random_effect=="Random_intercept"],Results_table_red$R2_c[Results_table_red$Random_effect=="Random_intercept"]-Results_table_red$R2_m[Results_table_red$Random_effect=="Random_intercept"])

dim(Results_table_red)
dim(module_hierarchy_RowRed)

boxplot(Results_table_red$Actinobacteriota_slope~module_hierarchy_RowRed$Domain)

