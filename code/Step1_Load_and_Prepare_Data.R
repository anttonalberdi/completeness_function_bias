library(tidyverse)

# Load and prepare the datasets 
# *****************************

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
Apodemus_quality=Apodemus_quality[Apodemus_quality$completeness>=70,]
Apodemus_quality=Apodemus_quality[Apodemus_quality$genome%in%rownames(Apodemus_dram),]
Apodemus_dram=Apodemus_dram[rownames(Apodemus_dram)%in%Apodemus_quality$genome,]
Apodemus_quality=Apodemus_quality[match(rownames(Apodemus_dram),Apodemus_quality$genome),]
mean(rownames(Apodemus_dram)==Apodemus_quality$genome)

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
Apodemus_taxmatrix=Apodemus_taxmatrix[match(rownames(Apodemus_dram),rownames(Apodemus_taxmatrix)),]

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
Crocidura_quality=Crocidura_quality[Crocidura_quality$completeness>=70,]
Crocidura_quality=Crocidura_quality[Crocidura_quality$genome%in%rownames(Crocidura_dram),]
Crocidura_dram=Crocidura_dram[rownames(Crocidura_dram)%in%Crocidura_quality$genome,]
Crocidura_quality=Crocidura_quality[match(rownames(Crocidura_dram),Crocidura_quality$genome),]
mean(rownames(Crocidura_dram)==Crocidura_quality$genome)

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

# Keep functions present in at least 5 MAGs in all hosts in dram table
dram_product_RowRed_ColRed=dram_product_RowRed[,as.logical((colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Apodemus",]>0)>4)*
                                                             (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Crocidura",]>0)>4)*
                                                             (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Felis",]>0)>4)*
                                                             (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Gallus",]>0)>4)*
                                                             (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Mus",]>0)>4)*
                                                             (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Actinobacteriota",]>0)>4)*
                                                             (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Bacteroidota",]>0)>4)*
                                                             (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Firmicutes",]>0)>4)*
                                                             (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Proteobacteria",]>0)>4))]


# Keep functions present in at least 5 MAGs in all hosts in n_steps table
n_steps_RowRed=n_steps[as.logical((colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Apodemus",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Crocidura",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Felis",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Gallus",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Mus",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Actinobacteriota",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Bacteroidota",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Firmicutes",]>0)>4)*
                                    (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Proteobacteria",]>0)>4)),]

# Keep functions present in at least 5 MAGs all hosts in hierarchy table
module_hierarchy_RowRed=module_hierarchy[as.logical((colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Apodemus",]>0)>4)*
                                                      (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Crocidura",]>0)>4)*
                                                      (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Felis",]>0)>4)*
                                                      (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Gallus",]>0)>4)*
                                                      (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Host=="Mus",]>0)>4)*
                                                      (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Actinobacteriota",]>0)>4)*
                                                      (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Bacteroidota",]>0)>4)*
                                                      (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Firmicutes",]>0)>4)*
                                                      (colSums(dram_product_RowRed[Explanatory_dataset_RowRed$Phylum=="Proteobacteria",]>0)>4)),]


# Remove functions with 100% fullness for all species in all hosts in n_steps table
n_steps_RowRed=n_steps_RowRed[!colMeans(dram_product_RowRed_ColRed)==1,]
# Remove functions with 100% fullness for all species in all hosts in n_steps table
module_hierarchy_RowRed=module_hierarchy_RowRed[!colMeans(dram_product_RowRed_ColRed)==1,]
# Remove functions with 100% fullness for all species in all hosts
dram_product_RowRed_ColRed=dram_product_RowRed_ColRed[,!colMeans(dram_product_RowRed_ColRed)==1]

mean(colnames(dram_product_RowRed_ColRed)==n_steps_RowRed$ModuleID)
mean(colnames(dram_product_RowRed_ColRed)==module_hierarchy_RowRed$Module)

## Save unified DRAM product
write.csv(dram_product_RowRed_ColRed,file="data/dram_product_fullness_data.csv")

## Save unified number of steps information for each module
write.csv(n_steps_RowRed,file="data/number_of_steps_in_modules.csv",row.names = FALSE)

## Save unified module hierarchy information
write.csv(module_hierarchy_RowRed,file="data/module_hierarchy.csv",row.names = FALSE)

## Save the explanatory variables table
write.csv(Explanatory_dataset_RowRed,file="data/explanatory_dataset.csv",row.names = FALSE)
