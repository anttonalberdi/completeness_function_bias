library(tidyverse)

# Load and prepare the datasets
# *****************************

DRAM_product <- read.delim(gzfile("data/mag_data/DRAM/product.tsv.gz"),sep="\t")
rownames(DRAM_product) <- DRAM_product$genome
DRAM_product <- DRAM_product[,-1]
# Convert True/False character to 1/0 binary numeric data
for(i in 1:ncol(DRAM_product)){
  if(is.character(DRAM_product[,i])){
    DRAM_product[,i][DRAM_product[,i]=="True"] <- 1
    DRAM_product[,i][DRAM_product[,i]=="False"] <- 0
    DRAM_product[,i] <- as.integer(DRAM_product[,i])
  }
}

str(DRAM_product)

n_steps <- read.csv("data/mag_data/DRAM/step_numbers.csv",sep=";",header = TRUE)
n_steps <- n_steps[match(colnames(DRAM_product),n_steps$ModuleID),]

mean(colnames(DRAM_product)==n_steps$ModuleID)

# Discard completely absent functions from further exploration
dram <- DRAM_product[,colSums(DRAM_product)>=0]
n_steps <- n_steps[colSums(DRAM_product)>=0,]

mean(colnames(dram)==n_steps$ModuleID)

## BIN quality
quality <- read.csv("data/genome_stats.csv")
quality$genome <- gsub("^.{0,3}", "",quality$ď.żAccession)
quality <- quality[match(rownames(dram),quality$genome),]
mean(rownames(dram)==quality$genome)

## Taxonomy
taxonomyclean <- quality %>%
  mutate(classification = str_replace_all(Taxonomy, ".__", "")) %>%
  # Split the classification string into columns for each taxonnomic rank
  separate(col = classification, sep = ";", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
# Remove first column
taxmatrix <- taxonomyclean[,c(6:12)]
# Set MAG name as row names of the dataframe
rownames(taxmatrix) <- taxonomyclean[,5]

# Pool Firmicutes A:H into Firmicutes
taxmatrix$Phylum[grepl("Firmi",taxmatrix$Phylum)] <- "Firmicutes"
taxmatrix <- taxmatrix[match(rownames(dram),rownames(taxmatrix)),]

mean(rownames(taxmatrix)==rownames(dram))
mean(rownames(taxmatrix)==quality$genome)

Explanatory_dataset=data.frame(Phylum=taxmatrix$Phylum,
                               Completeness=quality$Completeness)

# Keep KEGG modules only
module_hierarchy <- read.delim("data/mag_data/DRAM/module_hierarchy.tsv",sep = "\t")
dram_product <- dram[,colnames(dram)%in%module_hierarchy$Module]
n_steps <- n_steps[n_steps$ModuleID%in%module_hierarchy$Module,]
module_hierarchy <- module_hierarchy[module_hierarchy$Module%in%colnames(dram_product),]
module_hierarchy <- module_hierarchy[match(colnames(dram_product),module_hierarchy$Module),]
mean(module_hierarchy$Module==colnames(dram_product))==1
mean(module_hierarchy$Module==n_steps$ModuleID)==1

# Keep functions present in at least 5% of the MAGs in dram table
dram_product_ColRed <- dram_product[,colSums(dram_product>0)>592]
# Keep functions present in at least 10 MAGs in n_steps table
n_steps_RowRed <- n_steps[colSums(dram_product>0)>592,]
# Keep functions present in at least 10 MAGs in hierarchy table
module_hierarchy_RowRed <- module_hierarchy[colSums(dram_product>0)>592,]

mean(colnames(dram_product_ColRed)==n_steps_RowRed$ModuleID)
mean(colnames(dram_product_ColRed)==module_hierarchy_RowRed$Module)

## Save unified DRAM product
write.csv(dram_product_ColRed,file="data/dram_product_fullness_data.csv")

## Save unified number of steps information for each module
write.csv(n_steps_RowRed,file="data/number_of_steps_in_modules.csv",row.names = FALSE)

## Save unified module hierarchy information
write.csv(module_hierarchy_RowRed,file="data/module_hierarchy.csv",row.names = FALSE)

## Save the explanatory variables table
write.csv(Explanatory_dataset,file="data/explanatory_dataset.csv",row.names = FALSE)
