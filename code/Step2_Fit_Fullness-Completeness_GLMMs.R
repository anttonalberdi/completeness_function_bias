# Load datasets
# *************
dram_product <- read.csv("data/dram_product_fullness_data.csv",header = TRUE)
rownames(dram_product) <- dram_product$X
dram_product <- dram_product[,-1]
n_steps <- read.csv("data/number_of_steps_in_modules.csv",header = TRUE)
Explanatory_dataset <- read.csv("data/explanatory_dataset.csv",header = TRUE)

Results_table <- data.frame(matrix(ncol=21,nrow = ncol(dram_product)))
colnames(Results_table)=c("Function","Actinobacteriota_slope","Bacteroidota_slope",
                          "Firmicutes_slope","Proteobacteria_slope",
                          "Actinobacteriota_pred_70","Actinobacteriota_pred_80",
                          "Actinobacteriota_pred_90","Actinobacteriota_pred_100",
                          "Bacteroidota_pred_70","Bacteroidota_pred_80",
                          "Bacteroidota_pred_90","Bacteroidota_pred_100",
                          "Firmicutes_pred_70","Firmicutes_pred_80",
                          "Firmicutes_pred_90","Firmicutes_pred_100",
                          "Proteobacteria_pred_70","Proteobacteria_pred_80",
                          "Proteobacteria_pred_90","Proteobacteria_pred_100")
models_list=vector(mode="list",length = ncol(dram_product))

for(i in 1:ncol(dram_product)){
  response=dram_product[,i]
  weights=rep(n_steps[i,3],nrow(dram_product))
  M <- glm(response~Phylum*Completeness,data = Explanatory_dataset,
           weights = weights, family = binomial)
  models_list[[i]] <- M
  Results_table[i,1] <- n_steps$ModuleID[i]
  summ <- summary(M)
  Results_table[i,2] <- summ$coefficients[5,1]
  Results_table[i,3] <- summ$coefficients[5,1]+summ$coefficients[6,1]
  Results_table[i,4] <- summ$coefficients[5,1]+summ$coefficients[7,1]
  Results_table[i,5] <- summ$coefficients[5,1]+summ$coefficients[8,1]
  Results_table[i,6] <- predict(M,newdata = data.frame(Phylum="Actinobacteriota",
                                                          Completeness=70,
                                                          weights=weights[1]),
                             type = "response")
  Results_table[i,7] <- predict(M,newdata = data.frame(Phylum="Actinobacteriota",
                                                       Completeness=80,
                                                       weights=weights[1]),
                                type = "response")
  Results_table[i,8] <- predict(M,newdata = data.frame(Phylum="Actinobacteriota",
                                                       Completeness=90,
                                                       weights=weights[1]),
                                type = "response")
  Results_table[i,9] <- predict(M,newdata = data.frame(Phylum="Actinobacteriota",
                                                       Completeness=100,
                                                       weights=weights[1]),
                                type = "response")
  Results_table[i,10] <- predict(M,newdata = data.frame(Phylum="Bacteroidota",
                                                       Completeness=70,
                                                       weights=weights[1]),
                                type = "response")
  Results_table[i,11] <- predict(M,newdata = data.frame(Phylum="Bacteroidota",
                                                       Completeness=80,
                                                       weights=weights[1]),
                                type = "response")
  Results_table[i,12] <- predict(M,newdata = data.frame(Phylum="Bacteroidota",
                                                       Completeness=90,
                                                       weights=weights[1]),
                                type = "response")
  Results_table[i,13] <- predict(M,newdata = data.frame(Phylum="Bacteroidota",
                                                       Completeness=100,
                                                       weights=weights[1]),
                                type = "response")
  Results_table[i,14] <- predict(M,newdata = data.frame(Phylum="Firmicutes",
                                                        Completeness=70,
                                                        weights=weights[1]),
                                 type = "response")
  Results_table[i,15] <- predict(M,newdata = data.frame(Phylum="Firmicutes",
                                                        Completeness=80,
                                                        weights=weights[1]),
                                 type = "response")
  Results_table[i,16] <- predict(M,newdata = data.frame(Phylum="Firmicutes",
                                                        Completeness=90,
                                                        weights=weights[1]),
                                 type = "response")
  Results_table[i,17] <- predict(M,newdata = data.frame(Phylum="Firmicutes",
                                                        Completeness=100,
                                                        weights=weights[1]),
                                 type = "response")
  Results_table[i,18] <- predict(M,newdata = data.frame(Phylum="Proteobacteria",
                                                        Completeness=70,
                                                        weights=weights[1]),
                                 type = "response")
  Results_table[i,19] <- predict(M,newdata = data.frame(Phylum="Proteobacteria",
                                                        Completeness=80,
                                                        weights=weights[1]),
                                 type = "response")
  Results_table[i,20] <- predict(M,newdata = data.frame(Phylum="Proteobacteria",
                                                        Completeness=90,
                                                        weights=weights[1]),
                                 type = "response")
  Results_table[i,21] <- predict(M,newdata = data.frame(Phylum="Proteobacteria",
                                                        Completeness=100,
                                                        weights=weights[1]),
                                 type = "response")
}

saveRDS(models_list,file = "data/models_list.rds")
write.csv(Results_table,file = "data/Results_table.csv",row.names = FALSE)
