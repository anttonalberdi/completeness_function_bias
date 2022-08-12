library(glmmTMB)

# Load datasets
# *************
dram_product_RowRed_ColRed=read.csv("data/dram_product_fullness_data.csv",header = TRUE)
rownames(dram_product_RowRed_ColRed)=dram_product_RowRed_ColRed$X
dram_product_RowRed_ColRed=dram_product_RowRed_ColRed[,-1]
n_steps_RowRed=read.csv("data/number_of_steps_in_modules.csv",header = TRUE)
Explanatory_dataset_RowRed=read.csv("data/explanatory_dataset.csv",header = TRUE)

# Run Fullness vs. Completeness models 
# ************************************

Results_table=data.frame(matrix(ncol=22,nrow = ncol(dram_product_RowRed_ColRed)))
colnames(Results_table)=c("Function","Random_effect","Actinobacteriota_slope","Bacteroidota_slope",
                          "Firmicutes_slope","Proteobacteria_slope",
                          "Actinobacteriota_pred_70","Actinobacteriota_pred_80",
                          "Actinobacteriota_pred_90","Actinobacteriota_pred_100",
                          "Bacteroidota_pred_70","Bacteroidota_pred_80",
                          "Bacteroidota_pred_90","Bacteroidota_pred_100",
                          "Firmicutes_pred_70","Firmicutes_pred_80",
                          "Firmicutes_pred_90","Firmicutes_pred_100",
                          "Proteobacteria_pred_70","Proteobacteria_pred_80",
                          "Proteobacteria_pred_90","Proteobacteria_pred_100")
models_list=vector(mode="list",length = ncol(dram_product_RowRed_ColRed))
for(i in 1:ncol(dram_product_RowRed_ColRed)){
  response=dram_product_RowRed_ColRed[,i]
  weights=rep(n_steps_RowRed[i,3],nrow(dram_product_RowRed_ColRed))
  M_slope=glmmTMB(response~Phylum*Completeness+(Completeness|Host),data = Explanatory_dataset_RowRed,
                           weights = weights,family = binomial)
  
  M_intercept=glmmTMB(response~Phylum*Completeness+(1|Host),data = Explanatory_dataset_RowRed,
                               weights = weights,family = binomial)
  AIC_slope=AIC(M_slope)
  AIC_intercept=AIC(M_intercept)
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
  }else{
    models_list[[i]]=NA
    Results_table[i,1]=n_steps_RowRed$ModuleID[i]
    Results_table[i,2]="No_convergence"
    Results_table[i,3]=NA
    Results_table[i,4]=NA
    Results_table[i,5]=NA
    Results_table[i,6]=NA
    Results_table[i,7]=NA
    Results_table[i,8]=NA
    Results_table[i,9]=NA
    Results_table[i,10]=NA
    Results_table[i,11]=NA
    Results_table[i,12]=NA
    Results_table[i,13]=NA
    Results_table[i,14]=NA
    Results_table[i,15]=NA
    Results_table[i,16]=NA
    Results_table[i,17]=NA
    Results_table[i,18]=NA
    Results_table[i,19]=NA
    Results_table[i,20]=NA
    Results_table[i,21]=NA
    Results_table[i,22]=NA
  }
}
saveRDS(models_list,file = "data/models_list.rds")
write.csv(Results_table,file = "data/Results_table.csv",row.names = FALSE)