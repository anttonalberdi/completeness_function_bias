library(tidyverse)
library(lme4)

# Load datasets
# *************
Results_table=read.csv("data/Results_table.csv")
n_steps_RowRed=read.csv("data/number_of_steps_in_modules.csv",header = TRUE)
module_hierarchy_RowRed=read.csv("data/module_hierarchy.csv")

# Prepare data to fit the model
# *****************************
Slope_data=data.frame(Results_table$Actinobacteriota_slope,Results_table$Bacteroidota_slope,
                      Results_table$Firmicutes_slope,Results_table$Proteobacteria_slope) %>%
  pivot_longer(everything())
Phylum=rep(c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria"),nrow(Results_table))
Domain=rep(module_hierarchy_RowRed$Domain,each=4)
Steps=rep(n_steps_RowRed$Steps,each=4)
Function=rep(module_hierarchy_RowRed$Module,each=4)

data_modelling=data.frame(Function=Function,Domain=Domain,Steps=Steps,Phylum=Phylum,Slope=Slope_data$value)

# "Xenobiotics biodegradation" domain contains 1 function.
# "Biosynthesis of other secondary metabolites" domain contains 2 functions
# "Glycan metabolism" contains 4 functions
# The three of the will be discarded from the last analysis
data_modelling_red=data_modelling[data_modelling$Domain!="Xenobiotics biodegradation"&
                                    data_modelling$Domain!="Biosynthesis of other secondary metabolites"&
                                    data_modelling$Domain!="Glycan metabolism",]


# Fit the model and make predictions and inferences
# *************************************************

M=lmer(Slope~Phylum+Domain+Steps+(1|Function),data = data_modelling_red)
set.seed(1)
confint.merMod(M,method = "boot")

# Predictions of the slopes for different Phyla with bootstrap 95% CI

newdata_phylum=data.frame(Phylum=unique(data_modelling_red$Phylum),
                          Domain=rep("Amino acid metabolism",4),
                          Steps=rep(mean(data_modelling_red$Steps),4))
predFun_Phylum <- function(fit) {
  predict(fit,newdata_phylum,re.form=NA)
}

boot_CI_Phylum=bootMer(M,nsim = 999,seed = 1,FUN = predFun_Phylum)
mean_Actinobacteriota=mean(boot_CI_Phylum$t[,1])
upr_Actinobacteriota=quantile(boot_CI_Phylum$t[,1],probs = 0.975)
lwr_Actinobacteriota=quantile(boot_CI_Phylum$t[,1],probs = 0.025)
mean_Bacteroidota=mean(boot_CI_Phylum$t[,2])
upr_Bacteroidota=quantile(boot_CI_Phylum$t[,2],probs = 0.975)
lwr_Bacteroidota=quantile(boot_CI_Phylum$t[,2],probs = 0.025)
mean_Firmicutes=mean(boot_CI_Phylum$t[,3])
upr_Firmicutes=quantile(boot_CI_Phylum$t[,3],probs = 0.975)
lwr_Firmicutes=quantile(boot_CI_Phylum$t[,3],probs = 0.025)
mean_Proteobacteria=mean(boot_CI_Phylum$t[,4])
upr_Proteobacteria=quantile(boot_CI_Phylum$t[,4],probs = 0.975)
lwr_Proteobacteria=quantile(boot_CI_Phylum$t[,4],probs = 0.025)
boot_data_phylum=data.frame(mean=c(mean_Actinobacteriota,mean_Bacteroidota,
                                   mean_Firmicutes,mean_Proteobacteria),
                            upr=c(upr_Actinobacteriota,upr_Bacteroidota,
                                  upr_Firmicutes,upr_Proteobacteria),
                            lwr=c(lwr_Actinobacteriota,lwr_Bacteroidota,
                                  lwr_Firmicutes,lwr_Proteobacteria),
                            phylum=unique(data_modelling_red$Phylum))
boot_data_phylum$CI_phylum_labels=c("[0.044,0.062]","[0.041,0.058]","[0.031,0.05]","[0.057,0.075]")

# Figure 2B in the main manuscript
ggplot()+
  geom_point(data=data_modelling_red,mapping=aes(x=Phylum,y=Slope,color=Phylum),
             position = position_jitter(w = 0.2, h = 0),alpha=0.3,show.legend = FALSE)+
  scale_color_manual(values=phylum_palette)+
  geom_point(data = boot_data_phylum,mapping=aes(x=phylum,y=mean),size=2)+
  geom_errorbar(data = boot_data_phylum,mapping=aes(ymin=lwr, ymax=upr,x=phylum),size=1, width=.5)+
  geom_text(data=boot_data_phylum,mapping=aes(x=phylum,y=0.2,label=CI_phylum_labels),size=5, position=position_dodge(width=0.9))+
  ylim(c(0,0.25))+
  ggtitle("b)")+
  ylab(expression(Slope[Fullness/Completeness]))+
  theme_classic()+
  coord_flip()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=12,hjust=1),
        axis.text.y = element_text(size=12),
        title = element_text(size=18))

# Predictions of the slopes for different Domains with bootstrap 95% CI

newdata_domain=data.frame(Domain=unique(data_modelling_red$Domain),
                          Phylum=rep("Actinobacteriota",7),
                          Steps=rep(mean(data_modelling_red$Steps),7))
predFun_Domain <- function(fit) {
  predict(fit,newdata_domain,re.form=NA)
}
boot_CI_Domain=bootMer(M,nsim = 999,seed = 1,FUN = predFun_Domain)
mean_carbohydrate_metabolism=mean(boot_CI_Domain$t[,1])
upr_carbohydrate_metabolism=quantile(boot_CI_Domain$t[,1],probs = 0.975)
lwr_carbohydrate_metabolism=quantile(boot_CI_Domain$t[,1],probs = 0.025)
mean_amino_acid_metabolism=mean(boot_CI_Domain$t[,2])
upr_amino_acid_metabolism=quantile(boot_CI_Domain$t[,2],probs = 0.975)
lwr_amino_acid_metabolism=quantile(boot_CI_Domain$t[,2],probs = 0.025)
mean_nucleotide_metabolism=mean(boot_CI_Domain$t[,3])
upr_nucleotide_metabolism=quantile(boot_CI_Domain$t[,3],probs = 0.975)
lwr_nucleotide_metabolism=quantile(boot_CI_Domain$t[,3],probs = 0.025)
mean_lipid_metabolism=mean(boot_CI_Domain$t[,4])
upr_lipid_metabolism=quantile(boot_CI_Domain$t[,4],probs = 0.975)
lwr_lipid_metabolism=quantile(boot_CI_Domain$t[,4],probs = 0.025)
mean_biosynthesis_terpenoids=mean(boot_CI_Domain$t[,5])
upr_biosynthesis_terpenoids=quantile(boot_CI_Domain$t[,5],probs = 0.975)
lwr_biosynthesis_terpenoids=quantile(boot_CI_Domain$t[,5],probs = 0.025)
mean_metabolism_cofactors=mean(boot_CI_Domain$t[,6])
upr_metabolism_cofactors=quantile(boot_CI_Domain$t[,6],probs = 0.975)
lwr_metabolism_cofactors=quantile(boot_CI_Domain$t[,6],probs = 0.025)
mean_energy_metabolism=mean(boot_CI_Domain$t[,7])
upr_energy_metabolism=quantile(boot_CI_Domain$t[,7],probs = 0.975)
lwr_energy_metabolism=quantile(boot_CI_Domain$t[,7],probs = 0.025)

boot_data_domain=data.frame(mean=c(mean_carbohydrate_metabolism,mean_amino_acid_metabolism,
                                   mean_nucleotide_metabolism,mean_lipid_metabolism,
                                   mean_biosynthesis_terpenoids,mean_metabolism_cofactors,
                                   mean_energy_metabolism),
                            upr=c(upr_carbohydrate_metabolism,upr_amino_acid_metabolism,
                                  upr_nucleotide_metabolism,upr_lipid_metabolism,
                                  upr_biosynthesis_terpenoids,upr_metabolism_cofactors,
                                  upr_energy_metabolism),
                            lwr=c(lwr_carbohydrate_metabolism,lwr_amino_acid_metabolism,
                                  lwr_nucleotide_metabolism,lwr_lipid_metabolism,
                                  lwr_biosynthesis_terpenoids,lwr_metabolism_cofactors,
                                  lwr_energy_metabolism),
                            domain=unique(data_modelling_red$Domain))

boot_data_domain$CI_domain_labels=c("[0.049,0.067]","[0.044,0.062]","[0.067,0.1]",
                                    "[0.03,0.065]","[0.015,0.047]","[0.042,0.06]","[0.029,0.053]")

ggplot()+
  geom_point(data=data_modelling_red,mapping=aes(x=Domain,y=Slope,color=Domain),
             position = position_jitter(w = 0.2, h = 0),alpha=0.3,show.legend = FALSE)+
  scale_color_manual(values=domain_palette)+
  geom_point(data = boot_data_domain,mapping=aes(x=domain,y=mean),size=2)+
  geom_errorbar(data = boot_data_domain,mapping=aes(ymin=lwr, ymax=upr,x=domain), size=1,width=.5)+
  geom_text(data=boot_data_domain,mapping=aes(x=domain,y=0.2,label=CI_domain_labels),size=4, position=position_dodge(width=0.9))+
  ylim(c(0,0.25))+
  ggtitle("c)")+
  ylab(expression(Slope[Fullness/Completeness]))+
  coord_flip()+
  theme_classic()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=8,hjust=1),
        axis.text.y = element_text(size=12),
        title = element_text(size=18))

# Predicted association between the slopes and number of steps in modules with bootstrap 95% CI

newdata_steps=data.frame(Steps=seq(from=min(data_modelling_red$Steps),
                                   to=max(data_modelling_red$Steps),
                                   length=100),
                         Domain=rep("Amino acid metabolism",100),
                         Phylum=rep("Actinobacteriota",100))
predFun_Steps <- function(fit) {
  predict(fit,newdata_steps,re.form=NA)
}
boot_CI_steps=bootMer(M,nsim = 999,seed = 1,FUN = predFun_Steps)
boot_data_steps=data.frame(mean=apply(boot_CI_steps$t,2,mean),
                           upr=apply(boot_CI_steps$t,2,quantile,0.975),
                           lwr=apply(boot_CI_steps$t,2,quantile,0.025),
                           steps=newdata_steps$Steps)


ggplot()+
  geom_point(data=data_modelling_red,mapping=aes(x=Steps,y=Slope),
             alpha=0.3,show.legend = FALSE)+
  geom_line(data = boot_data_steps,mapping=aes(x=steps,y=mean),size=1,color="blue")+
  geom_ribbon(data = boot_data_steps,mapping=aes(ymin=lwr, ymax=upr,x=steps),alpha=0.3,fill="blue")+
  geom_text(aes(x=10,y=0.2,label=as.character(expression(beta~"="~"[-0.004,-0.001]"))),parse=TRUE,color="black",size=5)+
  ggtitle("d)")+
  ylab(expression(Slope[Fullness/Completeness]))+
  xlab("Number of steps")+
  theme_classic()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=12,hjust=1),
        axis.text.y = element_text(size=12),
        title = element_text(size=18))