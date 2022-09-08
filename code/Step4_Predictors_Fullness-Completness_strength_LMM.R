library(tidyverse)
library(lme4)

# Load datasets
# *************
Results_table <- read.csv("data/Results_table.csv")
n_steps <- read.csv("data/number_of_steps_in_modules.csv",header = TRUE)
module_hierarchy <- read.csv("data/module_hierarchy.csv")

# Prepare data to fit the model
# *****************************
Slope_data <- data.frame(Results_table$Actinobacteriota_slope,Results_table$Bacteroidota_slope,
                      Results_table$Firmicutes_slope,Results_table$Proteobacteria_slope) %>%
  pivot_longer(everything())
Phylum <- rep(c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria"),nrow(Results_table))
Domain <- rep(module_hierarchy$Domain,each=4)
Steps <- rep(n_steps$Steps,each=4)
Function <- rep(module_hierarchy$Module,each=4)

data_modelling <- data.frame(Function=Function,Domain=Domain,Steps=Steps,Phylum=Phylum,Slope=Slope_data$value)

# Fit the model and make predictions and inferences
# *************************************************
set.seed(1)

domain_palette <- c('Carbohydrate metabolism'="#DB93A9",'Amino acid metabolism'="#E89D74",
                 'Biosynthesis of other secondary metabolites'="#69AD86",
                 'Nucleotide metabolism'="#E9CA80",'Glycan metabolism'="#D7BFAF",
                 'Lipid metabolism'="#99796D",'Biosynthesis of terpenoids and polyketides'="#2D758C",
                 'Metabolism of cofactors and vitamins'="#A3DBDD",'Energy metabolism'="#ABA5D2",
                 'Xenobiotics biodegradation'="#8E8E8E")
phylum_palette=c('Actinobacteriota'="#F08E49",'Bacteroidota'="#7ACCC2",
                 'Firmicutes'="#DACC68",'Proteobacteria'="#E6695D")

M <- lmer(Slope~Phylum+Domain+Steps+(1|Function),data = data_modelling)
set.seed(1)
confint.merMod(M,method = "boot")

# Predictions of the slopes for different Phyla with bootstrap 95% CI

newdata_phylum <- data.frame(Phylum=unique(data_modelling$Phylum),
                          Domain=rep("Amino acid metabolism",4),
                          Steps=rep(mean(data_modelling$Steps),4))
predFun_Phylum <- function(fit) {
  predict(fit,newdata_phylum,re.form=NA)
}

boot_CI_Phylum <- bootMer(M,nsim = 999,seed = 1,FUN = predFun_Phylum)
mean_Actinobacteriota <- mean(boot_CI_Phylum$t[,1])
upr_Actinobacteriota <- quantile(boot_CI_Phylum$t[,1],probs = 0.975)
lwr_Actinobacteriota <- quantile(boot_CI_Phylum$t[,1],probs = 0.025)
mean_Bacteroidota <- mean(boot_CI_Phylum$t[,2])
upr_Bacteroidota <- quantile(boot_CI_Phylum$t[,2],probs = 0.975)
lwr_Bacteroidota <- quantile(boot_CI_Phylum$t[,2],probs = 0.025)
mean_Firmicutes <- mean(boot_CI_Phylum$t[,3])
upr_Firmicutes <- quantile(boot_CI_Phylum$t[,3],probs = 0.975)
lwr_Firmicutes <- quantile(boot_CI_Phylum$t[,3],probs = 0.025)
mean_Proteobacteria <- mean(boot_CI_Phylum$t[,4])
upr_Proteobacteria <- quantile(boot_CI_Phylum$t[,4],probs = 0.975)
lwr_Proteobacteria <- quantile(boot_CI_Phylum$t[,4],probs = 0.025)
boot_data_phylum <- data.frame(mean=c(mean_Actinobacteriota,mean_Bacteroidota,
                                   mean_Firmicutes,mean_Proteobacteria),
                            upr=c(upr_Actinobacteriota,upr_Bacteroidota,
                                  upr_Firmicutes,upr_Proteobacteria),
                            lwr=c(lwr_Actinobacteriota,lwr_Bacteroidota,
                                  lwr_Firmicutes,lwr_Proteobacteria),
                            phylum=unique(data_modelling$Phylum))
boot_data_phylum$CI_phylum_labels=c("[0.032,0.045]","[0.02,0.033]","[0.033,0.046]","[0.035,0.048]")
boot_data_phylum$phylum<-factor(boot_data_phylum$phylum,
                                levels = boot_data_phylum$phylum[order(boot_data_phylum$mean)])
data_modelling$Phylum<-factor(data_modelling$Phylum,
                              levels = boot_data_phylum$phylum[order(boot_data_phylum$mean)])

# Figure 2B in the main manuscript

ggplot()+
  geom_point(data=data_modelling,mapping=aes(x=Phylum,y=Slope,color=Phylum),
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

ggsave("figures/Fig_2B.pdf",width = 8,height = 8)

# Predictions of the slopes for different Domains with bootstrap 95% CI

newdata_domain <- data.frame(Domain=unique(data_modelling$Domain),
                          Phylum=rep("Actinobacteriota",10),
                          Steps=rep(mean(data_modelling$Steps),10))
predFun_Domain <- function(fit) {
  predict(fit,newdata_domain,re.form=NA)
}
boot_CI_Domain <- bootMer(M,nsim = 999,seed = 1,FUN = predFun_Domain)
mean_carbohydrate_metabolism <- mean(boot_CI_Domain$t[,1])
upr_carbohydrate_metabolism <- quantile(boot_CI_Domain$t[,1],probs = 0.975)
lwr_carbohydrate_metabolism <- quantile(boot_CI_Domain$t[,1],probs = 0.025)
mean_amino_acid_metabolism <- mean(boot_CI_Domain$t[,2])
upr_amino_acid_metabolism <- quantile(boot_CI_Domain$t[,2],probs = 0.975)
lwr_amino_acid_metabolism <- quantile(boot_CI_Domain$t[,2],probs = 0.025)
mean_biosynthesis_secondary_metabolites <- mean(boot_CI_Domain$t[,3])
upr_biosynthesis_secondary_metabolites <- quantile(boot_CI_Domain$t[,3],probs = 0.975)
lwr_biosynthesis_secondary_metabolites <- quantile(boot_CI_Domain$t[,3],probs = 0.025)
mean_nucleotide_metabolism <- mean(boot_CI_Domain$t[,4])
upr_nucleotide_metabolism <- quantile(boot_CI_Domain$t[,4],probs = 0.975)
lwr_nucleotide_metabolism <- quantile(boot_CI_Domain$t[,4],probs = 0.025)
mean_glycan_metabolism <- mean(boot_CI_Domain$t[,5])
upr_glycan_metabolism <- quantile(boot_CI_Domain$t[,5],probs = 0.975)
lwr_glycan_metabolism <- quantile(boot_CI_Domain$t[,5],probs = 0.025)
mean_lipid_metabolism <- mean(boot_CI_Domain$t[,6])
upr_lipid_metabolism <- quantile(boot_CI_Domain$t[,6],probs = 0.975)
lwr_lipid_metabolism <- quantile(boot_CI_Domain$t[,6],probs = 0.025)
mean_biosynthesis_terpenoids <- mean(boot_CI_Domain$t[,7])
upr_biosynthesis_terpenoids <- quantile(boot_CI_Domain$t[,7],probs = 0.975)
lwr_biosynthesis_terpenoids <- quantile(boot_CI_Domain$t[,7],probs = 0.025)
mean_metabolism_cofactors <- mean(boot_CI_Domain$t[,8])
upr_metabolism_cofactors <- quantile(boot_CI_Domain$t[,8],probs = 0.975)
lwr_metabolism_cofactors <- quantile(boot_CI_Domain$t[,8],probs = 0.025)
mean_energy_metabolism <- mean(boot_CI_Domain$t[,9])
upr_energy_metabolism <- quantile(boot_CI_Domain$t[,9],probs = 0.975)
lwr_energy_metabolism <- quantile(boot_CI_Domain$t[,9],probs = 0.025)
mean_xenobiotics_biodegradation <- mean(boot_CI_Domain$t[,10])
upr_xenobiotics_biodegradation <- quantile(boot_CI_Domain$t[,10],probs = 0.975)
lwr_xenobiotics_biodegradation <- quantile(boot_CI_Domain$t[,10],probs = 0.025)
boot_data_domain <- data.frame(mean=c(mean_carbohydrate_metabolism,mean_amino_acid_metabolism,
                                   mean_biosynthesis_secondary_metabolites,mean_nucleotide_metabolism,
                                   mean_glycan_metabolism,mean_lipid_metabolism,
                                   mean_biosynthesis_terpenoids,mean_metabolism_cofactors,
                                   mean_energy_metabolism,mean_xenobiotics_biodegradation),
                            upr=c(upr_carbohydrate_metabolism,upr_amino_acid_metabolism,
                                  upr_biosynthesis_secondary_metabolites,upr_nucleotide_metabolism,
                                  upr_glycan_metabolism,upr_lipid_metabolism,
                                  upr_biosynthesis_terpenoids,upr_metabolism_cofactors,
                                  upr_energy_metabolism,upr_xenobiotics_biodegradation),
                            lwr=c(lwr_carbohydrate_metabolism,lwr_amino_acid_metabolism,
                                  lwr_biosynthesis_secondary_metabolites,lwr_nucleotide_metabolism,
                                  lwr_glycan_metabolism,lwr_lipid_metabolism,
                                  lwr_biosynthesis_terpenoids,lwr_metabolism_cofactors,
                                  lwr_energy_metabolism,lwr_xenobiotics_biodegradation),
                            domain=unique(data_modelling$Domain))

boot_data_domain$CI_domain_labels=c("[0.037,0.05]","[0.032,0.045]","[0.03,0.06]",
                                    "[0.04,0.068]","[0.022,0.044]","[0.024,0.044]",
                                    "[0.021,0.041]","[0.031,0.047]","[0.02,0.034]",
                                    "[0.019,0.056]")
boot_data_domain$domain<-factor(boot_data_domain$domain,
                                levels = boot_data_domain$domain[order(boot_data_domain$mean)])
data_modelling$Domain<-factor(data_modelling$Domain,
                              levels = boot_data_domain$domain[order(boot_data_domain$mean)])

# Figure 2C in the main manuscript
ggplot()+
  geom_point(data=data_modelling,mapping=aes(x=Domain,y=Slope,color=Domain),
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

ggsave("figures/Fig_2C.pdf",width = 8,height = 8)

# Predicted association between the slopes and number of steps in modules with bootstrap 95% CI

newdata_steps <- data.frame(Steps=seq(from=min(data_modelling$Steps),
                                   to=max(data_modelling$Steps),
                                   length=100),
                         Domain=rep("Amino acid metabolism",100),
                         Phylum=rep("Actinobacteriota",100))
predFun_Steps <- function(fit) {
  predict(fit,newdata_steps,re.form=NA)
}
boot_CI_steps <- bootMer(M,nsim = 999,seed = 1,FUN = predFun_Steps)
boot_data_steps <- data.frame(mean=apply(boot_CI_steps$t,2,mean),
                           upr=apply(boot_CI_steps$t,2,quantile,0.975),
                           lwr=apply(boot_CI_steps$t,2,quantile,0.025),
                           steps=newdata_steps$Steps)


# Figure 2D in the main manuscript
ggplot()+
  geom_point(data=data_modelling,mapping=aes(x=Steps,y=Slope),
             alpha=0.3,show.legend = FALSE)+
  geom_line(data = boot_data_steps,mapping=aes(x=steps,y=mean),size=1,color="blue")+
  geom_ribbon(data = boot_data_steps,mapping=aes(ymin=lwr, ymax=upr,x=steps),alpha=0.3,fill="blue")+
  geom_text(aes(x=10,y=0.2,label=as.character(expression(beta~"="~"[-0.003,-0.001]"))),parse=TRUE,color="black",size=5)+
  ggtitle("d)")+
  ylab(expression(Slope[Fullness/Completeness]))+
  xlab("Number of steps")+
  theme_classic()+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=12,hjust=1),
        axis.text.y = element_text(size=12),
        title = element_text(size=18))

ggsave("figures/Fig_2D.pdf",width = 8,height = 8)
