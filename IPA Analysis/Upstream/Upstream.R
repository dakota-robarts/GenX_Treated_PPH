library(readr)
library(ggplot2)
library(extrafont)
library(dplyr)

setwd(".../IPA Analysis/Upstream")

UR0.1=read_delim("0.1um_Upstream.txt", 
                 "\t", escape_double = FALSE, trim_ws = TRUE, 
                 skip = 1)
UR10=read_delim("10um_Upstream.txt", 
                "\t", escape_double = FALSE, trim_ws = TRUE, 
                skip = 1)
UR100=read_delim("100um_Upstream.txt", 
                 "\t", escape_double = FALSE, trim_ws = TRUE, 
                 skip = 1)


for (i in 1:length(UR0.1$`Target Molecules in Dataset`)){
UR0.1$MoleculeCount[i]=length(unlist(strsplit(UR0.1$`Target Molecules in Dataset`[i],",")))
}

for (i in 1:length(UR10$`Target Molecules in Dataset`)){
  UR10$MoleculeCount[i]=length(unlist(strsplit(UR10$`Target Molecules in Dataset`[i],",")))
}

for (i in 1:length(UR100$`Target Molecules in Dataset`)){
  UR100$MoleculeCount[i]=length(unlist(strsplit(UR100$`Target Molecules in Dataset`[i],",")))
}
UR100$`p-value of overlap`
UR100 <- subset(UR100,`Molecule Type` == "transcription regulator" | `Molecule Type` == "translation regulator" |
                 `Molecule Type` == "ligand-dependent nuclear receptor" |`Molecule Type` =="growth factor" |
                  `Molecule Type` =="group" | `Molecule Type` == "kinase"| `Molecule Type` == "cytokine") %>% subset(`p-value of overlap` < 0.05)
UR10 <- subset(UR10, `Molecule Type` == "transcription regulator" | `Molecule Type` == "translation regulator" |
               `Molecule Type` == "ligand-dependent nuclear receptor"|`Molecule Type` =="growth factor" |
                 `Molecule Type` =="group"| `Molecule Type` == "kinase"| `Molecule Type` == "cytokine") %>% subset(`p-value of overlap` < 0.05)
UR0.1 <- subset(UR0.1, `Molecule Type` == "transcription regulator" | `Molecule Type` == "translation regulator" |
                 `Molecule Type` == "ligand-dependent nuclear receptor"|`Molecule Type` =="growth factor" |
                  `Molecule Type` =="group"| `Molecule Type` == "kinase"| `Molecule Type` == "cytokine") %>% subset(`p-value of overlap` < 0.05)

UR0.1_1=UR0.1[!is.na(UR0.1$`Predicted Activation State`),] #Filtering only the Ones with predicted
UR10_1=UR10[!is.na(UR10$`Predicted Activation State`),] #Filtering only the Ones with predicted
UR100_1=UR100[!is.na(UR100$`Predicted Activation State`),] #Filtering only the Ones with predicted

UR10_1$`Upstream Regulator`[grepl("miR",UR10_1$`Upstream Regulator`)]= "miR-16-5p"
UR10_1=UR10_1[!grepl("trichostatin",UR10_1$`Upstream Regulator`),]
UR10_1=UR10_1[!grepl("cyclohe",UR10_1$`Upstream Regulator`),]

UR100_1$`Upstream Regulator`[grepl("miR",UR100_1$`Upstream Regulator`)]= "miR-291a-3p"


tiff("0.1um_UpstreamPathways.tiff", width = 6, height = 4, units = 'in', res = 300, bg = "transparent")
ggplot(data=UR0.1_1,aes(x=UR0.1_1$`Activation z-score`, 
                                   y=UR0.1_1$`Upstream Regulator`, 
                                   colour=UR0.1_1$`p-value of overlap`, 
                                   size=UR0.1_1$MoleculeCount)) +
  geom_point() +
  expand_limits(x=0) +
  geom_vline(xintercept=0,linetype="dashed",color="red")+
  labs(x="Activation Z-Score", y="Upstream Regulator", colour="P-Value", size=" Gene Count", title=expression(paste("0.1",mu,"M Upstream Regulator Pathways")))+
  theme(plot.title=element_text(hjust = 0.5,size=16)) +
  scale_colour_gradient(low="#CB6600",high = "#FFC891")
dev.off()

tiff("10um_UpstreamPathways.tiff", width = 6, height = 4, units = 'in', res = 300, bg = "transparent")
ggplot(data=UR10_1,aes(x=UR10_1$`Activation z-score`, 
                        y=UR10_1$`Upstream Regulator`, 
                        colour=UR10_1$`p-value of overlap`, 
                        size=UR10_1$MoleculeCount)) +
  geom_point() +
  expand_limits(x=0) +
  geom_vline(xintercept=0,linetype="dashed",color="red")+
  labs(x="Activation Z-Score", y="Upstream Regulator", colour="P-Value", size=" Gene Count", title=expression(paste("10",mu,"M Upstream Regulator Pathways")))+
  theme(plot.title=element_text(hjust = 0.5,size=16)) +
  scale_colour_gradient(low="#CB6600",high = "#FFC891")
dev.off()


tiff("100um_UpstreamPathways.tiff", width = 6, height = 4, units = 'in', res = 600, bg = "transparent")
ggplot(data=UR100_1,aes(x=UR100_1$`Activation z-score`, 
                       y=UR100_1$`Upstream Regulator`, 
                       colour=UR100_1$`p-value of overlap`, 
                       size=UR100_1$MoleculeCount)) +
  geom_point() +
  expand_limits(x=0) +
  geom_vline(xintercept=0,linetype="dashed",color="red")+
  labs(x="Activation Z-Score", y="Upstream Regulator", colour="P-Value", size=" Gene Count", 
       title=expression(paste("100",mu,"M Upstream Regulator Pathways")))+
  theme(plot.title=element_text(hjust = 0.5,size=16))+
  scale_colour_gradient(low="#CB6600",high = "#FFC891")
dev.off()

