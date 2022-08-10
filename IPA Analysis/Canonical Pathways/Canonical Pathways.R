library(readr)
library(ggplot2)

setwd("C:/Users/codyr/OneDrive/Lab File/GenX/RNA Gene Array/Newest Analysis/IPA Analysis/Canonical Pathways")

CP0.1um=read_delim("0.1um_Canonical.txt", 
                   "\t", escape_double = FALSE, trim_ws = TRUE, 
                   skip = 1)
CP10um=read_delim("10um_Canonical.txt", 
                   "\t", escape_double = FALSE, trim_ws = TRUE, 
                   skip = 1)
CP100um=read_delim("100um_Canonical.txt", 
                   "\t", escape_double = FALSE, trim_ws = TRUE, 
                   skip = 1)

for (i in 1:length(CP0.1um$Molecules )){
  CP0.1um$MoleculeCount[i]=length(unlist(strsplit(CP0.1um$Molecules[i],",")))
}

for (i in 1:length(CP10um$Molecules )){
  CP10um$MoleculeCount[i]=length(unlist(strsplit(CP10um$Molecules[i],",")))
}

for (i in 1:length(CP100um$Molecules )){
  CP100um$MoleculeCount[i]=length(unlist(strsplit(CP100um$Molecules[i],",")))
}

CP0.1um$`p-value`= 2^-CP0.1um$`-log(p-value)`
CP10um$`p-value`= 2^-CP10um$`-log(p-value)`
CP100um$`p-value`= 2^-CP100um$`-log(p-value)`

CP0.1um=CP0.1um[order(CP0.1um$`p-value`),]
CP10um=CP10um[order(CP10um$`p-value`),]
CP100um=CP100um[order(CP100um$`p-value`),]

