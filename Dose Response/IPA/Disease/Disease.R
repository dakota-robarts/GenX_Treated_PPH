library(heatmap.plus)
library(gplots)
library(readxl)
library(RColorBrewer)
library(readr)
library(xlsx)

setwd(".../IPA/Diseases")

Diseases_100uM <- read_excel("Disease_100uM.xls", 
                              skip = 1)
Diseases_100uM=Diseases_100uM[order(Diseases_100uM$`Activation z-score`),]

Diseases_10uM <- read_excel("Disease_10uM.xls", 
                             skip = 1)
Diseases_10uM=Diseases_10uM[order(Diseases_10uM$`Activation z-score`),]

Diseases_0.1uM <- read_excel("Disease_0.1uM.xls", 
                             skip = 1)
Diseases_0.1uM=Diseases_0.1uM[order(Diseases_0.1uM$`Activation z-score`),]

Diseases_100uM=Diseases_100uM[!is.na(Diseases_100uM$`Activation z-score`),]
Diseases_10uM=Diseases_10uM[!is.na(Diseases_10uM$`Activation z-score`),]
Diseases_0.1uM=Diseases_0.1uM[!is.na(Diseases_0.1uM$`Activation z-score`),]

sharedpathways=intersect(Diseases_100uM$`Diseases or Functions Annotation`,
                         intersect(Diseases_10uM$`Diseases or Functions Annotation`,
                                   Diseases_0.1uM$`Diseases or Functions Annotation`))
which(duplicated(sharedpathways))#no overlapping pathways


-log(0.05,10)

Pathways=as.data.frame(sharedpathways)

#Adding z score for 0.1 uM
for (i in 1: length(Pathways$sharedpathways)){
  Pathways$`0.1uM`[i]=Diseases_0.1uM$`Activation z-score`[which(Diseases_0.1uM$`Diseases or Functions Annotation` %in%
                                                        Pathways$sharedpathways[i])]
}

#Adding z score for 10 uM
for (i in 1: length(Pathways$sharedpathways)){
  Pathways$`10uM`[i]=Diseases_10uM$`Activation z-score`[which(Diseases_10uM$`Diseases or Functions Annotation` %in%
                                                      Pathways$sharedpathways[i])]
}


#Adding z score for 100 uM
for (i in 1: length(Pathways$sharedpathways)){
  Pathways$`100uM`[i]=Diseases_100uM$`Activation z-score`[which(Diseases_100uM$`Diseases or Functions Annotation` %in%
                                                        Pathways$sharedpathways[i])]
}

numbers=c(0.1,10,100)
cor(numbers,as.numeric(Pathways[2,2:4]))


for( i in 1:length(Pathways$sharedpathways)){
  Pathways$Corr[i]=cor(numbers,as.numeric(Pathways[i,2:4]))
}

#Only 0.99 correlation============================================================
dose.dependent=Pathways[which(Pathways$Corr > 0.9),]
dose.dependent=dose.dependent[order(dose.dependent$`100uM`,decreasing = T),]
# write.xlsx(dose.dependent,file="Correlated(-0.99)AllGenes.xlsx")



reds=brewer.pal(6,"Oranges")[-1]
blues=rev(brewer.pal(6,"Blues"))[-6]
rb=colorRampPalette(c(blues,reds))(n=100)
rb


# tiff("Correlated(.99)Disease.tiff", width = 10, height = 5, units = 'in', res = 600, bg = "transparent")
heatmap.2(as.matrix(dose.dependent[c(1:14,34:40), 2:4]), trace="none", Colv = NA,
          col=rb,
          dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
                                         expression(paste("100 ", mu,"M"))), labRow = dose.dependent$sharedpathways[c(1:14,34:40)],
          distfun = function(x) dist(x,"manhattan"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          cexRow=.8,
          cexCol = 1.5,
          margins=c(5.1,35),
          lhei=c(.5,1.75),
          lwid=c(.2,1),
          key.par = list(cex=.5),
          key.xlab="z-score",
          main="Positive Correlated\nPhenotypes",
          key.title="Color key",
          density.info = "none")

# dev.off()

#Only -0.99 correlation============================================================
dose.dependent=Pathways[which(Pathways$Corr < -0.9),]
dose.dependent=dose.dependent[order(dose.dependent$`100uM`,decreasing = F),]
# write.xlsx(dose.dependent,file="Correlated(-0.99)AllGenes.xlsx")



reds=brewer.pal(6,"Oranges")[-1]
blues=rev(brewer.pal(6,"Blues"))[-6]
rb=colorRampPalette(c(blues,reds))(n=100)
rb


# tiff("Correlated(-.99)Disease.tiff", width = 10, height = 3, units = 'in', res = 600, bg = "transparent")
heatmap.2(as.matrix(dose.dependent[-5, 2:4]), trace="none", Colv = NA,
          col=rb,
          dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
                                         expression(paste("100 ", mu,"M"))), labRow = dose.dependent$sharedpathways[-5],
          distfun = function(x) dist(x,"manhattan"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          cexRow=1.25,
          cexCol = 1.5,
          margins=c(5.1,35),
          lhei=c(.5,.9),
          lwid=c(.2,1),
          key.par = list(cex=.5),
          key.xlab="z-score",
          main="Negative Correlated\nPhenotypes",
          key.title="Color key",
          density.info = "none")

# dev.off()


