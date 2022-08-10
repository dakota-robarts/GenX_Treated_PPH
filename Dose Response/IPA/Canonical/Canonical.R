library(heatmap.plus)
library(gplots)
library(readxl)
library(RColorBrewer)
library(readr)
library(xlsx)

setwd("D:/OneDrive/Lab File/GenX/RNA Gene Array/Newest Analysis/Dose Response/IPA/Canonical")

Canonical_100uM <- read_excel("Canonical_100uM.xls", 
                              skip = 1)
Canonical_100uM=Canonical_100uM[order(Canonical_100uM$`z-score`),]

Canonical_10uM <- read_excel("Canonical_10uM.xls", 
                              skip = 1)
Canonical_10uM=Canonical_10uM[order(Canonical_10uM$`z-score`),]

Canonical_0.1uM <- read_excel("Canonical_0.1uM.xls", 
                             skip = 1)
Canonical_0.1uM=Canonical_0.1uM[order(Canonical_0.1uM$`z-score`),]

Canonical_100uM=Canonical_100uM[!is.na(Canonical_100uM$`z-score`),]
Canonical_10uM=Canonical_10uM[!is.na(Canonical_10uM$`z-score`),]
Canonical_0.1uM=Canonical_0.1uM[!is.na(Canonical_0.1uM$`z-score`),]

sharedpathways=intersect(Canonical_100uM$`Ingenuity Canonical Pathways`,
                         intersect(Canonical_10uM$`Ingenuity Canonical Pathways`,
                                   Canonical_0.1uM$`Ingenuity Canonical Pathways`))
which(duplicated(sharedpathways))#no overlapping pathways


-log(0.05,10)

Pathways=as.data.frame(sharedpathways)

#Adding z score for 0.1 uM
for (i in 1: length(Pathways$sharedpathways)){
  Pathways$`0.1uM`[i]=Canonical_0.1uM$`z-score`[which(Canonical_0.1uM$`Ingenuity Canonical Pathways` %in%
                              Pathways$sharedpathways[i])]
}

#Adding z score for 10 uM
for (i in 1: length(Pathways$sharedpathways)){
  Pathways$`10uM`[i]=Canonical_10uM$`z-score`[which(Canonical_10uM$`Ingenuity Canonical Pathways` %in%
                                                        Pathways$sharedpathways[i])]
}


#Adding z score for 100 uM
for (i in 1: length(Pathways$sharedpathways)){
  Pathways$`100uM`[i]=Canonical_100uM$`z-score`[which(Canonical_100uM$`Ingenuity Canonical Pathways` %in%
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


# tiff("Correlated(.99)ConicalPath.tiff", width = 10, height = 5, units = 'in', res = 600, bg = "transparent")
# heatmap.2(as.matrix(dose.dependent[1:20, 2:4]), trace="none", Colv = NA,
#           col=rb,
#           dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
#                                          expression(paste("100 ", mu,"M"))), labRow = dose.dependent$sharedpathways[1:20],
#           distfun = function(x) dist(x,"manhattan"),
#           hclustfun = function(x) hclust(x, method="ward.D2"),
#           cexRow=.8,
#           cexCol = 1.5,
#           margins=c(5.1,35),
#           lhei=c(.5,1.75),
#           lwid=c(.2,1),
#           key.par = list(cex=.5),
#           key.xlab="z-score",
#           main="Positive Correlated\nCanonical Pathways",
#           key.title="Color key",
#           density.info = "none")
# 
# dev.off()

#Only -0.99 correlation============================================================
dose.dependent=Pathways[which(Pathways$Corr < -0.9),]
dose.dependent=dose.dependent[order(dose.dependent$`100uM`,decreasing = F),]
# write.xlsx(dose.dependent,file="Correlated(-0.99)AllGenes.xlsx")



reds=brewer.pal(6,"Oranges")[-1]
blues=rev(brewer.pal(6,"Blues"))[-6]
rb=colorRampPalette(c(blues,reds))(n=100)
rb

# 
tiff("Correlated(-.99)ConicalPath.tiff", width = 10, height = 4, units = 'in', res = 600, bg = "transparent")
heatmap.2(as.matrix(dose.dependent[-5, 2:4]), trace="none", Colv = NA,
          col=rb,
          dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
                                         expression(paste("100 ", mu,"M"))), labRow = dose.dependent$sharedpathways[-5],
          distfun = function(x) dist(x,"manhattan"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          cexRow=1.25,
          cexCol = 1.5,
          margins=c(5.1,35),
          lhei=c(.5,1.25),
          lwid=c(.2,1),
          key.par = list(cex=.5),
          key.xlab="z-score",
          main="Negative Correlated\nCanonical Pathways",
          key.title="Color key",
          density.info = "none")

dev.off()


