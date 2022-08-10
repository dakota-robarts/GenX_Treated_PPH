library(heatmap.plus)
library(gplots)
library(readxl)
library(RColorBrewer)
library(readr)
library(xlsx)

setwd("D:/OneDrive/Lab File/GenX/RNA Gene Array/Newest Analysis/Dose Response/IPA/Upstream Regs")

Upstream_100uM <- read_excel("Upstream_100uM.xls", 
                             skip = 1)
Upstream_100uM=Upstream_100uM[order(Upstream_100uM$`Activation z-score`),]

Upstream_10uM <- read_excel("Upstream_10uM.xls", 
                            skip = 1)
Upstream_10uM=Upstream_10uM[order(Upstream_10uM$`Activation z-score`),]

Upstream_0.1uM <- read_excel("Upstream_0.1uM.xls", 
                             skip = 1)
Upstream_0.1uM=Upstream_0.1uM[order(Upstream_0.1uM$`Activation z-score`),]

Upstream_100uM=Upstream_100uM[!is.na(Upstream_100uM$`Activation z-score`),]
Upstream_10uM=Upstream_10uM[!is.na(Upstream_10uM$`Activation z-score`),]
Upstream_0.1uM=Upstream_0.1uM[!is.na(Upstream_0.1uM$`Activation z-score`),]

sharedpathways=intersect(Upstream_100uM$`Upstream Regulator`,
                         intersect(Upstream_10uM$`Upstream Regulator`,
                                   Upstream_0.1uM$`Upstream Regulator`))
which(duplicated(sharedpathways))#no overlapping pathways


-log(0.05,10)

Pathways=as.data.frame(sharedpathways)

#Adding z score for 0.1 uM
for (i in 1: length(Pathways$sharedpathways)){
  Pathways$`0.1uM`[i]=Upstream_0.1uM$`Activation z-score`[which(Upstream_0.1uM$`Upstream Regulator` %in%
                                                                  Pathways$sharedpathways[i])]
}

#Adding z score for 10 uM
for (i in 1: length(Pathways$sharedpathways)){
  Pathways$`10uM`[i]=Upstream_10uM$`Activation z-score`[which(Upstream_10uM$`Upstream Regulator` %in%
                                                                Pathways$sharedpathways[i])]
}


#Adding z score for 100 uM
for (i in 1: length(Pathways$sharedpathways)){
  Pathways$`100uM`[i]=Upstream_100uM$`Activation z-score`[which(Upstream_100uM$`Upstream Regulator` %in%
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


tiff("Correlated(.99)Upstream.tiff", width = 10, height = 5, units = 'in', res = 600, bg = "transparent")
heatmap.2(as.matrix(dose.dependent[, 2:4]), trace="none", Colv = NA,
          col=rb,
          dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
                                         expression(paste("100 ", mu,"M"))), labRow = dose.dependent$sharedpathways[],
          distfun = function(x) dist(x,"manhattan"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          cexRow=.8,
          cexCol = 1.5,
          margins=c(5.1,35),
          lhei=c(.5,1.75),
          lwid=c(.2,1),
          key.par = list(cex=.5),
          key.xlab="z-score",
          main="Positive Correlated\nUpstream Regulators",
          key.title="Color key",
          density.info = "none")

dev.off()

#Only -0.99 correlation============================================================
dose.dependent=Pathways[which(Pathways$Corr < -0.9),]
dose.dependent=dose.dependent[order(dose.dependent$`100uM`,decreasing = F),]
# write.xlsx(dose.dependent,file="Correlated(-0.99)AllGenes.xlsx")



reds=brewer.pal(6,"Oranges")[-1]
blues=rev(brewer.pal(6,"Blues"))[-6]
rb=colorRampPalette(c(blues,reds))(n=100)
rb


tiff("Correlated(-.99)Upstream.tiff", width = 10, height = 3, units = 'in', res = 600, bg = "transparent")
heatmap.2(as.matrix(dose.dependent[, 2:4]), trace="none", Colv = NA,
          col=rb,
          dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
                                         expression(paste("100 ", mu,"M"))), labRow = dose.dependent$sharedpathways,
          distfun = function(x) dist(x,"manhattan"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          cexRow=1,
          cexCol = 1.5,
          margins=c(5.1,35),
          lhei=c(.5,.9),
          lwid=c(.2,1),
          key.par = list(cex=.5),
          key.xlab="z-score",
          main="Negative Correlated\nUpstream Regulators",
          key.title="Color key",
          density.info = "none")

dev.off()


