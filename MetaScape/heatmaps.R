library(chorddiag)
library(heatmap.plus)
library(gplots)
library(readxl)
library(readr)
library(RColorBrewer)
dat1 <- read_excel(".../MetaScape/2 or -2 cut off.xlsx")
significant_genes <- read_csv(".../significant_genes.csv")

loc0.1=which(significant_genes$`0.1uMFold Change`< 0)
significant_genes$`0.1uMFold Change`[loc0.1]=-1/significant_genes$`0.1uMFold Change`[loc0.1]

loc10=which(significant_genes$`10uMFold Change`< 0)
significant_genes$`10uMFold Change`[loc10]=-1/significant_genes$`10uMFold Change`[loc10]

loc100=which(significant_genes$`GenX_100uM vs. GenX_UT`< 0)
significant_genes$`GenX_100uM vs. GenX_UT`[loc100]=-1/significant_genes$`GenX_100uM vs. GenX_UT`[loc100]

significant_genes$`0.1uMFold Change`=log(significant_genes$`0.1uMFold Change`,2)
significant_genes$`10uMFold Change`=log(significant_genes$`10uMFold Change`,2)
significant_genes$`GenX_100uM vs. GenX_UT`=log(significant_genes$`GenX_100uM vs. GenX_UT`,2)

Cytochrome=allgenes[which(grepl("CYP",allgenes))][-13] #also getting rid of AHCYP4

cyp.df.2=significant_genes[which(significant_genes$`Gene Symbol` %in% Cytochrome),c(3,4,5,2)]
cyp.df=as.matrix(cyp.df.2[,-1])
rownames(cyp.df)=cyp.df.2$`Gene Symbol`

significant_genes[1489,]

reds=brewer.pal(6,"Oranges")
blues=rev(brewer.pal(6,"Blues"))
rb=colorRampPalette(c(blues,reds))(n=100)
rb

#==========================================================CYP450 Panel

# tiff("CYP450.tiff", width = 4, height = 5, units = 'in', res = 300, bg = "transparent")
# heatmap.2(cyp.df, trace="n", Colv = NA,
#           col=rb,
#           dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
#                                         expression(paste("100 ", mu,"M"))), labRow = row.names(cyp.df),
#           distfun = function(x) dist(x,"euclidean"),
#           hclustfun = function(x) hclust(x, method="ward.D2"),
#           cexRow=.8,
#           cexCol = 1.5,
#           margins=c(5.1,4.5),
#           lhei=c(.5,2),
#           # lwid=c(.2,1),
#           key.par = list(cex=.5),
#           key.xlab="Log2(FC)",
#           main="Cytochrome\nP450 Panel",
#           key.title="Color key")
# dev.off()


#===================================================Top 20 Enriched Clusters===========

setwd(".../MetaScape")

HeatmapSelectedGO <- read_csv("Enrichment_heatmap/HeatmapSelectedGO.csv")
all <- read_csv("Enrichment_GO/_FINAL_GO.csv")
mitotic=all[which(grepl("Mitotic Prometaphase",all$Description)),c(7,18,4:6)]
HeatmapSelectedGO=HeatmapSelectedGO[-20,]
HeatmapSelectedGO=rbind(HeatmapSelectedGO,mitotic)

x=paste(HeatmapSelectedGO$GO,":",HeatmapSelectedGO$Description)
HeatmapSelectedGO=HeatmapSelectedGO[, -c(1:2)]
HeatmapSelectedGO=as.matrix(abs(HeatmapSelectedGO))
rownames(HeatmapSelectedGO)=x

rb=c("#C0C0C0",colorRampPalette(reds)(n=99))
numbers=c(0,seq(0.01,10,length=100))


tiff("Top20.tiff", width = 10, height = 6, units = 'in', res = 500, bg = "transparent")


heatmap.2(HeatmapSelectedGO, trace="n", Colv = NA,
          col=rb,
          dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
                                        expression(paste("100 ", mu,"M"))), labRow = row.names(HeatmapSelectedGO),
          distfun = function(x) dist(x,"manhattan"),
          hclustfun = function(x) hclust(x, method="ward.D"),
          cexRow=1.1,
          cexCol = 1.5,
          margins=c(5.1,39),
          lhei=c(.48,2),
          lwid=c(.2,1.2),
          key.par = list(cex=.5),
          key.xlab=expression(paste('log'[10],"(p-value)")),
          main="Top 20\nEnriched\nClusters",
          key.title="Color key",
          breaks=numbers,
          density.info="none"
          )
dev.off()


HeatmapSelectedGOParent <- read_csv("Enrichment_heatmap/HeatmapSelectedGOParent.csv")
x=paste(HeatmapSelectedGOParent$GO,":",HeatmapSelectedGOParent$Description)
HeatmapSelectedGOParent=as.matrix(abs(HeatmapSelectedGOParent[, -c(1:2)]))
rownames(HeatmapSelectedGOParent)=x


rb=c("#C0C0C0",colorRampPalette(reds)(n=99))
numbers=c(0,seq(0.01,7,length=100))

tiff("Biological2.tiff", width = 10, height = 6, units = 'in', res = 500, bg = "transparent")


heatmap.2(HeatmapSelectedGOParent, trace="n", Colv = NA,
          col=rb,
          dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
                                         expression(paste("100 ", mu,"M"))), labRow = row.names(HeatmapSelectedGOParent),
          distfun = function(x) dist(x,"manhattan"),
          hclustfun = function(x) hclust(x, method="ward.D"),
          cexRow=1.1,
          cexCol = 1.5,
          margins=c(5.1,39),
          lhei=c(.48,2),
          lwid=c(.2,1.2),
          key.par = list(cex=.5),
          key.xlab=expression(paste('log'[10],"(p-value)")),
          main="Biological\nPathways",
          key.title="Color key",
          breaks=numbers,
          density.info="none"
)
dev.off()

