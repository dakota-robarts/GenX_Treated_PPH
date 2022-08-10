library(chorddiag)
library(heatmap.plus)
library(gplots)
library(readxl)
library(RColorBrewer)
library(readr)
library(xlsx)
library(ggplot2)


raw.df <- read_csv("D:/OneDrive/Lab File/GenX/RNA Gene Array/Newest Analysis/Genes_withValues_with2foldcutoff.txt")
raw.df=raw.df[,c(1,3,4,5,2)]
setwd("D:/OneDrive/Lab File/GenX/RNA Gene Array/Newest Analysis/Dose Response")

names=c(.1,10,100)

# loc=c()
# 
# for (i in 1:length(raw.df$`Gene Symbol`)){
# if (raw.df$`10uMFold Change`[i] == raw.df$`0.1uMFold Change`[i]){
#   loc=c(loc,i)
# }
#   print(i)
# }
# 
# raw.df=raw.df[-loc,]

for( i in 1:length(raw.df$`Transcript Cluster ID`)){
  raw.df$cor[i]=cor(names,as.numeric(raw.df[i,3:5]))
}


#looking for duplicates gene values from 10 and 0.1 concentrations 
# loc=c()
# for (i in 1:length(raw.df$`0.1uMFold Change`)){
#   if(raw.df$`0.1uMFold Change`[i] ==raw.df$`10uMFold Change`[i]){
#     loc=c(loc,i)
#   }
# }
# 
# raw.df=raw.df[-loc,]
#All 0.99 and -0.99 correlation============================================================
dose.dependent=raw.df[which(raw.df$cor >0.999 | raw.df$cor < -0.999),] # 951 genes

#Calculating Slope for each dose response
# x=c(0.1,10,100)
# for( i in 1:length(dose.dependent$`Gene Symbol`)){
#   model=summary(lm(as.numeric(as.data.frame(dose.dependent[i,3:5])) ~ x))
#   dose.dependent$Slope[i]=model$coefficients[2,1]
#   }
dose=c(0.1,10,100)

# Example plot of linear model 
# FC.gene=as.numeric(as.data.frame(dose.dependent[70,3:5]))
# df=as.data.frame(cbind(FC.gene,dose))
# col=brewer.pal(3,"Set2")
# tiff("Neg Correlation (ex.TBC1D1).tiff", width = 6, height = 4, units = 'in', res = 600, bg = "transparent")
# ggplot(df,aes(dose, FC.gene))+ 
#   geom_point(size=2.5,col=col[1])+
#   labs(y= "Gene Fold Change", x = expression(paste("Dose ",mu,"M")))+
#   ggtitle("Fitting Linear Model")+
#   geom_smooth(method = "lm",se = T, col=col[2],linetype = "dashed", size=.5 )+
#   theme(plot.title = element_text(hjust = 0.5,size = 18))
# dev.off()


# write.xlsx(dose.dependent,file="Correlation with slopes.xlsx")

# write.xlsx(dose.dependent,file="Correlated(All0.99)AllGenes.xlsx")
dose.dependent.up=raw.df[which(raw.df$cor >0.99),] # Only above .95 (848 genes)
dose.dependent.down=raw.df[which(raw.df$cor < -0.99),] #Only below -0.95 (375 genes)

dose.dependent$`0.1uMFold Change`[which(dose.dependent$`0.1uMFold Change` < 0)]=
  -1/dose.dependent$`0.1uMFold Change`[which(dose.dependent$`0.1uMFold Change` < 0)]
dose.dependent$`10uMFold Change`[which(dose.dependent$`10uMFold Change` < 0)]=
  -1/dose.dependent$`10uMFold Change`[which(dose.dependent$`10uMFold Change` < 0)]
dose.dependent$`GenX_100uM vs. GenX_UT`[which(dose.dependent$`GenX_100uM vs. GenX_UT` < 0)]=
  -1/dose.dependent$`GenX_100uM vs. GenX_UT`[which(dose.dependent$`GenX_100uM vs. GenX_UT` < 0)]

dose.dependent[,3:5]=log(dose.dependent[,3:5],2)

reds=brewer.pal(6,"Oranges")
blues=rev(brewer.pal(6,"Blues"))
rb=colorRampPalette(c(blues,reds))(n=99)
rb
numbers=c(seq(-1.75,1.75,length=100))

# tiff("Correlated(All.99)AllGenes.tiff", width = 4, height = 8, units = 'in', res = 600, bg = "transparent")
heatmap.2(as.matrix(dose.dependent[, -c(1,2,6)]), trace="none", Colv = NA,
          col=rb,
          dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
                                         expression(paste("100 ", mu,"M"))), labRow = "",
          distfun = function(x) dist(x,"manhattan"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          cexRow=.8,
          cexCol = 1.5,
          margins=c(5.1,4.5),
          lhei=c(.5,3.5),
          # lwid=c(.2,1),
          key.par = list(cex=.5),
          key.xlab="Log2(FC)",
          main="All Correlations",
          key.title="Color key",
          breaks=numbers,
          density.info = "none")

# dev.off()


  #Only 0.99 correlation============================================================
# dose.dependent=raw.df[which(raw.df$cor >0.99),]
# write.xlsx(dose.dependent,file="Correlated(+0.99)AllGenes.xlsx")

reds=brewer.pal(6,"Oranges")
blues=rev(brewer.pal(6,"Blues"))
rb=colorRampPalette(c(blues,reds))(n=99)
rb
numbers=c(seq(-1.75,1.75,length=100))

tiff("Correlated(.99)AllGenes.tiff", width = 4, height = 8, units = 'in', res = 600, bg = "transparent")
heatmap.2(as.matrix(dose.dependent[which(dose.dependent$cor > 0.99), -c(1,2,6)]), trace="none", Colv = NA,
          col=rb,
          dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
                                         expression(paste("100 ", mu,"M"))), labRow = "",
          distfun = function(x) dist(x,"manhattan"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          cexRow=.8,
          cexCol = 1.5,
          margins=c(5.1,4.5),
          lhei=c(.5,3.5),
          # lwid=c(.2,1),
          key.par = list(cex=.5),
          key.xlab="Log2(FC)",
          main="Positive\nCorrelations",
          key.title="Color key", 
          breaks=numbers,
          density.info = "none")

dev.off()

#Only -0.99 correlation============================================================
dose.dependent=raw.df[which(raw.df$cor < -0.99),]
# write.xlsx(dose.dependent,file="Correlated(-0.99)AllGenes.xlsx")

reds=brewer.pal(6,"Oranges")
blues=rev(brewer.pal(6,"Blues"))[-6]
rb=colorRampPalette(c(blues,reds))(n=99)
rb
numbers=c(seq(-1.75,2,length=100))

tiff("Correlated(-.99)AllGenes.tiff", width = 4, height = 8, units = 'in', res = 600, bg = "transparent")
heatmap.2(as.matrix(dose.dependent[which(dose.dependent$cor < -0.99), -c(1,2,6)]), trace="none", Colv = NA,
          col=rb,
          dendrogram = "row", labCol = c(expression(paste("0.1 ", mu,"M")),expression(paste("10 ", mu,"M")),
                                         expression(paste("100 ", mu,"M"))), labRow = "",
          distfun = function(x) dist(x,"manhattan"),
          hclustfun = function(x) hclust(x, method="ward.D2"),
          cexRow=.8,
          cexCol = 1.5,
          margins=c(5.1,4.5),
          lhei=c(.5,3.5),
          # lwid=c(.2,1),
          key.par = list(cex=.5),
          key.xlab="Log2(FC)",
          main="Negative\nCorrelations",
          key.title="Color key", 
          breaks=numbers,
          density.info = "none")

dev.off()
