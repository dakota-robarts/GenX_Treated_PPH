library(ggpubr)
library(readxl)
library(dplyr)
library(tidyr)
library(RColorBrewer)
setwd("D:/OneDrive/Lab File/GenX/WesternBlot/New/Densitometry")
dat <-read_excel("D:/OneDrive/Lab File/GenX/WesternBlot/New/Densitometry/values_long.xlsx", 
                            skip = 0)
dat$`Liver #`[which(dat$`Liver #` == 6)]=1
dat$`Liver #`[which(dat$`Liver #` == 7)]=2
dat$`Liver #`[which(dat$`Liver #` == 8)]=3
dat$`Liver #` <- factor(dat$`Liver #`,levels = c("1","2","3"))
dat$Treatment <- factor(dat$Treatment,levels=c("UT","0.1","10","100"))
table(dat$Value)
data <- dat %>% 
  filter(.$Value=="Norm to GAPDH/Paired UT")
colnames(data)

genes <- names(table(dat$Group))

for (gene in genes){
data %>% 
  filter(.$Group == gene)%>% 
  filter(!.$Treatment == "UT")%>% 
  ggplot()+
  geom_bar(aes(x=Treatment,y=Dens),shape=8,fill="gainsboro",stat = "summary", fun.y = "mean")+
  geom_point(aes(x=Treatment,y=Dens,color=`Liver #`,shape=`Liver #`,group=`Liver #`),size=2)+
  geom_line(aes(x=Treatment,y=Dens,color=`Liver #`,shape=`Liver #`,group=`Liver #`),size=.5)+ggtitle(paste(gene,"\nDensitometry"))+
  scale_color_brewer(palette="Dark2")+ theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"))+
  ylab("Relative Density to\n Paired UT")+
  xlab(expression(paste("Concentration (",mu,"M)")))
ggsave(paste(gene,"_Den.tiff",collapse = ""),device = "tiff",width = 3.5,height = 4,units = "in",dpi = 600)
}
