########################################################
# 2021-05-10
# DNA analysis: cluster, hotmap, and RNA degradtion
# dzr17723980497@gmail.com
########################################################
library('affyPLM')
library('affy')
library('affydata')
library('factoextra')
library('NbClust')
library('cowplot')

data(Dilution)
dil_stand <- gcrma(Dilution) ##Apply the robust multi-array average (RMA)
pset <- fitPLM(Dilution)  ##Do the PLM fitting
exp_dil<-exprs(dil_stand) ##Find the expression data for the chip

##kmeans
km<-kmeans(scale(c1),9)
factoextra::fviz_cluster(km, data = c1, 
                         ellipse.type = "confidence",
                         show.legend.text = FALSE)
##Hotmap for the DNA expression
pheatmap(exp_dil,show_rownames = F)

##RNA degredation map for chip
image(Dilution)


