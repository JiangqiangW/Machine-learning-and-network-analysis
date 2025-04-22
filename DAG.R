library("grDevices")
library("qgraph")
library("psych")
library("dplyr")
library("tidyverse")
library("showtext") 
library("NetworkComparisonTest")
library("readxl")
library("bootnet")

setwd("C:\\R YY\\AIGC")
###Longitudinal Trajectory DAG
## DAG Bayesian Network
#install.packages('huge')
library("huge")                 
library("bnlearn")
#BiocManager::install("Rgraphviz",force=TRUE)
library(Rgraphviz)
set.seed(123)
library("readxl")
data_all<-read_excel("Z.xlsx")

blacklist <- data.frame(
  from = c("CTQ1", "CTQ3", "CTQ4", "CTQ5","OCS","IS","Dep","Anx","AP","Mal","EI","GCAIESN","Gra","Onl","Fat","Per","UAI",
           "CTQ1", "CTQ3", "CTQ4", "CTQ5","OCS","IS","Dep","Anx","AP","Mal","EI","GCAIESN","Gen","Onl","Fat","Per","UAI",
           "CTQ1", "CTQ3", "CTQ4", "CTQ5","OCS","IS","Dep","Anx","AP","Mal","EI","GCAIESN","Gen","Gra","Fat","Per","UAI",
           "CTQ1", "CTQ3", "CTQ4", "CTQ5","OCS","IS","Dep","Anx","AP","Mal","EI","GCAIESN","Gen","Gra","Onl","Per","UAI",
           "GCAIESN","GCAIESN","GCAIESN","GCAIESN","GCAIESN","GCAIESN","GCAIESN","GCAIESN","GCAIESN","GCAIESN","GCAIESN",
           "CTQ1", "CTQ3", "CTQ4", "CTQ5","OCS","IS","Dep","Anx","AP","Mal","EI","GCAIESN","Gra","Onl","Fat","Gen","UAI","GCAIESN"
           ),  # 起点变量名
  to   = c("Gen", "Gen", "Gen","Gen", "Gen", "Gen","Gen", "Gen", "Gen","Gen",  "Gen","Gen", "Gen", "Gen","Gen", "Gen", "Gen",
           "Gra", "Gra", "Gra","Gra", "Gra", "Gra","Gra", "Gra", "Gra","Gra", "Gra",  "Gra","Gra", "Gra", "Gra","Gra", "Gra", 
           "Onl", "Onl", "Onl","Onl", "Onl", "Onl","Onl", "Onl", "Onl","Onl", "Onl", "Onl","Onl", "Onl","Onl", "Onl", "Onl",
           "Fat", "Fat", "Fat","Fat", "Fat", "Fat","Fat", "Fat", "Fat","Fat", "Fat", "Fat","Fat", "Fat","Fat", "Fat", "Fat",
           "CTQ1", "CTQ3", "CTQ4", "CTQ5","OCS","IS","Dep","Anx","AP","Mal","EI",
           "Per","Per","Per","Per","Per","Per","Per","Per","Per","Per","Per","Per","Per","Per","Per","Per","Per","UAI"
           )  # 终点变量名
)
# 50 random re-starts and 100 perturbations for each re-start
data_DAG_2<-as.data.frame(apply(data_all,2,as.numeric))
data_DAG_2<-data_DAG_2[complete.cases(data_DAG_2),]
DAG2 <- hc(data_DAG_2, restart = 50, perturb = 100,blacklist = blacklist)  # hc gives directed graph
bnlearn::score(DAG2,data_DAG_2) # Global Network Score
DAG2_astr<-arc.strength(DAG2,data_DAG_2,"bic-g") # Strength of Connection
DAG2_astr[order(DAG2_astr[,3]),]

set.seed(123)
bootnet2 <- boot.strength(data_DAG_2, R = 1000,algorithm = "hc", algorithm.args = list(restart = 5,
    perturb = 10,
    blacklist = blacklist),  
  debug = TRUE
)

bootnet2
save(bootnet2,file = "bootnet2.Rdata")
## strength: connection strength, e.g. 0.86 means that this connection appears in 86% of the fitted networks.
## direction: probability of the direction, e.g. 0.57 means that in 57% of the fitted networks the connection goes in 
## the direction depicted in the graph.

## filter the ones with a strength larger than 0.85 and a direction probability > 0.5
bootnet2[bootnet2$strength > 0.85 & bootnet2$direction > 0.5, ]
library(openxlsx)
# Filter the data
filtered_data <- bootnet2[bootnet2$strength > 0.85 & bootnet2$direction > 0.5, ]
# Export as an Excel file
write.xlsx(filtered_data, file = "filtered_data.xlsx")
###########################
#### net1: Build an average network of 10,000 times. (Sachs et al., 2005, Science)
avgnet1.1 <- averaged.network(bootnet2, threshold = 0.85)
avgnet1.1
bnlearn::score(avgnet1.1, data_DAG_2)
DAG2_astr <- arc.strength(avgnet1.1, data_DAG_2, "bic-g")   ## compute edge strengths
#picture
pdf(file = 'lope_DAG_mean.pdf', width = 10, height = 7)
par(cex=3)
strength.plot(avgnet1.1, DAG2_astr, shape = "circle",fontsize = 22)
dev.off()