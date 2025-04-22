#Load the toolkit you need for network analysis
library("igraph")
library("qgraph")
library("networktools")
library("ggplot2")
library("bootnet")
library("grDevices")
library("qgraph")
library("psych")
library("dplyr")
library("tidyverse")
library("showtext")
library("NetworkComparisonTest")
library("readxl")
library("bootnet")
library("boot")

setwd("C:\\R YY\\AIGC")
mydata<-read_excel("ZZZ.xlsx")
#Give each column a name based on your own data
myname<-c("Gen",	"Gra",	"Onl",	"Per",	"Fat",	"UAI"	,"GCAIESN",	"OCS","IS","Dep",	"Anx",	"AP",	"Mal",	"EI",
          "CTQ1",	"CTQ3",	"CTQ4",	"CTQ5")
colnames(mydata)<-myname
mydata.frame<-myname

#Grouping according to the attributes of each column scale, e.g. 1-4 columns belong to Dpress (the name is arbitrary)
feature_group<-list(A=c(1:6),B=c(7),C=c(8:14),D=c(15:18))
# Compute node predictability 
library(mgm) 
mydata <- as.matrix(mydata) 
p <- ncol(mydata)
fit_obj <- mgm(data = mydata,
               type = rep('g', p),
               level = rep(1, p),
               lambdaSel = 'CV',
               ruleReg = 'OR',
               pbar = FALSE)


pred_obj <- predict(fit_obj, mydata)
pred_obj
# Compute graph with tuning = 0.5 (EBIC) 
# That is, the network is constructed and the correlation matrix is sparsified by the Lasso algorithm, and the tuning parameter is set to 0.5

CorMat <- cor(mydata,method="spearman")
EBICgraph <- qgraph(CorMat, graph = "glasso", sampleSize = nrow(mydata), groups=feature_group,nodeNames=myname,threshold = TRUE,
                    tuning = 0.5, layout = "spring", details = TRUE, pie = pred_obj$error[,2])


#The method is the same as above, the basic relationship remains the same #partial correlation calculates the correlation of the depression scale, but this step must be run because MyNetwork needs to be used later
mynetwork <- estimateNetwork(mydata, default = "EBICglasso", tuning = 0.5,threshold = TRUE,corMethod="cor",corArgs=list(method="spearman",use="pairwise.complete.obs")) 
plot(mynetwork,layout = "spring",nodeNames=myname,groups=feature_group,pie = pred_obj$error[,2])

pdf(file = 'Network(0417).pdf', width = 18, height = 10) 
myplot<-plot(EBICgraph,layout = "spring",nodeNames=myname,pie = pred_obj$error[,2])   
dev.off()


# expected influence
dev.new()
myplot<-plot(mynetwork,layout = "spring",nodeNames=myname, weighted = TRUE, signed = TRUE)
pdf(file = 'EI.pdf', width = 4, height = 8)  
c1<-centralityPlot(myplot,include = "ExpectedInfluence",scale="z-scores",orderBy="ExpectedInfluence")
dev.off()

#CS
Result_all <- bootnet(mynetwork,statistics = c("ExpectedInfluence"),
                      nBoot=1000,nCore=8,type = "case")  #nboot一般要设置为1000

pdf(file='EI stability.pdf',width = 10,height = 7)
plot(Result_all, statistics = c("ExpectedInfluence"))
dev.off()
corStability(Result_all)

# Bootnet verifies the robustness of the network

b1<-bootnet(mynetwork,boots=1000,nCores=4,statistics=c("expectedInfluence","edge"))
#save bootstrapped files 保存的意义是免得下次要看结果的时候再跑，也可以选择不保存
save(b1,file = "b1.Rdata")

pdf(file='95CI.pdf',width = 15,height = 15)
plot(b1, lables = TRUE, order= "sample")
dev.off()

# bridge expected influence 
library(networktools)

#constructing a partial correlation matrix
myedges<-getWmat(mynetwork)
# write.csv(myedges,"MyNetcorkEdgess.csv")

#estimate bridge values for each node
#Here 1 represents scale group 1, and 4 1s are your data, and the first 4 columns are a group
mybridge<-bridge(myplot,communities = c('1','1','1','1','1','1',
                                        '2',
                                        '3','3','3','3','3','3','3',
                                        '4','4','4','4' ),
                 useCommunities="all", directed=NULL,nodes=NULL)

#create bridge expected influence graph
pdf(file = 'bridgeEI.pdf', width = 4, height = 8) 
plot(mybridge,include="Bridge Expected Influence (1-step)", order = "value", zscore = TRUE)
dev.off()

# 提取Bridge Expected Influence (1-step)
bridge_ei_1step <- mybridge$`Bridge Expected Influence (1-step)`
print(bridge_ei_1step)
# 计算z-score进行标准化
bridge_ei_1step_zscore <- scale(bridge_ei_1step)
print(bridge_ei_1step_zscore)

#sort and find bridge symptoms
a<-sort(mybridge[[4]])

#Bridge stability
caseDroppingBoot<-bootnet(mynetwork,boots=1000,type="case",
                          statistics = c("bridgeStrength","bridgeExpectedInfluence"),
                          communities = feature_group)

#get stability coefficients
corStability(caseDroppingBoot)

#plot centrality stability

pdf(file='Bridge EI stability.pdf',width = 10,height = 7)
#plot(Result_all, statistics = c("ExpectedInfluence"))
plot(caseDroppingBoot,statistics=c("bridgeExpectedInfluence"))
dev.off()