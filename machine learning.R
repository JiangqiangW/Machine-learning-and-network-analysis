#Set up a working directory
setwd("C:\\R YY\\ROC+\\ROC")
####ML####
data<-read_excel("ZZZ.xlsx")
train_sub= sample(nrow(data),7/10*nrow(data))

#Divide the dataset into a training set and a test set
train_data=data[train_sub,]
test_data=data[-train_sub,]
# Export the training set to a CSV file
write.csv(train_data, file = "train_data.csv", row.names = FALSE)

# Export the test set to a CSV file
write.csv(test_data, file = "test_data.csv", row.names = FALSE)
#安装包
# install.packages(c("seqinr", "plyr", "openxlsx", "randomForestSRC", "glmnet", "RColorBrewer"))
# install.packages(c("ade4", "plsRcox", "superpc", "gbm", "plsRglm", "BART", "snowfall"))
# install.packages(c("caret", "mboost", "e1071", "BART", "MASS", "pROC", "xgboost"))
# install.packages('maps')
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("mixOmics",force = TRUE)
# BiocManager::install("survcomp",force = TRUE)
# BiocManager::install("ComplexHeatmap",force = TRUE)
# install.packages('limma')
#加载包
library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
library(ggpubr)
library("readxl")


#Loads scripts for model training and model evaluation
source("ML.R")

#Read in the training set expression matrix and group information (outcome events)

Train_expr <- read.table("trains.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_class <- read.table("trainc.txt", header = T, sep ="\t", row.names = 1,check.names = F,stringsAsFactors = F)
#Read in the training set expression matrix and group information (outcome events)
Test_expr <- read.table('tests.txt', header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_class <- read.table('testc.txt', header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

Train_set <- Train_expr
Test_set <- Test_expr

Train_set <- as.matrix(Train_set)
Test_set <- as.matrix(Test_set)

#Get a list of model names
methods <- read.xlsx("methods.xlsx", startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)

# settings, column names for outcome events, and the minimum number of variables to be included

classVar = "outcome"
min.selected.var = 3

#Pre-training, which runs the algorithm for variable filtering first
preTrain.var <- list()
set.seed(seed = 777)

Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1]) 
preTrain.method = unique(unlist(preTrain.method)) 
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, 
                                 Train_set = Train_set, 
                                 Train_label = Train_class, 
                                 #Variable(筛选变量)和Model(获取模型)
                                 mode = "Variable", 
                                 classVar = classVar)
}
preTrain.var[["simple"]] <- colnames(Train_set)

#Model training
model <- list()
set.seed(seed = 777)
Train_set_bk = Train_set
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method 
  method <- strsplit(method, "\\+")[[1]] #各步骤算法名称
  if (length(method) == 1) method <- c("simple", method) 
  Variable = preTrain.var[[method[1]]]
  Train_set = Train_set_bk[, Variable]
  Train_label = Train_class 
  model[[method_name]] <- RunML(method = method[2],
                                Train_set = Train_set,
                                Train_label = Train_label,
                                #Variable and Model
                                mode = "Model",
                                classVar = classVar)
  
  #If the final model includes less than min.selected.var = 3, it is meaningless
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
Train_set = Train_set_bk

#save
saveRDS(model, "model.rds") 

#read
model <- readRDS("model.rds")
methodsValid <- names(model)

#Sample risk scores are calculated
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set))
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(data.frame(ID=rownames(RS_mat),RS_mat), "RS_mat.txt",sep = "\t", row.names = F, quote = F)

#Predictive classification
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
write.table(data.frame(ID=rownames(Class_mat),Class_mat),"Class_mat.txt", sep = "\t", row.names = F, quote = F)

#Get a list of variables in the model
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}
fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, "fea_df.txt",sep = "\t", row.names = F, col.names = T, quote = F)

#Calculate the C-index of each model
AUC_list <- list()

Test_class$Outcome <- Test_class$outcome
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],
                                Test_set = Test_set,
                                Test_label = Test_class,
                                Train_name = "Training", 
                               
                                Train_set = Train_set,
                                Train_label = Train_class,
                                # Train_set = NULL,
                                # Train_label = NULL,
                                cohortVar = "Cohort",
                                classVar = classVar)
}

AUC_mat <- do.call(rbind, AUC_list)
write.table(data.frame(ID=rownames(AUC_mat),AUC_mat), "AUC_mat.txt",sep = "\t", row.names = F, col.names = T, quote = F)

#Heatmap
AUC_mat <- read.table("AUC_mat.txt",sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean) 
avg_AUC <- sort(avg_AUC, decreasing = T) 
AUC_mat <- AUC_mat[names(avg_AUC), ] 
fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3)) 
#colour
if(ncol(AUC_mat) < 3) { 
  CohortCol <- c("red","blue")
} else { 
  CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired") 
}
names(CohortCol) <- colnames(AUC_mat)
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(AUC_mat, 
                    avg_AUC, 
                    CohortCol, "steelblue", 
                    cellwidth = cellwidth, cellheight = cellheight, 
                    cluster_columns = F, cluster_rows = F) 
pdf("AUC.pdf", width = cellwidth * ncol(AUC_mat) + 7, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
invisible(dev.off())

########Visualize the screening process of the model
library(ggplot2)
library(dplyr)
library(forcats)
library(gridExtra)
library(MuMIn)
library(ggplot2)

options(na.action = "na.fail")

Train_class <- as.factor(Train_class)
Train_set_df <- as.data.frame(Train_set, na.action = "na.omit") # 移除含有缺失值的行

fit <- glm(Train_class ~ ., data = Train_set_df, family = "binomial")
set.seed(777)
dredge_results <- dredge(fit)
coef_data <- data.frame(coef(dredge_results))
coef_data$Model <- rownames(coef_data)
coef_data_long <- reshape2::melt(coef_data, id.vars = "Model", variable.name = "Variable", value.name = "Coefficient")

# Extract the variable selection process of the step function
step_history <- step_fit$anova # Obtain the detailed step record of the step function

# Extract variable names and actions
step_df <- data.frame(
  Step = 1:nrow(step_history),
  Variable = sub("- ", "", step_history$Step), # Extract the variable name and remove the prefix
  Action = ifelse(is.na(step_history$Step), "None", "Removed"), # Determine the type of operation
  AIC = step_history$AIC
)

# Remove the NA value
step_df <- na.omit(step_df)

#Make sure that all variable types are consistent
step_df$Variable <- as.character(step_df$Variable)
step_df$Action <- as.character(step_df$Action)

# Sort the variables by step
step_df <- step_df %>%
  arrange(desc(Step)) %>%
  mutate(Variable = fct_reorder(Variable, Step)) 

# Draw a bar chart of the variable selection process
p1 <- ggplot(step_df, aes(x = Step, y = Variable, fill = Action)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Removed" = "red", "None" = "grey")) +
  labs(x = "Step", y = "Variable", fill = "Action", title = "Variable Selection Process") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5), # 居中标题
    panel.border = element_rect(colour = "black", fill=NA, size=1), # 添加图形边框
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # 设置图的边距
  ) +
  coord_fixed(ratio = 1) # 固定坐标轴比例
p1
# Draw a line plot of AIC changes
p2 <- ggplot(step_df, aes(x = Step, y = AIC)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  labs(x = "Step", y = "AIC", title = "AIC Changes Over Steps") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5), 
    panel.border = element_rect(colour = "black", fill=NA, size=1), 
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # Set the margins of the diagram
  ) +
  coord_fixed(ratio = 1)
p2

############################################计算最好模型的重要程度###########
library(randomForest)
library(glmnet)
library(caret)

# Read variables in the model
model.gene <- read.table("fea_df.txt", header = TRUE, sep = "\t", check.names = FALSE)
best.model <- "Stepglm[both]+RF"  # 最优算法的名字
best.model.gene <- subset(model.gene, algorithm == best.model)[, 1]

# Training set expression matrix
Train_expr <- read.table("trains.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
Train_class <- read.table("trainc.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# Extract the variables used in the model
Train_expr <- Train_expr[, best.model.gene]

# Convert the data to a matrix format
Train_matrix <- as.matrix(Train_expr)
Train_class_matrix <- as.matrix(Train_class)

# Train a random forest model
set.seed(777)
rf_model <- randomForest(as.factor(outcome) ~ ., data = data.frame(Train_class_matrix, Train_matrix), importance = TRUE)

# Extract variable importance
importance <- importance(rf_model)
importance_df <- data.frame(Gene = rownames(importance), Importance = importance[, 1])

# Sort by importance
importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]

# Create a lollipop chart
ggplot(importance_df, aes(x = Importance, y = reorder(Gene, Importance))) +
  geom_point(aes(color = Importance), size = 4) +
  geom_segment(aes(x = ifelse(Importance > 0, 0, Importance), xend = Importance, y = reorder(Gene, Importance), yend = reorder(Gene, Importance)), color = "grey", size = 1) +
  scale_color_gradientn(colors = c("#66c2a5", "#fc8d62"), 
                        values = scales::rescale(c(-10, quantile(importance_df$Importance, 0.25), 
                                                   quantile(importance_df$Importance, 0.5), 
                                                   quantile(importance_df$Importance, 0.75), 
                                                   max(importance_df$Importance))), 
                        limits = c(-10, max(importance_df$Importance))) +
  labs(title = "Variable Importance in Random Forest Model",
       x = "Importance",
       y = "Variable") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "grey"),
    plot.background = element_rect(fill = "white", color = "grey"),
    legend.position = "right"
  ) +
  guides(color = guide_colorbar(title = "Importance", title.position = "top", title.hjust = 0.5)) +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  coord_cartesian(xlim = c(-10, 40))
