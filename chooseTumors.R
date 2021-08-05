rm(list=ls())

#COMMON GENE FINDER

#Install all packages necessary to call these libraries
library(survminer)
library(survival)
library(MASS)
library(enrichR)
library(combinat)
library(gtools)
library(enrichR)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorRamps)
library(colorspace)
library(lubridate)
library(grid)
library(ggtext)
library(ggthemes)
library(ggrepel)

geneOfInterest = "COL1A1"

#exact name of the cancer in your data file USER_INPUT
cancerListNames = list(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  #"Cervical",
  # "Colon_and_Rectal",
  #"Colon",
  "Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  "Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  #"Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  "Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  #"Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid"
  #"Uterine"
)

cancerListNames1 = c(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  #"Cervical",
  # "Colon_and_Rectal",
  #"Colon",
  "Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  "Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  #"Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  "Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  #"Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid"
  #"Uterine"
)

cancerListNames = list(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  #"Cervical",
  # "Colon_and_Rectal",
  #"Colon",
  "Endometrioid",
  #"Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  "Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  #"Lung_Squamous_Cell",
  #"Melanoma",
  #"Mesothelioma",
  "Ocular_Melanomas",
  #"Ovarian",
  "Pancreatic",
  #"Pheochromocytoma_Paraganglioma",
  "Prostate",
  #"Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid"
  #"Uterine"
)

cancerListNames1 = c(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal",
  "Colon",
  "Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  #"Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  "Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  #"Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid",
  "Uterine"
)

cancerListNames = list(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal",
  "Colon",
  "Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  #"Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  "Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  #"Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid",
  "Uterine"
)

cancerListNames1 = c(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal",
  "Colon",
  "Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  #"Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  "Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  #"Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid",
  "Uterine"
)

cancerListNames = list(
  #"Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal",
  "Colon",
  "Endometrioid",
  "Esophageal",
  #"Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  #"Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  "Lung_Squamous_Cell",
  #"Melanoma",
  "Mesothelioma",
  #"Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  #"Thyroid",
  "Uterine"
)

cancerListNames1 = c(
  #"Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal",
  "Colon",
  "Endometrioid",
  "Esophageal",
  #"Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  #"Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  "Lung_Squamous_Cell",
  #"Melanoma",
  "Mesothelioma",
  #"Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  #"Thyroid",
  "Uterine"
)

#COL4A1

cancerListNames = list(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal",
  "Colon",
  #"Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  "Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  "Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  "Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  #"Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid",
  "Uterine"
)

cancerListNames1 = c(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal",
  "Colon",
  #"Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma",
  "Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung",
  "Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  "Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  #"Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid",
  "Uterine"
)

setwd("/Applications/UPenn/Summer2020/TCGA_Data/TCGA_data")
numberofgenes=data.frame(matrix(ncol = 1, nrow = 1))
totalgenelist=data.frame(matrix(ncol = 1, nrow = 1))
#Creates list of names for file read-in's, and a list of the trak1 dataframes once read in
#Modify as necessary if your trak1 files labeled in a different manner
trak1CancerListNames = paste("trak1_", cancerListNames, "_",geneOfInterest, ".txt", sep = "")
trak1CancerList <- vector(mode = "list", length = length(cancerListNames))
kmsig2SurvivalList <- vector(mode = "list", length = length(cancerListNames))
for (i in 1:length(cancerListNames)) {
  trak1CancerList[[i]] = na.omit(read.table(trak1CancerListNames[i], header = FALSE,sep=",",stringsAsFactors = FALSE))
}

kmsig2SurvivalListNames = paste("kmsig2",cancerListNames, "survfit", sep = "_") 
for (i in 1:length(cancerListNames)) {
  kmsig2SurvivalList[[i]] = na.omit(read.table(kmsig2SurvivalListNames[i], header = TRUE,sep=" ",stringsAsFactors = FALSE))
}

genesurv=data.frame(matrix(ncol = 6, nrow = 1))
for (i in 1:length(cancerListNames)) {
  genesurv[i,1:5] = kmsig2SurvivalList[[i]][kmsig2SurvivalList[[i]]$Gene=="LMNB1",]
  genesurv[i,6]=cancerListNames1[i]
}
colnames(genesurv)=c('pvalue','se','z','gene','hazard_ratio','cancer')

write.table(genesurv, file = "TOP2Acoxph",  sep = "\t",row.names = FALSE,
            col.names = TRUE )

cancerListNames1[which(genesurv$X1<0.05)]

cancerListNames1[which(genesurv$X5<1)]


a= na.omit(read.table("kmsig2_Liver_coxph", header = TRUE,sep=" ",stringsAsFactors = FALSE))
b = na.omit(read.table("trak1_Liver_LMNB1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE))

c=a[a$Gene %in% b$V5,]
d=c[c$p.value_Coeff<0.05 & c$Hazard_Ratio>1,]

Reduce(intersect, list(trak1CancerList[[1]]$V5,trak1CancerList[[2]]$V5,trak1CancerList[[3]]$V5,trak1CancerList[[4]]$V5,trak1CancerList[[5]]$V5,trak1CancerList[[6]]$V5,trak1CancerList[[7]]$V5,trak1CancerList[[8]]$V5,trak1CancerList[[9]]$V5,trak1CancerList[[10]]$V5,trak1CancerList[[11]]$V5,trak1CancerList[[12]]$V5,trak1CancerList[[13]]$V5,trak1CancerList[[14]]$V5,trak1CancerList[[15]]$V5,trak1CancerList[[16]]$V5,trak1CancerList[[17]]$V5,trak1CancerList[[18]]$V5,trak1CancerList[[19]]$V5,trak1CancerList[[20]]$V5,trak1CancerList[[22]]$V5,trak1CancerList[[23]]$V5,trak1CancerList[[24]]$V5,trak1CancerList[[25]]$V5,trak1CancerList[[26]]$V5,trak1CancerList[[27]]$V5,trak1CancerList[[28]]$V5,trak1CancerList[[29]]$V5,trak1CancerList[[30]]$V5,trak1CancerList[[31]]$V5,trak1CancerList[[32]]$V5))

Reduce(intersect, list(trak1CancerList[[1]]$V5,trak1CancerList[[2]]$V5,trak1CancerList[[3]]$V5,trak1CancerList[[4]]$V5,trak1CancerList[[5]]$V5,trak1CancerList[[6]]$V5,trak1CancerList[[7]]$V5,trak1CancerList[[8]]$V5,trak1CancerList[[9]]$V5,trak1CancerList[[10]]$V5,trak1CancerList[[11]]$V5,trak1CancerList[[12]]$V5,trak1CancerList[[13]]$V5,trak1CancerList[[14]]$V5,trak1CancerList[[15]]$V5,trak1CancerList[[16]]$V5,trak1CancerList[[17]]$V5,trak1CancerList[[18]]$V5,trak1CancerList[[19]]$V5,trak1CancerList[[20]]$V5,trak1CancerList[[21]]$V5,trak1CancerList[[22]]$V5,trak1CancerList[[23]]$V5,trak1CancerList[[24]]$V5,trak1CancerList[[25]]$V5))


Reduce(intersect, list(trak1CancerList[[5]]$V5,trak1CancerList[[3]]$V5,trak1CancerList[[16]]$V5,trak1CancerList[[17]]$V5))

cancerListNames=cancerListNames[b>5]

int=intersect(key,key2)

numberofgenes=data.frame(matrix(ncol = 1, nrow = 1))
totalgenelist=trak1CancerList[c(1:32)]
#plot(1:32,log(choose(32,1:32),base=10))

numberofgenes=dim(trak1CancerList[c(1:32)])[1]

for (i in 1:length(cancerListNames))
{
  numberofgenes[i]=dim(trak1CancerList[[i]])[1]
}

totalgenelist=trak1CancerList[[1]]$V5

for (i in 2:length(cancerListNames))
{
  totalgenelist=c(totalgenelist,trak1CancerList[[i]]$V5)
}

totalgenelist1=unique(totalgenelist)
mat=data.frame(matrix(ncol = length(cancerListNames), nrow = length(totalgenelist1)))
for (i in 1:length(totalgenelist1))
{
  for(j in 1:length(cancerListNames))
  {
    mat[i,j]=totalgenelist1[i] %in% trak1CancerList[[j]]$V5
  }
}

mat1=mat
mat1[mat]=1
mat1[!mat]=0
a=rowSums(mat1)
b=colSums(mat1)
mat2=mat[,b>2]
rownames(mat1)=totalgenelist1
pbmc <- CreateSeuratObject(mat1)

cancerListNames=cancerListNames[b>2]
trak1CancerList=trak1CancerList[[b>5]]

all.genes <- rownames(pbmc)
all.genes=CNN1$V5
pbmc <- RunUMAP(pbmc,features = rownames(pbmc))
pbmc[["cancer"]]=cancerListNames1
pbmc[["genes"]]=numberofgenes

DimPlot(pbmc, reduction = "umap",group.by = "cancer",label=TRUE,repel = TRUE,label.size = 5)+NoLegend()
DimPlot(pbmc, reduction = "umap",group.by = "genes",label=TRUE,repel = TRUE,label.size = 5)+NoLegend()



mat2=mat1[a>1,]
mat=mat[a>1,]
key=totalgenelist1[a>1]

for (r in 19:27)
{
  r=14
  test=combinations(20,r)
  #genes1=data.frame(matrix(ncol = 15, nrow = 1))
  maxtest1=data.frame(matrix(ncol = r, nrow = 1))
  
  maxsum=0
  maxtest=0
  maxsum1=data.frame(matrix(ncol = 1, nrow = 1))
  k=1
  for (i in 1:dim(test)[1])
  {
    combsum=sum(rowSums(mat2[,test[i,]])==r)
    
    if(combsum>maxsum)
    {
      #genes=key[which(rowSums(mat2[,test[i,]])==r)]
      maxsum=combsum
      maxtest=test[i,]
    }
    if(combsum==maxsum & maxsum>1)
    {
      #genes1[k,]=key[which(rowSums(mat2[,test[i,]])==r)]
      maxtest1[k,]=test[i,]
      maxsum1[k]=combsum
      k=k+1
    }
  }
  
  max9=cbind(t(maxsum1),maxtest1)
}

maxtest1=as.matrix(maxtest1)
key[which(rowSums(mat2[,maxtest1[1,]])==r)]


key2=totalgenelist1[a>1]

cancerListNames1[which(mat2[which(key=="FAP"),]==1)]


test=combinations(27,20)

for (r in 19:27)
{
 r=17
test=combinations(30,r)
maxtest1=data.frame(matrix(ncol = r, nrow = 1))
genes=data.frame(matrix(ncol = 4, nrow = 1))
maxsum=0
maxtest=0
maxsum1=data.frame(matrix(ncol = 1, nrow = 1))
k=1
for (i in 1:dim(test)[1])
{
    combsum=sum(rowSums(mat2[,test[i,]])==r)
  
  if(combsum>maxsum)
  {
    #genes=key[which(rowSums(mat2[,test[i,]])==r)]
    maxsum=combsum
    maxtest=test[i,]
  }
  if(combsum==maxsum & maxsum>3)
  {
  genes[k,]=key[which(rowSums(mat2[,test[i,]])==r)]
  maxtest1[k,]=test[i,]
  maxsum1[k]=combsum
  k=k+1
  }
}

max9=cbind(t(maxsum1),maxtest1)

#maxoverlap <- data.frame(matrix(ncol = 1, nrow = 1))
#maxoverlap[r,1]=maxsum
#fileNames = paste("MaxOverlap", r, sep = "")
write.table(max9, file = "MaxOverlap_LMNA_18choose3)b_atleast2", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(genes, file = "GenesMaxOverlap_COL1A1_30choose17pt2", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(key[which(rowSums(mat2[,maxtest])==r)], file = "9genesCOL1A117tumors", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

}
max20=maxtest

ab=max9[max9$`t(maxsum1)`==5,2:18]
ab=as.matrix(ab)
genes=data.frame(matrix(ncol = 1, nrow = 1))
cancers=data.frame(matrix(ncol = 1, nrow = 1))
for (i in 1:dim(ab)[1])
{
  ab_trak=trak1CancerList[ab[i,]]
  genes[i,]=ab_trak[[i]][ab_trak[[i]]$V5 %in% key[which(rowSums(mat2[,ab[i,]])==r)],5]
}

acta2 = read.table("MaxOverlap_COL4A1_30choose17", header = TRUE,sep=" ",stringsAsFactors = FALSE)
ab=acta2[acta2$t.maxsum1.==5,2:18]
ab=as.matrix(ab)
genes[2,]=key[which(rowSums(mat2[,ab[2,]])==r)]

write.table(genes, file = "GENESACTA217tumors", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

genes = read.table("GenesMaxOverlap_COL1A1_30choose17pt2", header = TRUE,sep=" ",stringsAsFactors = FALSE)

ab=acta2[acta2$t.maxsum1.==10,2:18]
ab=as.matrix(ab)
genes[2,]=key[which(rowSums(mat2[,ab[2,]])==r)]

n=1
ab_trak=trak1CancerList[ab[n,]]
genes[n,]=key[which(rowSums(mat2[,ab[n,]])==r)]
cancerListNames1[as.matrix(ab[n,])]

write.table(genes1, file = "GenesMaxOverlap_COL1A1_30choose17", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

Reduce(intersect, list(genes[1,],genes[2,]))

#
acta2 = read.table("MaxOverlap19", header = TRUE,sep=" ",stringsAsFactors = FALSE)
ab=acta2[acta2$t.maxsum1.==5,2:18]
ab=as.matrix(ab)
genes[2,]=key[which(rowSums(mat2[,ab[2,]])==r)]


#genes1=data.frame(matrix(ncol = 10, nrow = 1))
#names(ab_trak)=cancerListNames1[ab[2,]]

genes1=data.frame(matrix(ncol = 5, nrow = 1))
exp=data.frame(matrix(ncol = 10, nrow = 17))
for(i in 1:17)
{
  exp[i,]=ab_trak[[i]][ab_trak[[i]]$V5 %in% genes1[n,],1]
  genes2[i,]=ab_trak[[i]][ab_trak[[i]]$V5 %in% genes1[n,],5]
}

which(trak1CancerList[[c(ab[1,],5)]] == genes[1,])
geneExponent[k,l] = trak1CancerList[[c(cancerIndexes[j,l],1,rowID)]]

ab_trak=trak1CancerList[ab[2,]]
names(ab_trak)=cancerListNames1[ab[2,]]

exp=data.frame(matrix(ncol = 7, nrow = 17))
for(i in 1:17)
{
exp[i,]=ab_trak[[i]][ab_trak[[i]]$V5 %in% genes[2,],1]
}
cancer=cancerListNames1[ab[2,]]

n=6
ab_trak=trak1CancerList[as.matrix(maxtest1[n,])]
genes1[n,]=key[which(rowSums(mat2[,as.matrix(maxtest1[n,])])==r)]
cancerListNames1[as.matrix(maxtest1[n,])]

n=21
ab_trak=trak1CancerList[as.matrix(ab[n,])]
genes1[n,]=key[which(rowSums(mat2[,as.matrix(ab[n,])])==r)]
cancerListNames1[as.matrix(ab[n,])]
cancers=cancerListNames1[as.matrix(ab[1,])]

for (i in 2:21)
{
  cancers=c(cancerListNames1[as.matrix(ab[n,])],cancers)
}
c=unique(cancers)

exp=data.frame(matrix(ncol = 4, nrow = 18))
for(i in 1:18)
{
  ind=which(cancerListNames1==c[i])
  exp[i,1]=trak1CancerList[[ind]][trak1CancerList[[ind]]$V5== "COL4A2",1]
  exp[i,2]=trak1CancerList[[ind]][trak1CancerList[[ind]]$V5== "COL15A1",1]
  exp[i,3]=trak1CancerList[[ind]][trak1CancerList[[ind]]$V5== "NID2",1]
  exp[i,4]=trak1CancerList[[ind]][trak1CancerList[[ind]]$V5== "LAMA4",1]
}


#genes1=data.frame(matrix(ncol = 10, nrow = 1))
#names(ab_trak)=cancerListNames1[ab[2,]]

genes2=data.frame(matrix(ncol = 5, nrow = 1))
exp=data.frame(matrix(ncol = 5, nrow = 17))
for(i in 1:17)
{
  exp[i,]=ab_trak[[i]][ab_trak[[i]]$V5 %in% genes1[n,],1]
  genes2[i,]=ab_trak[[i]][ab_trak[[i]]$V5 %in% genes1[n,],5]
}

cancer=cancerListNames1[ab[2,]]
genes2[1,]=ab_trak[[1]][ab_trak[[1]]$V5 %in% genes1[1,],5]

which(ab_trak[[1]]$V5 %in% genes1[1,])

acta2 = read.table("MaxOverlap17", header = TRUE,sep=" ",stringsAsFactors = FALSE)
ab=acta2[acta2$t.maxsum1.==26,2:18]
genes[1,]=key[which(rowSums(mat2[,ab[1,1:17]])==r)]


ab=as.matrix(ab)

tumornames19=cancerListNames[max3]
tumornames=cancerListNames[max3]

tumornames_17=cancerListNames[max3]

tumornames_18=cancerListNames[max3]

b1=intersect(tumornames_17,tumornames_18)

barplot(as.matrix(numberofgenes),names.arg = cancerListNames,horiz = TRUE)

barplot(table(totalgenelist))

a=data.frame(cancerListNames1,(numberofgenes))
colnames(a)=c("Tumor","Number_of_Genes")
ggplot(data=a,aes(reorder(Tumor,Number_of_Genes),Number_of_Genes))+geom_bar(stat="identity")+coord_flip()

b=data.frame(table(totalgenelist))
b1=b[b$Freq>16,]
#b1=b1[order(b1$Freq,decreasing = TRUE),]
ggplot(data=b1,aes(x=reorder(totalgenelist,Freq),y=Freq))+geom_bar(stat="identity")+coord_flip()

hist(b$Freq)
c=data.frame(table(b$Freq))
ggplot(data=c[-1,],aes(x=reorder(Var1,Freq),y=Freq))+geom_bar(stat="identity")+coord_flip()

overlap20=Reduce(intersect, trak1CancerList[maxtest]$V5)

maxintersect1=0

mean(CNN1[CNN1$V1<1.25,3])

hist(CNN1[,1])
sd(CNN1[,1])



cancerListNames = list(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal", #IGNORE
  "Colon",
  "Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma", #IGNORE
  "Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung", #IGNORE
  "Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  "Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid",
  "Uterine"
)

cancerListNames1 = c(
  "Adrenocortical",
  "Bile_Duct",
  "Bladder",
  "Breast",
  "Cervical",
  # "Colon_and_Rectal", #IGNORE
  "Colon",
  "Endometrioid",
  "Esophageal",
  "Glioblastoma",
  "Head_and_Neck",
  "Kidney_Chromophobe",
  "Kidney_Clear_Cell",
  "Kidney_Papillary_Cell",
  "Large_B-cell_Lymphoma",
  "Liver",
  # "Lower_Grade_Glioma_and_Glioblastoma", #IGNORE
  "Lower_Grade_Glioma",
  "Lung_Adenocarcinoma",
  # "Lung", #IGNORE
  "Lung_Squamous_Cell",
  "Melanoma",
  "Mesothelioma",
  "Ocular_Melanomas",
  "Ovarian",
  "Pancreatic",
  "Pheochromocytoma_Paraganglioma",
  "Prostate",
  "Rectal",
  "Sarcoma",
  "Stomach",
  "Testicular",
  "Thymoma",
  "Thyroid",
  "Uterine"
)

#Creates list of names for file read-in's, and a list of the trak1 dataframes once read in
#Modify as necessary if your trak1 files labeled in a different manner
trak1CancerListNames = paste("trak1_", cancerListNames, "_",geneOfInterest, ".txt", sep = "")
trak1CancerList <- vector(mode = "list", length = length(cancerListNames))
kmsig2SurvivalList <- vector(mode = "list", length = length(cancerListNames))
for (i in 1:length(cancerListNames)) {
  trak1CancerList[[i]] = na.omit(read.table(trak1CancerListNames[i], header = FALSE,sep=",",stringsAsFactors = FALSE))
}













