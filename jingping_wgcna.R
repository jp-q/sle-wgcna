#Jingping Qiao
#
#
#Description:
#
#This R script performed two WGCNA to find new Signatures(Clusters) which related to SLE. The first WGCNA was 
#performed on stopsle.rsem.logfc with cell type CD16 Monocytes. After the first WGCNA, we got 31 clusters. We 
#combine these clusters with msigdb and built in signature list, got a new big signature list. Then we perform 
#score signaturing on logfc data using this big signature list and get z-scores. Perform t-test on these z-score 
#between health and sle, get p values, filter the signatures with p <0.05(which means significant), and get a
#new list of signatures. Now we perform the second WGCNA on new list of signatures.find small sized signature 
#(around 100 genes), check what other signature was clustered with them, and try to interpret the clusters. 
#
# --WGCNA on stopsle.rsem.logfc/signature list
# --Start date: 07-12-2017
# --End date:
# --Package: bioconducter:bioLite("WGCNA")
# --Reference: https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#                 FemaleLiver-02-networkConstr-blockwise.pdf;
#              https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html;
#              https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#                 Simulated-05-NetworkConstruction.pdf;
#              https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA
###########################################################################################

#################################################################################################

#####      #####  #     #   ###     ###    ####   #####   ##      ##   #   #####     #####
#####        #    ##   ##  #   #   #   #   #   #    #     # #     # #  #     #       #####
#####        #    # # # #  #   #  #     #  ####     #     ####    #  # #     #       #####
#####        #    #  #  #  ####    #   #   #  #     #     #   #   #   ##     #       #####
#####     #####   #     #  #        ###    #   #    #     #    #  #    #     #       #####

#If you want to see the result directly WITHOUT run any codes(the whole code takes about 1.5 hour to run)
#you could load these datasets to check the results, result presenting is at line 362
load("/home/x189259/GITCLONES/stop-sle-expression-matrix/JINGPING_WGCNA.Rdata")

################################################################################################




#setting workindg dir and source data
setwd("~/GITCLONES/stop-sle-expression-matrix/")
source(sprintf("%s/stop-sle-expression-matrix/load-rsem.R", Sys.getenv("GITCLONES")))
set.seed(1)


##Data shaping---------------------------------------------------------------------------------
#filter out baseline data and match logfc data with it
stopsle.rsem.metadata_baseline <- stopsle.rsem.metadata[stopsle.rsem.metadata$timepoint == "baseline",]
baseline_list <- stopsle.rsem.metadata_baseline$Sequences
stopsle.rsem.logfc_matched <- stopsle.rsem.logfc[,baseline_list]
stopsle.rsem.logtpm_matched <- stopsle.rsem.logtpm[,baseline_list]

#check hist and set up a cut-off "0" using log tpm data
log_tpm_data <- as.matrix(log(stopsle.rsem.logtpm_matched))
hist(log_tpm_data, breaks=200)
abline(v=0, col = "red")

#seperate data by cell types
cell_types <- split( stopsle.rsem.metadata_baseline, f = stopsle.rsem.metadata_baseline$cells)


#write a loop to seperate the list to data frame and clean out cut-offs (still log tpm)
name_list <- names(cell_types)
  for(i in 1:length(name_list)){
    #a<- gsub("\\s","",name_list[i])
    a <- data.frame(cell_types[i])[,3]
    c <- a %in% names(stopsle.rsem.logtpm_matched)
    #check if name is in database, if not, remove it.
    if("FALSE" %in% c){
      a <- as.data.frame(a)
      a <- a[c,]
      
    }
    #cut off compare to median
a <- stopsle.rsem.logtpm_matched[,as.character(a)]
a <- a[apply(a,1,median) > 0,]
a <- stopsle.rsem.logfc_matched[rownames(a),colnames(a)]
assign(gsub("\\s","",name_list[i]),a)


  }


msigdb <- scan("msigdb.txt", what="list", sep="\n")
names <- msigdb[c(TRUE,FALSE)]
a1 <- msigdb[c(FALSE,TRUE)]
# Separate elements by one or more whitepace
y1 <- strsplit(a1, "\\s+")
y1 <- y1[y1 != " "]
# Extract the first vector element and set it as the list element name
names(y1) <- names


#focusing on monocytes data only, transfer log tpm to log fc

#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------


# Load WGCNA and flashClust
library(WGCNA)
library(flashClust)
library(doMC)
library("parallel")
Cl <- makeCluster(detectCores()-1)
source(sprintf("%s/bioinformatics-signatures/load.R", Sys.getenv("GITCLONES")))
#save data
datExpr <- Neutrophils
soft <- NULL
firstWGCNA <- function(datExpr){
  registerDoMC(8)
  datExpr2 <-datExpr
##reverse dataframe so it fit feeds
row.names(datExpr) = datExpr$X
datExpr$X = NULL
datExpr = as.data.frame(t(datExpr)) # now samples are rows and genes are columns
dim(datExpr)

## Run this to check if there are gene outliers
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK 

#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}



# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
soft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
yd <- -sign(soft$fitIndices[,3])*soft$fitIndices[,2]
if(max(yd>0.9)){
index_yd <- which(yd>0.9)
# id <- diff(index_yd)
# if (all(id == c(rep(1,length(id))))  ){
  sft <- soft$fitIndices$Power[index_yd[1]]}else{
    index_yd <- which(yd>0.8)
    sft <- soft$fitIndices$Power[index_yd[1]]
  }
print(paste("The first WGCNA softPower we choose is :", sft))
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
# Mean connectivity as a function of the soft-thresholding power
plot(soft$fitIndices[,1], soft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(soft$fitIndices[,1], soft$fitIndices[,5], labels=powers, cex=cex1,col="red")


## Construct Networks------------------------------------------------------------------------------
#seems that 4 would be a nice sft power. build a adjacency "correlation" matrix with power=4
enableWGCNAThreads()
softPower = sft
adjacency = adjacency(datExpr, power = softPower, type = "signed") #specify network type
head(adjacency)

# Construct Networks
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM

# Generate Modules ------------------------------------------------------------------
# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
#This sets the minimum number of genes to cluster into a module
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=4, pamRespectsDendro= TRUE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = sft)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")


#plots tree showing how the eigengenes cluster together
png(file="clusterw_gene.png")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0.0
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
dev.off()

#plot dendrogram with module colors below it
png(file="cluster_gene.png")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
dev.off()



###########################################################################################
#use WGCNA again to cluster Jingping_cluster, mgsidb, and built in sig
#sort(as.character(stopsle.rtg[cluster_list[["Jingping_7"]]]))

#found that these are 11 clusters in CD16, so we label clusters with name, and change emsbl to symbles
names_clusters <- dynamicMods
genes <- as.vector(rownames(datExpr2))
names(genes) <- names_clusters
cluster_list <- as.list(split(genes,names(genes)))
cluster_list <- lapply(1:length(cluster_list),function(x) unname(cluster_list[[x]]))
cluster_list <- cluster_list[order(sapply(cluster_list,length),decreasing=T)]

#since the order was changed, manully change names
#check order with summary(cluster_list)

a <-  NULL
for (i in 1: length(cluster_list)){
  b <- paste("Jingping_",i,sep = "")
  a <- c(a,b)
}

names(cluster_list) <- a
for( i in 1:length(cluster_list)){
  cluster_list[[i]]  = as.character(stopsle.rtg[cluster_list[[i]]])
}
y2 <- cluster_list

#combine three lists and do score_sig test
whole_sig <- c(y1,y2,signature.list)
result_sig <- scoreSignatures(datExpr2,stopsle.rtg,whole_sig,cluster = Cl)


#remove NAs from result_sig, and seperate healthy and sle to two variables
na_r <- apply(result_sig,2,function(x) any(is.na(x)))
result_sig<- result_sig[,!na_r]
var_x <- result_sig[grep("HC",rownames(result_sig)),]
var_y <- result_sig[grep("SLE",rownames(result_sig)),]

#do a t-test, and fdr, save as new data frame
t_test_results <- lapply(1:ncol(result_sig),function(x,y) t.test(var_x[,x],var_y[,y],alternative="two.sided",mu=0)$p.value)
names(t_test_results) <- colnames(result_sig)
p_fdr <- as.data.frame(t_test_results)
p_fdr <- as.data.frame(t(p_fdr))
colnames(p_fdr) <- "p.value"
p_fdr$fdr <- p.adjust(p_fdr$p.value)
rownames(p_fdr) <- colnames(result_sig)
save(p_fdr,result_sig,file = "p_fdr.rdata")
# rule out sig with p <0.05
significant <- p_fdr[p_fdr$fdr < 0.05,]
significant_score <- t(result_sig[,as.vector(rownames(significant))])



#do WGCNA again
#rename and save one copy for CAMERA manipulation
dat <- as.data.frame(significant_score)
dim(dat)


##reverse dataframe so it fit feeds
row.names(dat) = dat$X
dat$X = NULL
dat = as.data.frame(t(dat)) # now samples are rows and genes are columns
dim(dat)

## Run this to check if there are outliers
gsg = goodSamplesGenes(dat, verbose = 3)
gsg$allOK 

# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
  printFlush(paste("Removing genes:", paste(names(dat)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dat)[!gsg$goodSamples], collapse=", ")))
  dat= dat[gsg$goodSamples, gsg$goodGenes]
}



# Choose a set of soft-thresholding powers
powers = c(seq(from = 1, to=50, by=2))
# Call the network topology analysis function
soft <- NULL
sft <- NULL
registerDoMC(6)
soft = pickSoftThreshold(dat, powerVector = powers, verbose = 5)

#If the lack of scale-free topology fit turns out to be caused by 
#an interesting biological variable that one does not want to remove 
#(i.e., adjust the data for), the appropriate soft-thresholding power can be 
#chosen based on the number of samples as in the table below. 
yd <- -sign(soft$fitIndices[,3])*soft$fitIndices[,2]
    index_yd <- which(yd>0.8)
    sft <- soft$fitIndices$Power[index_yd[1]]
    if(is.na(sft) == TRUE){
      sft = 6
    }
print(paste("The second WGCNA softPower we choose is :", sft))
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(soft$fitIndices[,1], soft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(soft$fitIndices[,1], soft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## Construct Networks------------------------------------------------------------------------------
#build a adjacency "correlation" matrix with power=15

enableWGCNAThreads()
softPower = sft
adjacency = adjacency(dat, power = softPower, type = "unsigned") #specify network type
head(adjacency)

# Construct Networks
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM


# Generate Modules ------------------------------------------------------------------
# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Sig Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
#This sets the minimum number of genes to cluster into a module
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= TRUE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(dat, colors= dynamicColors,softPower = 15)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")


#plots tree showing how the sig cluster together
png(file="clusterw_sig.png")
plot(METree, main= "Clustering of score.sig", xlab= "", sub= "")
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0.0
merge = mergeCloseModules(dat, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
dev.off()

#plot dendrogram with module colors below it
png(file="cluster_sig.png")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
dev.off()

#check out the results
cluster_analysis <- as.data.frame(significant_score)
cluster_analysis$cluster <- dynamicMods
Jingping_cluster_result <- cluster_analysis[grep("Jingping+", rownames(cluster_analysis), perl=TRUE, value=TRUE),]


save(cluster_list,Jingping_cluster_result,cluster_analysis,significant_score, file = "JINGPING_WGCNA.Rdata")
my_list <- c(cluster_list,Jingping_cluster_result,cluster_analysis,significant_score)
return(my_list)
}



##############################################################################################################################
##############################################################################################################################
#
#
#                                               RESULT PRESENTING SECTION
#
#
##############################################################################################################################
##############################################################################################################################



#load the results
setwd("~/Desktop/sle/CD16Monocytes_new/")
#These two command could check gene symbols in a cluster, or check other signatures within same clusters after the second WGCNA.
#for example: we use cluster Jingping_Monocyte16 as an example, this cluster was pre-cutted for gene expression and filtered by fdr <0.05.
#Jingping_Monocytes16 is within the cluster 20 after the second WGCNA.
summary(cluster_list)
#Jingping_cluster_result[Jingping_cluster_result$cluster != 0,41]
assigned <- as.data.frame(Jingping_cluster_result[,41])
colnames(assigned) <- "sig_cluster"
assigned$Jingping_cluster <- rownames(Jingping_cluster_result)
view(assigned)
#BTW, these numbers of Jingping_Monocytes clusters represents the order of size of genes, in a decreasing sorting.
#check sizes of Jingpingclusters
summary(cluster_list)

#cluster_analysis$cluster is the cluster number after second WGCNA

#This is a list of gene symbol in Jingping_Monocytes16
cluster_list$Jingping_11 #3
rownames(cluster_analysis[cluster_analysis$cluster == 
                            "3",])

cluster_list$Jingping_13 #7
rownames(cluster_analysis[cluster_analysis$cluster == "7",])

cluster_list$Jingping_16 #9
rownames(cluster_analysis[cluster_analysis$cluster == "9",])

cluster_list$Jingping_19 #7
rownames(cluster_analysis[cluster_analysis$cluster == "7",])

cluster_list$Jingping_20 #23
rownames(cluster_analysis[cluster_analysis$cluster == "23",])

cluster_list$Jingping_25 #7
rownames(cluster_analysis[cluster_analysis$cluster == "7",])

cluster_list$Jingping_26 #8
rownames(cluster_analysis[cluster_analysis$cluster == "8",])

cluster_list$Jingping_27 #3
rownames(cluster_analysis[cluster_analysis$cluster == "3",])

cluster_list$Jingping_28 #27
rownames(cluster_analysis[cluster_analysis$cluster == "27",])

cluster_list$Jingping_29 #27
rownames(cluster_analysis[cluster_analysis$cluster == "27",])

cluster_list$Jingping_31 #9
rownames(cluster_analysis[cluster_analysis$cluster == "9",])



#T8
cluster_list$Jingping_33 #7
rownames(cluster_analysis[cluster_analysis$cluster == "7",])

cluster_list$Jingping_22 #5
rownames(cluster_analysis[cluster_analysis$cluster == "5",])


cluster_list$Jingping_24  #3
rownames(cluster_analysis[cluster_analysis$cluster == "3",])

cluster_list$Jingping_25  #5
rownames(cluster_analysis[cluster_analysis$cluster == "5",])

###################MONOCYTES
t_data <- t(significant_score)
# idx <- grep("Jingping+", colnames(Jingping24),perl = TRUE)
# #move the cluster to the bottom
# Jingping24 <- Jingping24[, c((1:ncol(Jingping24))[-idx],idx)]
# save_names <- colnames(Jingping24)
# colnames(Jingping24) <- c(1:(ncol(Jingping24)))
# colnames(Jingping24)[269] <- "J24"
# Jingping24[,269]
#REACTOME_GPCR_LIGAND_BINDING (drug target)
#GO_CALCIUM_ION_TRANSPORT
#GO_CYTOSOLIC_CALCIUM_ION_TRANSPORT
#OKUMURA_INFLAMMATORY_RESPONSE_LPS
#REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL
#TRANSMEMBRANE
#GO_C_C_CHEMOKINE_RECEPTOR_ACTIVITY
#RENAL_ABSORPTION
Jingping11 <- as.data.frame(t_data[,rownames(cluster_analysis[cluster_analysis$cluster == "3",])])
c_Jingping11 <- cor(Jingping11)
c_Jingping11 <- data.frame(c_Jingping11[,grep("Jingping+", colnames(c_Jingping11),perl = TRUE,value = TRUE)])
as.vector(rownames(c_Jingping11))


#APOPTOSIS
# [97] "GSE10325_BCELL_VS_LUPUS_BCELL_UP :"                                                    
# "GSE22886_NEUTROPHIL_VS_MONOCYTE_UP :"  
#GO_REGULATION_OF_MONONUCLEAR_CELL_MIGRATION 
#GSE34515_CD16_NEG_MONOCYTE_VS_DC_UP :"
Jingping13 <- as.data.frame(t_data[,rownames(cluster_analysis[cluster_analysis$cluster == "7",])])
c_Jingping13 <- cor(Jingping13)
c_Jingping13 <- data.frame(c_Jingping13[,grep("Jingping+", colnames(c_Jingping13),perl = TRUE,value = TRUE)])
as.vector(rownames(c_Jingping13))

#GSE22886_DC_VS_MONOCYTE_UP 
#GO_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE
#GSE21927_HEALTHY_VS_TUMOROUS_BALBC_MOUSE_MONOCYTE_DN
#GO_REGULATION_OF_MONOCYTE_DIFFERENTIATION
#GO_REGULATION_OF_INFLAMMATORY_RESPONSE
Jingping16 <- as.data.frame(t_data[,rownames(cluster_analysis[cluster_analysis$cluster == "9",])])
c_Jingping16 <- cor(Jingping16)
c_Jingping16 <- data.frame(c_Jingping16[,grep("Jingping+", colnames(c_Jingping16.31),perl = TRUE,value = TRUE)])
as.vector(rownames(c_Jingping16))


#mitotic, NUMA, dna damage
Jingping26 <- as.data.frame(t_data[,rownames(cluster_analysis[cluster_analysis$cluster == "8",])])
c_Jingping26 <- cor(Jingping26)
c_Jingping26 <- data.frame(c_Jingping26[,grep("Jingping+", colnames(c_Jingping26),perl = TRUE,value = TRUE)])
colnames(c_Jingping26) <- "Jingping_26"
as.vector(rownames(c_Jingping26))



#GSE10325_MYELOID_VS_LUPUS_MYELOID_UP
##?no
Jingping18 <- as.data.frame(t_data[,rownames(cluster_analysis[cluster_analysis$cluster == "0",])])
c_Jingping18 <- cor(Jingping18)
c_Jingping18 <- data.frame(c_Jingping18[,grep("Jingping+", colnames(c_Jingping18),perl = TRUE,value = TRUE)])
as.vector(rownames(c_Jingping18))

#B T MONO differe
Jingping20 <- as.data.frame(t_data[,rownames(cluster_analysis[cluster_analysis$cluster == "23",])])
c_Jingping20 <- cor(Jingping20)
c_Jingping20 <- data.frame(c_Jingping20[,grep("Jingping+", colnames(c_Jingping20),perl = TRUE,value = TRUE)])
as.vector(rownames(c_Jingping20))

#BOHN_PRIMARY_IMMUNODEFICIENCY_SYNDR
Jingping28 <- as.data.frame(t_data[,rownames(cluster_analysis[cluster_analysis$cluster == "27",])])
c_Jingping28 <- cor(Jingping28)
c_Jingping28 <- data.frame(c_Jingping28[,grep("Jingping+", colnames(c_Jingping28),perl = TRUE,value = TRUE)])
as.vector(rownames(c_Jingping28))