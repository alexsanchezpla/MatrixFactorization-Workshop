## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------------------
require(knitr)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, comment = NA, prompt = TRUE, tidy = FALSE, fig.width = 7, fig.height = 7, fig_caption = TRUE,cache=FALSE)
Sys.setlocale("LC_TIME", "C")


## ----packages, include=FALSE------------------------------------------------------------------------------------------------------------------------------
# require(devtools)
# if(!require(installifnot)) install_github("uebvhir/installifnot")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
library(dbplyr)
outline <- c("Workshop presentation and Overview", "10m",
"Exploring single omics datasets", "20m",
"Building a SummarizedExperiment from scratch",	"10m",
"TCGA multi-assay dataset", "10m",
"Multiple Factor Analysis and iCluster", "30m",
"Enrichment Analysis and Factor interpretation",	"20m")

timeTable<- matrix(outline, ncol=2, byrow = TRUE)
colnames(timeTable) <- c("Activity", "Time")

kableExtra::kable(timeTable) %>% kableExtra::kable_styling()


## ----prepareLeukemiaSubset, eval=FALSE--------------------------------------------------------------------------------------------------------------------
## library(stringi)
## myDir <- getwd()
## myDir <-unlist(stri_split_regex(myDir, "/"))
## currDir <- myDir[length(myDir)]
## if (currDir=="202107-GRBio-Integration-Workshop") setwd("lab-Matrix_factorization")
## 
## if (!exists("datasets/leukemiaExpressionSubset.rds")){
##   library(leukemiasEset)
##   data("leukemiasEset")
##   library(SummarizedExperiment)
##   leukemiaSummExp <-as(leukemiasEset, "SummarizedExperiment")
##   mat0=assay(leukemiaSummExp)
##   mat1=(limma::normalizeQuantiles(mat0))
##   mat2=mat1[order(matrixStats::rowSds(mat1)/rowMeans(mat1)),][1:1000,]
##   saveRDS(mat2,"datasets/leukemiaExpressionSubset.rds")
##   saveRDS(mat1,"datasets/leukemiaExpression.rds")
## }else{
##   mat=readRDS("datasets/leukemiaExpressionSubset.rds")
## }


## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------
## expFile=system.file("extdata","leukemiaExpressionSubset.rds",
##                     package="compGenomRData")
## mat=readRDS(expFile)


## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------
## leukemiaSummExp
## dim(assay(leukemiaSummExp))
## metadata(leukemiaSummExp)
## colData(leukemiaSummExp)
## dim(colData(leukemiaSummExp))


## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------
## leukemiaSub <-leukemiaSummExp[,1:10]
## dim(assay(leukemiaSub))
## dim(colData(leukemiaSub))


## ----exploreLeukemia--------------------------------------------------------------------------------------------------------------------------------------
library(stringi)
myDir <- getwd() 
myDir <-unlist(stri_split_regex(myDir, "/"))
currDir <- myDir[length(myDir)]
if (currDir=="202107-GRBio-Integration-Workshop") setwd("lab-Matrix_factorization")


library(pheatmap)
expFile="datasets/leukemiaExpressionSubset.rds"
mat=readRDS(expFile)

# set the leukemia type annotation for each sample
# this has to be a dta frame
annotation_col = data.frame(
                    LeukemiaType =substr(colnames(mat),1,3))
rownames(annotation_col)=colnames(mat)
pheatmap(mat,show_rownames=FALSE,show_colnames=FALSE,
         annotation_col=annotation_col,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean")


## ----scatterb4PCA-----------------------------------------------------------------------------------------------------------------------------------------
plot(mat[rownames(mat)=="ENSG00000100504",],
     mat[rownames(mat)=="ENSG00000105383",],pch=19,
     ylab="CD33 (ENSG00000105383)",
     xlab="PYGL (ENSG00000100504)")


## ----PCArot-----------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
# create the subset of the data with two genes only
#
sub.mat=t(mat[rownames(mat) %in% c("ENSG00000100504","ENSG00000105383"),])
# ploting our genes of interest as scatter plots
plot(scale(mat[rownames(mat)=="ENSG00000100504",]),
     scale(mat[rownames(mat)=="ENSG00000105383",]),
     pch=19,
     ylab="CD33 (ENSG00000105383)",
     xlab="PYGL (ENSG00000100504)",
     col=as.factor(annotation_col$LeukemiaType),
     xlim=c(-2,2),ylim=c(-2,2))
# create the legend for the Leukemia types
legend("bottomright",
       legend=unique(annotation_col$LeukemiaType),
       fill =palette("default"),
       border=NA,box.col=NA)
# calculate the PCA only for our genes and all the samples
pr=princomp(scale(sub.mat))
# plot the direction of eigenvectors
# pr$loadings returned by princomp has the eigenvectors
arrows(x0=0, y0=0, x1 = pr$loadings[1,1], 
         y1 = pr$loadings[2,1],col="pink",lwd=3)
arrows(x0=0, y0=0, x1 = pr$loadings[1,2], 
         y1 = pr$loadings[2,2],col="gray",lwd=3)
# plot the samples in the new coordinate system
plot(-pr$scores,pch=19,
     col=as.factor(annotation_col$LeukemiaType),
     ylim=c(-2,2),xlim=c(-4,4))
# plot the new coordinate basis vectors
arrows(x0=0, y0=0, x1 =-2, 
         y1 = 0,col="pink",lwd=3)
arrows(x0=0, y0=0, x1 = 0, 
         y1 = -1,col="gray",lwd=3)
legend("bottomright",
       legend=unique(annotation_col$LeukemiaType),
       fill =palette("default"),
       border=NA,box.col=NA)



## ----eigenOnCovMat----------------------------------------------------------------------------------------------------------------------------------------
## ----eigenOnCovMat,eval=FALSE--------------------------------------------------------------------------------
## cov.mat=cov(sub.mat) # calculate covariance matrix
## cov.mat
## eigen(cov.mat) # obtain eigen decomposition for eigen values and vectors


## ----SVDcartoon, echo=FALSE, fig.align='center',out.width='60%',fig.cap="Singular value decomposition (SVD) explained in a diagram"-----------------------
knitr::include_graphics("images/SVDcartoon.png")


## ----svd,out.width='65%',fig.width=8.5,fig.cap="SVD on the matrix and its transpose"----------------------------------------------------------------------
par(mfrow=c(1,2))
d=svd(scale(mat)) # apply SVD
assays=t(d$u) %*% scale(mat) # projection on eigenassays
plot(assays[1,],assays[2,],pch=19,
     col=as.factor(annotation_col$LeukemiaType))
#plot(d$v[,1],d$v[,2],pch=19,
#     col=annotation_col$LeukemiaType)
pr=prcomp(t(mat),center=TRUE,scale=TRUE) # apply PCA on transposed matrix
# plot new coordinates from PCA, projections on eigenvectors
# since the matrix is transposed eigenvectors represent 
plot(pr$x[,1],pr$x[,2],col=as.factor(annotation_col$LeukemiaType))



## ----out.width='70%',fig.cap= "Singular value decomposition (SVD) reorganized as multiplication of m-by-n weights matrix and eigenvectors"----------------
knitr::include_graphics("images/SVDasWeights.png")


## ---- fig.cap="Gene expression of a gene can be regarded as a linear combination of eigenvectors. "-------------------------------------------------------
knitr::include_graphics("images/SVDlatentExample.png")


## ---- fig.cap="Independent Component Analysis (ICA)"------------------------------------------------------------------------------------------------------
knitr::include_graphics("images/ICAcartoon.png")


## ---- out.width='50%',fig.width=5,fig.cap="Leukemia gene expression values per patient on reduced dimensions by ICA."-------------------------------------
## ----fastICAex, ----
library(fastICA)
ica.res=fastICA(t(mat),n.comp=2) # apply ICA
# plot reduced dimensions
plot(ica.res$S[,1],ica.res$S[,2],col=as.factor(annotation_col$LeukemiaType))


## ---- fig.cap="Non-negative matrix factorization summary",out.width='70%'---------------------------------------------------------------------------------
knitr::include_graphics("images/NMFcartoon.png")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
library(NMF)
res=NMF::nmf(mat,rank=3,seed="nndsvd") # nmf with 3 components/factors
w <- basis(res) # get W
h <- coef(res)  # get H
# plot 1st factor against 3rd factor
plot(h[1,],h[3,],col=as.factor(annotation_col$LeukemiaType),pch=19)


## __Want to know more ?__

## The NMF package vignette has extensive information on how to run NMF to get stable results and an estimate of components: https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf



## ----MDS--------------------------------------------------------------------------------------------------------------------------------------------------
mds=cmdscale(dist(t(mat)))
isomds=MASS::isoMDS(dist(t(mat)))
# plot the patients in the 2D space
par(mfrow=c(1,2))
plot(mds,pch=19,col=as.factor(annotation_col$LeukemiaType),
     main="classical MDS")
plot(isomds$points,pch=19,col=as.factor(annotation_col$LeukemiaType),
     main="isotonic MDS")


## ----tSNE-------------------------------------------------------------------------------------------------------------------------------------------------
## ----tsne,eval=TRUE, out.width='60%',fig.width=5, fig.cap="t-SNE of leukemia expression dataset"-------------
library("Rtsne")
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(t(mat),perplexity = 10) # Run TSNE
 #image(t(as.matrix(dist(tsne_out$Y))))
# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col=as.factor(annotation_col$LeukemiaType),
     pch=19)
# create the legend for the Leukemia types
legend("bottomleft",
       legend=unique(annotation_col$LeukemiaType),
       fill =palette("default"),
       border=NA,box.col=NA)


## __Want to know more ?__

## - How perplexity affects t-sne, interactive examples:  https://distill.pub/2016/misread-tsne/

## - More on perplexity: https://blog.paperspace.com/dimension-reduction-with-t-sne/

## - Intro to t-SNE: https://www.oreilly.com/learning/an-illustrated-introduction-to-the-t-sne-algorithm



## ---------------------------------------------------------------------------------------------------------------------------------------------------------
library(stringi)
myDir <- getwd() 
myDir <-unlist(stri_split_regex(myDir, "/"))
currDir <- myDir[length(myDir)]
if (currDir=="202107-GRBio-Integration-Workshop") setwd("lab-Matrix_factorization")

library(magrittr)
# read in the csv from the multi-omics-folder
csvfile <- "datasets/multi-omics/COREAD_CMS13_gex.csv" 
x1 <- read.csv(csvfile, row.names=1)
# Fix the gene names in the data frame
rownames(x1) <- sapply(strsplit(rownames(x1), "\\|"), function(x) x[1])
# Output a table
knitr::kable(head(t(head(x1))), caption="Example gene expression data (head)") %>% kableExtra::kable_styling() 


## ------------------------------------------------------------------------------------------------------------
csvfile <- "datasets/multi-omics/COREAD_CMS13_muts.csv"
x2 <- read.csv(csvfile, row.names=1)
# Set mutation data to be binary (so if a gene has more than 1 mutation,
# we only count one)
x2[x2>0]=1
# output a table
knitr::kable(head(t(head(x2))), caption="Example mutation data (head)") %>% kableExtra::kable_styling()


## ------------------------------------------------------------------------------------------------------------
# read in the csv from the companion package as a data frame
csvfile <- "datasets/multi-omics/COREAD_CMS13_cnv.csv"
x3 <- read.csv(csvfile, row.names=1)
# output a table
knitr::kable(head(t(head(x3))), 
             caption="Example copy number data for CRC samples") %>% kableExtra::kable_styling()


## ------------------------------------------------------------------------------------------------------------
# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_subtypes.csv",
                       package="compGenomRData")
covariates <- read.csv(csvfile, row.names=1)
# Fix the TCGA identifiers so they match up with the omics data
rownames(covariates) <- gsub(pattern = '-', replacement = '\\.',
                             rownames(covariates))
covariates <- covariates[colnames(x1),]
# create a dataframe which will be used to annotate later graphs
anno_col <- data.frame(cms=as.factor(covariates$cms_label))
rownames(anno_col) <- rownames(covariates)
# output a table
knitr::kable(head(anno_col), 
             caption="Clinical information (covariates)")



## ---------------------------------------------------------------------------------------------------------------------------------------------------------
pheatmap::pheatmap(x1,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Gene expression data")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------
pheatmap::pheatmap(x2,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Mutation data")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------
pheatmap::pheatmap(x3,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Copy number data")



## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------
# run the MFA function from the FactoMineR package
r.mfa <- FactoMineR::MFA(
  t(rbind(x1,x2,x3)), # binding the omics types together
  c(dim(x1)[1], dim(x2)[1], dim(x3)[1]), # specifying the dimensions of each
  graph=FALSE)


## ------------------------------------------------------------------------------------------------------------
# first, extract the H and W matrices from the MFA run result
mfa.h <- r.mfa$global.pca$ind$coord
mfa.w <- r.mfa$quanti.var$coord

# create a dataframe with the H matrix and the CMS label
mfa_df <- as.data.frame(mfa.h)
mfa_df$subtype <- factor(covariates[rownames(mfa_df),]$cms_label)

# create the plot
ggplot2::ggplot(mfa_df, ggplot2::aes(x=Dim.1, y=Dim.2, color=subtype)) +
ggplot2::geom_point() + ggplot2::ggtitle("Scatter plot of MFA")


## ------------------------------------------------------------------------------------------------------------
pheatmap::pheatmap(t(mfa.h)[1:2,], annotation_col = anno_col,
                  show_colnames = FALSE,
                  main="MFA for multi-omics integration")




## ----jointNMF---------------------------------------------------------------------------------------------------------------------------------------------
## ----warning=FALSE-------------------------------------------------------------------------------------------
# Feature-normalize the data
x1.featnorm <- x1 / rowSums(x1)
x2.featnorm <- x2 / rowSums(x2)
x3.featnorm <- x3 / rowSums(x3)

# Normalize by each omics type's frobenius norm
matExpr<-x1.featnorm.frobnorm <- x1.featnorm / norm(as.matrix(x1.featnorm), type="F")
matSNP<-x2.featnorm.frobnorm <- x2.featnorm / norm(as.matrix(x2.featnorm), type="F")
matCNV<-x3.featnorm.frobnorm <- x3.featnorm / norm(as.matrix(x3.featnorm), type="F")

# Split the features of the CNV matrix into two non-negative features each
split_neg_columns <- function(df) {
  n <- dim(df)[1]
  k <- dim(df)[2]
  df2 <- matrix(rep(0, n*2*k), ncol=2*k)
  for (i in 1:k){
    df2[,2*i-1] <- pmax(df[,i],0)
    df2[,2*i]   <- pmax(-df[,i], 0)
  }
  as.data.frame(df2)
}
matCNVpos<- x3.featnorm.frobnorm.nonneg <- t(split_neg_columns(t(x3.featnorm.frobnorm)))
colnames(matCNVpos) <- colnames(x3.featnorm.frobnorm.nonneg) <- colnames(x3.featnorm.frobnorm)

# run the nmf function from the NMF package
require(NMF)
jointNames <-data.frame(exprNames = colnames(matExpr),
                        snpNames = colnames(matSNP),
                        cnvNames = colnames(matCNVpos)
                        )
jointMat <- rbind(matExpr, matSNP, matCNVpos)
                     
r.nmf <- nmf(t(jointMat),
             2,
             method='Frobenius')


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
nmf.h <- NMF::basis(r.nmf)
nmf.w <- NMF::coef(r.nmf)
nmfw <- t(nmf.w)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
nmf_df <- as.data.frame(nmf.h)
colnames(nmf_df) <- c("dim1", "dim2")
nmf_df$subtype <- factor(covariates[rownames(nmf_df),]$cms_label)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
ggplot2::ggplot(nmf_df, ggplot2::aes(x=dim1, y=dim2, color=subtype)) +
ggplot2::geom_point() +
ggplot2::ggtitle("Scatter plot of 2-component NMF for multi-omics integration")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
pheatmap::pheatmap(t(nmf_df[,1:2]),
                   annotation_col = anno_col,
                   show_colnames=FALSE,
                   main="Heatmap of 2-component NMF")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::include_graphics("images/icluster.png" )


## ----integrativeClustering--------------------------------------------------------------------------------------------------------------------------------
library(iClusterPlus)
# run the iClusterPlus function
r.icluster <- iClusterPlus::iClusterPlus(
  t(x1), # Providing each omics type
  t(x2),
  t(x3),
  type=c("gaussian", "binomial", "multinomial"), # Providing the distributions
  K=2, # provide the number of factors to learn
  alpha=c(1,1,1), # as well as other model parameters
  lambda=c(.03,.03,.03)
  )
# extract the H and W matrices from the run result
# here, we refer to H as z, to keep with iCluster terminology
icluster.z <- r.icluster$meanZ
rownames(icluster.z) <- rownames(covariates) # fix the row names
icluster.ws <- r.icluster$beta
# construct a dataframe with the H matrix (z) and the cancer subtypes
# for later plotting
icp_df <- as.data.frame(icluster.z)
colnames(icp_df) <- c("dim1", "dim2")
rownames(icp_df) <- colnames(x1)
icp_df$subtype <- factor(covariates[rownames(icp_df),]$cms_label)


## ---- moiclusterplusscatter,fig.cap="iCluster+ learns factors which allow tumor sub-types CMS1 and CMS3 to be discriminated.", echo=FALSE-----------------

ggplot2::ggplot(icp_df, ggplot2::aes(x=dim1, y=dim2, color=subtype)) +
ggplot2::geom_point() +
ggplot2::ggtitle("Scatter plot of iCluster+ factors")


## ----iclusterFactors--------------------------------------------------------------------------------------------------------------------------------------
pheatmap::pheatmap(t(icp_df[,1:2]),
                   annotation_col = anno_col, 
                   show_colnames = FALSE,border_color = NA,
                   main="Heatmap of iCluster+ factors")



## ---------------------------------------------------------------------------------------------------------------------------------------------------------
nmf.clusters <- max.col(nmf.h)
names(nmf.clusters) <- rownames(nmf.h)

# create an annotation data frame indicating the NMF one-hot clusters
# as well as the cancer subtypes, for the heatmap plot below
anno_nmf_cl <- data.frame(
  nmf.cluster=factor(nmf.clusters),
  cms.subtype=factor(covariates[rownames(nmf.h),]$cms_label)
)

# generate the plot
pheatmap::pheatmap(t(nmf.h[order(nmf.clusters),]),
  cluster_cols=FALSE, cluster_rows=FALSE,
  annotation_col = anno_nmf_cl,
  show_colnames = FALSE,border_color=NA,
  main="Joint NMF factors\nwith clusters and molecular subtypes")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
icluster.clusters <- kmeans(icluster.z, 2)$cluster
names(icluster.clusters) <- rownames(icluster.z)

# create an annotation dataframe for the heatmap plot
# containing the kmeans cluster assignments and the cancer subtypes
anno_icluster_cl <- data.frame(
  iCluster=factor(icluster.clusters),
  cms.subtype=factor(covariates$cms_label))

# generate the figure
pheatmap::pheatmap(
  t(icluster.z[order(icluster.clusters),]), # order z by the kmeans clusters
  cluster_cols=FALSE, # use cluster_cols and cluster_rows=FALSE
  cluster_rows=FALSE, # as we want the ordering by k-means clusters to hold
  show_colnames = FALSE,border_color=NA,
  annotation_col = anno_icluster_cl,
  main="iCluster factors\nwith clusters and molecular subtypes")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
# create an annotation dataframe for the heatmap
# for each feature, indicating its omics-type
data_anno <- data.frame(
  omics=c(rep('expression',dim(x1)[1]),
          rep('mut',dim(x2)[1]),
          rep('cnv',dim(x3.featnorm.frobnorm.nonneg)[1])))
rownames(data_anno) <- c(rownames(x1),
                         paste0("mut:", rownames(x2)),
                         rownames(x3.featnorm.frobnorm.nonneg))
rownames(nmfw) <- rownames(data_anno)
# generate the heat map
pheatmap::pheatmap(nmfw,
                   cluster_cols = FALSE,
                   annotation_row = data_anno,
                   main="NMF coefficients",
                   clustering_distance_rows = "manhattan",
                   fontsize_row = 1)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
# select genes associated preferentially with each factor
# by their relative loading in the W matrix
library(enrichR)
genes.factor.1 <- names(which(nmfw[1:dim(x1)[1],1] > nmfw[1:dim(x1)[1],2]))
genes.factor.2 <- names(which(nmfw[1:dim(x1)[1],1] < nmfw[1:dim(x1)[1],2]))
# call the enrichr function to find gene sets enriched
# in each latent factor in the GO Biological Processes 2018 library
go.factor.1 <- enrichr(genes.factor.1,
                                databases = c("GO_Biological_Process_2018")
                                )$GO_Biological_Process_2018
go.factor.2 <- enrichr(genes.factor.2,
                                databases = c("GO_Biological_Process_2018")
                                )$GO_Biological_Process_2018



## ---------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
go.factor.2$Genes <- gsub(";", "; ", go.factor.2$Genes)
the.table <- knitr::kable(head(go.factor.2, 3)[,c("Term", "Adjusted.P.value", "Combined.Score")],
                 caption="GO-terms associated with NMF factor 2",
                 format="latex")
#the.table <- kableExtra::column_spec(the.table, 1, width="10em")
the.table <- kableExtra::kable_styling(the.table ,latex_options = c( "scale_down"))
#the.table <- kableExtra::column_spec(the.table, 4, width="10em")
the.table


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
# create a data frame holding covariates (age, gender, MSI status)
a <- data.frame(age=covariates$age,
                gender=as.numeric(covariates$gender),
                msi=covariates$msi)
b <- nmf.h
colnames(b) <- c('factor1', 'factor2')
# concatenate the covariate dataframe with the H matrix
cov_factor <- cbind(a,b)
# generate the figure
ggplot2::ggplot(cov_factor, ggplot2::aes(x=msi, y=factor1, group=msi)) +
  ggplot2::geom_boxplot() +
  ggplot2::ggtitle("NMF factor 1 microsatellite instability")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
ggplot2::ggplot(cov_factor, ggplot2::aes(x=msi, y=factor2, group=msi)) +
  ggplot2::geom_boxplot() +
  ggplot2::ggtitle("NMF factor 2 and microsatellite instability")


