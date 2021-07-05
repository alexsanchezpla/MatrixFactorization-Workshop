## ----setup, include=FALSE------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------------------------------------------
# install.packages("fastICA", "NMF", "Rtsne", dep=TRUE)


## ------------------------------------------------------------------------------------------------------------
library (compGenomRData)
library(pheatmap)
expFile=system.file("extdata","leukemiaExpressionSubset.rds",
                    package="compGenomRData")
mat=readRDS(expFile)
dim(mat)
# canType<- substr(colnames(mat),1,3)
# sampleNum <- substr(colnames(mat),11,13)
# colnames(mat)<- paste(canType,sampleNum, sep="_")
# table(canType)
# head(mat)
# # set the leukemia type annotation for each sample
annotation_col = data.frame(
                    LeukemiaType =substr(colnames(mat),1,3))
rownames(annotation_col)=colnames(mat)


## ----scatterb4PCA,out.width='60%',fig.width=5.5, fig.cap="Gene expression values of CD33 and PYGL genes across leukemia patients."----
plot(mat[rownames(mat)=="ENSG00000100504",],
     mat[rownames(mat)=="ENSG00000105383",],pch=19,
     ylab="CD33 (ENSG00000105383)",
     xlab="PYGL (ENSG00000100504)")


## ----pcaRot,out.width='60%',fig.width=8.5,fig.cap="Geometric interpretation of PCA finding eigenvectors that point to the direction of highest variance. Eigenvectors can be used as a new coordinate system."----
par(mfrow=c(1,2))
# create the subset of the data with two genes only
# notice that we transpose the matrix so samples are 
# on the columns
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


## ----eigenOnCovMat,eval=FALSE--------------------------------------------------------------------------------
## cov.mat=cov(sub.mat) # calculate covariance matrix
## cov.mat
## eigen(cov.mat) # obtain eigen decomposition for eigen values and vectors


## ----SVDcartoon,echo=FALSE,fig.align='center',out.width='60%',fig.cap="Singular value decomposition (SVD) explained in a diagram. "----
knitr::include_graphics("images/SVDcartoon.png")


## ----svd,out.width='65%',fig.width=8.5,fig.cap="SVD on the matrix and its transpose"-------------------------
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


## ----SVDasWeigths,echo=FALSE,out.width='70%',fig.cap="Singular value decomposition (SVD) reorganized as multiplication of m-by-n weights matrix and eigenvectors "----
knitr::include_graphics("images/SVDasWeights.png")


## ----SVDlatentExample,echo=FALSE,fig.cap="Gene expression of a gene can be regarded as a linear combination of eigenvectors. "----
knitr::include_graphics("images/SVDlatentExample.png")


## ----ICAcartoon,echo=FALSE,fig.cap="Independent Component Analysis (ICA)"------------------------------------
knitr::include_graphics("images/ICAcartoon.png")


## ----fastICAex, out.width='50%',fig.width=5,fig.cap="Leukemia gene expression values per patient on reduced dimensions by ICA."----
library(fastICA)
ica.res=fastICA(t(mat),n.comp=2) # apply ICA
# plot reduced dimensions
plot(ica.res$S[,1],ica.res$S[,2],col=as.factor(annotation_col$LeukemiaType))


## ----NMFcartoon,echo=FALSE,fig.cap="Non-negative matrix factorization summary",out.width='70%'---------------
knitr::include_graphics("images/NMFcartoon.png")


## ----nmfCode,out.width='60%',fig.width=5,fig.cap="Leukemia gene expression values per patient on reduced dimensions by NMF. Components 1 and 3 is used for the plot."----
library(NMF)
res=NMF::nmf(mat,rank=3,seed="nndsvd") # nmf with 3 components/factors
w <- basis(res) # get W
h <- coef(res)  # get H
# plot 1st factor against 3rd factor
plot(h[1,],h[3,],col=as.factor(annotation_col$LeukemiaType),pch=19)


## __Want to know more ?__

## The NMF package vignette has extensive information on how to run NMF to get stable results and an estimate of components: https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf


## ----mds2,out.width='60%',fig.width=8.5,fig.cap="Leukemia gene expression values per patient on reduced dimensions by classical MDS and isometric MDS."----
mds=cmdscale(dist(t(mat)))
isomds=MASS::isoMDS(dist(t(mat)))
# plot the patients in the 2D space
par(mfrow=c(1,2))
plot(mds,pch=19,col=as.factor(annotation_col$LeukemiaType),
     main="classical MDS")
plot(isomds$points,pch=19,col=as.factor(annotation_col$LeukemiaType),
     main="isotonic MDS")


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


## ----dataLoadClu,eval=FALSE----------------------------------------------------------------------------------
## expFile=system.file("extdata",
##                     "leukemiaExpressionSubset.rds",
##                     package="compGenomRData")
## mat=readRDS(expFile)


## ------------------------------------------------------------------------------------------------------------
if (!(require(compGenomRData))) {
devtools::install_github("compgenomr/compGenomRData")
}
# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_gex.csv", 
                       package="compGenomRData")
x1 <- read.csv(csvfile, row.names=1)
# Fix the gene names in the data frame
rownames(x1) <- sapply(strsplit(rownames(x1), "\\|"), function(x) x[1])
# Output a table
knitr::kable(head(t(head(x1))), caption="Example gene expression data (head)")


## ------------------------------------------------------------------------------------------------------------
# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_muts.csv",
                       package="compGenomRData")
x2 <- read.csv(csvfile, row.names=1)
# Set mutation data to be binary (so if a gene has more than 1 mutation,
# we only count one)
x2[x2>0]=1
# output a table
knitr::kable(head(t(head(x2))), caption="Example mutation data (head)")


## ------------------------------------------------------------------------------------------------------------
# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_cnv.csv", 
                       package="compGenomRData")
x3 <- read.csv(csvfile, row.names=1)
# output a table
knitr::kable(head(t(head(x3))), 
             caption="Example copy number data for CRC samples")


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


## ------------------------------------------------------------------------------------------------------------
pheatmap::pheatmap(x1,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Gene expression data")


## ------------------------------------------------------------------------------------------------------------
pheatmap::pheatmap(x2,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Mutation data")


## ------------------------------------------------------------------------------------------------------------
pheatmap::pheatmap(x3,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Copy number data")


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

# exctract the H and W matrices from the nmf run result
nmf.h <- NMF::basis(r.nmf)
nmf.w <- NMF::coef(r.nmf)
nmfw <- t(nmf.w)


## ------------------------------------------------------------------------------------------------------------
# create a dataframe with the H matrix and the CMS label (subtype)
nmf_df <- as.data.frame(nmf.h)
colnames(nmf_df) <- c("dim1", "dim2")
nmf_df$subtype <- factor(covariates[rownames(nmf_df),]$cms_label)

# create the scatter plot
ggplot2::ggplot(nmf_df, ggplot2::aes(x=dim1, y=dim2, color=subtype)) +
ggplot2::geom_point() +
ggplot2::ggtitle("Scatter plot of 2-component NMF for multi-omics integration")


## ------------------------------------------------------------------------------------------------------------
pheatmap::pheatmap(t(nmf_df[,1:2]),
                   annotation_col = anno_col,
                   show_colnames=FALSE,
                   main="Heatmap of 2-component NMF")


## ----moiCluster,fig.cap="Sketch of iCluster model. Each omics datatype is decomposed to a coefficient matrix and a shared latent variable matrix, plus noise.",fig.align = 'center',out.width='75%',echo=FALSE----
knitr::include_graphics("images/icluster.png" )


## ------------------------------------------------------------------------------------------------------------
if(!require(iClusterPlus)) BiocManager::install("iClusterPlus")


## ----momultiOmicsiclusterplus--------------------------------------------------------------------------------
# run the iClusterPlus function
r.icluster <- iClusterPlus::iClusterPlus(
  t(x1), # Providing each omics type
  t(x2),
  t(x3),
  type=c("gaussian", "binomial", "multinomial"), # Providing the distributions
  K=2, # provide the number of factors to learn
  alpha=c(1,1,1), # as well as other model parameters
  lambda=c(.03,.03,.03))
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


## ----moiclusterplusscatter,fig.cap="iCluster+ learns factors which allow tumor sub-types CMS1 and CMS3 to be discriminated.", echo=FALSE----
ggplot2::ggplot(icp_df, ggplot2::aes(x=dim1, y=dim2, color=subtype)) +
ggplot2::geom_point() +
ggplot2::ggtitle("Scatter plot of iCluster+ factors")


## ----moiclusterplusheatmap,fig.cap="iCluster+ factors, shown in a heatmap, separate tumors into their subtypes well.", echo=FALSE,fig.height=3----
pheatmap::pheatmap(t(icp_df[,1:2]),
                   annotation_col = anno_col, 
                   show_colnames = FALSE,border_color = NA,
                   main="Heatmap of iCluster+ factors")


## __Want to know more ?__

## - Read the original iCluster paper: Shen R., Olshen A. B., Ladanyi M. (2009). Integrative clustering of multiple genomic data types using a joint latent variable model with application to breast and lung cancer subtype analysis. _Bioinformatics_ 25, 2906–2912. 10.1093/bioinformatics/btp543 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2800366/

## - Read the original iClusterPlus paper: an extension of iCluster: Shen R., Mo Q., Schultz N., Seshan V. E., Olshen A. B., Huse J., et al. (2012). Integrative subtype discovery in glioblastoma using iCluster. _PLoS ONE_ 7:e35236. 10.1371/journal.pone.0035236 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3335101/

## - Learn more about the LASSO for model regularization: Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. _J. Royal. Statist. Soc B._, Vol. 58, No. 1, pages 267-288: http://www-stat.stanford.edu/%7Etibs/lasso/lasso.pdf

## - Learn more about the EM algorithm: Dempster, A. P., et al. Maximum likelihood from incomplete data via the EM algorithm. _Journal of the Royal Statistical Society. Series B (Methodological)_, vol. 39, no. 1, 1977, pp. 1–38. JSTOR, JSTOR: http://www.jstor.org/stable/2984875

## - Read about MCMC algorithms: Hastings, W.K. (1970). Monte Carlo sampling methods using Markov chains and their applications. _Biometrika._ 57 (1): 97–109. doi:10.1093/biomet/57.1.97: https://www.jstor.org/stable/2334940


## ------------------------------------------------------------------------------------------------------------
# one-hot clustering in one line of code:
# assign each sample the cluster according to its dominant NMF factor
# easily accessible using the max.col function
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


## ------------------------------------------------------------------------------------------------------------
# use the kmeans function to cluster the iCluster H matrix (here, z)
# using 2 as the number of clusters.
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


## ----moNMFHeatmap,fig.cap="Heatmap showing the association of input features from multi-omics data (gene expression, copy number variation, and mutations), with JNMF factors. Gene expression features dominate both factors, but copy numbers and mutations mostly affect only one factor each."----
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


## ----moenrichr,message=FALSE,warning=FALSE,results='hide',error=FALSE----------------------------------------
if(!require(enrichR)) BiocManager::install("enrichR")
# select genes associated preferentially with each factor
# by their relative loading in the W matrix
genes.factor.1 <- names(which(nmfw[1:dim(x1)[1],1] > nmfw[1:dim(x1)[1],2]))
genes.factor.2 <- names(which(nmfw[1:dim(x1)[1],1] < nmfw[1:dim(x1)[1],2]))
# call the enrichr function to find gene sets enriched
# in each latent factor in the GO Biological Processes 2018 library
go.factor.1 <- enrichR::enrichr(genes.factor.1,
                                databases = c("GO_Biological_Process_2018")
                                )$GO_Biological_Process_2018
go.factor.2 <- enrichR::enrichr(genes.factor.2,
                                databases = c("GO_Biological_Process_2018")
                                )$GO_Biological_Process_2018


## ----moNMFGOTerms,caption="Top GO-terms associated with NMF factor 2.",echo=FALSE,results='as.is'------------
library(kableExtra)
go.factor.2$Genes <- gsub(";", "; ", go.factor.2$Genes)
the.table <- knitr::kable(head(go.factor.2, 3)[,c("Term", "Adjusted.P.value", "Combined.Score")],
                 caption="GO-terms associated with NMF factor 2",
                 format="latex")
#the.table <- kableExtra::column_spec(the.table, 1, width="10em")
the.table <- kableExtra::kable_styling(the.table ,latex_options = c( "scale_down"))
#the.table <- kableExtra::column_spec(the.table, 4, width="10em")
the.table


## ----moNMFClinicalCovariates,fig.cap="Box plot showing MSI/MSS status distribution and NMF factor 1 values."----
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

## ----moNMFClinicalCovariates2,fig.cap="Box plot showing MSI/MSS status distribution and NMF factor 2 values."----
ggplot2::ggplot(cov_factor, ggplot2::aes(x=msi, y=factor2, group=msi)) +
  ggplot2::geom_boxplot() +
  ggplot2::ggtitle("NMF factor 2 and microsatellite instability")


## ----moNMFExerciseColumnSplitting, eval=FALSE, echo=TRUE-----------------------------------------------------
## # Implement this function
## split_neg_columns <- function(x) {
##     # your code here
## }
## # a test that shows the function above works
## test_split_neg_columns <- function() {
##     input <- as.data.frame(cbind(c(1,2,1),c(0,1,-2)))
##     output <- as.data.frame(cbind(c(1,2,1), c(0,0,0), c(0,1,0), c(0,0,2)))
##     stopifnot(all(output == split_neg_columns(input)))
## }
## # run the test to verify your solution
## test_split_neg_columns()


## ----moNMFCIMP,echo=FALSE, eval=FALSE------------------------------------------------------------------------
## a <- data.frame(age=covariates$age, gender=as.factor(covariates$gender), msi=covariates$msi, cimp=as.factor(covariates$cimp))
## b <- nmf.h
## colnames(b) <- c('factor1', 'factor2')
## cov_factor <- cbind(a,b)
## ggplot2::ggplot(cov_factor, ggplot2::aes(x=cimp, y=factor1, group=cimp)) + ggplot2::geom_boxplot() + ggplot2::ggtitle("NMF factor 1 and CIMP status")
## ggplot2::ggplot(cov_factor, ggplot2::aes(x=cimp, y=factor2, group=cimp)) + ggplot2::geom_boxplot() + ggplot2::ggtitle("NMF factor 2 and CIMP status")

