
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(UpSetR)

(patient_data <- data.frame(sex = c("M", "F", "M", "F", "M"),
                            age = 38:42,
                            row.names = c("Alex", "Lupe", "Klaus", 
                                          "Mireia", "Guille")))

# expression

(arraydat <- matrix(seq(101, 110),
                    ncol = 5,
                    dimnames = list(c("ENST00000294241", 
                                      "ENST00000355076"),
                                    c("array1", "array2", 
                                      "array3", "array4",
                                      "array5"))))


coldat <- data.frame(slope53 = rnorm(5),
                     row.names = c("array1", "array2", 
                                   "array3", "array4",
                                   "array5"))

exprdat <- SummarizedExperiment(arraydat, colData = coldat)
exprdat

(exprmap <- data.frame(primary = rownames(patient_data)[c(1, 2, 4, 3, 5)],
                       colname = c("array1", "array2", "array3", 
                                   "array4", "array5"),
                       stringsAsFactors = FALSE))

# methylation

(methyldat <- matrix(1:10, 
                     ncol = 5,
                     dimnames = list(c("ENST00000355076", "ENST00000383706"),
                                     c("methyl1", "methyl2", "methyl3",
                                       "methyl4", "methyl5"))))


(methylmap <- data.frame(primary = c("Alex", "Lupe", "Mireia", 
                                     "Klaus", "Guille"),
                         colname = c("methyl1", "methyl2", "methyl3", 
                                     "methyl4", "methyl5"),
                         stringsAsFactors = FALSE))

# microRNA

(microdat <- matrix(201:212, 
                    ncol = 3,
                    dimnames = list(c("hsa-miR-21", "hsa-miR-191",
                                      "hsa-miR-148a", "hsa-miR148b"),
                                    c("micro1", "micro2", "micro3"))))

(micromap <- data.frame(primary = c("Alex", "Mireia", "Klaus"),
                        colname = c("micro1", "micro2", "micro3"), 
                        stringsAsFactors = FALSE))

# metabolomics

(metaboldat <- matrix(round(rnorm(16, 1000, 500), 2), 
                      ncol = 4,
                      dimnames = list(c("alanine", "histidine", 
                                        "phloretin", "resveratrol"),
                                      c("met004", "met005",
                                        "met006", "met007"))))


(metabolmap <- data.frame(primary = c("Klaus", "Alex"),
                          colname = c("met004", "met005"),
                          stringsAsFactors = FALSE))

# create MultiAssayExperiment

listmap <- list(exprmap, methylmap, micromap, metabolmap)
names(listmap) <- c("Expression", "Methylation", "microRNA", "Metabolomics")
listmap

dfmap <- listToMap(listmap)
dfmap

objlist <- list("Expression" = exprdat, 
                "Methylation" = methyldat,
                "microRNA" = microdat,
                "Metabolomics" = metaboldat)

myMultiAssay <- MultiAssayExperiment(objlist, patient_data, dfmap)
myMultiAssay

# extract data

experiments(myMultiAssay)
colData(myMultiAssay)
sampleMap(myMultiAssay)
metadata(myMultiAssay)

# visualize sample intersection

upsetSamples(myMultiAssay)

