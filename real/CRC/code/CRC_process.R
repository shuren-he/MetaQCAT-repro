rm(list=ls())

library(dirmult)
library(HMP)
library(stringr)
library(zoo)
library(qvalue)
library(grid)
library(MetaQCAT)

# Locate the project from this script (Rscript/source) or the working directory.
# This keeps paths portable if the project folder is moved or renamed.
project.root <- local({
  source.files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
  }, character(1))
  script.args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script.files <- c(rev(source.files[nzchar(source.files)]),
                    sub("^--file=", "", script.args))
  start.dirs <- unique(c(dirname(script.files), getwd()))
  for (start.dir in start.dirs) {
    if (!dir.exists(start.dir)) next
    candidate <- normalizePath(start.dir, winslash = "/", mustWork = TRUE)
    repeat {
      if (all(dir.exists(file.path(candidate, c("real/CRC/data"))))) {
        return(candidate)
      }
      parent <- dirname(candidate)
      if (identical(parent, candidate)) break
      candidate <- parent
    }
  }
  stop("Cannot locate the MetaQCAT project. Run from the project directory or source this script by its file path.")
})

data.dir <- file.path(project.root, "real", "CRC", "data")

load(file.path(data.dir, "count.genus.RData"))
load(file.path(data.dir, "tax.RData"))
load(file.path(data.dir, "meta.RData"))

count.genus <- count.genus[-19,] # remove extremely small  observation
meta <- meta[-19,]

otu.tab.lst <- split.data.frame(count.genus, meta$Country)

rename.mat <- function(x,col,names){
  a <- matrix(x,ncol=col,byrow = F)
  colnames(a) = names
  return(a)
}

## For each matrix in the list, its columns names must be the scientific name of the genus (genus level)
otu.tab.lst <- lapply(otu.tab.lst, rename.mat, col = ncol(count.genus), names=colnames(count.genus))

## Find the column indexes for each otu table which has more than 90% zero observations
.Remove.ind <- function(otu.tab.lst, ratio){
  n = length(otu.tab.lst)
  rev.ind.lst = NULL
  for(i in 1:n)
  {
    Y <- otu.tab.lst[[i]]
    Y0 = Y
    Y0[Y>0] = 0
    Y0[Y==0] = 1
    rev.ind.lst[[i]] <-which(colSums(Y0)/nrow(Y0)>ratio)
  }
  return (rev.ind.lst)
}

## Denote the genus that has more than 90% zero observations
remove.index.lst  <- .Remove.ind(otu.tab.lst, 0.9)

## Remove the indexes for each OTU table in the list
.Remove.col  <- function(otu.tab.lst,rev.ind.lst){
  if(length(rev.ind.lst) !=  length(otu.tab.lst))
  {
    stop("the number of elements of otu table list and remove index list should be the same")
  }
  else
  {
    n = length(otu.tab.lst)
    for(i in 1:n){
      rev.ind = rev.ind.lst[[i]]
      if(length(rev.ind) == 0){
        next
      }
      otu.tab.lst[[i]] = otu.tab.lst[[i]][ ,-rev.ind, drop = FALSE]
    }
  }
  return(otu.tab.lst)
}

counts <- .Remove.col(otu.tab.lst, remove.index.lst)

## Rename the columns for Tax
## Kingdom -> Rank1, Phylum -> Rank5, Class -> Rank4, ..., Genus -> Rank1
tax = tax[6:7] # Extract the Tax information needed
tax <- tax[!duplicated(tax$genus),]
colnames(tax) <- paste("Rank",c(1:ncol(tax)),sep="")
rownames(tax) <- tax$Rank2

## Construct the design matrix (Case group:1; Control Group: 0). Can add more covariates of interest like age, gender, etc...
disease <- rep(1,nrow(count.genus))
disease[meta$Group == 'CTR'] <- 0
disease <- split(disease, meta$Country)

FE_Time = Sys.time()
CRC_FE <- MetaQCAT(counts, disease, 1, Tax = tax, n.perm = 1000, meta.method = "FE-MetaQCAT-O", n.threads = 8)
FE_Time = difftime(Sys.time(), FE_Time, units = "secs")

RE_Time = Sys.time()
CRC_RE <- MetaQCAT(counts, disease, 1, Tax = tax, n.perm = 1000, meta.method = "RE-MetaQCAT-O", n.threads = 8)
RE_Time = difftime(Sys.time(), RE_Time, units = "secs")


CRC_FE_pval <- CRC_FE$lineage.pval[2,]

# q-values for the paper's discovery plots (qvalue package).
CRC_FE_qval <- qvalue::qvalue(p = CRC_FE_pval)$qvalues


CRC_RE_pval <- CRC_RE$lineage.pval[1,]

# q-values for the paper's discovery plots (qvalue package).
CRC_RE_qval <- qvalue::qvalue(p = CRC_RE_pval)$qvalues
