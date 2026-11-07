\name{QCAT_GEE}
\alias{QCAT_GEE}

\title{
QCAT+GEE Two-Part Test
}
\description{
This function performs Quasi-Conditional Association Test (QCAT) for the positive taxon counts and (Generalized Estimating Equation) GEE-based association test for the zero counts. Then, the two-part test is from the combination of those two tests.  
This function allows users to (a) perform QCAT+GEE on multivariate taxon counts; (b) perform QCAT+GEE on the taxonomic tree to localize the covariate-associated lineages; and (c) assess the overall association of the microbial community with the covariate of interest.
}
\usage{
QCAT_GEE(OTU, X, X.index, Z, Z.index, Tax = NULL, 
min.depth = 0, n.perm = NULL, fdr.alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OTU}{
a numeric matrix contains counts with each row corresponds to a sample and each column corresponds to an OTU or a taxon. Column name is mandatory. No missing values are allowed. 
}
  \item{X}{
a numeric matrix contains covariates for the positive-part test with each column pertains to one variable. Samples in the OTU and X matrices should be identical and in the same order. No missing values are allowed. 
}
  \item{X.index}{
a vector indicates the columns in X for the covariate(s) of interest in the positive-part test. The remaining columns in X will be treated as confounders in modeling the positive abundance. If X.index is not specified, all the columns in X will be treated as the covariates of interest for the positive-part test.
}
  \item{Z}{
a numeric matrix contains covariates for the zero-part test with each column pertains to one variable. Samples in the OTU and Z matrices should be identical and in the same order. No missing values are allowed. If Z is not specified, Z will be set to X.   
}
  \item{Z.index}{
a vector indicates the columns in Z for the covariate(s) of interest in the zero-part test. The remaining columns in Z will be treated as confounders in modeling the proportion of zero. If Z.index is not specified, all the columns in Z will be treated as the covariates of interest for the zero-part test.
}
  \item{Tax}{
a matrix define the taxonomic ranks with each row corresponds to an OTU or a taxon and each column corresponds to a rank (start from the higher taxonomic level, e.g., from kingdom to genus). Row name is mandatory and should be consistent with the column name of the OTU table,  Column name should be formatted as "Rank1", "Rank2", etc, representing the taxonomic levels from highest to lowest. 
}
  \item{min.depth}{
lower bound for sample read depth. Samples with read depth less than min.depth will be removed before the analysis.
}
  \item{n.perm}{
number of permutations. Under default (n.perm=NULL), only asymptotic test will be performed.
}
  \item{fdr.alpha}{
desired false discovery rate for multiple tests on the lineages.
}
}
\value{
If Tax=NULL (Default), tests are performed using all the OTUs/taxa.

If Tax is provided, tests are performed for lineages derived from the taxonomic hierarchy. The output is a list that contains 3 components
\item{lineage.pval}{a list of zero-part, positive-part and two-part (combined) p-values for all lineages. By default (n.perm=NULL), only the asymptotic test will be performed. }
\item{sig.lineage}{a list of significant lineages based on the zero-part, positive-part and two-part p-values}
\item{global.pval}{p-values of the global tests}


}
\references{
Tang, Z.-Z., Chen, G., Alekseyenko, A. V., and Li, H. (2017). A General Framework for Association Analysis of Microbial Communities on a Taxonomic Tree. Bioinformatics, 33(9), 1278-1285.

}
\author{
Zheng-Zheng Tang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
data(data.toy)
OTU.toy = data.toy$OTU.toy
Tax.toy = data.toy$Tax.toy
case = data.toy$meta.toy
# the OTUs should be consistent between the OTU table and the taxonomy table
OTU.toy.reorder = OTU.toy[,match(rownames(Tax.toy), colnames(OTU.toy))]
# perform QCAT+GEE test
QCAT_GEE(OTU=OTU.toy.reorder, X=case, Tax=Tax.toy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution-free association test }
