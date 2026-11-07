\name{QCAT}
\alias{QCAT}

\title{
Quasi-Conditional Association Test (QCAT)
}
\description{
This function allows users to (a) perform QCAT on multivariate taxon counts; (b) perform QCAT on the taxonomic tree to localize the covariate-associated lineages; and (c) assess the overall association of the microbial community with the covariate of interest.
}
\usage{
QCAT(OTU, X, X.index, Tax = NULL, 
min.depth = 0, n.perm = NULL, fdr.alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OTU}{
a numeric matrix contains counts with each row corresponds to a sample and each column corresponds to an OTU or a taxon. Column name is mandatory. No missing values are allowed.
}
  \item{X}{
a numeric matrix contains covariates with each column pertains to one variable. Samples in the OTU and X matrices should be identical and in the same order. No missing values are allowed. 
}
  \item{X.index}{
a vector indicates the columns in X for the covariate(s) of interest. The remaining columns in X will be treated as confounders. If X.index is not specified, all the columns in X will be treated as the covariates of interest.
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
If Tax=NULL (Default), a test is performed using all the OTUs/taxa.

If Tax is provided, tests are performed for lineages derived from the taxonomic hierarchy. The output is a list that contains 3 components
\item{lineage.pval}{p-values for all lineages. By default (n.perm=NULL), only the asymptotic test will be performed.}
\item{sig.lineage}{a vector of significant lineages}
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
# perform QCAT
QCAT(OTU=OTU.toy.reorder, X=case, Tax=Tax.toy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution-free association test }
