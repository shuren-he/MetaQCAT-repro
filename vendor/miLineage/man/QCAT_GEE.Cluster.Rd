\name{QCAT_GEE.Cluster}
\alias{QCAT_GEE.Cluster}

\title{
QCAT+GEE Two-Part Test for Clustered Microbiome Data (QCAT-C: two-part version)
}
\description{
This function allows users to (a) perform two-part QCAT-C on multivariate taxon counts; (b) perform two-part QCAT-C on the taxonomic tree to localize the covariate-associated lineages; and (c) assess the overall association of the microbial community with the covariate of interest.
}
\usage{
QCAT_GEE.Cluster(ID, OTU, OTU.base=NULL, X, X.index, Z, Z.index, Tax=NULL,
min.depth=0, perm.type=NULL, n.perm=NULL, fdr.alpha=0.05, test="chisq")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
  a vector contains cluster IDs
  }
  \item{OTU}{
a numeric matrix contains counts with each row corresponds to a sample and each column corresponds to an OTU or a taxon. Column name is mandatory. No missing values are allowed. 
}
  \item{OTU.base}{
If the clustered data is generated from a longituidnal study where microbiome samples are repeatedly measured. The OTU.base is a numeric matrix contains the baseline abundance if it needs to be adjusted for in the model. The data format should be the same as OTU. 
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
  \item{perm.type}{
type of permutation. Under default (perm.type=NULL), simple permutation (i.e. randomly shuffle covariates of interest across all samples) will be performed. If perm.type="WTH", within-cluster permutation will be performed. If perm.type="BTW", cluster-level permutation will be performed. 
}
  \item{n.perm}{
number of permutations. Under default (n.perm=NULL), only asymptotic test will be performed.
}
  \item{fdr.alpha}{
desired false discovery rate for multiple tests on the lineages.
}
  \item{test}{
specify "chisq" (default) for the multivariate chi-square test for both zero and positive parts or "mix" for the multivariate chi-square test for the zero part and the variance-component test for the positive part. 
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
Tang, Z.-Z. and Chen, G. (2019). Robust and Powerful Differential Composition Tests for
Clustered Microbiome Data. Statistics in Biosciences, 1-17.

}
\author{
Zheng-Zheng Tang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
data(data.toy.clustered)
OTU.toy = data.toy.clustered$OTU.toy
Tax.toy = data.toy.clustered$Tax.toy
meta = data.toy.clustered$meta.toy
case = meta[,1,drop=FALSE]
id = meta[,2]
# the OTUs should be consistent between the OTU table and the taxonomy table
OTU.toy.reorder = OTU.toy[,match(rownames(Tax.toy), colnames(OTU.toy))]
# perform QCAT+GEE test
QCAT_GEE.Cluster(ID=id, OTU=OTU.toy.reorder, X=case, Tax=Tax.toy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution-free association test }
