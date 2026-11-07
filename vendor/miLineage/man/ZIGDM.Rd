\name{ZIGDM}
\alias{ZIGDM}

\title{
Zero-Inflated Generalized Dirichlet Multinomial (ZIGDM) Tests
}
\description{
Different from the distribution-free QCAT and QCAT+GEE, the ZIGDM tests are based on a parametric model (ZIGDM) for multivariate taxon counts. 
The ZIGDM tests can not only detect differential mean but also differential dispersion level or presence-absence frequency in microbial compositions.
This function allows users to (a) perform ZIGDM tests on multivariate taxon counts; (b) perform ZIGDM tests on the taxonomic tree to localize the covariate-associated lineages; and (c) assess the overall association of the microbial community with the covariate of interest.
}
\usage{
ZIGDM(OTU, X4mean=NULL, X4disp=NULL, X4zero=NULL, test.type = "Mean", 
X.index, ZI.LB = 10, Tax = NULL, min.depth = 0, 
n.perm = NULL, fdr.alpha = 0.05, details=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OTU}{
a numeric matrix contains counts with each row corresponds to a sample and each column corresponds to an OTU or a taxon. Column name is mandatory. No missing values are allowed. 
}
  \item{X4mean}{
a numeric matrix contains covariates that link to mean abundance in microbial compositions. Each column pertains to one variable. Set X4mean=NULL if only intercept term is needed. Samples in the OTU and X4mean matrices should be identical and in the same order. No missing values are allowed.  
}
  \item{X4disp}{
a numeric matrix contains covariates that link to dispersion level in microbial compositions. Each column pertains to one variable. Set X4disp=NULL if only intercept term is needed. Samples in the OTU and X4disp matrices should be identical and in the same order. No missing values are allowed.   
}
  \item{X4zero}{
a numeric matrix contains covariates that link to presence-absence frequency in microbial compositions. Each column pertains to one variable. Set X4zero=NULL if only intercept term is needed. Samples in the OTU and X4zero matrices should be identical and in the same order. No missing values are allowed. 
}
  \item{test.type}{
  If test.type = "Mean", the function will test for differential mean (Default). If test.type = "Disp", the function will test for differential dispersion. If test.type = "Zero", the function will test for differential presence-absence frequency. 
}
  \item{X.index}{
If test.type = "Mean", X.index is a vector indicates the columns in X4mean for the covariate(s) of interest. The remaining columns in X4mean will be treated as confounders in modeling the abundance of the present taxa.
If test.type = "Disp", X.index is a vector indicates the columns in X4disp for the covariate(s) of interest. The remaining columns in X4disp will be treated as confounders in modeling the dispersion level of the present taxa.
If test.type = "Zero", X.index is a vector indicates the columns in X4zero for the covariate(s) of interest. The remaining columns in X4zero will be treated as confounders in modeling the presence-absence frequency.
If X.index is not specified, all the columns in X4* matrix will be treated as the covariates of interest. 
}
  \item{ZI.LB}{
lower bound of zero counts for the taxon that needs zero-inflated model. 
The counts for a taxon is assumed to be zero-inflated if the number of zero observations for the taxon is greater than ZI.LB (Default is 10). If ZI.LB=NULL, GDM model (i.e. non-zero-inflated version of ZIGDM) will be applied.
}
  \item{Tax}{
a matrix define the taxonomic ranks with each row corresponds to an OTU or a taxon and each column corresponds to a rank (start from the higher taxonomic level, e.g., from kingdom to genus). Row name is mandatory and should be consistent with the column name of the OTU table,  Column name should be formatted as "Rank1", "Rank2", etc, representing the taxonomic levels from highest to lowest. 
}
  \item{min.depth}{
lower bound of sample read depth. Samples with read depth less than min.depth will be removed before the analysis.
}
  \item{n.perm}{
number of permutations. By default (n.perm=NULL), only the asymptotic test will be performed.
}
  \item{fdr.alpha}{
desired false discovery rate for multiple tests on the lineages.
}
  \item{details}{
If TRUE, the function will output parameter estimates (under the null model) and score test statistics.
}
}
\value{
If Tax=NULL (Default), a test is performed using all the OTUs/taxa.

If Tax is provided, tests are performed for lineages derived from the taxonomic hierarchy. The function will output a list that contains 3 components
\item{lineage.pval}{p-values for all lineages. By default (n.perm=NULL), only the asymptotic test will be performed. }
\item{sig.lineage}{a vector of significant lineages}
\item{global.pval}{p-values of the global tests}

If details=TRUE, The function will output another list that contains parameter estimates (under the null model) and the score test statistic for each lineage.
}
\references{
Tang, Z.-Z. and Chen, G. (2019). Zero-Inflated Generalized Dirichlet Multinomial Regression Model for Microbiome Compositional Data Analysis. Biostatistics, 20(4), 698-713.

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
# perform ZIGDM test for detecting differential dispersion
ZIGDM(OTU=OTU.toy.reorder, X4disp=case, test.type="Disp", Tax=Tax.toy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ parametric association test }
