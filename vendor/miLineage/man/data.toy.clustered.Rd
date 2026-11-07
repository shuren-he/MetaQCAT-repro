\name{data.toy.clustered}
\alias{data.toy.clustered}
\docType{data}
\title{
Simulated data mimicking a matched case-control study. 
}
\description{
The data contain an OTU table, a taxonomy table, and a meta-data table. The data include 50 samples (25 matched case-control pairs) and 7 OTUs. 
}
\usage{data("data.toy.clustered")}
\format{
data.toy contains the following objects:
\itemize{
\item OTU.toy: a matrix of OTU counts for 50 samples and 7 OTUs
\item Tax.toy: a matrix of taxonomy table with 3 ranks  
\item meta.toy: a matrix of two variables. One variable indicates case/control status, the other variable is the case-control pair IDs. 
}
}
\source{
This is a simulated data based on a mock taxonomic tree. 
}
\examples{
data(data.toy.clustered)
}
\keyword{datasets}
