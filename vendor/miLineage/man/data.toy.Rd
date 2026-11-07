\name{data.toy}
\alias{data.toy}
\docType{data}
\title{
Simulated data that include an OTU table, a taxonomy table, and a meta-data table
}
\description{
This data is simulated based on a mock taxonomic tree. The data contain an OTU table, a taxonomy table, and a meta-data table. The data include 50 samples (25 cases and 25 controls) and 7 OTUs.   
}
\usage{data("data.toy")}
\format{
data.toy contains the following objects:
\itemize{
\item OTU.toy: a matrix of OTU counts for 50 samples and 7 OTUs
\item Tax.toy: a matrix of taxonomy table with 3 ranks  
\item meta.toy: a matrix contains a variable indicating case/control status
}
}
\source{
This is a simulated data based on a mock taxonomic tree. 
}
\examples{
data(data.toy)
}
\keyword{datasets}
