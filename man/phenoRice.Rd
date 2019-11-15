\name{phenorice}

\alias{phenorice}

\title{
R implementation of the phenorice model
}

\description{
Run the phenorice model
}

\usage{
phenorice(evi, ndfi, lst, p, checkLST=TRUE)
}

\arguments{
 \item{evi}{vector with evi values} 
 \item{ndfi}{vector with ndfi values} 
 \item{lst}{vector with lst values} 
 \item{p}{list with named parameters}
 \item{checkLST}{boolean. If \code{TRUE} the LST is used}
}


\details{
\tabular{lll}{
	\tab \bold{Parameters} \tab \bold{Explanation} \cr
	\tab \code{x} \tab bla bla  \cr
	\tab \code{y} \tab bla bla  \cr
}
}


\value{
vector
}

\references{
  Boschetti, M., L. Busettoa, G. Manfrona, A. Laborte, S. Asilo, S. Pazhanivelan & A Nelson, 2017. PhenoRice: A method for automatic extraction of spatio-temporal information on rice crops using satellite data time series. Remote Sensing of Environment 194: 347-365;

  Busetto, L., Zwart, S.J. and Boschetti, M., 2019.  Analysing spatial?temporal changes in rice cultivation practices in the Senegal River Valley using MODIS time-series and the PhenoRice algorithm. International Journal of Applied Earth Observation and Geoinformation, 75, pp.15-28.
}



\examples{

# Get index for one pixel
evi <- readRDS(system.file("evi.rds", package="phenorice"))
ndfi <- readRDS(system.file("ndfi.rds", package="phenorice"))
lst <- readRDS(system.file("lst.rds", package="phenorice"))

# Smooth each index
evi <- filter_VI(evi)
ndfi <- fill_VI(ndfi)
lst <- fill_VI(lst)

# Get parameters
p <- getPars()

# Run and get result (planting date, date of max EVI value, flowering date, heading date, harvest date) for tha pixel 
re <- phenorice(evi, ndfi, lst, p)

# Change parameters
p['evi_minth'] = .2
re <- phenorice(evi,ndfi,lst,p)


f <- system.file("ricevi.rds", package="phenorice")
d <- readRDS(f)
evi <- unlist(d[1,2:47])
ndfi <- unlist(d[1,48:93])
lst <- unlist(d[1,94:139])
p <- getPars()

phenorice(evi, ndfi, lst, p)


x <- matrix(NA, 25, 5)
for (i in 1:25) {
	evi <- unlist(d[i,2:47])
	ndfi <- unlist(d[i,48:93])
	lst <- unlist(d[i,94:139])
	x[i,] <- phenorice(evi, ndfi, lst, p)
}
x

}
