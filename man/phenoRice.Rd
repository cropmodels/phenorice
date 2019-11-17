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
 \item{evi}{numeric, rasterstack or SpatRaster with evi values} 
 \item{ndfi}{numeric, rasterstack or SpatRaster with ndfi values} 
 \item{lst}{numeric, rasterstack or SpatRaster with lst values} 
 \item{p}{list with named parameters based on phenorice method}
 \item{checkLST}{boolean. If \code{TRUE} the LST is used}
}


\details{
\tabular{lll}{
	\tab \bold{Parameters} \tab \bold{Explanation} \cr
	\tab \code{evi_meanth} \tab threshold for annual mean EVI \cr
	\tab \code{evi_maxth} \tab threshold for maximum EVI  \cr
	\tab \code{evi_minth} \tab threshold for minimum EVI \cr
	\tab \code{pos_start} \tab start of heading \cr
	\tab \code{pos_end} \tab end of heading \cr
	\tab \code{vl1} \tab shortest vegetative growth length \cr
	\tab \code{vl2} \tab longest vegetative growth length \cr
	\tab \code{winfl} \tab period for flooding  \cr
	\tab \code{minndfi} \tab threshold for NDFI  \cr
	\tab \code{windecr} \tab period after EVI maximum  \cr
	\tab \code{decr} \tab percent decrease of EVI after EVI maximum  \cr
	\tab \code{tl1} \tab shortest total growing length  \cr
	\tab \code{tl2} \tab longest total growing length  \cr
	\tab \code{lst_th} \tab minmum land surface temperature for planting, needs \code{checkLST=TRUE}  \cr
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
