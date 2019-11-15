\name{getPars}

\alias{getPars}

\title{
Default parameters
}

\description{
Get default parameters for the phenorice model
}

\usage{
  getPars()
}


\details{
  \tabular{lll}{
    \tab \bold{Parameters} \tab \bold{Explanation} \cr
    \tab \code{evi_meanth} \tab threshold for annual mean EVI  \cr
    \tab \code{evi_maxth} \tab threshold for maximum evi   \cr
    \tab \code{evi_minth} \tab threshold for minimum evi   \cr
    \tab \code{pos_start} \tab start of heading  \cr
    \tab \code{pos_end} \tab end of heading  \cr
    \tab \code{vl1} \tab the shortest vegetative growth length  \cr
    \tab \code{vl2} \tab the longest vegetative growth length  \cr
    \tab \code{winfl} \tab half period of flooding \cr
    \tab \code{minndfi} \tab threshold for ndfi  \cr
    \tab \code{windecr} \tab period after maximum EVI  \cr
    \tab \code{decr} \tab percent decrease of EVI after EVI maximum  \cr
    \tab \code{tl1} \tab the shortest total growing length \cr
    \tab \code{tl2} \tab the longest total growing length  \cr
    \tab \code{lst_th} \tab the minmum temperature for planting   \cr
  }
}


\value{
  list
}


\examples{
p <- getPars()
str(p)
}
