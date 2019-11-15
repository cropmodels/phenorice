# Functions for PhenoRice methods

# check an increase (decrease) in EVI before (after) the EVI maximum
# at least 3 increasing (decreasing) observations within 5 consecutive periods around maximum EVI points
.checkChangeRate <- function(maxevidate, evi){
  startdate <- maxevidate - 5
  enddate <- maxevidate + 5
  evibefore <- evi[startdate:maxevidate]
  eviafter <- evi[maxevidate:enddate]
  inc <- sum(sign(diff(evibefore)))
  dec <- sum(sign(diff(eviafter)))
  if (inc >= 1 & dec <= -1){
    return(maxevidate)
  } else {
    return(0)
  }
}


# check an increase in EVI after the first minimum or start of season
# at least 3 increasing observations within 5 consecutive periods after minimim
.checkChangeRate_min <- function(mind, evi){
	if (length(mind) == 0) return(mind)
	enddate <- mind + 5
	enddate[enddate > length(evi)] <- length(evi)
	for (i in 1:length(mind)) {
		eviafter <- evi[mind[i]:enddate[i]]
		inc <- sum(sign(diff(eviafter)))
  # 3 (-)ive and 2 (+) will result in 1 that is the minimum threhold
		if (inc < 1) {
			mind[i] <- NA
		}
	}
    mind
}


# Check for NDFI signal >= minndfi around minimum EVI for a period (by setting winfl) before and after min
.ndficheck <- function(mind, ndfi, evi, winfl, minndfi){
  startdate <- mind - winfl 
  enddate <- mind + winfl 
  startdate <- max(1, startdate)
  enddate <- min(enddate, length(evi))
  ndfisub <- ndfi[startdate:enddate]
  k <- (ndfisub >= minndfi) 
  if (sum(k) >= 1){
    return(mind)
  } else {
    return(NA)
  } 
}


# Check the period between min EVI date and max EVI date
.evirangeCheck2 <- function(min_evi_dates, max_evi_date, vl1, vl2) {
	dif <- max_evi_date - min_evi_dates
	min_evi_dates[dif >= vl1 & dif <= vl2]
}

	
# Check for land surface temperature
.lstcheck2 <- function(min_evi_dates, lst, lst_th) {
    i <- lst[min_evi_dates] >= lst_th
	min_evi_dates[stats::na.omit(i)]
}	


# Check if EVI decrease after max following PhenoRice definition of sharp decrease
# The algorithm checks if EVI decreases by more than decr of the amplitude of the min-max range 
# decrease in EVI after a period (windecr)
.checkEVIaftermax2 <- function(dates, evi, windecr, decr){
 
  good <- rep(TRUE, nrow(dates))
  for (i in 1:nrow(dates)) {
	maxdate <- dates[i,2]
	enddate <- maxdate + windecr
	mindate <- dates[i,1]
	enddate <- min(enddate, length(evi))
	dtrng <- maxdate:enddate 
	mv <- min(evi[dtrng])
 
	test = evi[maxdate]- decr*(evi[maxdate] - evi[mindate])
	good[i] <- mv <= test	
  }
  dates[good, ,drop=FALSE]
}


# main function
phenorice <- function(evi, ndfi, p, lst,...){
  
  # making lst optional
  checkLST <- TRUE
  if(missing(lst)) checkLST = FALSE
  
  # output for single season rice
	rice <- cbind(NA, NA, NA, NA, NA)
	if (any(is.na(evi))) return(rice)
	#if (any(is.na(ndfi))) return(rice)
	if(checkLST){
	  if (any(is.na(lst))) return(rice)
	}

	rice <- cbind(0, 0, 0, 0, 0)
	
	# check if mean evi is not too high
	if (mean(evi) > p$evi_meanth) return(rice)
	
    # Identify all local EVI maxima
    max_evi_dates <- which(diff(sign(diff(evi))) == -2)+1
    # Identify large EVI maxima
    max_evi_dates  <- max_evi_dates[evi[max_evi_dates] >= p$evi_maxth] # threshold can be changed, for MOD09A1 some values are pretty low
    # we are interested in max_evi_dates within certain date range within the year. 
    max_evi_dates <- max_evi_dates[max_evi_dates >= p$pos_start & max_evi_dates <= p$pos_end]

    if (length(max_evi_dates) == 0) return(rice)
    max_evi_dates <- sapply(max_evi_dates, .checkChangeRate, evi)
    max_evi_dates <- max_evi_dates[max_evi_dates > 0]
    if (length(max_evi_dates) == 0) return(rice)
    max_evi_dates <- max_evi_dates[which.max(evi[max_evi_dates])]
	
	# Find dates for EVI minima
    min_evi_dates <- which(diff(sign(diff(c(evi)))) == 2) + 1
    names(min_evi_dates) <- NULL
    if (length(min_evi_dates) == 0) return(rice)
    min_evi_dates <- min_evi_dates[evi[min_evi_dates] <= p$evi_minth]
    
    if(checkLST){
    min_evi_dates <- .lstcheck2(min_evi_dates, lst, p$lst_th)
    }
    
    min_evi_dates <- .evirangeCheck2(min_evi_dates, max_evi_dates, p$vl1, p$vl2)
    min_evi_dates <- sapply(min_evi_dates, .ndficheck, ndfi, evi, p$winfl, p$minndfi)
    min_evi_dates <- min_evi_dates[!is.na(min_evi_dates)]
	  min_evi_dates <- sapply(min_evi_dates, .checkChangeRate_min, evi)
    min_evi_dates <- min_evi_dates[!is.na(min_evi_dates)]
    
    i <- which(sapply(min_evi_dates, length) == 0)
    
    if (length(i) == length(min_evi_dates)) {
      return(rice)
      } else if (length(i) > 0) {
        max_evi_dates <- max_evi_dates[-i]	
        min_evi_dates <- min_evi_dates[-i]
        }
		
    min_evi_dates <- max(min_evi_dates)
    dates <- cbind(min_evi_dates, max_evi_dates)
    dates <- .checkEVIaftermax2(dates, evi, p$windecr, p$decr)
    i <- which.max(evi[dates[,2]])
    dates <- dates[i, ,drop=FALSE]
	
	if (nrow(dates) == 0) {
		return(rice)
	} else {
	  # Flowering date 
	  flowering <- stats::median((dates[1]:length(evi))[evi[dates[1]:length(evi)] >= (evi[dates[2]] - (0.1*(evi[dates[2]]-evi[dates[1]])))])
	  # Heading date
	  heading <- max((dates[2]:length(evi))[evi[dates[2]:length(evi)] >= (evi[dates[2]] - (0.1*(evi[dates[2]]-evi[dates[1]])))])
	  # Harvest date
	  eos <- min((dates[2]:length(evi))[evi[dates[2]:length(evi)] <= (evi[dates[2]] - (p$decr*(evi[dates[2]]-evi[dates[1]])))])
	  if((eos-dates[1]>=p$tl1) & (eos-dates[1]<=p$tl2)){
	    return(c(dates, flowering, heading, eos))
	  } else {
	    return(rice)
	  }
	}
}


# Parameters for paddy rice in CA
getPars <- function() {
	pars <- list()
	pars$evi_meanth <- 0.4 # threshold for annual mean EVI
	pars$evi_maxth <- 0.5 # threshold for maximum evi 
	pars$evi_minth <- 0.25 # threshold for minimum evi 
	pars$pos_start <- 23 # start of heading
	pars$pos_end <- 32 # end of heading
	pars$vl1 <- 5 # the shortest vegetative growth length
	pars$vl2 <- 14 # the longest vegetative growth length
	pars$winfl <- 2 # period for flooding
	pars$minndfi <- 0 # threshold for ndfi
	pars$windecr <- 15 # period after EVI maximum 
	pars$decr <- 0.7 # percent decrease of EVI after EVI maximum
	pars$tl1 <- 13 # the shortest total growing length
	pars$tl2 <- 23 # the longest total growing length
	pars$lst_th <- 15 # the minmum temperature for planting 
	pars
}

