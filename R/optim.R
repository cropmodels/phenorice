

.kappa <- function(obs, prd) {
	conmat <- table(obs, prd)
	n <- sum(conmat) 
	# observed (true) cases per class
	rs <- rowSums(conmat) 
	p <- rs / n
	# predicted cases per class
	cs <- colSums(conmat) 
	q <- cs / n 
	expAccuracy <- sum(p*q)
	OA <- sum(diag(conmat))/ n
	(OA - expAccuracy) / (1 - expAccuracy)
}

.optim_pheno <- function(pars, d) {
	n <- nrow(d)
	x <- matrix(NA, n, 2)
	for (i in 1:n) {
		evi <- unlist(d[i,2:47])
		ndfi <- unlist(d[i,48:93])
		lst <- unlist(d[i,94:139])
		x[i,] <- phenorice(evi, ndfi, lst, as.list(pars))
	}
	predrice <- x[,1] > 0
	obsrice  <- d[,1] > 0
	1 - .kappa(obsrice, predrice)
	#x <- table(predrice, obsrice)
	#1 - sum(diag(x)) / sum(x)
}



