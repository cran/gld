plotgld <- function(lambda1 = 0, lambda2 = 1, lambda3, lambda4, param =
"fmkl", new.plot = TRUE, truncate = 0, bnw = FALSE, col.or.type = 1, granularity = 4000, 
xlab = NULL, ylab=NULL, quant.probs = seq(0,1,.25),...)
{
	if (is.null(xlab)){xlab <- "x"}
	u <- seq(from = 0, to = 1, by = 1/granularity)
	# Only difference across parameterisations is calculating the 
	# quantiles and density
	quantiles <- qgl(u, lambda1, lambda2, lambda3, lambda4,param)
	density <- qdgl(u,lambda1,lambda2,lambda3,lambda4,param)
	if(truncate > 0) {
		if(new.plot) {
			if (is.null(ylab)){
				ylab <- paste( "probability density (values below", deparse(substitute(truncate)), "not shown)")
			}
			if(bnw) {
				plot(quantiles[density > truncate], 
				density[density > truncate], type = "l", 
				xlab = xlab, ylab = ylab, lty = col.or.type,...)
			}
			else {
				plot(quantiles[density > truncate], 
				density[density > truncate], type = "l", 
				xlab = xlab, ylab = ylab, col = col.or.type,...)
			}

		}
		else {
			if(bnw) {
				lines(quantiles[density > truncate], density[
				  density > truncate], lty = col.or.type)
			}
			else {
				lines(quantiles[density > truncate], density[
				  density > truncate], col = col.or.type)
			}
		}
	}
	else {
		if(new.plot) {
			if (is.null(ylab)){
				ylab <- "probability density"
			}
			if(bnw) {
				plot(quantiles, density, type = "l", xlab = xlab,ylab = ylab, lty=col.or.type, ...) 	
			}
			else {
				plot(quantiles, density, type = "l", xlab =
xlab,ylab = ylab, col=col.or.type, ...) 	
			}
		}
		else {
			if(bnw) {
				lines(quantiles, density, lty = col.or.type)
			}
			else {
				lines(quantiles, density, col = col.or.type)
			}
		}
	}
	if (!is.null(quant.probs)){quantile(quantiles,quant.probs) } 
}

plotglc <- function(lambda1 = 0, lambda2 = 1, lambda3, lambda4,param="fmkl", granularity=4000,
xlab="x",ylab="cumulative probability",...)
{
	u <- seq(from = 1/granularity, to = 1 - 1/granularity, length = 
		granularity - 1)
	x <- qgl(u,lambda1,lambda2,lambda3,lambda4,param)
	plot(x, u, pch = ".",xlab=xlab,ylab=ylab,...)
}
