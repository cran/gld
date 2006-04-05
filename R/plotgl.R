plotgld <- function(lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,
  param="fmkl",lambda5=NULL, new.plot = TRUE, truncate = 0, bnw = FALSE,
  col.or.type = 1, granularity = 4000,xlab = "x", ylab=NULL, 
  quant.probs = seq(0,1,.25), ...)
{
# standard parameter fixin - copied directly from dgl, but we want the 
# warnings to happen in this function.
# Tidy the parameters so gl.check.lambda will work
lambdas <- .gl.parameter.tidy(lambda1,lambda2,lambda3,lambda4,param,lambda5)
# Check the parameters
if(!gl.check.lambda(lambdas,param=param,vect=TRUE)) {
        stop(paste("The parameter values", lambdas,
"\ndo not produce a proper distribution with the",param,
"parameterisation - see \ndocumentation for gl.check.lambda"))
        }
u <- seq(from = 0, to = 1, by = 1/granularity)
# Only difference across parameterisations is calculating the 
# quantiles and density
quantiles <- qgl(u,lambda1=lambdas,param=param)
density <- qdgl(u,lambda1=lambdas,param=param)
# Check for pathological densities - include check for all,
dots.required <- FALSE
nearzero <- .Machine$double.eps*c(0,1,3)
check.for.jump <- qgl(nearzero,lambda1=lambdas,param=param)
check.for.jump2 <- qgl(1-nearzero,lambda1=lambdas,param=param)
if (is.finite(check.for.jump[1]))
	{
	first.space <- check.for.jump[2] - check.for.jump[1]
	from1to3 <- check.for.jump[3] - check.for.jump[2]
	if (first.space > (from1to3 * 1e10)) {
		warning(paste("These parameter values give a pathological density.  \n",
		"The theoretical minimum: F^{-1}(0)=",signif(check.for.jump[1],4),"\n is much less than F^{-1}("
		,signif(nearzero[2],4),")=",signif(check.for.jump[2],4),"The density is undefined between these points"))
		dots.required <- TRUE
		}
	}
if (is.finite(check.for.jump2[1]))
	{
	last.space <- check.for.jump2[1] - check.for.jump2[2]
	from1to3 <- check.for.jump2[2] - check.for.jump2[3]
	print(c(first.space,from1to3))
	if (first.space > (from1to3 * 1e10)) {
		warning(paste("These parameter values give a pathological density.  \n",
		"The theoretical maximum: F^{-1}(1)=",signif(check.for.jump2[1],4),
		"\n is much more than F^{-1}(1-",signif(nearzero[2],4),")=",signif(check.for.jump2[2],4),".\n","The density is undefined between these points"))
		dots.required <- TRUE
		}
	}
# warning - or different plot
# Fix this for infinity 
if(truncate > 0) { 
	# If truncated, not RS pathological problem, because density at the endpoint is zero
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
			if (dots.required) {
				plot(quantiles, density, type = "n", xlab = xlab,ylab = ylab, lty=col.or.type, ...)
				points(quantiles[1],density[1])
				points(quantiles[granularity+1],density[granularity+1])
				lines(quantiles[-c(1,granularity+1)], density[-c(1,granularity+1)],  lty=col.or.type, ...)
				}
			else {
				plot(quantiles, density, type = "l", xlab = xlab,ylab = ylab, lty=col.or.type, ...) 
				}
			}
		else {
			if (dots.required) {
				plot(quantiles, density, type = "n", xlab = xlab,ylab = ylab, lty=col.or.type, ...)
				points(quantiles[1],density[1],...)
				points(quantiles[granularity+1],density[granularity+1],...)
				lines(quantiles[-c(1,granularity+1)], density[-c(1,granularity+1)], col=col.or.type, ...)
				print("hello")
				}
			else {
				plot(quantiles, density, type = "l", xlab = xlab, ylab = ylab, col=col.or.type, ...)
				}
			}
		}
	else {
		# Not a new plot - so the initial plots with type "n" arent needed in the
		# dots required case
		if(bnw) {
			if (dots.required) {
				points(quantiles[1],density[1])
				points(quantiles[granularity+1],density[granularity+1])
				lines(quantiles[-c(1,granularity+1)], density[-c(1,granularity+1)], lty=col.or.type, ...)
				}
			else {
				lines(quantiles, density, lty = col.or.type)
				}
			}
		else {
			if (dots.required) {
				points(quantiles[1],density[1],...)
				points(quantiles[granularity+1],density[granularity+1],...)
				lines(quantiles[-c(1,granularity+1)], density[-c(1,granularity+1)], col=col.or.type, ...)
				}
			else {		
				lines(quantiles, density, col = col.or.type)
				}
			}
		}
	}
if (!is.null(quant.probs)){quantile(quantiles,quant.probs) } 
}

plotglc <- function(lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,
  param="fmkl",lambda5=NULL, granularity=4000, xlab="x",
  ylab="cumulative probability",...)
{
# standard parameter fixin' - copied directly from dgl, but we want the 
# warnings to happen in this function.
# Tidy the parameters so gl.check.lambda will work
lambdas <- .gl.parameter.tidy(lambda1,lambda2,lambda3,lambda4,param,lambda5)
# Check the parameters
if(!gl.check.lambda(lambdas,param=param,vect=TRUE)) {
        stop(paste("The parameter values", lambdas,
"\ndo not produce a proper distribution with the",param,
"parameterisation - see \ndocumentation for gl.check.lambda"))
        }
	u <- seq(from = 1/granularity, to = 1 - 1/granularity, length = 
		granularity - 1)
	x <- qgl(u,lambda1=lambdas,param=param)
	plot(x, u, pch = ".",xlab=xlab,ylab=ylab,...)
}
