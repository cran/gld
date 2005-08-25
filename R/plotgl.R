plotgld <- function(lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,
  param="fmkl",lambda5=NULL, new.plot = TRUE, truncate = 0, bnw = FALSE,
  col.or.type = 1, granularity = 4000,xlab = "x", ylab=NULL, 
  quant.probs = seq(0,1,.25), ...)
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
u <- seq(from = 0, to = 1, by = 1/granularity)
# Only difference across parameterisations is calculating the 
# quantiles and density
quantiles <- qgl(u,lambda1=lambdas,param=param)
density <- qdgl(u,lambda1=lambdas,param=param)
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
