qqgl <- function (y, lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,
  param="fmkl",lambda5=NULL,abline=TRUE,...)
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
n <- length(y)
u <- seq(from = 1/(n + 1), by = 1/(n + 1), length = n)
q <- qgl(u, lambda1=lambdas, param=param)
if(abline) { 
	ret <- qqplot(q,y,...)
	abline(0,1)
	}
else {ret <- qqplot(q,y,...)}
invisible(ret)
}
