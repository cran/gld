rgl <- function(n,lambda1=0,lambda2=1,lambda3,lambda4,param="fmkl")
{
# Check the parameters
if(!gl.check.lambda(lambda1,lambda2,lambda3,lambda4,param)) {
        stop(paste("The parameter values", lambda1, lambda2, lambda3, lambda4,
"\ndo not produce a proper distribution with the",param,
"parameterisation - see \ndocumentation for gl.check.lambda"))
        }
# Produce the uniform data
p <- runif(n)
# convert to gl
res <- qgl(p,lambda1,lambda2,lambda3,lambda4,param)
res
}
