rgl <- function(n,lambda1=0,lambda2=1,lambda3,lambda4,parameterisation="fmkl")
{
# Check the parameters
if(!gl.check.lambda(lambda1,lambda2,lambda3,lambda4,parameterisation)) {
        stop("illegal value for one of the parameters - see documentation for gl.check.lambda")
        }
# Produce the uniform data
u _ runif(n)
# convert to gl
res _ qgl(u,lambda1,lambda2,lambda3,lambda4,parameterisation)
res
}
