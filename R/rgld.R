rgl <- function(n,lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,param="fkml",lambda5=NULL)
{
# Check the parameters
lambdas <- .gl.parameter.tidy(lambda1,lambda2,lambda3,lambda4,param,lambda5)
# Check the parameters
if(!gl.check.lambda(lambdas,param=param,vect=TRUE)) {
  stop(paste("The parameter values (", paste(lambdas,collapse=", "),
             ")\ndo not produce a proper distribution with the ",param,
             " parameterisation - see \ndocumentation for gl.check.lambda",sep=""))
}
# Produce the uniform data
p <- runif(n)
# convert to gl
res <- qgl(p,lambda1=lambdas,param=param)
res
}
