dgl <- function(x,lambda1=0,lambda2=1,lambda3,lambda4,param="fmkl",
  inverse.eps=1e-8,max.iterations=500)
{
# Check the parameters
if(!gl.check.lambda(lambda1,lambda2,lambda3,lambda4,param)) {
        stop(paste("The parameter values", lambda1, lambda2, lambda3, lambda4,
"\ndo not produce a proper distribution with the",param,
"parameterisation - see \ndocumentation for gl.check.lambda"))
        }
# calculate u=F(x) numerically, then use qdgl
u <- pgl(x,lambda1,lambda2,lambda3,lambda4,param,inverse.eps,
  max.iterations)
dens <- qdgl(u,lambda1,lambda2,lambda3,lambda4,param)
dens
}


pgl <- function(q,lambda1=0,lambda2=1,lambda3,lambda4,param="fmkl",
  inverse.eps=1e-8,max.iterations=500)
{
# Check the parameters
if(!gl.check.lambda(lambda1,lambda2,lambda3,lambda4,param)) {
        stop(paste("The parameter values", lambda1, lambda2, lambda3, lambda4,
"\ndo not produce a proper distribution with the",param,
"parameterisation - see \ndocumentation for gl.check.lambda"))
        }
# F(x) needs to be calculated numerically
# We need to tell the C program some information
length.of.vector <- length(q)
# Check to see that the given quantiles are possible with these values of the 
# parameters
extremes <- qgl(c(0,inverse.eps,1-inverse.eps,1),lambda1,lambda2,lambda3,lambda4,param)
# These should really be warnings. It should allow values out of range, and 
# give them density=0 and prob = 0 or 1 respectively.  Problem is, how
# do I get the right places kept in the vector
if ( min(q) < extremes[1] ) {stop(paste("Smallest 'x' value: ",min(q),
  " Theoretical minimum possible value for gen lambda with parameters",
  lambda1,lambda2,lambda3,lambda4," is ",extremes[1]))}
if ( min(q) < extremes[2] ) {stop(paste("Smallest 'x' value: ",min(q),
  " Calculated quantile for u=",inverse.eps,
  "\n(the requested accuracy for numerical calculation of F(x))",
  "\nfor gen lambda with parameters\n",lambda1,lambda2,lambda3,lambda4," is ",
  extremes[2]))}
if ( max(q) > extremes[4] ) {stop(paste("Largest 'x' value: ",max(q),
  " Theoretical possible value for gen lambda with parameters",
  lambda1,lambda2,lambda3,lambda4," is ",extremes[4]))}
if ( max(q) > extremes[3] ) {stop(paste("Largest 'x' value: ",min(q),
  " Calculated quantile for u=",1-inverse.eps,
  "\n(1 minus the requested accuracy for numerical calculation of F(x))",
  "\nfor gen lambda with parameters\n",lambda1,lambda2,lambda3,lambda4," is ",
  extremes[3]))}
# Need a blank u to send them.
u <- q*0
result <- switch(param,
        freimer=,  # allows for alternate expressions
        frm=,  # allows for alternate expressions
        FMKL=,
        fmkl=.C("gl_fmkl_distfunc",lambda1,lambda2,lambda3,lambda4,
		as.double(0),as.double(1),inverse.eps,as.integer(max.iterations),
		as.double(q),as.double(u),as.integer(length.of.vector)),
        ramberg=, # Ramberg & Schmeiser
        ram=,
        RS=,
        rs=.C("gl_rs_distfunc",lambda1,lambda2,lambda3,lambda4,
		as.double(0),as.double(1),inverse.eps,max.iterations,
		as.double(q),as.double(u),as.integer(length.of.vector)),
        stop("Error: Parameterisation must be either fmkl or rs")
        ) # closes "switch"
if (!(is.numeric(result[[1]]))){ 
	stop("Values for quantiles outside range.  This shouldn't happen") 
	}
u <- result[[10]]
u
}
