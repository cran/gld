dgl <- function(x,lambda1=0,lambda2=1,lambda3,lambda4,param="fmkl",
  inverse.eps=1e-8,max.iterations=500)
{
# Check the parameters
if(!gl.check.lambda(lambda1,lambda2,lambda3,lambda4,param)) {
        stop("illegal value for one of the parameters - see documentation for gl.check.lambda")
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
        stop("illegal value for one of the parameters - see documentation for gl.check.lambda")
        }
# F(x) needs to be calculated numerically
# We need to tell the C program some information
length.of.vector <- length(q)
# Check to see that the given quantiles are possible with these values of the 
# parameters
extremes <- qgl(c(0,1),lambda1,lambda2,lambda3,lambda4,param)
# These should really be warnings. It should allow values out of range, and 
# give them density=0 and prob = 0 or 1 respectively.  Problem is, how
# do I get the right places kept in the vector
if ( min(q) < extremes[1] ) {stop(paste("Smallest data value: ",min(q),
  " Minimum possible value for gen lambda with parameters",
  lambda1,lambda2,lambda3,lambda4," is ",extremes[1]))}
if ( max(q) > extremes[2] ) {stop(paste("Largest data value: ",max(q),
  " Maximum possible value for gen lambda with parameters",
  lambda1,lambda2,lambda3,lambda4," is ",extremes[2]))}
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
