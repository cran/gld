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
# Unless x is outside the range, then density should be zero
extreme<-qgl(c(0,1),lambda1,lambda2,lambda3,lambda4,param)
# It may be better to change this to simply  
# (x <= extreme[2])*(x >= extreme[1])
outside.range <- !as.logical((x<=extreme[2])*(x>=extreme[1]))
u <- pgl(x,lambda1,lambda2,lambda3,lambda4,param,inverse.eps,
  max.iterations)
dens <- qdgl(u,lambda1,lambda2,lambda3,lambda4,param)
dens[outside.range] <- 0
dens
}

pgl <- function(q,lambda1=0,lambda2=1,lambda3,lambda4,param="fmkl",
    inverse.eps=1e-8,max.iterations=500)
{
# Thanks to Steve Su, <s.su@qut.edu.au>, for improvements to this code
# Check the parameters
if(!gl.check.lambda(lambda1,lambda2,lambda3,lambda4,param)) {
	stop(paste("The parameter values", lambda1, lambda2, lambda3,lambda4,
    	"\ndo not produce a proper distribution with the",param,
    	"parameterisation - see \ndocumentation for gl.check.lambda"))
    	} 
jr <- q; jr[sort.list(q)] <- seq(along=q) 
order.x<-order(q) 
xx<-sort(q) 
# xx has the sorted data, and jr & order.x the information to get back to the
# original order.
extreme<-qgl(c(inverse.eps,1-inverse.eps),lambda1,lambda2,lambda3,lambda4,param)
max.e<-extreme[2]
min.e<-extreme[1]
ind.min<-xx<=min.e
ind.max<-xx>=max.e 
# This simpler comparison works here because we are using inverse.eps as our
# tolerance
q<-xx[as.logical((xx<max.e)*(xx>min.e))] 
# We only want to calculate the probabilities for q values inside the support
length.of.vector <- length(q) 
# Need a blank u to send to C
u <- 0*q 
result <- switch(param, 
	freimer=, # allows for alternate expressions 
	frm=, # allows for alternate expressions 
	FMKL=, 
	fmkl=.C("gl_fmkl_distfunc",lambda1,lambda2,lambda3,lambda4, 
		as.double(0),as.double(1),inverse.eps,
		as.integer(max.iterations),as.double(q),as.double(u),
		as.integer(length.of.vector),PACKAGE="gld"), 
    	ramberg=, # Ramberg & Schmeiser 
    	ram=, 
    	RS=, 
    	rs=.C("gl_rs_distfunc",lambda1,lambda2,lambda3,lambda4, 
    		as.double(0),as.double(1),inverse.eps,max.iterations, 
    		as.double(q),as.double(u),as.integer(length.of.vector),
		PACKAGE="gld"), 
    	stop("Error: Parameterisation must be either fmkl or rs") 
    	) # closes "switch" 
if (!(is.numeric(result[[1]]))){ 
	stop("Values for quantiles outside range. This shouldn't happen") 
} 
u <- result[[10]] 
xx[as.logical((xx<max.e)*(xx>min.e))]<-u 
xx[ind.min]<-0 
xx[ind.max]<-1 
# Revert to the original order of the dataset: 
xx[jr] 
} 
