gl.check.lambda _ function(lambda1,lambda2,lambda3,lambda4,param="fmkl")
# Checks to see that the lambda values given are allowed.
{
# Check all the parameters are finite
if (sum(is.finite(c(lambda1,lambda2,lambda3,lambda4)))<4) 
	{ return(FALSE)
	}
param <- switch(param,  
# Different tests apply for each parameterisation
	freimer=,  # allows for alternate expressions
	frm=,  # allows for alternate expressions
	FMKL=,
	fmkl={
	if (lambda2<=0) {ret _ FALSE}
	else {ret _ TRUE}
	},
	ramberg=, # Ramberg & Schmeiser
	ram=,
	RS=,
	rs={
	if (lambda3*lambda4>0) { # regions 3 and 4 
				 # all values of lambda 3 and lambda 4 OK
				 # check lambda 2
		if ((lambda3>0)&(lambda4>0)) { # region 3 - l2 >0
			if (lambda2<=0) {ret _ FALSE}
			else {ret _ TRUE}
			}
		if ((lambda3<0)&(lambda4<0)) { # region 4 - l2 <0
			if (lambda2>=0) {ret _ FALSE}
			else {ret _ TRUE}
			}
		}	
	else { 	# other quadrants - lambda 2 must be negative, and lambda3 
		# lambda 4 values need checking.
		if (lambda2>=0) {return(FALSE)}
		# Rectangular regions where RS is not defined 
		if ((lambda3>0)&(lambda3<1)&(lambda4<0)) {return(FALSE)}
		if ((lambda4>0)&(lambda4<1)&(lambda3<0)) {return(FALSE)}
		# Different here because there are a 
		# number of ways in which the parameters can fail.
		# 
		# Curved regions where RS is not defined
		# change to shorter var names
		lc _ lambda3
		ld _ lambda4
		if ((lambda3>-1)&(lambda3<0)&(lambda4>1)) {  # region 5 or not?
			if ( ((1-lc)^(1-lc)*(ld-1)^(ld-1))/((ld-lc)^(ld-lc)) > -lc/ld )	
				{return(FALSE)}
			else 	{return(TRUE)}
			}
		# Second curved region 
		if ((lambda4>-1)&(lambda4<0)&(lambda3>1)) {  # region 6 or not?
			if ( ((1-ld)^(1-ld)*(lc-1)^(lc-1))/((lc-ld)^(lc-ld)) > -ld/lc )
				{return(FALSE)}
			else 	{return(TRUE)}
			}
		# There may be some limit results that mean these are not correct, but
		# I'll check that later
		if (lambda3 == 0) {
			warning('lambda 3 = 0 - could be a problem')
			return(FALSE)
			}
		if (lambda4 == 0) {
			warning('lambda 4 = 0 - could be a problem')
			return(FALSE)
			}
		# If we get here, then the parameters are OK.
		ret _ TRUE
		}
	},
	stop("Error when checking validity of parameters.\n Parameterisation must be either fmkl or rs")
	) # closes "switch"
ret
}


qgl.fmkl _ function(p,lambda1,lambda2,lambda3,lambda4)
{
u <- p
# Check the values are OK)
if(!gl.check.lambda(lambda1,lambda2,lambda3,lambda4,param="fmkl")) {
	stop("illegal value for one of the parameters - see documentation for gl.check.lambda")
	}
# If OK, determine special cases
if (lambda3 == 0) { 
	if (lambda4 == 0) { # both log
		quants <- lambda1 + (log(u) - log(1 - u))/lambda2
		}
	else	{ # l3 zero, l4 non-zero
		quants <- lambda1 + 
			(log(u) - ((1 - u)^lambda4 - 1)/lambda4)/lambda2
		}
	}
else 	{ # lambda3 non-zero
	if (lambda4 == 0) { # non-zero, l4 zero
		quants <- lambda1 + 
			((u^lambda3 - 1)/lambda3 - log(1 - u))/lambda2
		}
	else	{ # both non-zero - use usual formula
		quants _ lambda1 + ( ( u ^ lambda3 - 1)	/ lambda3 
			- ( (1-u)^lambda4 - 1) /lambda4 ) / lambda2
		}
	}
quants
}

qgl.rs _ function(p,lambda1,lambda2,lambda3,lambda4)
{
u <- p
# Check the values are OK)
if(!gl.check.lambda(lambda1,lambda2,lambda3,lambda4,param="rs")) {
	stop("illegal value for one of the parameters - see documentation for gl.check.lambda")
	}
# At present, I'm rejecting zero values for l3 and l4, though I think there 
# are limit results, so one functional form.
quants _ lambda1 + ( u ^ lambda3 - (1-u)^lambda4 ) / lambda2
quants
}

qgl _ function(p,lambda1,lambda2,lambda3,lambda4,param="fmkl")
{
u <- p
result <- switch(param,  
# Different tests apply for each parameterisation
	freimer=,  # allows for alternate expressions
	frm=,  # allows for alternate expressions
	FMKL=,
	fmkl=qgl.fmkl(u,lambda1,lambda2,lambda3,lambda4),
	ramberg=, # Ramberg & Schmeiser
	ram=,
	RS=,
	rs=qgl.rs(u,lambda1,lambda2,lambda3,lambda4),
	stop("Error: Parameterisation must be either fmkl or rs")
	) # closes "switch"
result
}

qdgl _ function(p,lambda1,lambda2,lambda3,lambda4,param="fmkl")
{
u <- p
result <- switch(param,  
# Different tests apply for each parameterisation
	freimer=,  # allows for alternate expressions
	frm=,  # allows for alternate expressions
	FMKL=,
	fmkl=qdgl.fmkl(u,lambda1,lambda2,lambda3,lambda4),
	ramberg=, # Ramberg & Schmeiser
	ram=,
	RS=,
	rs=qdgl.fmkl(u,lambda1,lambda2,lambda3,lambda4),
	stop("Error: Parameterisation must be either fmkl or rs")
	) # closes "switch"
result
}


qdgl.rs _ function(p,lambda1=0,lambda2=1,lambda3,lambda4)
{
u <- p
# Check the values are OK)
if(!gl.check.lambda(lambda1,lambda2,lambda3,lambda4,param="rs")) {
	stop("illegal value for one of the parameters - see documentation for gl.check.lambda")
	}
dens _  lambda1/(lambda3 * (u^(lambda3 -1)) + lambda4 * ((1 - u)^(lambda4 -1)))
dens
}


qdgl.fmkl _ function(p,lambda1,lambda2,lambda3,lambda4)
{
u <- p
# Check the values are OK)
if(!gl.check.lambda(lambda1,lambda2,lambda3,lambda4,param="fmkl")) {
	stop("illegal value for one of the parameters - see documentation for gl.check.lambda")
	}
# The density is given by 1/Q'(u)
dens _ lambda2/(u^(lambda3 - 1) + (1 - u)^(lambda4 - 1))
dens
}
