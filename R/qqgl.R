qqgl <- function (y, lambda1, lambda2, lambda3, lambda4, param = "fmkl",abline=TRUE,...) 
{
	n <- length(y)
	u <- seq(from = 1/(n + 1), by = 1/(n + 1), length = n)
    	q <- switch(param, freimer = , frm = , FMKL = , fmkl = qgl.fmkl(u, 
       	 lambda1, lambda2, lambda3, lambda4), ramberg = , ram = , 
        RS = , rs = qgl.rs(u, lambda1, lambda2, lambda3, lambda4), 
        stop("Error: Parameterisation must be either fmkl or rs"))
	if(abline) { 
		ret <- qqplot(q,y,...)
		abline(0,1)
		}
	else {ret <- qqplot(q,y,...)}
invisible(ret)
}
