# $Id: starship.R,v 1.3 2003-05-13 16:28:27+10 king Exp $
starship <- function(data,optim.method="Nelder-Mead",initgrid=NULL,
	param="FMKL",optim.control=NULL) {
# call adaptive grid first to find a first minimum
if (is.null(initgrid) ) {gridmin <- starship.adaptivegrid(data,param=param) }
else 	{ 
	warning("No checks for grids implemented")
	gridmin <- starship.adaptivegrid(data,initgrid$lcvect,initgrid$ldvect, 
		param=param)
	}
# If they haven't sent any control parameters, scale by max(lambda1,1),
# lambda2 (can't be <= 0), don't scale for lambda3, lambda4
if (is.null(optim.control) ) {
	optim.control <- list(
		parscale=c(max(1,abs(gridmin$lambda[1])),abs(gridmin$lambda[2]),
		1,1))
	}
# else use what they sent - this should allow them to change other stuff in the
# control while keeping parscale 
# call optimiser
optimmin <- optim(par=gridmin$lambda,fn=starship.obj,method=optim.method,
	control=optim.control,data=data,param=param)
list(lambda=optimmin$par,grid.results=gridmin,optim.results=optimmin)
}

starship.adaptivegrid <- function(data,lcvect=c(-1.5,-1,-.5,-.1,0,.1,.2,.4,
	.8,1,1.5), ldvect=c(-1.5,-1,-.5,-.1,0,.1,.2,.4,.8,1,1.5),param="FMKL") 
{
data <- sort(data)
quarts <- quantile(data)
nombo <- length(data)
minresponse <- 1000
minlambda <- c(NA,NA,NA,NA)

for (ld in ldvect) {
	for (lc in lcvect) {
		# calculate expected lambda2 on basis of IQR
		# IQR for 0,1,lc,ld
		iqr1 <- qgl(.75,0,1,lc,ld,param=param) -
qgl(.25,0,1,lc,ld,param=param)
		# actual IQR
		iqr <- quarts[4] - quarts[2]
		# so estimated lambda 2 from IQR
		lbguess <- iqr1/iqr
		for (lb in (c(0.5,0.7,1,1.5,2)*lbguess) )
			{
			# calculate expected lambda1 on the basis of median
			lavect <- quarts[3] -
c(qgl(0.65,0,lb,lc,ld,param=param),
			qgl(0.55,0,lb,lc,ld,param=param),qgl(0.5,0,lb,lc,ld,param=param), 
			qgl(0.45,0,lb,lc,ld,param=param),qgl(0.35,0,lb,lc,ld,param=param) )
			for (la in lavect) {
				# calculate uniform g-o-f
				response <- starship.obj(c(la,lb,lc,ld),data,param)
				if (response < minresponse) {
					minresponse <- response
					minlambda <- c(la,lb,lc,ld)
				} # new minimum
				# otherwise try the next
			} #lavect
		} #lbvect
	} # lcvct
} # ldvect
list(response=minresponse,lambda=minlambda)
}

starship.obj <- function(par,data,param="fmkl") 
{
l1 <- par[1]; l2 <- par[2];
l3 <- par[3]; l4 <- par[4];

# Check that these are legitimate parameter values.  If not, give a very
# large internal g-o-f measure.  We do this instead of NAs to make the
# optimistations easy to code.  Should investigate using NAs instead
if (!gl.check.lambda(l1, l2, l3, l4, param)) {
	return(54321)
	}

x <- sort(data)
# defining other variables
nombo <- length(x)

# MAIN FUNCTION This was the original guts of the C starship program (King
# and MacGillivray 1999)

xacc <- 1e-8
# values sent to pgl
# lower & upper limit on u
x1 <- xacc 
x2 <- 1.0 - xacc

u <- pgl(x,l1,l2,l3,l4,param=param);
# write.table(matrix(c(u),nrow=1),"interim-output/u1.txt",append=T,sep=",",quote=F,col.names=F)
response <- .anddarl(u,nombo);
# write.table(matrix(c(response),nrow=1),"interim-output/response.txt",append=T,sep=",",quote=F,col.names=F)
return(response)
}

# ANDERSON-DARLING TEST
.anddarl <- function(u,nombo)
{
ad <- c(dim(nombo),dim(1));
for(j in 1:nombo) {ad[j] <- ((2*j-1)/nombo)*(log(u[j])+log(1-u[nombo+1-j]))}
# This is useful in opt to illustrate what its doing
# hist(u)
asum <- -nombo-sum(ad);
asum;
}
# END OF ANDERSON-DARLING TEST
