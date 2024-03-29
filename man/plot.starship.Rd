\name{plot.starship}
\alias{plot.starship}
\title{Plots to compare a fitted generalised lambda distribution to data}
\description{Plots to compare a fitted Generalised Lambda Distribution to data.  This works with fitted gld objects from \code{\link{starship}} or  \code{\link{fit.fkml}}.} 
\usage{
\method{plot}{starship}(x, data, ask = NULL, one.page = FALSE, 
    breaks = "Sturges", plot.title = "default",...)
}
\arguments{
  \item{x}{An object of class \code{\link{starship}}.  These are produced by the fitting functions \code{\link{fit.fkml}} and \code{\link{starship}}.}
  \item{data}{Data to which the gld was fitted.  Leave this as NULL if the return.data argument was TRUE in the call that created (\code{\link{starship}} or \code{\link{fit.fkml}})}
  \item{ask}{Ask for user input before next plot.  The default of NULL changes to TRUE if one.page is FALSE and \code{plot} is called interactively, otherwise it changes to FALSE this is then passed to \code{\link{par}(ask)}.  Does not permanently change this setting. The argument is ignored if \code{one.page} is TRUE}
  \item{one.page}{If TRUE, put the two plots on one page using \code{\link{par}(mfrow=c(2,1))}. Does not permanently change this setting.}
  \item{breaks}{Control the number of histogram bins --- passed to \code{\link{hist}}.}
  \item{plot.title}{Main title for histogram and QQ --- passed 
      to \code{\link{hist}(main=)} and \code{\link{qqgl}(main=)}.  
      If you set this to \code{"default"}, it will include the fitting
      method and gld type, for example 
      ``Starship fit of FMKL type GLD''.}
  \item{\dots}{arguments passed to \code{\link{plot}} AND \code{\link{hist}}}
}
\details{
This function produces a Q-Q plot of the data against the fitted generalised lambda distribution and a histogram of the data overlaid with the fitted probability density function.
}
\references{ Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988), 
	\emph{A study of the generalized tukey lambda family}, Communications 
	in Statistics - Theory and Methods \bold{17}, 3547--3567.

Ramberg, J. S. & Schmeiser, B. W. (1974), \emph{An approximate method for
generating asymmetric random variables}, Communications of the ACM \bold{17}, 
78--82.
 
King, R.A.R. & MacGillivray, H. L. (1999), \emph{A starship method for
fitting the generalised \eqn{\lambda}{lambda} distributions},
Australian and New Zealand Journal of
Statistics \bold{41}, 353--374

\url{https://github.com/newystats/gld/}
}
\author{Robert King, \email{robert.king.newcastle@gmail.com},
\url{https://github.com/newystats/}
}
\seealso{
\code{\link{starship}},
}
\examples{
data <- rgl(100,0,1,.2,.2)
starship.result <- starship(data,optim.method="Nelder-Mead",initgrid=list(lcvect=(0:4)/10,
ldvect=(0:4)/10),return.data=TRUE)
plot(starship.result)
}
\keyword{distribution}


