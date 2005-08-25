\name{qqgl}
\alias{qqgl}
\title{Quantile-Quantile plot against the generalised lambda distribution}
\description{\code{qqgl} produces a Quantile-Quantile plot of data against the
generalised lambda distribution.  It does for the generalised lambda
distribution what \code{\link{qqnorm}} does for the normal.
}
\usage{
qqgl(y, lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,
  param="fmkl",lambda5=NULL, abline = TRUE, ...)
}
\arguments{
  \item{y}{The data sample}
  \item{lambda1}{This can be either a single numeric value or a vector.  
  
  If it is a vector, it must be of length 4 for parameterisations 
  \code{fmkl} or \code{rs} and of length 5 for parameterisation \code{fm5}.  
  If it is a vector, it gives all the parameters of the generalised lambda 
  distribution (see below for details) and the other \code{lambda} arguments
  must be left as NULL.
        
 If it is a a single value, it is \eqn{\lambda_1}{lambda 1}, the location 
 parameter of the distribution and the other parameters are given by the
 following arguments
 
 \emph{Note that the numbering of the \eqn{\lambda}{lambda} parameters for
 the \emph{fmkl} parameterisation is different to that used by Freimer, 
 Mudholkar, Kollia and Lin.}
	}
  \item{lambda2}{\eqn{\lambda_2}{lambda 2} - scale parameter}
  \item{lambda3}{\eqn{\lambda_3}{lambda 3} - first shape parameter}
  \item{lambda4}{\eqn{\lambda_4}{lambda 4} - second shape parameter}
  \item{lambda5}{\eqn{\lambda_5}{lambda 5} - a skewing parameter, in the
    fm5 parameterisation}
  \item{param}{choose parameterisation: 
  \code{fmkl} uses \emph{Freimer, Mudholkar, Kollia and Lin (1988)} (default).
  \code{rs} uses \emph{Ramberg and Schmeiser (1974)}
  \code{fm5} uses the 5 parameter version of the FMKL parameterisation 
  (paper to appear)}
  \item{abline}{A logical value, TRUE adds a line through the origian
with a slope of 1 to the plot}
  \item{\dots}{graphical parameters}
}
\details{
See \code{\link{gld}} for more details on the Generalised Lambda
Distribution.  A Q-Q plot provides a way to visually assess the
correspondence between a dataset and a particular distribution. 
}
\value{
  A list of the same form as that returned by \code{\link{qqline}}
  \item{x}{The x coordinates of the points that were/would be plotted,
corresponding to a generalised lambda distibution with parameters 
\eqn{\lambda_1}{lambda 1}, \eqn{\lambda_2}{lambda 2},
\eqn{\lambda_3}{lambda 3}, \eqn{\lambda_4}{lambda 4}. }
  \item{y}{The original \code{y} vector, i.e., the corresponding y
    coordinates.}
}
\references{
King, R.A.R. & MacGillivray, H. L. (1999), \emph{A starship method for
fitting the generalised \eqn{\lambda}{lambda} distributions},
Australian and New Zealand Journal of
Statistics \bold{41}, 353--374

\url{http://tolstoy.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{http://tolstoy.newcastle.edu.au/~rking/}}

\seealso{\code{\link{gld}},\code{\link{starship}} }
\examples{
qqgl(rgl(100,0,1,0,-.1),0,1,0,-.1)
}
\keyword{distribution}
\keyword{hplot}
\keyword{aplot}
