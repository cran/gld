\name{qqgl}
\alias{qqgl}
\title{Quantile-Quantile plot against the generalised lambda distribution}
\description{\code{qqgl} produces a Quantile-Quantile plot against the
generalised lambda distribution.  It does for the generalised lambda
distribution what \code{\link{qqnorm}} does for the normal.
}
\usage{
qqgl(y, lambda1, lambda2, lambda3, lambda4, param = "fmkl", abline = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The data sample}
\item{lambda1}{\eqn{\lambda_1}{lambda 1} - location parameter
% This note needed to be put inside an item, otherwise R cmd check can't find
% the documentation for the other arguments
  \emph{Note that the numbering of the \eqn{\lambda}{lambda} parameters 
	is different to that used by Freimer, Mudholkar, Kollia and Lin.}
	}
  \item{lambda2}{\eqn{\lambda_2}{lambda 2} - scale parameter}
  \item{lambda3}{\eqn{\lambda_3}{lambda 3} - first shape parameter}
  \item{lambda4}{\eqn{\lambda_4}{lambda 4} - second shape parameter}
  \item{param}{choose parameterisation: 
  \code{fmkl} uses \emph{Freimer, Mudholkar, Kollia and Lin (1988)} (default).
  \code{rs} uses \emph{Ramberg and Schmeiser (1974)}}
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

\url{http://maths.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au},
\url{http://maths.newcastle.edu.au/~rking/}}

\seealso{\code{\link{gld}},\code{\link{starship}} }
\examples{
qqgl(rgl(100,0,1,0,-.1),0,1,0,-.1)
}
\keyword{distribution}
\keyword{hplot}
\keyword{aplot}
