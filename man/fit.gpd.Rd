\name{fit.gpd}
\alias{fit.gpd}
\alias{fit.gpd.lmom}
\alias{fit.gpd.lmom.given} 
\title{Estimate parameters of the GPD type generalised lambda distribution}
\description{Estimates parameters of the GPD type generalised \eqn{\lambda}{lambda} Distribution.  Currently, only estimation via method of L moments is implemented.

The Method of L-Moments estimates for the GPD type are the only estimates for any
generalised lambda distribution type with closed form expressions.
}
\usage{
fit.gpd(x, method = "LM", na.rm = TRUE, record.cpu.time = FALSE, return.data = FALSE)
fit.gpd.lmom(data, na.rm = TRUE)
fit.gpd.lmom.given(lmoms, n = NULL)
}
\arguments{
  \item{x}{Data to be fitted, as a vector}
  \item{method}{A character string, to select the estimation method.  Only Method of
  L-Moments "LM" is implemented.}
  \item{na.rm}{Logical: Should missing values be removed?}
  \item{record.cpu.time}{Logical: should the CPU time used in fitting be recorded in the fitted model object?}
  \item{return.data}{Logical: Should the function return the data (from the argument \code{x})?}
  \item{data}{Data to be fitted, as a vector}
  \item{lmoms}{A numeric vector containing two L-moments and two L-moment ratios, in the order \eqn{l_1,l_2,t_3,t_4}{l_1, l_2, t_3, t_4}.}
  \item{n}{the sample size, defaults to NULL}
}
\details{The method of L-Moments equates sample L-Moments with expressions 
for the L-Moments of the GPD type GLD.  Closed form expressions exist to give these estimates.

For many values there are two possible estimates for the same L Moment
values, one in each of two regions of the GPD GLD parameter space, 
denoted region A and region B in van Staden (2013).  More details on these 
regions can be found on page 154 of van Staden (2013).

If the 4th L-Moment ratio, \eqn{\tau_4}{tau4} is less than the minimum 
value that \eqn{\tau_4}{tau4} can obtain for the GPD generalised lambda
distribution;

\deqn{\tau_4^{(min)}=\frac{12-5\sqrt{6}}{12+5\sqrt{6}} \approx -0.0102051,}{(12-5*sqrt(6))/(12+5*sqrt(6)) = approx -0.0102051,}

there is no possible estimate (from either region A or B), and this 
function returns NA for the estimates.
}
\value{
These funciton return an object of \code{\link[base]{class}} 
"\code{GldGPDFit}".  It is a list, containing these components (optional components noted here);
  \item{estA}{The estimate in region A.  This will be NA if there is no estimate in region A}
  \item{estB}{The estimate in region B.  This will be NA if there is no estimate in region B}
  \item{warn}{(only if estA and estB are both NA), the reason there are no
estimates.  If this is the case, the function also issues a \code{\link{warning}}.}
  \item{cpu}{A vector showing the computing time used, returned if 
    \code{record.cpu.time} is TRUE (only for \code{fit.gpd}).}     
  \item{data}{The data, if \code{return.data} is TRUE (only for \code{fit.gpd}).}
  \item{param}{The character \code{"gpd"}, indicating the GPD type of the
  generalised lambda distribution.}

  Each of the estimate elements (if they are not NA) are a vector of length 4,
  giving the estimated parameters, in order;

  \eqn{\alpha}{alpha} - location parameter
  \eqn{\beta}{beta} - scale parameter
  \eqn{\delta}{delta} - skewness parameter
  \eqn{\lambda}{lambda} - kurtosis parameter
}
\references{
Van Staden, Paul J., & M.T. Loots. (2009), \emph{Method of L-moment Estimation for the Generalized Lambda Distribution}. In Proceedings of the Third Annual ASEARC Conference. Callaghan, NSW 2308 Australia: School of Mathematical and Physical Sciences, University of Newcastle.

}
\author{Robert King, \email{robert.king@newcastle.edu.au}, \url{http://tolstoy.newcastle.edu.au/~rking/}

Paul van Staden}
\seealso{
\code{\link{GeneralisedLambdaDistribution}}
}
\examples{
fit.gpd.lmom.given(c(1,.3,.6,.8))
}
\keyword{distribution}