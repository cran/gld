\name{GeneralisedLambdaDistribution}
\alias{dgl}
\alias{qdgl}
\alias{qdgl.fmkl}
\alias{qdgl.rs}
\alias{pgl}
\alias{qgl}
\alias{qgl.fmkl}
\alias{qgl.rs}
\alias{rgl}
\alias{GeneralizedLambdaDistribution}
\alias{GeneralisedLambdaDistribution}
\alias{gld}
\title{The Generalised Lambda Distribution}
\description{Density, quantile density, distribution function, quantile 
function and random generation for the generalised lambda distribution 
(also known as the asymmetric lambda, or Tukey lambda).  Works for 
both the \code{fmkl} (recommended) and \code{rs} parameterisations.
}
\usage{
dgl(x, lambda1, lambda2, lambda3, lambda4,
param="fmkl",inverse.eps=1e-8,max.iterations=500)
qdgl(p, lambda1, lambda2, lambda3, lambda4, param="fmkl")
qdgl.fmkl(p, lambda1, lambda2, lambda3, lambda4)
qdgl.rs(p, lambda1, lambda2, lambda3, lambda4)
pgl(q, lambda1, lambda2, lambda3, lambda4,
param="fmkl",inverse.eps=1e-8,max.iterations=500)
qgl(p, lambda1, lambda2, lambda3, lambda4, param="fmkl")
qgl.fmkl(p, lambda1, lambda2, lambda3, lambda4)
qgl.rs(p, lambda1, lambda2, lambda3, lambda4)
rgl(n, lambda1, lambda2, lambda3, lambda4, param="fmkl")
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
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
  \item{inverse.eps}{Accuracy of calculation for the numerical determination of 
	\eqn{F(x)}, defaults to \eqn{10^{-8}}{1e-8}}
  \item{max.iterations}{Maximum number of iterations in the numerical 
	determination of \eqn{F(x)}, defaults to 500}
}
\details{
The generalised lambda distribution, also known as the asymmetric lambda, 
or Tukey lambda distribution, is a distribution with a wide range of shapes.  
The distribution is defined by its quantile function, the inverse of the
distribution function.  There are two parameterisations of the distribution.  
The default parameterisation is that due to \emph{Freimer Mudholkar, Kollia 
and Lin (1988)} (see references below), with a quantile function:
\deqn{F^{-1}(u)= \lambda_1 + { { \frac{u^{\lambda_3}-1}{\lambda_3} - 	%
\frac{(1-u)^{\lambda_4}-1}{\lambda_4} } \over \lambda_2 } 		%
}{F inverse (u) = lambda1 + ( (u^lambda3 -1)/lambda3 - ((1-u)^lambda4 	%
-1)/lambda4 ) / lambda 2}
for \eqn{\lambda_2 > 0}{lambda2 >0}.

The alternative parameterisation, chosen by setting
\code{param="rs"} is that due to \emph{Ramberg and Schmeiser
(1974)}, with the quantile function:
\deqn{F^{-1}(u)= \lambda_1 + \frac{u^{\lambda_3} - (1-u)^{\lambda_4}} 	%
{\lambda_2} }{F inverse (u) = lambda1 + ( u^lambda3 - (1-u)^lambda4 ) / % 
lambda 2}

This parameterisation has a complex series of rules determining which values 
of the parameters produce valid statistical distributions.  See 
\link{gl.check.lambda} for details.

The distribution is defined by its quantile function and its distribution and
density functions do not exist in closed form.  Accordingly, the results 
from \code{pgl} and \code{dgl} are the result of a numerical solutions to
equations, using the Newton-Raphson method.  Since the quantile density
function, \eqn{f(F^{-1}(u))}, does exist, an additional function, \code{qdgl},
computes this.

\code{qdgl.fmkl} and \code{qdgl.rs} are versions of \code{qdgl} that assume the 
FMKL and RS parameterisations, respectively.

\code{qgl.fmkl} and \code{qgl.rs} are versions of \code{qgl} that assume the 
FMKL and RS parameterisations, respectively.
}
\value{
  \code{dgl} gives the density (based on the quantile density and a 
  numerical solution to \eqn{F^{-1}(u)=x}{F inv (u)=x}), 

  \code{qdgl} gives the quantile density,
  
  \code{pgl} gives the distribution function (based on a numerical
	solution to \eqn{F^{-1}(u)=x}{F inv (u)=x}),

  \code{qgl} gives the quantile function, and

  \code{rgl} generates random deviates.
}
\references{ Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988), 
	\emph{A study of the generalized tukey lambda family}, Communications 
	in Statistics - Theory and Methods \bold{17}, 3547--3567.

Karian, Z.E., Dudewicz, E.J., and McDonald, P. (1996), \emph{The extended 
    generalized lambda distribution system for fitting distributions to data: 
    history, completion of theory, tables, applications, the ``Final Word'' 
    on Moment fits}, Communications in Statistics - Simulation and Computation 
    \bold{25}, 611--642.

Ramberg, J. S. & Schmeiser, B. W. (1974), \emph{An approximate method for
generating asymmetric random variables}, Communications of the ACM \bold{17}, 
78--82.

\url{http://maths.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au}, \url{http://maths.newcastle.edu.au/~rking/}}

\examples{
qgl(seq(0,1,0.02),0,1,0.123,-4.3)
pgl(seq(-2,2,0.2),0,1,-.1,-.2,param="fmkl",inverse.eps=.Machine$double.eps)
# calculate the probabilities more accurately than normal
}
\keyword{distribution}
