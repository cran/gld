\name{gl.check.lambda}
\alias{gl.check.lambda}
\title{Function to check the validity of parameters of the generalized lambda 
distribution}
\description{
  Checks the validity of parameters of the generalized lambda.  In the case of 
  the FMKL parameterisation, this is just a case of checking that 
  \eqn{\lambda_2}{lambda 2} is positive.  The RS parameterisation is much more 
  complex.
}
\usage{
gl.check.lambda(lambda1, lambda2, lambda3, lambda4, param="fmkl")
}
\arguments{
  \item{lambda1}{\eqn{\lambda_1}{lambda 1} - location parameter
% This note needed to be put inside an item, other R cmd check can't find
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
}
\details{
  See \code{\link{GeneralisedLambdaDistribution}} for details on the 
  generalised lambda distribution.  This function determines the validity of 
  parameters of the distribution.
  
  The FMKL parameterisation gives a valid 
  statistical distribution for any real values of \eqn{\lambda_1}{lambda 1},
  \eqn{\lambda_3}{lambda 3},\eqn{\lambda_4}{lambda 4} and any positive real 
  values of \eqn{\lambda_2}{lambda 2}.

  For the RS parameterisation, the combinations of parameters value that give 
  valid distributions are the following (the region numbers in the table 
  correspond to the labelling of the regions in \emph{Ramberg and 
  Schmeiser (1974)} and \emph{Karian, Dudewicz and McDonald (1996)}):

\tabular{cccccc}{
region \tab \eqn{\lambda_1}{lambda 1} \tab \eqn{\lambda_2}{lambda 2}	
\tab \eqn{\lambda_3}{lambda 3} \tab \eqn{\lambda_4}{lambda 4} \tab note \cr
1 \tab	all 	\tab \eqn{<0}	\tab \eqn{< -1}		\tab \eqn{> 1} \tab \cr
2 \tab	all 	\tab \eqn{<0}	\tab \eqn{> 1}		\tab \eqn{< -1} \tab \cr
3 \tab	all 	\tab \eqn{>0}	\tab \eqn{\ge 0}	\tab \eqn{\ge 0} \tab 
one of \eqn{\lambda_3}{lambda 3} and \eqn{\lambda_4}{lambda 4} must be non-zero \cr
4 \tab  all     \tab \eqn{<0} 	\tab \eqn{\le 0}	\tab \eqn{\le 0} \tab 
one of \eqn{\lambda_3}{lambda 3} and \eqn{\lambda_4}{lambda 4} must be non-zero \cr
5 \tab 	all	\tab \eqn{<0}	\tab \eqn{> -1} and \eqn{< 0} 	\tab \eqn{>1}
\tab equation 1 below must also be satisfied \cr
6 \tab 	all	\tab \eqn{<0}	\tab \eqn{>1} 	\tab \eqn{> -1} and \eqn{< 0}
\tab equation 2 below must also be satisfied \cr
}

Equation 1

\deqn{  \frac{(1-\lambda_3) ^{1-\lambda_3}(\lambda_4-1)^{\lambda_4-1}}
{(\lambda_4-\lambda_3)^{\lambda_4-\lambda_3}} <
- \frac{\lambda_3}{\lambda_4} }{% non TeX version
( (1-lambda3) ^ ( 1 - lambda3) * (lambda4 -1) ^ (lambda4 -1) ) /
( (lambda4 - lambda3) ^ (lambda4 - lambda3) ) <
- lambda3 / lambda 4 }

Equation 2

\deqn{  \frac{(1-\lambda_4) ^{1-\lambda_4}(\lambda_3-1)^{\lambda_3-1}}
{(\lambda_3-\lambda_4)^{\lambda_3-\lambda_4}} <
- \frac{\lambda_4}{\lambda_3} }{%
( (1-lambda4) ^ ( 1 - lambda4) * (lambda3 -1) ^ (lambda3 -1) ) /
( (lambda3 - lambda4) ^ (lambda3 - lambda4) ) <
- lambda4 / lambda 3
}

}
\value{This logical function takes on a value of TRUE if the parameter values 
given produce a valid statistical distribution and FALSE if they don't}
\references{ Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988),
        \emph{A study of the generalized tukey lambda family}, Communications
       in Statistics - Theory and Methods \bold{17}, 3547--3567.

Karian, Z.E., Dudewicz, E.J., and McDonald, P. (1996), \emph{The extended 
    generalized lambda distribution system for fitting distributions to data:
    history, completion of theory, tables, applications, the ``Final Word'' 
    on Moment fits}, Communications in Statistics - Simulation and
Computation 
    \bold{25}, 611--642.

Ramberg, J. S. & Schmeiser, B. W. (1974), \emph{An approximate method for
generating asymmetric random variables}, Communications of the ACM \bold{17},
78--82.

\url{http://maths.newcastle.edu.au/~rking/gld/}
}
\author{Robert King, \email{robert.king@newcastle.edu.au}, \url{http://maths.newcastle.edu.au/~rking/}}
\note{The complex nature of the rules in this function for the RS
parameterisation are the reason for the invention of the FMKL
parameterisation and its status as the default parameterisation in the other
generalized lambda functions.}

\seealso{The generalized lambda functions \code{\link{GeneralisedLambdaDistribution}}  }

\examples{
gl.check.lambda(0,1,.23,4.5) ## TRUE
gl.check.lambda(0,-1,.23,4.5) ## FALSE 
gl.check.lambda(0,1,0.5,-0.5,"rs") ## FALSE
}
\keyword{distribution}
