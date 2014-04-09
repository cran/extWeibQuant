\name{simWbMix}
\alias{simWbMix}
\title{Simulate data from a mixture of two Weibull distributions}
\description{Same as title.
}
\usage{
simWbMix(n, mixParm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size}
  \item{mixParm}{A two-by-three matrix describing the Weibull mixture
  
  p, a_1, b_1;
  
  1-p, a_2, b_2
}
}
\value{A vector of length n}

\seealso{
\code{\link{rweibull}}, \code{\link{quanWbMix}}, \code{\link{emCenWbMix.T1}}, \code{\link{emCenWbMix.T2}}
}
\examples{
mmix = rbind(c(0.7, 5, 7),c(0.3, 15, 6))
set.seed(1)
y <- sort(simWbMix(300, mmix)) #Generate the data
hist(y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Weibull}
\keyword{Mixture}% __ONLY ONE__ keyword per line
