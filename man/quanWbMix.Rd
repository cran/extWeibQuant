\name{quanWbMix}
\alias{quanWbMix}
\title{
Calculate the quantiles of a mixture of two Weibull distributions.
}
\description{This is a function based on the \code{\link{uniroot}} to calculate the quantiles of a Weibull mixture model. 
}
\usage{
quanWbMix(intProb, mixParm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intProb}{A vector of the probability below the quantile of interest, e.g. 5\%}
  \item{mixParm}{A two-by-three matrix describing the Weibull mixture. See \code{\link{simWbMix}}.}
}

\value{
A two-row matrix: the first row is \code{intProb} and the second row is the corresponding quantiles.
}

\seealso{
\code{\link{qweibull}}, \code{\link{simWbMix}}, \code{\link{emCenWbMix.T1}}, \code{\link{emCenWbMix.T2}}
}
\examples{
mmix = rbind(c(0.7, 5, 7),c(0.3, 15, 6))
quanWbMix(c(0.1, 0.5, 0.7), mmix)
quanWbMix(0.1, mmix)
}
\keyword{Quantile}
\keyword{Weibull mixture}
