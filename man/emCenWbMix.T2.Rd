\name{emCenWbMix.T2}
\alias{emCenWbMix.T2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EM algorithm for the Type-II right censored Weibull Mixture with two populations.
}
\description{
This is the interface for the EM algorithm of censored Weibull mixture for Type II censored data. The description of the Weibull mixture model can be found in \code{\link{emCenWbMix.T1}}.
}
\usage{
emCenWbMix.T2(dat, n, iniParam = NULL, useC = FALSE, conCr = 1e-06, nIter = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A vector of the observations after censoring. None NA or negative values allowed.
}
  \item{n}{
The original sample size, including both the censored and uncensored observations. See \code{\link{cenWbMLE.T2}}. 
}
  \item{iniParam}{See \code{\link{emCenWbMix.T1}}}
  \item{useC}{See \code{\link{emCenWbMix.T1}}}
  \item{conCr}{See \code{\link{emCenWbMix.T1}}}
  \item{nIter}{See \code{\link{emCenWbMix.T1}}}
}
\value{See \code{\link{emCenWbMix.T1}}}
\references{See \code{\link{emCenWbMix.T1}}}


\seealso{
\code{\link{rweibull}}, \code{\link{quanWbMix}}, \code{\link{simWbMix}}, \code{\link{cenWbMLE.T2}}, \code{\link{emCenWbMix.T1}}
}
\examples{
mmix = rbind(c(0.7, 5, 7),c(0.3, 15, 6))
vmix = c(mmix) #A vector version of the paramters
set.seed(1)
y <- sort(simWbMix(300, mmix)) #Generate the data
#The uncensored EM
emCenWbMix.T2(y, 300, iniParam=vmix, useC=TRUE)
#Type II censoring with the largest half of the data censored.
emCenWbMix.T2(y[1:150], 300, iniParam=vmix, useC=TRUE)
}
\keyword{Weibull}
\keyword{mixture}
\keyword{censor}