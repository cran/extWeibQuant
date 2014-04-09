\name{cenWbMLE.T2}
\alias{cenWbMLE.T2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
censored Weibull MLE for Type II right-censored data
}
\description{
Calculate the censored Weibull MLE for the Type II censored data with the algorithm described in ASTM 5457(2004). 
Return the estimates of the shape and scale parameters. A comprehesive description of this algorithm can be found in Liu (2012).
}
\usage{
cenWbMLE.T2(dat, n, useC = FALSE, conCr = 1e-09, nIter = 1000)
}
\arguments{
  \item{dat}{
A vector of the observations after censoring. None NA or negative values allowed. 
}
\item{n}{
The original sample size, including both the censored and uncensored observations.
}
  \item{useC}{
See \code{\link{cenWbMLE.T1}}
}
  \item{conCr}{
See \code{\link{cenWbMLE.T1}}
}
  \item{nIter}{
See \code{\link{cenWbMLE.T1}}
}
}
\value{
See \code{\link{cenWbMLE.T1}}
}
\references{
See \code{\link{cenWbMLE.T1}}
}
\author{
Yang (Seagle) Liu <yang.liu@stat.ubc.ca>
}
\note{
Please report the numerical problems and inconvenience when using this function to the author.

Please notice that the ways of inputing data in \code{cenWbMLE.T1} and \code{cenWbMLE.T2} are differnt. For \code{cenWbMLE.T1}, the algorithm require a full "orginal" data set (with the uncensored observations as NA or a arbitary value larger than the threshold) and the original sample size is decided as the length of the input data, while \code{cenWbMLE.T2} requires the observed data points and the original sample size.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rweibull}}, \code{\link{cenWbMLE.T1}}, \code{\link{emCenWbMix.T2}}
}
\examples{
set.seed(1)
y <- sort(rweibull(100, 7, 7)) ##Generate the data
cenWbMLE.T2(y, 100) #The MLE for the complete data
cenWbMLE.T2(y[1:10], 100) #Censor the largerst 90\% of the data.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{censor}
\keyword{Weibull}
\keyword{subjective censoring}
