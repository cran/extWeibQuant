\name{emCenWbMix.T1}
\alias{emCenWbMix}
\alias{emCenWbMix.T1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EM algorithm for the Type-I right censored Weibull Mixture with two populations. 
}
\description{
EM algorithm to estimate the parameters of a mixture of two two-parameter Weibull distributions for possibly Type-I right censored data. The PDF of this mixture distribution is 

\deqn{f(x) = p (a_1/b_1)(x/b_1)^(a_1-1)exp(-(x/b_1)^a_1) + (1-p)(a_2/b_2)(x/b_2)^(a_2-1)exp(-(x/b_2)^a_2),}

p is the proportion of the first sub-population. a_1, a_2 are the shape parameters for the two sub-populations. b_1, b_2 are the scale parameters. More details of the mixture model and this algorithm can be found in Liu (2012).
}
\usage{
emCenWbMix.T1(dat, Cx=NULL, iniParam = NULL, useC = FALSE, conCr = 1e-06, nIter = 10000)
}
\arguments{
  \item{dat}{
The data vector, should not contain any negative observations, but NA is allowed. 
The NA's and observations larger than the censoring threshold \code{Cx} will be censored in the calculation. The algorithm will decide the original sample size (before censoring) from the length of this vector.
}
  \item{Cx}{
The censoring threshold of Type I right-censoring. If NULL, the uncensored mixture will be estimated.
}
  \item{iniParam}{
A vector of length six (p, 1-p, a_1, a_2, b_1, b_2). All the elements must be postive. If not provided, the algorithm will generate them randomly. See \code{Details}.
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
\details{
It is well known that the EM algorithm is highly influenced by the initial value of the parameters. We strongly recommend the users provide some proper/reasonable initial values for their data sets. It is usually very difficult for EM algorithm to converge from the random initials generated in this function. Also, it is better to try starting this function from different initial values to obtain a "global" maximum of the likelihood. 
}
\value{
\item{convergence}{An integer indicating why the algorithm terminated
\itemize{
\item 0, successfully converged;
\item 1-4, Error in the maximization step. Same error code as in \code{\link{cenWbMLE.T1}};
\item 5, Numerically under/over flow, maybe the current parameters do not fit the data;
\item 6, EM Iteration limit reached.
}
}
\item{nllh}{The negative log-likelihood}
\item{estimates}{A 2 by 3 matrix of the parameter estimates:

p, a_1, b_1;

1-p, a_2, b_2;
}
\item{iniParam}{A copy of the iniParam}
}
\references{
Liu Y. (2012). Lower Quantile Estimation of Wood Strength Data. \emph{Master Thesis, Department of Statistics, UBC.} \href{https://circle.ubc.ca/bitstream/handle/2429/43067/ubc_2012_fall_liu_yang.pdf?sequence=10}{Downloadable here.}
}
\author{
Yang (Seagle) Liu <yang.liu@stat.ubc.ca>
}
\note{
Please report the numerical problems and inconvenience when using this function to the author.
}

\seealso{
\code{\link{rweibull}}, \code{\link{quanWbMix}}, \code{\link{simWbMix}}, \code{\link{cenWbMLE.T1}}, \code{\link{emCenWbMix.T2}}
}
\examples{
mmix = rbind(c(0.7, 5, 7),c(0.3, 15, 6))
vmix = c(mmix) #A vector version of the paramters
set.seed(1)
y <- sort(simWbMix(300, mmix)) #Generate the data
#The uncensored mixture
emCenWbMix.T1(y, iniParam=vmix, useC=TRUE)
#The mixture if we censor the data around 9.
emCenWbMix.T1(y, Cx=9, iniParam=vmix, useC=TRUE) 
}
\keyword{Weibull}
\keyword{mixture}
\keyword{censor}
