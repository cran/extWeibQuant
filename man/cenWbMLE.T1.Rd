\name{cenWbMLE.T1}
\alias{cenWbMLE}
\alias{cenWbMLE.T1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
censored Weibull MLE for Type I right-censored data 
}
\description{
Calculate the censored Weibull (two parameter, shape and scale) MLE for the Type I censored data with the algorithm described in ASTM 5457(2004).
Return the estimates of the shape and scale parameters. A comprehesive description of this algorithm can be found in Liu (2012).
}
\usage{
cenWbMLE.T1(dat, Cx=NULL, useC = FALSE, conCr = 1e-09, nIter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dat}{
A vector of the data, should not contain any negative observations, but NA is allowed. 
The NA's and observations larger than the censoring threshold \code{Cx} will be censored in the calculation. The algorithm will decide the original sample size (before censoring) from the length of this vector.
}
  \item{Cx}{
The censoring threshold of Type I right-censoring. If NULL, the complete (uncensored MLE) will be calculated. 
}
  \item{useC}{
Default to be false and use the R routine to calculate estimates. 
If true, the function will use the C routine, which is much faster than the R routine, 
but harder for the user to identify the numerical issues (if there is any).
}
  \item{conCr}{
In terms of the relative change in the negative log-likelihood. 
The algorithm is viewed as converged if the relative change is smaller than \code{conCr}
}
  \item{nIter}{
The maximum numer of iterations allowed in the function.  
}
}
\value{
\item{convergence}{an integer indicating why the algorithm terminated
\itemize{
\item 0, successfully converged;
\item 1, Numerical over/under-flow, maybe the shape parameter diverges to infinity;
\item 2, Iteration limit reached;
\item 3, Shape parameter becomes negative;
\item 4, Shape parameter diverges to infinity.
}
}
\item{estimates}{Shape and scale parameter estimates}
}
\references{
ASTM (2004). Standard specfication for computing reference resistance of wood-based materials and structural connections for load and resistance factor design D5457. \emph{American Society for Testing Materials,}
Philadephia, Pa.

Liu Y. (2012). Lower Quantile Estimation of Wood Strength Data. \emph{Master Thesis, Department of Statistics, UBC.} \href{https://circle.ubc.ca/bitstream/handle/2429/43067/ubc_2012_fall_liu_yang.pdf?sequence=10}{Downloadable here.}
}
\author{
Yang (Seagle) Liu <yang.liu@stat.ubc.ca>
}
\note{
Please report the numerical problems and inconvenience when using this function to the author.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rweibull}}, \code{\link{cenWbMLE.T2}}, \code{\link{emCenWbMix.T1}}
}
\examples{
set.seed(1)
y <- sort(rweibull(100, 7, 7)) ##Generate the data
cenWbMLE.T1(y) #The MLE for the complete data
cenWbMLE.T1(y, 5) #Censor the data at 5 and calculate the censored MLE
##Or
newy <- rep(NA, 100)
newy[y<=5] <- y[y<=5] #Censor the data at 5
fit <- cenWbMLE.T1(newy, 5) #Calculate the censored MLE
qweibull(0.05, fit$estimates[1], fit$estimates[2]) 
#Calculate the 5\% quantile of the fitted distribution.
}
\keyword{censor}
\keyword{Weibull}
\keyword{subjective censoring}