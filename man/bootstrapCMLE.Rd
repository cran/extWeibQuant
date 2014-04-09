\name{bootstrapCMLE}
\alias{bootstrapCMLE}
\title{
Bootstrap Censored Weibull MLE for censoring threshold selection and standard error of the quantile estimates
}
\description{Proposed in Chapter 5 of Liu (2012). The original data set will be bootstrapped to obtain an estimate of the mean squared error of the quantile estimates, such that quantile estimates under different subjective censoring thresholds could be compared. We could also obtain the standard error of quantile estimates via bootstrap.}
\usage{
bootstrapCMLE(dat, qInt = 0.05, canSet = seq(0.1, 0.5, by = 0.1), 
			B = 5000, randSeed = NULL, conCr = 1e-09, nIter = 1000)
}
\arguments{
  \item{dat}{A complete data set without censoring. See \code{Details}}
  \item{qInt}{The quantile of interest, e.g. 5\% quantile or 10\% quantile.}
  \item{canSet}{A vector of the candidate subjective censoring thresholds, expressed as the proportion of data smaller it, e.g. 10\%, 20\%, ..., 100\% (non-censoring)}
  \item{B}{Number of bootstrap replicate data sets}
  \item{randSeed}{The seed for random number generation. If NULL, the random seed will not be set}
  \item{conCr}{Convergence criterion for the algorithm to calculate the censored Weibull MLE. See \code{\link{cenWbMLE.T1}}}
  \item{nIter}{See \code{\link{cenWbMLE.T1}}}
}
\details{
This function is designed to only work for a complete data set where every observation is fully observed (non-censoring). We could decide the best threshold (proportion) of subjective censoring based on bootstrap. For the advantage of subjective censoring, please see Liu (2012).

This function will call C to do all calculations. So it is recommended that the user should make sure the \code{\link{cenWbMLE.T2}} could work for their original data set.
}
\value{
\item{results}{A matrix of length(\code{canSet}) by six. The first column is the candidate threshold (proportion). The second and third column is the parameter estimates of the Weibull model for the original data set. 
The fourth column is the quantile estimate under this censoring threshold. The fifth column and sixth column are the bootstrap estimate of the standard error (SE) and root mean squared error (RMSE) of this quantile estimate.}
\item{bQEst}{A matrix of B-by-length(\code{canSet}). The quantile estimates under each censoring threshold for each bootstrap replicate.}
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

\seealso{\code{\link{bootstrapCenWbMix}}}

\examples{
set.seed(1)
y <- sort(rweibull(100, 7, 7))
tlist <- bootstrapCMLE(y, B=1000, canSet=c(0.1, 0.5, 1), randSeed=1)
tlist$results #Usually, we only need to look at the results part.
}
\keyword{Bootstrap}
\keyword{Threshold selection}% __ONLY ONE__ keyword per line
\keyword{Weibull}
\keyword{Censor}