\name{bootstrapCenWbMix}
\alias{bootstrapCenWbMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap Censored Weibull Mixture for censoring threshold selection and standard error of the quantile estimates
}
\description{
Similarly as in \code{\link{bootstrapCMLE}}, We will bootstrap the data set to obtain the MSE and SE of the quantile estimate under a certain subjective censoring threshold. Warning: The computation burden of this function is extremely heavy.
}
\usage{
bootstrapCenWbMix(dat, qInt = 0.05, canSet = c(0.5, 0.7, 1), B = 1000, 
			iniVec = NULL, randSeed = NULL, conCr = 1e-06, nIter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{The complect data set. See \code{\link{bootstrapCMLE}}}
  \item{qInt}{Quantile of interest. See \code{\link{bootstrapCMLE}}}
  \item{canSet}{Candidate set of thresholds. See \code{\link{bootstrapCMLE}}}
  \item{B}{Number of bootstrap replicates}
  \item{iniVec}{A vector of length 6*length(\code{canSet}), supplying the initial values of the EM algorithm under each censoring threshold for the original data set. Should be orginized in the following way (p, 1-p, a_1, a_2, b_1, b_2, p, 1-p, ...). If null, the initial values will be generated randomly, which is not recommended.}
  \item{randSeed}{See \code{\link{bootstrapCMLE}}}
  \item{conCr}{Same as in \code{\link{emCenWbMix.T1}}}
  \item{nIter}{Same as in \code{\link{emCenWbMix.T1}}}
}
\details{The EM algorithm for Weibull mixture is extremely senstive to the initial value, so the user should always try to find proper initial values for each threshold and data set. 

The initial values in the argument are for the original data set. The EM algorithm for the original data set will be first carried out from these initial values. If they succesfully converged, the EM algorithm for the bootstrap data sets will start from the estimates of the original data set in the previous step. Otherwise, the EM algorithm will start from the initial value in the arguements. Although we could not gurantee the "global" maximum can be reached for each bootstrap replicate under this scheme, our simulation shown that it works reasonably well.

Similarly as \code{\link{bootstrapCMLE}}, the main calculation of this function is done in C.
}
\value{
\item{results}{A matrix of length(\code{canSet}) by four. The first column is the candidate threshold (proportion). The second column is the quantile estimate under this censoring threshold. The third column and fourth column are the bootstrap estimate of the standard error (SE) and root mean squared error (RMSE) of this quantile estimate.}
\item{parameters}{A matrix of length(\code{canSet}) by six, recording the parameter estimates of the Weibull mixture under each threshold.}
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

\seealso{\code{\link{bootstrapCMLE}}}

\examples{
###
mmix = rbind(c(0.7, 5, 7),c(0.3, 15, 6))
vmix = c(mmix) #A vector version of the paramters
set.seed(1)
y <- sort(simWbMix(100, mmix)) #Generate Data
#Run the EM for 70\% Type II censoring
mix70 <- emCenWbMix.T2(y[1:70], 100, iniParam=vmix, useC=TRUE)
#Run the EM for the complete data set
mix100<- emCenWbMix.T2(y, 100, iniParam=vmix, useC=TRUE)
#Use bootstrap with the previous estimates as initials
bootstrapCenWbMix(y, canSet=c(0.7, 1), 
		iniVec=c(c(mix70$estimates), c(mix100$estimates)), B=5,randSeed=1)
#Or use the "true" value as the initials.
bootstrapCenWbMix(y, canSet=c(0.7, 1), iniVec=c(vmix, vmix), B=5,randSeed=1)
#B=5 is just for illustration. It should be at least 1000.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Bootstrap}
\keyword{Threshold selection}% __ONLY ONE__ keyword per line
\keyword{Weibull Mixture}
\keyword{Censor}
