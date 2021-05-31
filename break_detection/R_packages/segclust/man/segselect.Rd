\name{segselect}
\alias{segselect}
\title{segselect}
\description{
	model selection for segmentation 
}
\usage{
	out <- segselect(J,Kmax,S=0.5)
}
\details{}
\arguments{
	\item{J}{Residual Sum of Squares for segmentation models up to Kmax segments}
 	\item{Kmax}{maximum number of segments}
 	\item{S}{threshold for the adaptive method, set to 0.5 (default)}
}
\value{
	\item{Kselect}{Selected number of segments}
}
\examples{
	x1      <- rnorm(20,0,1)
	x2      <- rnorm(30,2,1)
	x       <- c(x1,x2)
	Kmax    <- 20
	vh      <- TRUE
	out     <- segmean(x,Kmax)
	Kselect <- segselect(out$J.est,Kmax)
	output  <- segout(x,K=Kselect,th = out$t.est,draw=TRUE)
}
\references{Picard, F., Robin, S., Lavielle, M., Vaisse, C., & Daudin, J. -J. (2005). A statistical approach for array CGH data analysis. BMC Bioinformatics, 6(1), 1-14.}
\keyword{cluster}









