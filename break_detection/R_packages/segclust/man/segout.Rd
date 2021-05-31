\name{segout}
\alias{segout}
\title{segout}
\description{
	Extraction of parameters for  segmentation model
}
\usage{
	out <- segclustout(x,K,th,draw,vh=TRUE)
}
\arguments{
	\item{x}{data vector (without missing values)},
	\item{K}{number of segments}
	\item{th}{estimated positions for breakpoints for segmentation models up to Kmax segments (t.est from segmean)}
	\item{draw}{TRUE for plotting}
	\item{vh}{Variance homogeneity, default = TRUE}

}
\value{
       \item{output}{dataframe containing results of the estimation procedure}
       \item{output$signal}{input signal x}
       \item{output$mean}{estimated mean according to the model, for each position}
       \item{output$sd}{estimated standard deviation according to the model, for each position}
       \item{output$bp}{breakpoint coordinates, equals 1 for a breakpoint (corresponding to the end of the segments)}
}
\examples{
	
	x1         <- rnorm(20,0,1)
	x2         <- rnorm(30,2,1)
	x3         <- rnorm(10,0,1)
	x4         <- rnorm(40,2,1)
	x          <- c(x1,x2,x3,x4)
	Kmax       <- 20
	out        <- segmean(x,Kmax)
	Kselect    <- segselect(out$J.est,Kmax)
	output     <- segout(x,K=Kselect,th = out$t.est,draw=TRUE)

}
\references{Picard, F., Robin, S., Lebarbier, E., & Daudin, J. -J. (2007). 
A segmentation/clustering model for the analysis of array CGH data. Biometrics, 63(3) 758-766}
\author{F. Picard, M. Hoebecke}
\keyword{cluster}


