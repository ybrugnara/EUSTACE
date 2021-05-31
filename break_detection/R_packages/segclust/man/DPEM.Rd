\name{DPEM}
\alias{DPEM}
\title{Dynamic Programing / EM  segmentation/clustering}
\description{
	Estimation of segmentation/clustering parameters in the Gaussian case, using dynamic programming and the EM algorithm. This 
	function performs a global analysis with estimation and model selection. It uses functions hybrid(), segclustselect(), 
	and seclustout(). Pmin and Pmax must be different. 
}

\usage{
	output <- DPEM(x,Pmin,Pmax,Kmax,method,draw,lmin=1,lmax=length(x),vh=TRUE,S=0.5)
}
\arguments{
	\item{x}{data vector (without missing values)},
	\item{Pmin}{minimum number of clusters}
	\item{Pmax}{maximum number of clusters}
	\item{Kmax}{max number of segments. Kmax must be greater than Pmax}
	\item{method}{model selection method. Equals "sequential" or "BIC"}
	\item{draw}{equals TRUE for a graphical display}
	\item{lmin}{minimum segment length, default value lmin = 1}
	\item{lmax}{maximum segment length, default value lmax = length(x)}
	\item{vh}{TRUE for homogeneous variances (default), FALSE otherwise}
	\item{S}{Threshold for model selection, set at 0.5}
}
\value{
      \item{output}{dataframe containing results of the estimation procedure}
       \item{output$signal}{input signal x}
       \item{output$mean}{estimated mean according to the model, for each position}
       \item{output$sd}{estimated standard deviation according to the model, for each position}
       \item{output$cluster}{cluster for each point}
       \item{output$bp}{breakpoint coordinates, equals 1 for a breakpoint (corresponding to the end of the segments)}
}
\examples{
	x1         <- rnorm(20,0,1)
	x2         <- rnorm(30,2,1)
	x3         <- rnorm(10,0,1)
	x4         <- rnorm(40,2,1)
	x          <- c(x1,x2,x3,x4)
	Pmin       <- 1
	Pmax       <- 4
	Kmax       <- 20
	output	   <- DPEM(x,Pmin,Pmax,Kmax,method="BIC",draw=TRUE)
	
}
\references{Picard, F., Robin, S., Lebarbier, E., & Daudin, J. -J. (2007). 
A segmentation/clustering model for the analysis of array CGH data. Biometrics, 63(3) 758-766}
\author{F. Picard, M. Hoebecke}
\keyword{cluster}


