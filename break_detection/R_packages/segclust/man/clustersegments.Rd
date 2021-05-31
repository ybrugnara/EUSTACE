\name{clustersegments}
\alias{clustersegments}
\title{clustersegments}
\description{
	Cluster segments using the EM algorithm when the number of clusters and the segmentation are given. The segmentation is not revised after clustering (contrary to hybrid)
}
\usage{
	out <- clustersegments(x,P,bp,vh=TRUE)
}
\arguments{
	\item{x}{data vector (without missing values)},
	\item{P}{number of clusters},
	\item{bp}{vector (size length(x)), such that bp[t]=1 if t is a breakpoint and 0 otherwise. (t corresponds to the end of a segment). bp[length(x)] is always equal to 1.}
	\item{vh}{TRUE for homogeneous variances (default), FALSE otherwise}
}
\value{
       \item{out}{dataframe}
       \item{output$signal}{input signal x}
       \item{output$mean}{estimated mean using a mixture model with P cluster, AFTER a segmentation}
       \item{output$sd}{estimated standard deviation using a mixture model with P cluster, AFTER a segmentation}
       \item{output$cluster}{cluster for each point AFTER segmentation}
       \item{output$bp}{breakpoint coordinates, equals 1 for a breakpoint (corresponding to the end of the segments)}
}
\examples{
	x1           <- rnorm(20,0,1)
	x2           <- rnorm(30,2,1)
	x            <- c(x1,x2)
	bp           <- rep(0,length(x))
	bp[c(20,50)] <- 1
	P            <- 2
	out          <- clustersegments(x,P,bp)
}
\references{Picard, F., Robin, S., Lebarbier, E., & Daudin, J. -J. (2007). 
A segmentation/clustering model for the analysis of array CGH data. Biometrics, 63(3) 758-766}
\author{F. Picard, M. Hoebecke}
\keyword{cluster}

