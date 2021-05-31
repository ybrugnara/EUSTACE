\name{segclustout}
\alias{segclustout}
\title{segclustout}
\description{
	Extraction of parameters for  a segmentation/clustering model
}
\usage{
	out <- segclustout(x,param,P,K,draw)
}
\arguments{
	\item{x}{data vector (without missing values)},
	\item{param}{list of parameters estimated by hybrid for a given P},
	\item{P}{number of clusters}
	\item{K}{number of segments (must be smaller than P)}
	\item{draw}{TRUE for plotting}
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
	Linc       <- matrix(-Inf, ncol=Pmax,nrow= Kmax)
	param.list <- list()

	for (P in (Pmin:Pmax)){   
	    out.hybrid      <- hybrid(x,P,Kmax)
	    param.list[[P]] <- out.hybrid$param
	    Linc[,P]        <- out.hybrid$Linc    
  	}
	out.select <- segclustselect(x,param,Pmin,Pmax,Kmax,Linc, method = "sequential")
	output     <- segclustout(x,param.list[[out.select$Pselect]],out.select$Pselect,out.select$Kselect,draw=TRUE)

}
\references{Picard, F., Robin, S., Lebarbier, E., & Daudin, J. -J. (2007). 
A segmentation/clustering model for the analysis of array CGH data. Biometrics, 63(3) 758-766}
\author{F. Picard, M. Hoebecke}
\keyword{cluster}


