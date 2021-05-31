\name{hybrid}
\alias{hybrid}
\title{hybrid algorithm for segmentation/clustering}
\description{
	Estimation of segmentation/clustering parameters in the Gaussian case, using dynamic programming and the EM algorithm. This function estimates the 
	parameters for a model with P clusters and for a number of segments from P to Kmax (Kmin=P).
}

\usage{
	out.hybrid <- hybrid(x,P,Kmax,lmin=1,lmax=length(x),vh=TRUE)
}
\arguments{
	\item{x}{data vector (without missing values)},
	\item{P}{number of clusters}
	\item{Kmax}{max number of segments. The minimum number of segments is P}
	\item{lmin}{minimum segment length, default value lmin = 1}
	\item{lmax}{maximum segment length, default value lmax = length(x)}
	\item{vh}{TRUE for homogeneous variances (default), FALSE otherwise}
}
\value{
       \item{Linc}{incomplete-data log-likelihood for a model with P clusters and P to Kmax segments}
       \item{param[[K]]}{list of estimated parameters for a segmentation/clustering model with P clusters and K segments}
       \item{param[[K]]$phi}{mixture parameters of size 3*P : P means, P standard deviations, P mixture proportions}
       \item{param[[K]]$rupt}{breakpoints position : matrix with 2 columns (begin and end of segments) and K rows}
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
	out.select <- segclustselect(x,param=param.list,Pmin,Pmax,Kmax,Linc,method = "BIC")
	output     <- segclustout(x,param.list[[out.select$Pselect]],out.select$Pselect,out.select$Kselect,draw=TRUE)
}
\references{Picard, F., Robin, S., Lebarbier, E., & Daudin, J. -J. (2007). 
A segmentation/clustering model for the analysis of array CGH data. Biometrics, 63(3) 758-766}
\author{F. Picard, M. Hoebecke}
\keyword{cluster}


