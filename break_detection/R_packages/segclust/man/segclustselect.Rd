\name{segclustselect}
\alias{segclustselect}
\title{segclustselect}
\description{
	Model selection for segmentation/clustering
}
\usage{
	out <- segclustselect(x,param,Pmin,Pmax,Kmax,Linc,method,S=0.5,lmin=1,lmax=length(x),vh=TRUE)
}
\details{ This function is used to select simulteaneously Pselect and Kselect, the number of clusters and the number of segments
in a segmentation/clustering model. It is based on the penalization of the incomplete-data log-likelhood Linc. Two 
methods are implemented. The first one is based on a sequential choice of Pselect and Kselect as described in Picard et al. (2007). The
second one is based on a modified BIC criterion as described in Zhang et al. (2007). The function uses the Stirling approximation of the Gamma function,
such that :
\deqn{\log \Gamma(x) \sim (x-1/2)\times \log (x) - x + 1/2 \times \log (2 \pi)}{log Gamma(x) ~ (x-1/2)* log(x) - x + 1/2 * log(2 \pi)} 
}
\arguments{
 	\item{x}{data vector},
 	\item{param}{parameters estimated by hybrid()}
	\item{Pmin}{minimum number of clusters}
	\item{Pmax}{maximum number of clusters}
 	\item{Kmax}{maximum number of segments}
 	\item{Linc}{incomplete-data log-likelihood calculated by hybrid()}
	\item{method}{Method used of the selection. Equals "sequential" or "BIC"}
 	\item{S}{threshold for the adaptive method, default value S = 0.5}
	\item{lmin}{minimal segment length, default value lmin = 1}
	\item{lmax}{maximal segment length, default value lmax = length(x)}
	\item{vh}{TRUE for homogeneous variances (default), FALSE otherwise}
	
}
\value{
\item{Pselect}{Selected number of clusters}
\item{Kselect}{Selected number of segments}
}
\examples{
	x1      <- rnorm(20,0,1)
	x2      <- rnorm(30,2,1)
	x3      <- rnorm(10,0,1)
	x4      <- rnorm(40,2,1)
	x       <- c(x1,x2,x3,x4)
	Pmin    <- 1
	Pmax    <- 4
	Kmax    <- 20
	Linc    <- matrix(-Inf, ncol=Pmax,nrow= Kmax)
	param.list   <- list()
	for (P in (Pmin:Pmax)){   
	    out.hybrid      <- hybrid(x,P,Kmax)
	    param.list[[P]] <- out.hybrid$param
	    Linc[,P]        <- out.hybrid$Linc    
  	}
	out.select <- segclustselect(x,param=param.list,Pmin,Pmax,Kmax,Linc,method = "BIC")
}
\references{
Picard, F., Robin, S., Lebarbier, E., & Daudin, J. -J. (2007). 
 A segmentation/clustering model for the analysis of array CGH data. Biometrics, 63(3) 758-766 \\
Zhang NR, Siegmund DO. A modified Bayes information criterion with applications to the analysis of comparative genomic hybridization data.
Biometrics. 2007 63(1):22-32.
}
\keyword{cluster}









