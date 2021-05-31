 \name{segmixt}
 \alias{segmixt}
 \title{segmixt}
 \description{
 	segmentation of a signal when considering changes in a mixture of Gaussian vectors
 }

 \usage{
 	out <- segmixt(x,P,Kmax,phi,lmin=1,lmax=length(x))
 }
 \details{ 

 This function can be used for a segmentation/clustering model with P clusters and up to Kmax segments. phi gives 
 the parameters of the mixture (P means, P standard deviations, P mixture proportions). It uses dynamic programming to find 
 the best breakpoints, and is based on the calculus of the incomplete-data log-likelihood:
 \deqn{J.est_{K,P} = \sum_{k=1}^{K} \log \left\{ \sum_{p=1}^{P} \pi_p f(y^k;\theta_p)\right\}}{J.est_{K,P} = \sum_{k=1}^{K} log [ \sum_{p=1}^{P} (pi_p f(y^k;theta_p)) ] }
 where \eqn{f(y^k;\theta_p)}{f(y^k;theta_p)} is the density of a Gaussian vector of length \eqn{n_k}{n_k}.
 }
 \arguments{
 	\item{x}{data vector},
	\item{P}{number of clusters}
	\item{Kmax}{maximum number of segments}
	\item{phi}{parameters of the mixture}
 	\item{lmin}{minimal segment length, default value lmin = 1}
	\item{lmax}{maximal segment length, default value lmax = length(x)} 		
 }
 \value{ 
 	\item{J.est}{incomplete-data log-likelihood for a segmentation/clustering model with P clusters and up to Kmax segments}
 	\item{t.est}{estimated positions of breakpoints for a segmentation/clustering model up to Kmax segments}
 }
 \references{Picard, F., Robin, S., Lebarbier, E., & Daudin, J. -J. (2007). 
 A segmentation/clustering model for the analysis of array CGH data. Biometrics, 63(3) 758-766}
 \author{F. Picard, M. Hoebecke}
 \keyword{ts}
