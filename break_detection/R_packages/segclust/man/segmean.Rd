 \name{segmean}
 \alias{segmean}
 \title{segmean}
 \description{
 	segmentation of a signal when considering changes in the mean
 }
 \usage{
 	out <- segmean(x,Kmax,lmin=1,lmax=length(x),vh=TRUE)
 }
 \details{ This function can be used for a segmentation model such that:
 \deqn{\forall t \in I_k \,\, Y_t = \mu_k + \varepsilon_t}{for t in I_k, Y_t = \mu_k + epsilon_t,} and the variance of \eqn{\epsilon_t} can be either
homoscedastic (vh=TRUE) or heteroscedastic (vh=FALSE). It uses dynamic programming to find the best breakpoints, and is based on the 
calculus of the Residual Sum of Squares: \deqn{J.est_K = \sum_{k=1}^{K} \sum_{t \in I_k} (y_t-\hat{\mu}_k)^2.}{J.est_K = \sum_{k=1}^{K} \sum_{t \in I_k} (y_t-\hat{\mu}_k)^2.}
 }
 \arguments{
 	\item{x}{data vector},
 	\item{Kmax}{maximum number of segments}
 	\item{lmin}{minimal segment length, default value lmin = 1}
	\item{lmax}{maximal segment length, default value lmax = length(x)}
 	\item{vh}{TRUE for an homogeneous variance (default), FALSE otherwise}
 }
 \value{
 \item{J.est}{Residual Sum of Squares for segmentation models up to Kmax segments}
 \item{t.est}{estimated positions for breakpoints for segmentation models up to Kmax segments}
 }
\examples{
	x1      <- rnorm(20,0,1)
	x2      <- rnorm(30,2,1)
	x       <- c(x1,x2)
	Kmax    <- 20
	out     <- segmean(x,Kmax)
	Kselect <- segselect(out$J.est,Kmax)
	output  <- segout(x,K=Kselect,th = out$t.est,draw=TRUE)

}
 \references{Picard, F., Robin, S., Lavielle, M., Vaisse, C., & Daudin, J. -J. (2005). 
 A statistical approach for array CGH data analysis. BMC Bioinformatics, 6(1), 1-14.}
 \keyword{ts}


 
