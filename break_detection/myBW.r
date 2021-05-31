## Code provided by Andrea Toreti

myBaumWelch<-function (object, control = bwcontrol(), ...) 
{
    x <- object$x
    Pi <- object$Pi
    delta <- object$delta
    distn <- object$distn
    pm <- object$pm
    tol <- control$tol
    if (distn[1] != "glm") {
        Mstep <- parse(text = paste("Mstep.", distn, "(x, cond, pm, object$pn)", 
            sep = ""))
    }
    else {
        Mstep <- parse(text = paste("Mstep.glm", "(x, cond, pm, object$pn, distn[2], distn[3])", 
            sep = ""))
    }
    m <- nrow(Pi)
    n <- length(x)
    oldLL <- -Inf
    condorel<-0

    for (iter in 1:control$maxiter) {
#browser()     
### code added by YB ###
    if(sum(pm$sd<0.001)>0) return(NULL)
########################
        cond <- Estep(x, Pi, delta, distn, pm, object$pn)
        diff <- cond$LL - oldLL
        if (control$prt) {
            cat("iter =", iter, "\n")
            cat("LL =", formatC(cond$LL, digits = log10(1/tol) + 
                2, format = "f"), "\n")
            cat("diff =", diff, "\n\n")
        }
        if (diff < 0 & control$posdiff) 
            stop("Worse log-likelihood on last iteration")
        if (eval(control$converge)) 
            break
        Pi <- diag(1/apply(cond$v, MARGIN = 2, FUN = sum)) %*% 
            apply(cond$v, MARGIN = c(2, 3), FUN = sum)
        if (object$nonstat) 
            delta <- cond$u[1, ]
        else delta <- compdelta(Pi)
        pm <- eval(Mstep)
        oldLL <- cond$LL
#controllo se ottengo una matrice di transizione con un elemento sulla diagonale 
#uguale a zero, in questo caso ho solo k-1 stati...quindi il mio modello a 
#k stati non è ragionevole.
        ladPi<-diag(Pi) 
       for (pidia in 1:length(ladPi)) {
### code added by YB ###
         if(is.na(ladPi[pidia])) return(NULL)
########################
         if(ladPi[pidia]<=0.5) {condorel<-1}
        }
        if(condorel==1) {break}
    }
    object$delta <- delta
    object$Pi <- Pi
    object$u <- cond$u
    object$v <- cond$v
    object$pm <- pm
    object$LL <- cond$LL
    object$iter <- iter
    object$diff <- diff
    return(object)
}

