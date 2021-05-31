segselect <- function(J,Kmax,S=0.5){
Jtild = (Kmax-1)*((J[Kmax]-J)/(J[Kmax]-J[1]))+1
D2    = diff(diff(Jtild))
if (length(which(D2>=S))>0)  Kselect = max(which(D2>=S))+1
else Kselect = 1
invisible(Kselect)
}
