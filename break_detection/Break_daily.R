###########################################################
###########################################################
###                                                     ###
###  EUSTACE Break-detection software                   ###
###  Version 1.1 for daily data                         ###
###  04-February-2021                                   ###
###  Author: Yuri Brugnara, University of Bern          ###
###  yuri.brugnara@giub.unibe.ch                        ###
###                                                     ###
###########################################################
###########################################################
###
##########################
###  FIRST READ THIS:  ###
##########################
###################################################################################################
###
# This software is designed to run on a Unix system.
###
# Following R packages must be installed:
# -----------
# segclust, HiddenMarkov, vegan, genalg, maps, permute, lattice
# NB: This software was tested using R version 3.1.2, it might not work with previous versions.
###
# External R scripts called by this script:
# -------------------
# penlik.r (Caussinus & Mestre, 2004);
# RHtestV5.r (Wang et al., 2007);
# gahmm.r, myBW.r, myeval_fortranized.r, myevalh.r, myrgba.r, myvit.r (Toreti et al., 2012)
###
# External Fortran scripts:
# -------------------------
# myeval.f, read_daily.f (same version as this program)
# NB: The Fortran scripts must be compiled using the command 'R CMD SHLIB scriptname.f'.
###
# Conventions for input data files:
# --------------------------
# All input data files must be in the same directory, one file can contain only one variable;
# filenames must begin with 'n' for Tmin, with 'x' for Tmax, with 'a' for Tavg ([min+max]/2), 
#   with 'd' for DTR (max-min), with 't' for Tmean;
# input files must be ASCII files;
# input files format must be strictly as follows (including spaces):
# - 1st line: '# Anytext' (usually the source)
# - 2nd line: '# coordinates: xxx.xxxxxxN, xxxx.xxxxxxE,   xxxx.xm'
# - 3rd line: '# station of: stationname' (station name cannot be more than 50 characters, 
#                                          spaces are allowed)
# - 4th line: '# Anytext' (usually the name of the variable)
# - 5th line: '# Anytext' (usually info on QC)
# - from 6th line: ' yyyy mm dd    xxx.xx' (Fortran format: i5,2i3,f10.2)
# missing values must be omitted or indicated by the value -99.99
###
###################################################################################################
#
#
#
########################################
###  THEN SET THE INPUT PARAMETERS:  ###
########################################
rm(list=objects()) ##clean R workspace
###################################################################################################

### Full path of the directory containing the input files, omit final '/' 
inpath <- "/home/yb/Desktop/MO/validation/Wyoming_world_4/formatted"

### Full path of the directory where you want to save the output, omit final '/'
outpath <- "/home/yb/Desktop/MO/validation/Wyoming_world_4/formatted/bd"

### First letter of the input filenames (see above) (NB: can also be a longer prefix, in order to
###                                                      read only a certain group of files)
fl <- "WY"

### [Optional] Full path of a ASCII list of input files (comment if not used)
#inlist <- "/home/users/ybrugnara/scripts/list_bd_5.dtr"

### [Optional] Filename of the candidate series (to be used to analyse only one specific series, 
###                                              comment otherwise)
#candidate <- "Ta_bath.txt"

### Maximum number of missing days to use a certain month
daysmax <- 5

### Minimum number of available months to use a certain series
ymin <- 10*12

### Maximum number of reference series to use (NB: directly proportional to computation time and
###                                                false positives!)
nrefmax <- 8

### Period to analyse
yrange <- c(1970,2011)

### Minimum correlation of the reference series
cormin <- 0.6

### Maximum distance in km of the reference series from the candidate (search starts using
### 'distmax'/10, then progressively increases the radius up to distmax if necessary)
distmax <- 1000

## Minimum fraction of overlapping non-missing months
overlapmin <- 0.8

### Additional semi-annual break detection (TRUE/FALSE)
semiy <- TRUE

### Graphical output (TRUE/FALSE)
graphout <- FALSE

### Recovery mode (TRUE/FALSE) - True to restart from the point where a run crashed 
### (e.g. because of segmentation fault)
recmode <- FALSE 

## Do monthly t-test?
monthly <- TRUE

## Window size for monthly t-test (years, must be an odd number)
w <- 5

## Likelihood threshold (the lower this number is, the more false alarms will result)
thr <- 2

###################################################################################################
#
#
#
############################
###  OUTPUT DESCRIPTION  ###
############################
###################################################################################################
#
# The software produces the following outputs:
# - info on reference stations available for each candidate series (standard output);
# - plot of standardized differences with breaks and map of the stations (one pdf file for each
#   series in outpath - only if relative test is performed);
# - R workspace image (saved in outpath/<fl>_breaks.RData) containing:
#   - list 'breaks' --> detected breakpoints for each candidate series
#   - list 'bdscore' --> detection score for each year and each season for each candidate
#                        series (i.e., sum of first-difference correlation coefficients of the
#                                reference series with at least 6 years of common data in a 
#                                11-year window centred on the target year)
#
###################################################################################################



#####################################
### START (do not edit from here) ###
#####################################


## Load R packages
require(segclust)
require(HiddenMarkov)
require(vegan)
require(genalg)
require(maps)

## Load external R scripts
source("penlik.r") #caume
assignInNamespace("segmean",segmean,ns="segclust")
source("RHtestsV5.r") #wang
source("gahmm.r") #gahmdi

## Load external Fortran scripts
dyn.load("myeval.so") #speeds up gahmdi

## Initialize variables
v <- "1.1"
workpath <- getwd()
if (exists("inlist")) system(paste0("cp ",inlist," ",outpath,"/list"))
if (!exists("inlist")) {
  setwd(inpath)
  system(paste0("ls ",fl,"* > ",outpath,"/list"))
}
filenames <- t(read.table(paste0(outpath,"/list")))
setwd(workpath)
if (length(unique(substr(c(filenames),1,1))) > 1) {
  stop("The first letter of the filenames should be identical among all input files")
}
if (exists("candidate")) {
  ic <- which(filenames==candidate)
  filenames <- candidate
}
seasons <- "y"

# Set up recovery mode
if (recmode) {
  load(paste0(outpath,"/",fl,"_breaks.RData"))
  ifile1 <- 1+which(filenames==names(breaks)[length(names(breaks))])
} else {
  breaks <- list()
  bdscore <- list()
  ifile1 <- 1
}


## Loop over candidate stations
for (candiname in filenames[ifile1:length(filenames)]) {

  ## Skip empty files
  if (substr(try(read.table(paste0(inpath,"/",candiname),nrows=1),silent=TRUE),1,5)[1]
      != "Error") {
    
    ### Read daily series and transform to monthly by calling external Fortran script
    ### (27 is the maximum number of data series that the foreing function interface will pass to R,
    ### using a larger limit might cause an error)
    if (!exists("candidate")) ic <- which(filenames==candiname)
    dyn.load("read_daily.so")
    mydat <- .Fortran("rddaily", 
                      files = as.character("list"), 
                      nf = as.integer(27), 
                      x = array(data=as.numeric(-99.99),dim=(yrange[2]-yrange[1]+1)*27*12),
                      rx = as.numeric(rep(0,27-1)),
                      lat = as.numeric(rep(0,27)),
                      lon = as.numeric(rep(0,27)),
                      ele = as.numeric(rep(0,27)),
                      yr1 = as.integer(yrange[1]),
                      yr2 = as.integer(yrange[2]),
                      ddmax = as.integer(daysmax),
                      dismax = as.numeric(distmax),
                      nmin = as.integer(ymin),
                      rmin = as.numeric(cormin),
                      nmax = as.integer(nrefmax),
                      overlap = as.numeric(overlapmin),
                      ic = as.integer(ic),
                      sname = character(1),
                      j = as.integer(1),
                      inpath = as.character(paste0(inpath,"//")),
                      outpath = as.character(paste0(outpath,"//"))) 
    dyn.unload("read_daily.so")
    nref <- mydat$j
    h <- nref
    cname <- rev(strsplit(candiname,"/")[[1]])[1]
    years <- mydat$yr1:mydat$yr2
    nyears <- length(years)
    
    ## Set up list 'breaks'
    breaks[[cname]] <- list()
    breaks[[cname]]$range <- c(min(years),max(years))
    if (mydat$ele[1] < (-418)) mydat$ele[1] <- NA
    breaks[[cname]]$coord <- c(mydat$lon[1],mydat$lat[1],mydat$ele[1])
    breaks[[cname]]$bd <- FALSE
    breaks[[cname]]$abs <- FALSE
    bdscore[[cname]] <- list()

    ## Format station and file names
    if (nref >= 3) {
      sname <- gsub("  ","",substr(mydat$sname,1,50))
      sname <- gsub("_"," ",sname)
      refnames <- t(read.table(paste0(outpath,"/list_ref")))
      for (i in 1:nref) refnames[i] <- rev(strsplit(refnames[i],"/")[[1]])[1]
    } else refnames <- c()
    
    ## Organize monthly data into a list of matrices
    mondata <- list()
    mondata[[candiname]] <- t(matrix(mydat$x[((years[1]-yrange[1])*12+1):
                                               ((years[1]-yrange[1]+nyears)*12)],
                                     nrow=12,
                                     ncol=nyears,
                                     dimnames=list(1:12,years)))
    mondata[[candiname]][mondata[[candiname]]<=(-90)] <- NA
    
    ## Calculate climatology of candidate
    clima <- array(dim=12)
    for (im in 1:12) {
      if (sum(!is.na(mondata[[1]][,im])) >= (ymin/12)) {
        clima[im] <- mean(mondata[[1]][,im],na.rm=TRUE)
      }
    }

    ## Check if there are enough data in the candidate series
    if (sum(!is.na(mondata[[candiname]])) >= ymin) {

    ## Populate list of monthly data with reference series
      if (nref >= 3) {
        for (i in 1:nref) {
          mondata[[refnames[i]]] <- t(matrix(mydat$x[((yrange[2]-yrange[1]+1)*12*i+
                                                        (years[1]-yrange[1])*12+1):
                                                       ((yrange[2]-yrange[1]+1)*12*i+
                                                          (years[1]-yrange[1]+nyears)*12)],
                                             nrow=12,
                                             ncol=nyears,
                                             dimnames=list(1:12,years)))
          mondata[[refnames[i]]][mondata[[refnames[i]]]<=(-90)] <- NA
        }
      }
      mydat$x <- NULL
      
      ## Define lists for difference series and correlations
      z <- list()
      rref <- list()
      
      
      ## Loop over seasons
      if (semiy) seasons <- c("s","w","y")
      for (season in seasons) {
        
        ################################
        ### RELATIVE BREAK DETECTION ###
        ################################
        
        ## Calculate annual/seasonal means and write them in the matrix 'mytx'
        mytx <- matrix(nrow=nyears,ncol=nref+1)
        if (!is.null(refnames)) dimnames(mytx) <- list(years,c(candiname,refnames))
        if (season=="y") mytx[,1] <- rowMeans(mondata[[candiname]])
        if (season=="s") mytx[,1] <- rowMeans(mondata[[candiname]][,4:9])
        if (season=="w") mytx[,1] <- rowMeans(cbind(mondata[[candiname]][,1:3],
                                                    rbind(rep(NA,3),mondata[[candiname]]
                                                          [1:(dim(mytx)[1]-1),10:12])))
        if (nref >= 3) {
          for (i in 1:nref) {
            if (season == "y") mytx[,i+1] <- rowMeans(mondata[[refnames[i]]])
            if (season == "s") mytx[,i+1] <- rowMeans(mondata[[refnames[i]]][,4:9])
            if (season == "w") mytx[,i+1] <- rowMeans(cbind(mondata[[refnames[i]]][,1:3],
                                                            rbind(rep(NA,3),mondata[[refnames[i]]]
                                                                  [1:(dim(mytx)[1]-1),10:12])))
          }
        }
        mytx <- round(mytx,2)
        
        ## Read correlations and coordinates
        rref[[season]] <- mydat$rx[mydat$rx>=cormin]
        lat <- mydat$lat[1:(nref+1)]
        lon <- mydat$lon[1:(nref+1)]
        ele <- mydat$ele[1:(nref+1)]
               
        ## Define matrix 'candi' containing the years and the candidate series
        if (is.null(dim(mytx))) candi <- cbind(years,mytx)
        else candi <- cbind(years,mytx[,1])

        # Check if there are enough reference series, if not skip relative break detection
        if (h >= 3) { 
          breaks[[cname]]$bd <- TRUE
          
          # Define array that will contain the maximum number of detectable segments for each
          # difference series based on their lengths (does not apply to RHTests)
          mykmax <- array(dim=h)
          
          ## Remove reference series with less than 'ymin'/12 common years
          ## Vector with the number of common years for each series: 'ncyrs'
          names(rref[[season]]) <- refnames
          h <- dim(mytx)[2]-1
          ko <- c()
          ncyrs <- array(dim=h)
          names(ncyrs) <- refnames
          for (i in 2:dim(mytx)[2]) {
            if (sum(!is.na(mytx[,1]-mytx[,i])) < ymin/12) {
              ko <- append(ko,i)
            } else {
              ncyrs[i-1] <- sum(!is.na(mytx[,1]-mytx[,i]))
            }
          }
          if (!is.null(ko)) {
            h <- h-length(ko)
            mytx <- mytx[,-ko]
            rref[[season]] <- rref[[season]][-(ko-1)]
            ncyrs <- ncyrs[which(!is.na(ncyrs))]
            lat <- lat[-ko]
            lon <- lon[-ko]
            ele <- ele[-ko]
          }
          
          ## Order reference series from best to worst based on overlapping and correlation
          ## with candidate
          neworder <- order(ncyrs,rref[[season]],decreasing=TRUE,na.last=NA)
          rref[[season]] <- rref[[season]][neworder]
          ncyrs <- ncyrs[neworder]
          lat <- lat[c(1,neworder+1)]
          lon <- lon[c(1,neworder+1)]
          ele <- ele[c(1,neworder+1)]
          
          ## Limit the number of reference series to nrefmax
          if (h > nrefmax) {
            h <- nrefmax
            rref[[season]] <- rref[[season]][1:h]
            ncyrs <- ncyrs[1:h]
          }
          
          ## Define matrix 'ref' containing the reference series
          ref <- matrix(nrow=nyears,ncol=h)
          for (i in 1:h) ref[,i] <- mytx[,which(colnames(mytx)==names(rref[[season]])[i])]
          colnames(ref) <- names(rref[[season]])[1:h]
          rownames(ref) <- years
          
          ## Calculate standardized difference series
          z[[season]] <- matrix(nrow=nyears,ncol=(h+1))
          colnames(z[[season]]) <- c("Year",colnames(ref))
          z[[season]][,1] <- years
          for (i in 2:(h+1)) {
            z[[season]][,i] <- candi[,2]-ref[,(i-1)]
            z[[season]][,i] <- z[[season]][,i]-mean(z[[season]][,i],na.rm=TRUE)
            z[[season]][,i] <- z[[season]][,i]/sd(z[[season]][,i],na.rm=TRUE)
            mykmax[i-1] <- round(sum(!is.na(z[[season]][,i]))/10,0)+1
          }
          
          ## Calculate detection score for each year, based on 11-year windows
          bdscore[[cname]][[season]] <- data.frame(years,rep(NA,length(years)))
          names(bdscore[[cname]][[season]]) <- c("Year","Score")
          if (!is.null(dim(z[[season]])[2])) {
            for (i in 6:(nyears-5)) {
              bdscore[[cname]][[season]]$Score[i] <- 0
              for (j in 1:(dim(z[[season]])[2]-1)) {
                ## At least 5 non-missing values are required within the 11 years
                if (sum(!is.na(z[[season]][(i-5):(i+5),j+1])) > 5) {
                  bdscore[[cname]][[season]]$Score[i] <- bdscore[[cname]][[season]]$Score[i]+
                                                         rref[[season]][j]
                }
              }
            }
            bdscore[[cname]][[season]]$Score[1:5] <- bdscore[[cname]][[season]]$Score[6]
            bdscore[[cname]][[season]]$Score[(nyears-4):nyears] <-
              bdscore[[cname]][[season]]$Score[nyears-5]
          }
          
          ## Caussinus and Mestre test
          ## Breakpoints are written into matrix 'inho_cm'
          inho_cm <- matrix(nrow=max(mykmax,na.rm=TRUE),ncol=h)
          colnames(inho_cm) <- colnames(ref)
          for (i in 1:h) {
            y <- z[[season]][,c(1,i+1)]
            y <- y[apply(y,1,function(x)all(!is.na(x))),]
            iso <- pen.lik(y[,2],kmax=mykmax[i])
            if (iso[1] != 0) {
              for (j in 1:length(iso)) {
                inho_cm[j,i] <- y[iso[j],1]
              }
            }
          }
          
          ## Wang and Feng test (RHtests)
          ## For this test a few temporary files are written on disk
          ## Breakpoints are written into matrix 'inho_wf'
          inho_wf <- matrix(nrow=30,ncol=h)
          colnames(inho_wf) <- colnames(ref)
          write.table(cbind(candi[,1],0,0,candi[,2]),
                      file=paste0(outpath,"/tmp"),
                      col.names=FALSE,
                      row.names=FALSE,
                      quote=FALSE)
          for (i in 1:h) {
            write.table(cbind(candi[,1],0,0,ref[,i]),
                        file=paste0(outpath,"/tmp_ref"),
                        col.names=FALSE,
                        row.names=FALSE,
                        quote=FALSE)
            ## The variable 'status' catches possible errors that will then trigger warnings
            status <- FindU.wRef(Bseries=paste0(outpath,"/tmp"),
                                 MissingValueCode="NA",
                                 p.lev=0.95,
                                 Iadj=10000,
                                 Mq=10,
                                 Ny4a=0,
                                 output=paste0(outpath,"/tmp"),
                                 Rseries=paste0(outpath,"/tmp_ref"))
            if (!is.null(status)) {
              if (length(list.files(outpath,pattern="tmp_mCs.txt")) > 0) {
                status <- StepSize.wRef(Bseries=paste0(outpath,"/tmp"),
                                        MissingValueCode="NA",
                                        p.lev=0.95,
                                        Iadj=10000,
                                        Mq=10,
                                        Ny4a=0,
                                        InCs=paste0(outpath,"/tmp_mCs.txt"),
                                        output=paste0(outpath,"/tmp"),
                                        Rseries=paste0(outpath,"/tmp_ref"))
                if(!is.null(status)) {
                  rawfile <- read.table(paste0(outpath,"/tmp_mCs.txt"),fill=TRUE)
                  if (dim(rawfile)[1] > 1) {
                    rawfile <- read.table(paste0(outpath,"/tmp_mCs.txt"),skip=1)
                    j <- which(rawfile[,1]!=1 | rawfile[,2]!="Yes")
                    nb <- dim(rawfile)[1]-length(j)
                    while(length(j) > 0) {
                      rawfile <- read.table(paste0(outpath,"/tmp_mCs.txt"),sep="\t")
                      write.table(rawfile[-(j+1),],
                                  file=paste0(outpath,"/tmp_mCs.txt"),
                                  quote=FALSE,
                                  row.names=FALSE,
                                  col.names=FALSE,
                                  append=FALSE)
                      if(length(rawfile[-(j+1),])[1] > 1) {
                        status <- StepSize.wRef(Bseries=paste0(outpath,"/tmp"),
                                                MissingValueCode="NA",
                                                p.lev=0.95,
                                                Iadj=10000,
                                                Mq=10,
                                                Ny4a=0,
                                                InCs=paste0(outpath,"/tmp_mCs.txt"),
                                                output=paste0(outpath,"/tmp"),
                                                Rseries=paste0(outpath,"/tmp_ref"))
                        if (!is.null(status)) {
                          rawfile <- read.table(paste0(outpath,"/tmp_mCs.txt"),skip=1)
                          j <- which(rawfile[,1]!=1 | rawfile[,2]!="Yes")
                          nb <- dim(rawfile)[1]-length(j)
                        } else {
                          warning(paste0(candiname,
                                         ": impossible to apply WANG with reference series ",
                                         colnames(ref)[i]))
                          j <- c()
                          nb <- 0
                        }
                      } else {
                        j <- c()
                        nb <- 0
                      }
                    }
                    if (nb > 0) inho_wf[1:nb,i] <- as.numeric(substr(rawfile[,3],1,4))
                  }
                } else warning(paste0(candiname,
                                      ": impossible to apply WANG with reference series ",
                                      colnames(ref)[i]))
              }
              system(paste0("rm ",outpath,"/tmp_*"))
            } else warning(paste0(candiname,
                                  ": impossible to apply WANG with reference series ",
                                  colnames(ref)[i]))
          }
          system(paste0("rm ",outpath,"/tmp"))
          
          ## GAHMDI test
          ## Breakpoints are written into matrix 'inho_ga'
          inho_ga <- matrix(nrow=max(mykmax,na.rm=TRUE),ncol=h)
          colnames(inho_ga) <- colnames(ref)
          for (i in 1:h) {
            y <- z[[season]][,c(1,i+1)]
            y <- y[apply(y,1,function(x)all(!is.na(x))),]
            ## The variable 'isos' catches possible errors that will then trigger warnings
            isos <- gahmm(y[,2],kmax=mykmax[i])
            if(!is.null(isos)) {
              iso <- isos$ICLChan
              if (iso[1] != 0) {
                for (j in 1:length(iso)) {
                  inho_ga[j,i] <- y[iso[j],1]
                }
              }
            } else warning(paste0(candiname,
                                  ": impossible to apply GAHMDI with reference series ",
                                  colnames(ref)[i]))
          }
          
        }
        
        
        ################################
        ### ABSOLUTE BREAK DETECTION ###
        ################################
        
        ## Absolute test on monthly anomalies (RHtests)
        ## For this test a few temporary files are written on disk
        ## Breakpoints are written into vector 'inho_ab'
        inho_ab <- c()
        yseq <- c()
        seasclima <- clima
        if (season == "s") seasclima[c(1:3,10:12)] <- NA
        if (season == "w") seasclima[4:9] <- NA
        if (sum(!is.na(seasclima)) > 0) {
          breaks[[cname]]$abs <- TRUE
          anoms <- c(t(mondata[[1]]))-rep(seasclima,dim(mondata[[1]])[1])
          for (iy in years) yseq <- append(yseq,rep(iy,12))
          write.table(cbind(yseq,rep(1:12,length(years)),0,anoms),
                      file=paste0(outpath,"/tmp"),
                      col.names=FALSE,
                      row.names=FALSE,
                      quote=FALSE)
          if(substr(try(FindU(InSeries=paste0(outpath,"/tmp"),
                              MissingValueCode="NA",
                              p.lev=0.95,
                              Iadj=10000,
                              Mq=10, Ny4a=0,
                              output=paste0(outpath,"/tmp")),
                        silent=TRUE),1,5)[1] != "Error") {
            ## The variable 'status' catches possible errors that will then trigger warnings
            status <- FindU(InSeries=paste0(outpath,"/tmp"),
                            MissingValueCode="NA",
                            p.lev=0.95,
                            Iadj=10000,
                            Mq=10,
                            Ny4a=0,
                            output=paste0(outpath,"/tmp"))
            if (!is.null(status)) {
              if (length(list.files(outpath,pattern="tmp_mCs.txt")) > 0) {
                status <- StepSize(InSeries=paste0(outpath,"/tmp"),
                                   MissingValueCode="NA",
                                   p.lev=0.95,
                                   Iadj=10000,
                                   Mq=10,
                                   Ny4a=0,
                                   InCs=paste0(outpath,"/tmp_mCs.txt"),
                                   output=paste0(outpath,"/tmp"))
                if(!is.null(status)) {
                  rawfile <- read.table(paste0(outpath,"/tmp_mCs.txt"),fill=TRUE)
                  if (dim(rawfile)[1]>1) {
                    rawfile <- read.table(paste0(outpath,"/tmp_mCs.txt"),skip=1)
                    j <- which(rawfile[,1]!=1 | rawfile[,2]!="Yes")
                    nb <- dim(rawfile)[1]-length(j)
                    while(length(j) > 0) {
                      rawfile <- read.table(paste0(outpath,"/tmp_mCs.txt"),sep="\t")
                      write.table(rawfile[-(j+1),],
                                  file=paste0(outpath,"/tmp_mCs.txt"),
                                  quote=FALSE, row.names=FALSE,
                                  col.names=FALSE,
                                  append=FALSE)
                      if(length(rawfile[-(j+1),])[1] > 1) {
                        status <- StepSize(InSeries=paste0(outpath,"/tmp"),
                                           MissingValueCode="NA",
                                           p.lev=0.95,
                                           Iadj=10000,
                                           Mq=10,
                                           Ny4a=0,
                                           InCs=paste0(outpath,"/tmp_mCs.txt"),
                                           output=paste0(outpath,"/tmp"))
                        if (!is.null(status)) {
                          rawfile <- read.table(paste0(outpath,"/tmp_mCs.txt"),skip=1)
                          j <- which(rawfile[,1]!=1 | rawfile[,2]!="Yes")
                          nb <- dim(rawfile)[1]-length(j)
                        } else {
                          warning(paste0(candiname,": impossible to apply absolute test"))
                          breaks[[cname]]$abs <- FALSE
                          j <- c()
                          nb <- 0
                        }
                      } else {
                        j <- c()
                        nb <- 0
                      }
                    }
                    if (nb>0) inho_ab <- as.numeric(substr(rawfile[,3],1,4))
                  }
                } else {
                  warning(paste0(candiname,": impossible to apply absolute test"))
                  breaks[[cname]]$abs <- FALSE
                }
              }
              system(paste0("rm ",outpath,"/tmp_*"))
            } else {
              warning(paste0(candiname,": impossible to apply absolute test"))
              breaks[[cname]]$abs <- FALSE
            }
          } else {
            warning(paste0(candiname,": impossible to apply absolute test"))
            breaks[[cname]]$abs <- FALSE
          }
          system(paste0("rm ",outpath,"/tmp"))
        }
        
        
        ## Sort out breakpoints
        ## The number of detections for each year are written into the data frame 'breaks_cnt' 
        breaks_cnt <- data.frame(years,rep(0,length(years)))
        names(breaks_cnt) <- c("Year","Count")
        if (h >= 3) {
          breaks_raw <- list(caume=inho_cm,wang=inho_wf,gahmdi=inho_ga)
          brm <- list()
          for (meth in names(breaks_raw)) {
            brm[[meth]] <- c()
            if (!is.null(breaks_raw[[meth]])) {
              br <- unique(breaks_raw[[meth]][!is.na(breaks_raw[[meth]])])
              for (ib in br) {
                if (sum(breaks_raw[[meth]]%in%c(ib-1,ib,ib+1)) >= 3) {
                  brm[[meth]] <- append(brm[[meth]],ib)
                }
              }
              if (length(brm[[meth]] > 1)) brm[[meth]] <- sort(brm[[meth]])
              breaks_cnt$Count[which(breaks_cnt$Year%in%brm[[meth]])] <-
                breaks_cnt$Count[which(breaks_cnt$Year%in%brm[[meth]])]+1
            }
          }
        }
        breaks_cnt$Abs <- rep(0,length(breaks_cnt$Count))
        if(!is.null(inho_ab)) breaks_cnt$Abs[which(breaks_cnt$Year%in%inho_ab)] <- 1
        breaks[[cname]][[season]] <- breaks_cnt
        
        
      } # end loop over seasons
      
      
      ## Merge yearly and seasonal breaks
      breaks[[cname]]$all <- breaks[[cname]]$y
      for (i in 1:dim(breaks[[cname]]$all)[1]) {
        breaks[[cname]]$all$Count[i] <- max(breaks[[cname]]$all$Count[i],
                                            breaks[[cname]]$w$Count[i],
                                            breaks[[cname]]$s$Count[i])
        breaks[[cname]]$all$Abs[i] <- max(breaks[[cname]]$all$Abs[i],
                                          breaks[[cname]]$w$Abs[i],
                                          breaks[[cname]]$s$Abs[i])
      }
      
      
      ######################################################
      ### BREAK DETECTION ON MONTHLY DATA (t-test) ###
      #######################################################
      
      if (!is.null(breaks[[cname]]$all)) {
        
        ## Extract breakpoint years (local maxima in detections) - Relative tests only
        ## Each breakpoint must have been detected at least twice within +- 1 year
        fdiff <- diff(c(0,breaks[[cname]]$all$Count))
        counts <- filter(breaks[[cname]]$all$Count, rep(1,3))
        counts[is.na(counts)] <- 0
        bpyears <- breaks[[cname]]$all$Year[c(diff(fdiff>=0),0)<0 & counts>=thr]
        breaks[[cname]]$breakpoints <- data.frame(Year=bpyears, Month=rep(NA,length(bpyears)), 
                                                  Size=rep(NA,length(bpyears)), p=rep(NA,length(bpyears)))
        
        if (monthly) {
          
          bnew <- data.frame(y=array(dim=length(bpyears)),m=array(dim=length(bpyears)),
                             p=array(dim=length(bpyears)),j=array(dim=length(bpyears)))
          ## Loop on breakpoints
          for (b in bpyears) {
            ## Define window (narrower near the edges and near other breakpoints)
            ib <- max(b-years[1]+1-as.integer(w/2),
                      ifelse(is.na(bnew$y[1]),1,max(bnew$y,na.rm=TRUE)-years[1]+2),
                      na.rm=TRUE):min(b-years[1]+1+as.integer(w/2),
                                      ifelse(length(which(bpyears>b))>0,min(bpyears[which(bpyears>b)])-years[1],dim(mondata[[1]])[1]),
                                      na.rm=TRUE)
            refcors <- c()
            #if (b%in%breaks[[cname]]$y) { #check which season is involved
            seamonths <- 1:12
            #} else if (!is.null(breaks[[cname]]$w) & !is.null(breaks[[cname]]$s)) {
            #if (min(abs(b-breaks[[cname]]$w))<=2) {
            #if (min(abs(b-breaks[[cname]]$s))<=2) {
            #seamonths <- 1:12
            #} else {
            #seamonths <- c(1:3,10:12)
            #}
            #} else if (min(abs(b-breaks[[cname]]$s))<=2) {
            #seamonths <- 4:9
            #}
            #} else if (!is.null(breaks[[cname]]$w) & is.null(breaks[[cname]]$s)) {
            #seamonths <- c(1:3,10:12)
            #} else if (is.null(breaks[[cname]]$w) & !is.null(breaks[[cname]]$s)) {
            #seamonths <- 4:9
            #} else {
            #stop(paste("ERROR: Impossible to obtain season for breakpoint in",b,"for station",cname))
            #}
            nm <- length(seamonths)
            refseries <- array(data=0,dim=length(ib)*nm)
            for (i in 1:length(refnames)) {
              if (sum(is.na(mondata[[refnames[i]]][ib,seamonths]))<0.05*length(ib)*nm) { #check that enough data are available in the references
                refseries <- refseries + 
                  (c(t(mondata[[refnames[i]]][ib,seamonths])) -
                     rep(colMeans(mondata[[refnames[i]]][ib,seamonths],na.rm=TRUE),length(ib))) * 
                  (mydat$rx[i]**2)
                refcors <- append(refcors, mydat$rx[i]**2)
              }
            }
            if (length(refcors)>0) refseries <- refseries / sum(refcors)
            q <- c(t(mondata[[candiname]][ib,seamonths])) - 
              rep(colMeans(mondata[[candiname]][ib,seamonths],na.rm=TRUE),length(ib)) -
              refseries
            p <- array(dim=length(q))
            for (i in nm:(nm*(length(ib)-1))) { #look for the break inside the defined window
              if (!is.na(q[i])) {
                if (sum(!is.na(q[(i-(nm-1)):i]))<2) {
                  p[i-(nm-1)] <- 0 #put the break where there is a large gap
                } else if (sum(!is.na(q[(i+1):(i+nm)]))<2) {
                  p[i+1] <- 0 #put the break where there is a large gap
                } else p[i] <- t.test(x=q[1:i],y=q[(i+1):(length(ib)*nm)],na.action="na.omit")$p.value
              }
            }
            bnew$y[which(bpyears==b)] <- as.integer(ib[1]+years[1]-1+(which.min(p)-1)/nm)
            bnew$m[which(bpyears==b)] <- seamonths[which.min(p)-as.integer((which.min(p)-1)/nm)*nm]
            bnew$p[which(bpyears==b)] <- min(p,na.rm=TRUE)
            if (sum(!is.na(q[1:which.min(p)]))>=12 & sum(!is.na(q[(which.min(p)+1):length(q)]))>=12) { #estimate break size (at least 12 values per side required)
              j_est <- t.test(x=q[1:which.min(p)],y=q[(which.min(p)+1):length(q)],na.action="na.omit")$estimate
              bnew$j[which(bpyears==b)] <- j_est[2]-j_est[1]
            }
          }
          breaks[[cname]]$breakpoints <- data.frame(Year=bnew$y, Month=bnew$m, Size=bnew$j, p=bnew$p)
          
        }
        
      }
      
      
      ##############
      ### OUTPUT ###
      ##############

      ## Graphical output produced only if there are at least 3 reference series for yearly data
      if (graphout & h>=3) {
        
        ## Create pdf
        setwd(outpath)
        pdf(paste0(cname,".pdf"),width=7.0,height=10.0,title=cname,paper="a4")
        
        ## Page 1
        oldpar <- par()
        par(mfrow=c(3,1), mar=c(5,4,5,2)+0.1)
        if (semiy) seasons <- c("y","w","s")
        for (season in seasons) {
          plot(years,rep(NA,length(years)),xlab="",ylab="",ylim=c(-4,4),lab=c(10,10,7))
          axis(4,-4:4)
          grid(lty=1,lwd=0.5)
          color <- rainbow(length(rref[[season]]))
          maintitle <- "Yearly standardized differences with reference series ("
          if (season=="y") {
            if(substr(cname,1,1) == "n") title(main=paste0(maintitle,"T Min)"),line=1,cex.main=1)
            if(substr(cname,1,1) == "x") title(main=paste0(maintitle,"T Max)"),line=1,cex.main=1)
            if(substr(cname,1,1) == "t") title(main=paste0(maintitle,"T Mean)"),line=1,cex.main=1)
            if(substr(cname,1,1) == "d") title(main=paste0(maintitle,"DTR)"),line=1,cex.main=1)
            if(substr(cname,1,1) == "a") title(main=paste0(maintitle,"T Avg)"),line=1,cex.main=1)
            text(mean(years),6.5,sname,cex=1.5,xpd=TRUE)
            abline(h=6,lwd=2,xpd=TRUE)
          }
          if (season == "w") title(main="ONDJFM semester",line=1,cex.main=1)
          if (season == "s") title(main="AMJJAS semester",line=1,cex.main=1)
          for (i in (which(bdscore[[cname]][[season]]$Score<cormin*3 &
                           !is.na(mytx[,1]))+years[1]-1)) {
            rect(i-0.5,-4.2,i+0.5,4.2,density=20,angle=45,col="red",border=NA,lwd=0.3)
          }
          for (i in (which(bdscore[[cname]][[season]]$Score>=cormin*3 &
                           bdscore[[cname]][[season]]$Score<5 & !is.na(mytx[,1]))+years[1]-1)) {
            rect(i-0.5,-4.2,i+0.5,4.2,density=20,angle=45,col="yellow",border=NA,lwd=0.3)
          }
          if (!is.null(dim(z[[season]])[2])) {
            for (i in 1:(dim(z[[season]])[2]-1)) {
              lines(years,z[[season]][,i+1],col=color[i])
            }
            lines(years,rowMeans(z[[season]][,-1],na.rm=TRUE),lwd=2,col="#4C4646")
            breaks_cnt <- breaks[[cname]][[season]][which(breaks[[cname]][[season]]$Count>0),]
            breaks_abs <- breaks[[cname]][[season]][which(breaks[[cname]][[season]]$Abs>0),]
            abline(v=breaks_cnt$Year,lwd=breaks_cnt$Count,lty=2)
            abline(v=breaks_abs$Year,lwd=1,lty=3,col="blue")
            legend(years[1],-5.5,
                   c(paste0(names(rref[[season]]),", r=",round(rref[[season]],2)),
                     "arithmetic mean"),
                   col=c(color,"#4C4646"),
                   lwd=c(rep(1,dim(z[[season]])[2]-1),2),
                   cex=0.8,
                   inset=c(0,-0.5),
                   bty="n",
                   xpd=TRUE,
                   ncol=4)
          }
        }
        
        ## Page 2
        par(oldpar)
        map(xlim=c(lon[1]-5,lon[1]+5), ylim=c(max(lat[1]-5,-90),min(lat[1]+5,90)), 
            fill=TRUE, col="#E6E6E6", lty=1)
        #map.axes()
        points(c(lon[2:(h+1)],lon[1]),c(lat[2:(h+1)],lat[1]),
               bg=c(color,"black"), cex=2, pch=c(rep(21,h),8))
        text(lon[1],min(lat[1]+10,90)+3,
             paste("Longitude:",lon[1],"   Latitude:",lat[1],"   Elevation:",ele[1],"m"),
             xpd=TRUE)
        abline(h=max(lat[1]-10,-90)-15.03,xpd=TRUE)
        text(lon[1],max(lat[1]-10,-90)-15.4,
             paste("Document created by the EUSTACE break-detection software version",v),
             pos=2, cex=0.5, xpd=TRUE)
        
        dev.off()
        
      }
      
      
    } # end if on data availability
    
  } else { 
    header <- read.table(paste0(inpath,"/",candiname),
                         comment.char="",
                         quote="",
                         sep="\t",
                         skip=1,
                         nrows=1,
                         stringsAsFactors=FALSE)[,1]
    cname <- rev(strsplit(candiname,"/")[[1]])[1]
    breaks[[cname]] <- list()
    breaks[[cname]]$range <- c(NA,NA)
    breaks[[cname]]$coord <- c(as.numeric(substr(header,29,39)),
                               as.numeric(substr(header,16,25)),
                               as.numeric(substr(header,44,50)))
    breaks[[cname]]$abs <- FALSE
    breaks[[cname]]$bd <- FALSE
    bdscore[[cname]] <- list()
  } #end if on empty files

  
  ## Save workspace image
  setwd(outpath)
  save(list=c("breaks","bdscore"),file=paste0(fl,"_breaks.RData"))
  setwd(workpath)

  
} # end loop over candidate stations


setwd(outpath)


## Write tabular output
out <- list()
for (nam in names(breaks)) {
  if (!is.null(breaks[[nam]]$breakpoints)) {
    if(nrow(breaks[[nam]]$breakpoints) > 0) {
      out[[nam]] <- cbind(nam, breaks[[nam]]$breakpoints)
      names(out[[nam]]) <- c("File", "Year", "Month", "Size", "p")
    }
  }
}
out <- Reduce(rbind, out)
out$File <- sub(paste0(fl,"_"), "", out$File)
write.table(out, file=paste0(fl,"_breakpoints.csv"), row.names=FALSE, quote=FALSE, sep=";")


## Create png images
if (!exists('candidate') & graphout) {
  
#  setwd(outpath)
  png(paste0(fl,'_breaks.png'),res=200,width=1920,height=1200)
  par(mar=c(5,4,4,4)+0.1)
  years <- yrange[1]:yrange[2]
  nyears <- length(years)
  allbreaks <- c()
  nstations <- array(dim=nyears,data=0)
  for (st in names(breaks)) {
    if (breaks[[st]]$bd) {
      if (!is.null(breaks[[st]]$breakpoints$Year)) {
        allbreaks <- append(allbreaks, breaks[[st]]$breakpoints$Year)
      }
      nstations[breaks[[st]]$all$Year-yrange[1]+1] <- nstations[breaks[[st]]$all$Year-yrange[1]+1] + 1
    }
  }
  nbreaks <- vector(length=nyears)
  names(nbreaks) <- years
  for (i in years) nbreaks[i-years[1]+1] <- sum(i==allbreaks)
  x <- barplot(height=nbreaks/nstations, names.arg=years, ylab='Breakpoints per series', cex.axis=1.5, 
               cex.lab=1.5, cex.names=1.5, cex.main=1.5, ylim=c(0,0.25), 
               main=paste('Total number of breaks:',length(allbreaks),'- Average homogeneous station-years:',
                          round(sum(nstations)/length(allbreaks),1)),
               xaxt="n", add=FALSE)
  grid(NA,NULL,lty=1)
  barplot(height=nbreaks/nstations, names.arg=years, ylab='Breakpoints per series', cex.axis=1.5, 
          cex.lab=1.5, cex.names=1.5, cex.main=1.5, ylim=c(0,0.25), 
          main=paste('Total number of breaks:',length(allbreaks),'- Average homogeneous station-years:',
                     round(sum(nstations)/length(allbreaks),1)),
          xaxt="n", add=TRUE)
  lines(x, 0.25*nstations/pretty(max(nstations))[2], lwd=2)
  axis(1, at=x[seq(1,length(years),20),1], labels=seq(yrange[1],yrange[2],20), cex.axis=1.5)
  axis(4, at=seq(0,0.25, by=0.05),
       labels=seq(0,pretty(max(nstations))[2],by=pretty(max(nstations))[2]/5), cex.axis=1.5, cex.lab=1.5)
  mtext("Number of stations", 4, 3, cex=1.5)
  box()
  dev.off()
  
  if (monthly) {
    png(paste0(fl,'_jumps.png'),res=200,width=1920,height=1200)
    alljumps <- c()
    for (st in names(breaks)) {
      if (!is.null(breaks[[st]]$breakpoints$Size)) alljumps <- append(alljumps,breaks[[st]]$breakpoints$Size)
    }
    hist(alljumps, pretty(c(min(alljumps,na.rm=TRUE),max(alljumps,na.rm=TRUE)),50), freq=FALSE, col='grey', 
         main=paste0('Average: ',round(mean(alljumps,na.rm=TRUE),3),'°C'), xlab='[°C]', 
         cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
    grid(NA,NULL,lty=1)
    hist(alljumps, pretty(c(min(alljumps,na.rm=TRUE),max(alljumps,na.rm=TRUE)),50), freq=FALSE, col='grey', 
         main=paste0('Average: ',round(mean(alljumps,na.rm=TRUE),3),'°C'), 
         cex.axis=1.5, cex.lab=1.5, cex.main=1.5, add=TRUE)
    box()
    dev.off()
  }
  
  print(paste("PNG files saved in",outpath),quote=FALSE)
  
}
system("rm list list_ref")
setwd(workpath)


###############
### THE END ###
###############
