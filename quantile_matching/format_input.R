## Transform the input/output of the breakpoint detection into the input for the homogenisation
## NB: RData files from the breakpoint detection must be first copied to the working directory

workdir <- "/home/yb/Desktop/MO/homogenisation/Input" # where the files will be written (.RData files from b.d. must be in this folder)
datapath <- "/home/yb/Desktop/MO/data/formatted" # where the daily data in the b.d. format are
codes <- c("Tx", "Tn", "Ta", "Tw") # prefix used in the filenames to distinguish variables (at least 2 characters, the 2nd must be unique)
filter_bp <- TRUE # TRUE to remove small breakpoints
merge_ta_tw <- TRUE # TRUE to merge breakpoints of Ta and Tw
inventory <- "/home/yb/Desktop/MO/data/station_list.csv"

## Work on lists of breakpoints
setwd(workdir)
filenames <- list()
inv <- read.csv(inventory)
inv$tw_id <- inv$ta_id <- inv$tn_id <- inv$tx_id <- NA
for (f in list.files(pattern="RData")) {
  load(f)
  prefix <- substr(names(breaks)[1], 1, 2)
  filenames[[prefix]] <- names(breaks)
  prefix <- paste0(substr(prefix,2,2), toupper(prefix), "_")
  for (i in 1:length(breaks)) {
    if (!is.null(breaks[[names(breaks[i])]]$breakpoints)) {
      if (nrow(breaks[[names(breaks[i])]]$breakpoints) > 0) {
        if (filter_bp) {
          ko <- which(abs(breaks[[names(breaks[i])]]$breakpoints$Size)<0.25 | 
                        breaks[[names(breaks[i])]]$breakpoints$p>0.01)
        } else {
          ko <- nrow(breaks[[names(breaks[i])]]$breakpoints) + 1
        }
        if (length(ko) == 0) ko <- nrow(breaks[[names(breaks[i])]]$breakpoints) + 1
        breaks[[names(breaks[i])]]$y <- breaks[[names(breaks[i])]]$breakpoints$Year[-ko]
        breaks[[names(breaks[i])]]$m <- breaks[[names(breaks[i])]]$breakpoints$Month[-ko]
        breaks[[names(breaks[i])]]$m[which(is.na(breaks[[names(breaks[i])]]$m))] <- 12
      } else {
        breaks[[names(breaks[i])]]$y <- c()
        breaks[[names(breaks[i])]]$m <- c()
      }
    }
    names(breaks)[i] <- paste0(prefix, 1000000+i)
  }
  save(breaks, file=sub("T","",f))
}

## Work on daily data files
file.remove(list.files(pattern="ref"))
for (cod in codes) {
  prefix <- paste0(substr(cod,2,2), toupper(cod), "_")
  newnames <- paste0(prefix, 1000001:(1000000+length(filenames[[cod]])), ".ref")
  file.copy(paste0(datapath,"/",cod,"/",filenames[[cod]]), newnames)
  for (i in 1:length(filenames[[cod]])) {
    nam1 <- filenames[[cod]][i]
    st <- substr(nam1, 4, nchar(nam1)-4)
    nam2 <- newnames[i]
    inv[[paste0(tolower(cod),"_id")]][inv$station_file_name==st] <- substr(nam2, 5, 11)
  }
}

## Write inventory
write.csv(inv, "station_list.csv", row.names=FALSE, quote=FALSE)

## Merge Ta and Tw
if (merge_ta_tw) {
  load("a_breaks.RData")
  a_breaks <- breaks
  load("w_breaks.RData")
  for (i in 1:length(breaks)) {
    filename_ta <- names(a_breaks)[i]
    id <- substr(filename_ta, 5, nchar(filename_ta))
    filename_tw <- paste0("wTW_", inv$tw_id[which(inv$ta_id==id)])
    if (!is.null(breaks[[filename_tw]]$y) & 
        !is.null(a_breaks[[filename_ta]]$y)) {
      breaks[[filename_tw]]$y <- append(breaks[[filename_tw]]$y, 
                                             a_breaks[[filename_ta]]$y)
      breaks[[filename_tw]]$m <- append(breaks[[filename_tw]]$m, 
                                             a_breaks[[filename_ta]]$m)
      if (length(breaks[[filename_tw]]$y) > 1) {
        neworder <- order(breaks[[filename_tw]]$y, breaks[[filename_tw]]$m)
        breaks[[filename_tw]]$y <- breaks[[filename_tw]]$y[neworder]
        breaks[[filename_tw]]$m <- breaks[[filename_tw]]$m[neworder]
        dupl <- diff(breaks[[filename_tw]]$y)
        dupl <- c(FALSE, dupl<=1)
        breaks[[filename_tw]]$y <- breaks[[filename_tw]]$y[!dupl]
        breaks[[filename_tw]]$m <- breaks[[filename_tw]]$m[!dupl]
      }
      a_breaks[[filename_ta]]$y <- breaks[[filename_tw]]$y
      a_breaks[[filename_ta]]$m <- breaks[[filename_tw]]$m
    } else if (!is.null(a_breaks[[names(breaks[i])]]$y)) {
      breaks[[filename_tw]]$y <- a_breaks[[filename_ta]]$y
      breaks[[filename_tw]]$m <- a_breaks[[filename_ta]]$m
    }
  }
  save(breaks, file="w_breaks.RData")
  breaks <- a_breaks
  save(breaks, file="a_breaks.RData")
}