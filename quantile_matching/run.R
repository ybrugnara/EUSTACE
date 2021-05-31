## Run homogenization on all stations for a certain variable

mainpath <- "/home/yb/Desktop/MO/homogenisation" # full path of folder containing the R scripts, including this one
datapath <- paste0(mainpath,"/Input") # folder with data to be homogenized
inpath <- paste0(mainpath,"/WorkFolder") # folder where pre-processed input data will be written - NB: any file in here will be deleted
outpath <- paste0(mainpath,"/Output") # folder where homogenized data will be written
pattern <- "aTA" # initial characters of input filenames in datapath (first character must be variable-specific)

setwd(mainpath)

## Copy files to the work directory
file.remove(list.files(inpath, full.names=TRUE))
files <- list.files(datapath, pattern)
file.copy(paste0(datapath,"/",files), paste0(inpath,"/",files))
ele <- substr(pattern, 1, 1)
breaksfile <- paste0(ele, "_breaks.RData")
file.copy(paste0(datapath,"/",breaksfile), paste0(inpath,"/",breaksfile))

## Remove empty files
for (f in list.files(inpath,pattern,full.names=TRUE)) {
  x <- tryCatch(
    {
      x <- read.table(f, nrows=1)
    },
    error = function(e){
      file.remove(f)
    }
  )
}

files <- list.files(inpath, pattern)
for (target in files) {
  ## Rename target file
  previous_target <- list.files(inpath, "ser", full.names=TRUE)
  if (length(previous_target) > 0) {
    file.rename(previous_target, sub("ser","ref",previous_target))
  }
  file.rename(paste0(inpath,"/",target), paste0(inpath,"/",sub("ref","ser",target)))
  ## Run homogenisation software
  source("AA_quant_match_MO.R")
  setwd(mainpath)
}
file.copy(paste0(inpath,"/report_hom.txt"), paste0(ele,"_report_hom.txt"))
