### SCRIPT FOR THE CALCULATION OF QUANTILE MATCHING ADJUSTMENTS
# Candidate and reference series must be in separate files in the inpath with format 'yyyy mm dd value'.

# Name of the file is fundamental:
# paste0(ele,'_',ser_id+1000000,'.ser') for the candidate
# paste0(ele,'_',ser_id+1000000,'.ref') for the references
# where ele is 'TN', 'TX', 'TA' or 'TW', and ser_id is a positive integer lower than 999999
# no other files can end with '.ser' or '.ref'

# breaks must be written in the format of a list 
# names of the elements have to be in the format xTX_1000000.ser
# (for all the series, including references the extension has to be '.ser')
# the fields of the element have to include $range, $y (year of breakpoint) and $m (month of breakpoint)
# the list of breaks must be stored in inpath as an R workspace
# with the name 'n_breaks.RData', 'x_breaks.RData', 'a_breaks.RData', or 'w_breaks.RData'

# version edited by Yuri Brugnara for the homogenisation of UKMO data (March 2021)
# NB: this script works for one station only; use script run.R to homogenise a whole data set


### PACKAGES
suppressMessages(library(gdata))
suppressMessages(library(trend))

### DIRECTORIES (defined in run.R)
#mainpath="/home/yb/Desktop/MO/homogenisation/"                       #where are the functions?
#inpath=paste0(mainpath,'WorkFolder')                                 #where are the series?
#outpath=paste0(mainpath,'Output')                               #where to put the temporary results?
#finalpath=outpath                                                  #where to put the final results?
rprt='report_hom.txt'

##PARAMETERS
#ele='n'
refyon=1    #do you want to use references for the correction calculation? Yes(1) or No(0)
anyon=0     #do you want to calculate correction factors through anomalies calculation? Yes(1) or No(0) 
thre_y=3    #minimum years of overlapping data needed to calculate the percentiles
ampl_bin=5  #amplitude bins of percentiles
error_report=1 #do you want errors to be reported in a text file? Yes(1) or No(0)
thre_corr=0.7 #threshold on correlation between candidate and reference series
cns_then_sm=0 #check negative slopes before smoothing out (1 for yes, 0 for no)


fq=ampl_bin
lq=100-fq
fb=(ampl_bin*3/2)          #first boundary
lb=100-fb                       #last boundary
nb=length(seq(fb,lb,ampl_bin))+2  #number of boundaries
nq=nb-1                         #number of intervals
endmonth=c(31,28,31,30,31,30,31,31,30,31,30,31)
length_ref=20
if (nq%%2==0)
{
  fb=fb-ampl_bin/2
  lb=lb+ampl_bin/2
  nb=length(seq(fb,lb,ampl_bin))  #number of boundaries
  nq=nb-1
}


### FUNCTIONS
setwd(mainpath)
source('AA_functions_hom.R')

setwd(inpath)

write('#####################################################',rprt,append=TRUE)
write(as.character(Sys.time()), rprt, append=TRUE)

##GET THE SERIES


setwd(inpath)
files.ser <- list.files(pattern=glob2rx(paste0(ele,"*ser"))) #list of targets
files.wrk = files.ser
files.ref <- list.files(pattern=glob2rx(paste0(ele,"*ref")))  #list of references
files=c(files.ser,files.ref)            #general list of files


series.ori=lapply(files.wrk,read.table) 
names(series.ori)=files.wrk

nox <- substr(names(series.ori)[1], 1, 1)
if (nox == "n") {
  load("n_breaks.RData")
  TNoTX <- "TN"
} else if (nox == "x") {
  load("x_breaks.RData")
  TNoTX <- "TX"
} else if (nox == "a") {
  load("a_breaks.RData")
  TNoTX <- "TA"
} else if (nox == "w") {
  load("w_breaks.RData")
  TNoTX <- "TW"
} 


if (length(files.ref)<3)
{
  print('Not enough references (maybe not downloaded since the series of this station have no breaks)')
  write('Not enough references (maybe not downloaded since the series of this station have no breaks). No corrections will be applied',rprt,append=TRUE)
  for (s in 1:length(series.ori)) #loop on the series of the same station
  {
    series=series.ori[[s]]
    colnames(series)=c('year','month','day','series')
    series_name=names(series.ori[s])
    ser_id=extract_ser_id(series_name)
    aux=substr(series_name,1,10)
    output_s=paste0(aux,'.hom')
    setwd(outpath)
    fs=as.Date(paste(series[1,1],series[1,2],series[1,3],sep='-'))
    ls=as.Date(paste(series[dim(series)[1],1],series[dim(series)[1],2],series[dim(series)[1],3],sep='-'))
    series_hom=series
    colnames(series_hom)=c('year','month','day','hom')
    series_hom=complete(series_hom,fs,ls)
    series_hom$qc=8
    series_hom$qca=8
    series_hom[series_hom$hom<(-90),c('qc','qca')]=9
    setwd(outpath)
    hom4out=series_hom[,1:4]
    write.fwf(hom4out,output_s,colnames=FALSE,sep=', ',eol='\n',append=FALSE)
    setwd(inpath)
  }
}else
{
  series.ref=lapply(files.ref,read.table)
  names(series.ori)=files.ser
  names(series.ref)=files.ref
  
  
  ## DIVIDE REFERENCE SERIES INTO HOMOGENEOUS SUBSERIES
  ref_names=names(series.ref)
  list_rm_ref=NULL
  for (j in 1:length(series.ref))
  {
    ref_id=extract_ser_id(ref_names[j])
    colnames(series.ref[[j]])=c('year','month','day','refer')
    if (ncol(series.ref[[j]])>4)
    {
      series.ref[[j]]=series.ref[[j]][,1:4]
    }
    if (dim(series.ref[[j]])[1]<thre_y*365)
    {
      list_rm_ref=c(list_rm_ref,j)
    }
  }
  if(length(list_rm_ref)>0)
  {
    series.ref=series.ref[-list_rm_ref]
  }
  refers=list()
    i=0
  for (j in 1:length(series.ref))
  {
  
    ref_ser=series.ref[[j]]
    colnames(ref_ser)=c('year','month','day','refer')
    ref_name=names(series.ref[j])
    ref_id=extract_ser_id(ref_name,base=1)
    chp_ref=read.breaks.simple(ref_id,nox,breaks=breaks)
    chp_y_ref=chp_ref$year
    chp_m_ref=chp_ref$month   

    if (length(chp_y_ref)==0)
    { 
      i=i+1
      refers[[i]]=ref_ser
      names(refers)[i]=paste0(ref_name,'-00')
    } else
    {
      ref_rest=ref_ser
      colnames(ref_rest)=c('year','month','day','refer')
      for (k in 1:length(chp_y_ref))
      {
        ref_fin=ref_rest[ref_rest$year<chp_y_ref[k] | (ref_rest$year==chp_y_ref[k] & ref_rest$month<chp_m_ref[k]),]
        ref_rest=ref_rest[ref_rest$year>chp_y_ref[k] | (ref_rest$year==chp_y_ref[k] & ref_rest$month>=chp_m_ref[k]),]
        if (dim(ref_fin)[1]>0)
        {
          i=i+1
          refers[[i]]=ref_fin
          names(refers)[i]=paste0(ref_name,'-0',k)
        }
      }
      if (dim(ref_rest)[1]>0)
      {
        i=i+1
        refers[[i]]=ref_rest
        names(refers)[i]=paste0(ref_name,'-0',k+1)
      }
    }
  }
  refers_all=refers
  
  
  
 for (s in 1:length(series.ori)) #loop on the series of the same station
 
   {
    series=series.ori[[s]]
    colnames(series)=c('year','month','day','series')
   
    if(ncol(series)>4)
    {
      series=series[,1:4]
    }
    series_name=names(series.ori[s])
    write(paste('Working on series',series_name),rprt,append=TRUE)
    ser_id=extract_ser_id(series_name)
  
    fs=as.Date(paste(series[1,1],series[1,2],series[1,3],sep='-'))
    ls=as.Date(paste(series[dim(series)[1],1],series[dim(series)[1],2],series[dim(series)[1],3],sep='-'))
    
 
    aux=substr(series_name,1,11)
    output_a=paste0(aux,'.adj')
    output_s=paste0(aux,'.qm')
 

    ## DIVIDE CANDIDATE SERIES INTO HOMOGENEOUS SUBSERIES
    
    fys=series[1,'year']
    lys=series[dim(series)[1],'year']
    chp=read.breaks.simple(ser_id,nox,breaks=breaks)
    chp_y=chp$year
    chp_m=chp$month

    chp_y_old=chp_y
    write( 'Breaks:',rprt,append=TRUE)
    chp_m=chp_m[which(chp_y>fys & chp_y<lys)]
    chp_y=chp_y[which(chp_y>fys & chp_y<lys)]
    write(paste(chp_y,chp_m),rprt,append=TRUE)
    
    
    if (length(chp_y)>0) #proceed only if there is at least one break
    {
   
      donats=list()
      i=0
      ser_rest=series
      for (k in 1:length(chp_y))
      {
        i=i+1
        ser_fin=ser_rest[ser_rest$year<chp_y[k] | (ser_rest$year==chp_y[k] & ser_rest$month<chp_m[k]),]
        ser_rest=ser_rest[ser_rest$year>chp_y[k] | (ser_rest$year==chp_y[k] & ser_rest$month>=chp_m[k]),]
        donats[[i]]=ser_fin
        colnames(donats[[i]])=c('year','month','day','donat')
        names(donats)[i]=paste0(series_name,'-0',k)
      }
      basis=ser_rest
      if (dim(basis)[1]<thre_y*365)
      {
        colnames(donats[[length(chp_y)]])=c('year','month','day','series')
        basis=rbind(donats[[length(chp_y)]],basis)
        donats[[length(chp_y)]]=NULL
      }
      basis_name=paste0(series_name,'-0',k+1)
      colnames(basis)=c('year','month','day','basis')
      basis_ser_id=extract_ser_id(series_name)
      benchmark=basis
      
      overlap_dr=list()
      qntl_drs=list()
      bndrs_ds=list()
      corrs=list()
      corrs01=list()
      corrs02=list()
      basis=benchmark
      
      
      setwd(outpath)
      intro=paste0("#",TNoTX," 2 iter adjustments for original series: ",ser_id)
      hdr=paste0(format('#ser_id,',width=12,justify='right'),
                 format('break_start,',width=12,justify='right'),
                 format('break_stop,',width=12,justify='right'),
                 format('ref_id,',width=12, justify='right'),
                 format('month,',width=12, justify='right'),
                 format('percentile,',width=12, justify='right'),
                 format('lower_bnd,',width=12, justify='right'),
                 format('upper_bnd,',width=12, justify='right'),
                 format('adj,',width=12, justify='right'))
      #write(intro,file=output_a)
      #write(hdr,file=output_a,append=TRUE)
      setwd(inpath)
      
      if (length(donats)>=1)
      {
        
        
        for (d in length(donats):1) #loop on donats
        {
          donat=donats[[d]]
          donat_ser_id=extract_ser_id(names(donats)[d])
          if (dim(donat)[1]>=thre_y*365)
          {
            # print("ciccia")
            
            
            ##CHOOSE A SET OF REFERENCE SERIES REMOVING NON OVERLAPPING SERIES
            
            refers=refers_all
            refers[[paste0(files.wrk[s],'-01')]]=NULL
            refers[[paste0(files.wrk[s],'-02')]]=NULL
            refers[[paste0(files.wrk[s],'-03')]]=NULL
            refers[[paste0(files.wrk[s],'-04')]]=NULL
            fyb=basis[1,'year']
            fmb=basis[1,'month']
            lyb=basis[dim(basis)[1],1]
            fyd=donat[1,'year']
            fmd=donat[1,'month']
            lyd=donat[dim(donat)[1],'year']
            lmd=donat[dim(donat)[1],'month']
         
            nyr=lapply(refers, function(x) {x[dim(x)[1],1]-x[1,1]})
            fyr=lapply(refers, function(x) {x[1,1]})
            lyr=lapply(refers, function(x) {x[dim(x)[1],1]})
            
            
            
            i=1
           
            
            i=1
            while (i<=length(refers))
            {
              if (fyr[[i]]>lyb-thre_y+1 | lyr[[i]]<fyb+thre_y-1)
              { 
                refers[[i]]=NULL
                fyr[[i]]=NULL
                lyr[[i]]=NULL
                nyr[[i]]=NULL
             
              }else
              {
                if (fyr[[i]]>lyd-thre_y+1 | lyr[[i]]<fyd+thre_y-1)
                { 
                  refers[[i]]=NULL
                  fyr[[i]]=NULL
                  lyr[[i]]=NULL
                  nyr[[i]]=NULL
            
                }else
                {
                  i=i+1
                }
              }
            } #check if there is enough overlap with basis and reference
            
            refers=lapply(refers,function(x) cut.series(series=x,y=fyb,m=fmb,ny=length_ref))
            nyr=lapply(refers, function(x) {x[dim(x)[1],1]-x[1,1]})
            fyr=lapply(refers, function(x) {x[1,1]})
            lyr=lapply(refers, function(x) {x[dim(x)[1],1]})
            overlaps=lapply(refers,function(x) overlapize_ref(x,basis))
            correlations=sapply(overlaps,function(x) round(cor(x$refer,x$basis,use='pairwise.complete.obs'),digits=3))
            
            value=ser_id+1000000
            filename=paste0('details_ref_',TNoTX,'_',value,'_',chp_y[d],'.txt')
            setwd(outpath)
            
            if (length(correlations)>0)
            {
              write(paste(names(correlations)[1],fyr[[1]],lyr[[1]],round(correlations[1],digits = 4)),filename,append=FALSE)
              if (length(correlations)>1)
              {
                for (i in 2:length(correlations)){write(paste(names(correlations)[i],fyr[[i]],lyr[[i]],round(correlations[i],digits = 4)),filename,append=TRUE)}
              }
            }else
            {write('No references',filename,append=FALSE)}
            setwd(inpath)
            
            NAcorr=correlations[is.na(correlations)]
            NAcorr0=NAcorr
            while(length(NAcorr)>0)
            {
              inac=which(is.na(correlations))[1]
              correlations=correlations[-inac]
              refers[[inac]]=NULL
              fyr[[inac]]=NULL
              lyr[[inac]]=NULL
              nyr[[inac]]=NULL
             
              NAcorr=correlations[is.na(correlations)]
            }
            
            correls=correlations
            i=1
            j=1
            low_correl_ref=list()
            fylc=list()
            lylc=list()
            nylc=list()
            latlc=list()
            lonlc=list()
            low_correls=NULL
            
            while (i<=length(refers))
            {
              if (correls[i]<thre_corr)
              { 
                low_correl_ref[[j]]=refers[[i]]
                low_correls[j]=correls[i]
                names(low_correl_ref)[j]=names(refers)[i]
                names(low_correls)[j]=names(refers)[i]
                fylc[[j]]=fyr[[i]]
                names(fylc)[j]=names(refers)[i]
                lylc[[j]]=lyr[[i]]
                names(lylc)[j]=names(refers)[i]
                nylc[[j]]=nyr[[i]]
                names(nylc)[j]=names(refers)[i]
                
                refers[[i]]=NULL
                fyr[[i]]=NULL
                lyr[[i]]=NULL
                nyr[[i]]=NULL
           
                correls=correls[-i]
                j=j+1
              }else
              {
                i=i+1
              }
            }
            
            i=1
            if (length(refers)>18) #if there are more than 18 reference, select the most correlated with the basis
            { 

              index_cr=order(correls,decreasing=TRUE)
              mincr=correls[[index_cr[19]]]
              setwd(outpath)
              
              write(paste(''),filename,append=TRUE)
              write(paste('threshold corr      ','    ','    ',round(mincr,digits = 4)),filename,append=TRUE)
              setwd(inpath)
              i=1
              while (i<=length(refers))
              {
                if (correls[i]<=mincr)
                {
                  refers[[i]]=NULL
                  fyr[[i]]=NULL
                  lyr[[i]]=NULL
                  nyr[[i]]=NULL
            
                  correls=correls[-i]
                }else
                {
                  i=i+1
                }
              }
            }
          
            
            ref_ser_id=lapply(names(refers),function(x) extract_ser_id(x))
          
            
            if (length(refers)<=5) ### add inhomogeneous series if the subseries are not enough
            {
              write('Not enough homogeneous reference series. Some non homogeneous series will be selected.',rprt,append=TRUE)
              print('Not enough homogeneous reference series. Some non homogeneous series will be selected.')
              ref_ser_id_ori=list()

              ref_ser_id_ori=lapply(names(series.ref),function(x) extract_ser_id(x))
              ref_recover=series.ref[c(!ref_ser_id_ori%in%ref_ser_id)]
              nyrr=lapply(ref_recover, function(x) {x[dim(x)[1],1]-x[1,1]})
              fyrr=lapply(ref_recover, function(x) {x[1,1]})
              lyrr=lapply(ref_recover, function(x) {x[dim(x)[1],1]})
        
             
              i=1
              while (i<=length(ref_recover))
              {
                if (fyrr[[i]]>lyb-thre_y | lyrr[[i]]<fyb+thre_y)
                { 
                  ref_recover[[i]]=NULL
                  fyrr[[i]]=NULL
                  lyrr[[i]]=NULL
                  nyrr[[i]]=NULL
               
                }else
                {
                  if (fyrr[[i]]>lyd-thre_y | lyrr[[i]]<fyd+thre_y)
                  { 
                    ref_recover[[i]]=NULL
                    fyrr[[i]]=NULL
                    lyrr[[i]]=NULL
                    nyrr[[i]]=NULL
                   
                  }else
                  {
                    i=i+1
                  }
                }
              }
              
              ref_recover=lapply(ref_recover,function(x) cut.series(series=x,y=fyb,m=fmb,ny=length_ref))
              
              nyrr=lapply(ref_recover, function(x) {x[dim(x)[1],1]-x[1,1]})
              fyrr=lapply(ref_recover, function(x) {x[1,1]})
              lyrr=lapply(ref_recover, function(x) {x[dim(x)[1],1]})
          
              i=1
              while (i<=length(ref_recover))
              {
                if (is.integer(fyrr) && is.integer(lyrr))
                {
                if (!is.na(fyrr) && !is.na(lyrr))
                {
                  if (fyrr[[i]]>lyb-thre_y | lyrr[[i]]<fyb+thre_y)
                  { 
                    ref_recover[[i]]=NULL
                    fyrr[[i]]=NULL
                    lyrr[[i]]=NULL
                    nyrr[[i]]=NULL
                 
                  }else
                  {
                    if (fyrr[[i]]>lyd-thre_y | lyrr[[i]]<fyd+thre_y)
                    { 
                      ref_recover[[i]]=NULL
                      fyrr[[i]]=NULL
                      lyrr[[i]]=NULL
                      nyrr[[i]]=NULL
   
                    }else
                    {
                      i=i+1
                    }
                  }
                }
                else
                {
                  ref_recover[[i]]=NULL
                  fyrr[[i]]=NULL
                  lyrr[[i]]=NULL
                  nyrr[[i]]=NULL
                
                }
                }
                else
                {
                  ref_recover[[i]]=NULL
                  fyrr[[i]]=NULL
                  lyrr[[i]]=NULL
                  nyrr[[i]]=NULL
                
                }
              }
              
              overlaps_rec=lapply(ref_recover,function(x) overlapize_ref(x,basis))
              i=1
              while (i<=length(ref_recover))
              {
                if (nrow(overlaps_rec[[i]])<=thre_y*365)
                { 
                  overlaps_rec[[i]]=NULL
                  ref_recover[[i]]=NULL
                  fyrr[[i]]=NULL
                  lyrr[[i]]=NULL
                  nyrr[[i]]=NULL
               
                }else
                {
                  i=i+1
                }
              }
              correlations_rec=sapply(overlaps_rec,function(x) cor(x$refer,x$basis,use='pairwise.complete.obs'))
              
              value=ser_id+1000000
              setwd(outpath)
             
              if (length(correlations_rec)>0)
              {
                filename=paste0('details_ref_',TNoTX,'_',value,'_',chp_y[d],'.txt')
                write(paste(''),filename,append=TRUE)
                write(paste('Recovery of non-homogenized stations:'),filename,append=TRUE)
                write(paste(names(correlations_rec)[1],fyrr[[1]],lyrr[[1]],round(correlations_rec[1],digits = 4)),filename,append=TRUE)
                if (length(correlations_rec)>1)
                {
                  for (i in 2:length(correlations_rec)){write(paste(names(correlations_rec)[i],fyrr[[i]],lyrr[[i]],round(correlations_rec[i],digits = 4)),filename,append=TRUE)}
                }
              }else
              {write('No more references series to recover',filename,append=TRUE)}
              setwd(inpath)
              
              NAcorr=correlations_rec[is.na(correlations)]
              NAcorr0=NAcorr
              while(length(NAcorr)>0)
              {
                inac=which(is.na(correlations_rec))[1]
                correlations_rec=correlations_rec[-inac]
                ref_recover[[inac]]=NULL
                fyrr[[inac]]=NULL
                lyrr[[inac]]=NULL
                nyrr[[inac]]=NULL
              
                NAcorr=correlations_rec[is.na(correlations_rec)]
              }
              
              correls=correlations_rec
              i=1
              j=1
              low_correl_ref=list()
              fylc=list()
              lylc=list()
              nylc=list()
            
              low_correls=NULL
              
              while (i<=length(ref_recover))
              {
                if (correls[i]<thre_corr)
                { 
                  low_correl_ref[[j]]=ref_recover[[i]]
                  low_correls[j]=correls[i]
                  names(low_correl_ref)[j]=names(ref_recover)[i]
                  names(low_correls)[j]=names(ref_recover)[i]
                  fylc[[j]]=fyrr[[i]]
                  names(fylc)[j]=names(ref_recover)[i]
                  lylc[[j]]=lyrr[[i]]
                  names(lylc)[j]=names(ref_recover)[i]
                  nylc[[j]]=nyrr[[i]]
                  names(nylc)[j]=names(ref_recover)[i]
                
                  ref_recover[[i]]=NULL
                  fyrr[[i]]=NULL
                  lyrr[[i]]=NULL
                  nyrr[[i]]=NULL
            
                  correls=correls[-i]
                  j=j+1
                }else
                {
                  i=i+1
                }
              }
              
              
              i=1
              num_ref_rec_need=5-length(refers)
              if (length(ref_recover)>num_ref_rec_need)
              {
                index_cr=order(correls,decreasing=TRUE)
                mincr=correls[[index_cr[num_ref_rec_need+1]]]
                setwd(output)
                
                write(paste(''),filename,append=TRUE)
                write(paste('thresh corr rec  ','    ','    ',round(mincr,digits = 4)),filename,append=TRUE)
                setwd(inpath)
                i=1
                while (i<=length(ref_recover))
                {
                  if (correls[i]<=mincr)
                  {
                    ref_recover[[i]]=NULL
                    fyrr[[i]]=NULL
                    lyrr[[i]]=NULL
                    nyrr[[i]]=NULL
                  
                    correls=correls[-i]
                  }else
                  {
                    i=i+1
                  }
                }
                # lrr=array(unlist(lapply(ref_recover,function(x) dim(x)[1] )))
                # index_lrr=order(lrr,decreasing=TRUE)
                # minlrr=lrr[[index_lrr[num_ref_rec_need+1]]]
                # i=1
                # while (i<=length(ref_recover))
                # {
                #   if (lrr[i]<=minlrr)
                #   {
                #     ref_recover[[i]]=NULL
                #     fyrr[[i]]=NULL
                #     lyrr[[i]]=NULL
                #     nyrr[[i]]=NULL
                #     lrr=lrr[-i]
                #   }else
                #   {
                #     i=i+1
                #   }
                # }
              }
              
              
              
              
              
              
              
              
              
              refers=c(refers,ref_recover)
              ref_ser_id=list()
              if (length(refers)>0)
              {
                for (j in 1:length(refers))
                {
                  aux=substr(names(refers)[[j]],5,11)
                  aux=as.numeric(aux)
                  ref_ser_id[[j]]=aux-1000000
                }
              }
            }
            
            if (length(refers)>=3)
            {
              
              write('List of references:',rprt,append=TRUE)
              write(names(refers),rprt,append=TRUE)
              # for (i in 1:length(refers))
              # {
              #   refer=refers[[i]]
              # }
              
              ##IDENTIFY OVERLAPPING PERIODS AND CALCULATE QUANTILES
              
              bndrs_ds[[d]]=list()
              corrs[[d]]=list()
              
              if (refyon==1)
              {
                if (error_report==1)
                {
                  tryCatch(  {overlap_br=lapply(refers, function(x) overlapize_ref(x,basis));
                  qntl_brs=lapply(overlap_br, function(x) corr_qntl_ref_3m(x,nb,nq,fq,lq,ampl_bin,thre_y))},
                  warning = function(w) {write(c(basis_name,"warning in qntl calculation of base series"), rprt, append=TRUE)},
                  error = function(e) {print('error in qntl calculation of base series');
                    write(c(sprintf('%-4s %7d %7d', nox, blend_sta_id, basis_ser_id),"error in qntl calculation of base series"), rprt, append=TRUE)})
                }else
                {
                  overlap_br=lapply(refers, function(x) overlapize_ref(x,basis))
                  qntl_brs=lapply(overlap_br, function(x) corr_qntl_ref_3m(x,nb,nq,fq,lq,ampl_bin,thre_y))
                }
              }
              
              donat_ser_ids=ser_id
              overlap_dr[[d]]=list()
              qntl_drs[[d]]=list()
              bndrs_ds[[d]]=list()
              corrs[[d]]=list()
              donat_name=names(donats)[d]
              donat_ser_id=ser_id
              print(paste('Corrs calculation of',donat_ser_id,'from',fyd,'to',lyd))
              write(paste('Corrs calculation of',donat_ser_id,'from',fyd,'to',lyd),rprt,append=TRUE)
              overlap=merge(basis,donat)
              error=0
              if (dim(overlap)[1]!=0)
              {
                error=1
                print('There is an overlap! Error somewhere')
              }
              if (error==0)
              {
                if (error_report==1)
                {
                  tryCatch(
                    {overlap_dr[[d]]=lapply(refers, function(x) overlapize_ref(x,donat)); #we need the corrections for each quantile
                    qntl_drs[[d]]=lapply(overlap_dr[[d]], function(x) corr_qntl_ref_3m(x,nb,nq,fq,lq,ampl_bin,thre_y));                   #we need also the values of the boundaries, necessary to assign each data to a quantile
                    bndrs_ds[[d]]=lapply(overlap_dr[[d]], function(x) bndrs_d_3m(x,nq,fb,lb,ampl_bin,thre_y))},
                    warning = function(w) {write(c(sprintf('%-4s %7d %7d', nox, sta_id, ser_id), "warning in qntl calculation of donating series"), rprt, append=TRUE)},
                    error = function(e) {print('error in qntl calculation of donating series');
                      write(c(sprintf('%-4s %7d %7d', nox, sta_id, ser_id),"error in qnt calculation of donating series"), rprt, append=TRUE)})
                } else
                {
                  overlap_dr[[d]]=lapply(refers, function(x) overlapize_ref(x,donat)) #we need the corrections for each quantile
                  qntl_drs[[d]]=lapply(overlap_dr[[d]], function(x) corr_qntl_ref_3m(x,nb,nq,fq,lq,ampl_bin,thre_y)) #we need also the values of the boundaries, necessary to assign each data to a quantile
                  bndrs_ds[[d]]=lapply(overlap_dr[[d]][[]], function(x) bndrs_d_3m(x,nq,fb,lb,ampl_bin,thre_y))
                }
                
                corrs_raw=list()
                corrs0=matrix(NA,nq,12)
                corrs02[[d]]=matrix(NA,nq,12)
                
                corrs_raw=list()
                corrs_cns=list()
                corrs_sm=list()
                for (j in 1:length(refers))
                {
                  corrs_raw[[j]]=round(qntl_brs[[j]]-qntl_drs[[d]][[j]],digits=1)
                  
                  if(length(corrs_raw[[j]][!is.na(corrs_raw[[j]])])>=nq)
                  {
                    if (cns_then_sm==1)
                    {
                      corrs_cns[[j]]=matrix(NA,nq,12)
                      overlap=overlap_br[[j]]
                      for (m in 1:12)
                      { 
                        corr_raw=corrs_raw[[j]][,m]
                        if (length(corr_raw[!is.na(corr_raw)])>0)
                        {
                          if (m==1) 
                          {vec=overlap[overlap$month %in% c(1,2,12) & !is.na(overlap$basis),'basis']
                          }else if (m==12)
                          {vec=overlap[overlap$month %in% c(11,12,1) & !is.na(overlap$basis),'basis']
                          }else
                          {vec=overlap[overlap$month %in% seq(m-1,m+1) & !is.na(overlap$basis),'basis']}
                          qntl=qntls_3m(vec,nb,nq,ampl_bin,thre_y)
                          if (length(qntl[!is.na(qntl)])>nq)
                          {
                            corrs_cns[[j]][,m]=check_neg_slope(corr_raw,qntl)
                          }else
                          {
                            corrs_cns[[j]][,m]=rep(NA,19)
                          }
                          #corrs[[d]][[j]][,m]=monotonize(corrs[[d]][[j]][,m],nq)
                        }
                        else
                        {
                          corrs_cns[[j]][,m]=rep(NA,19)
                        }
                        
                      }
                      corrs[[d]][[j]]=matrix(NA,19,12)
                      corrs[[d]][[j]]=round(smoothen_nqx12(corrs_cns[[j]]),digits=1)
                    }
                    else if (cns_then_sm==0)
                    {
                      corrs_sm[[j]]=round(smoothen_nqx12(corrs_raw[[j]]),digits=1)
                      overlap=overlap_br[[j]]
                      corrs[[d]][[j]]=matrix(NA,19,12)
                      for (m in 1:12)
                      {
                        corr_sm=corrs_sm[[j]][,m]
                        if (length(corr_sm[!is.na(corr_sm)])>0)
                        {
                          if (m==1)
                          {vec=overlap[overlap$month %in% c(1,2,12) & !is.na(overlap$basis),'basis']
                          }else if (m==12)
                          {vec=overlap[overlap$month %in% c(11,12,1) & !is.na(overlap$basis),'basis']
                          }else
                          {vec=overlap[overlap$month %in% seq(m-1,m+1) & !is.na(overlap$basis),'basis']}
                          qntl=qntls_3m(vec,nb,nq,ampl_bin,thre_y)
                          corrs[[d]][[j]][,m]=check_neg_slope(corr_sm,qntl)
                          #corrs[[d]][[j]][,m]=monotonize(corrs[[d]][[j]][,m],nq)
                        }else
                        {
                          corrs[[d]][[j]][,m]=rep(NA,19)
                        }
                      }
                      
                    }
                  }else
                  {
                    corrs[[d]][[j]]=matrix(NA,nq,12)
                  }
                  colnames(corrs[[d]][[j]])=seq(1,12)
                  rownames(corrs[[d]][[j]])=seq(fq,lq,ampl_bin)
                  
                  setwd(outpath)
                  #WRITE IN THE STATION FILE
                  if (lyd%%4==0 & (lyd%%100!=0 | lyd%%400==0)) {endmonth[2]=29}
                  for (m in 1:12)
                  {
                    for (q in 1:nq)
                    {
                      if(q<10)
                      {
                        if(!is.na(corrs[[d]][[j]][q,m]))
                        {
                          row=paste(format(ser_id,width=12),
                                    format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                    format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                    format(ref_ser_id[[j]],width=12),
                                    format(m,width=12),
                                    format(as.numeric(rownames(corrs[[d]][[j]])[q]),width=12),
                                    format(as.numeric(bndrs_ds[[d]][[j]][q,m]*10),width=12),
                                    format((as.numeric(bndrs_ds[[d]][[j]][q+1,m]-0.1)*10),width=12),
                                    paste0(format(corrs[[d]][[j]][q,m]*10,width=12)),sep=',')
                        }else
                        {
                          if(abs(bndrs_ds[[d]][[j]][q,m])>100 & abs(bndrs_ds[[d]][[j]][q+1,m])>100)
                          {
                            row=paste(format(ser_id,width=12),
                                      format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                      format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                      format(ref_ser_id[[j]],width=12),
                                      format(m,width=12),
                                      format(as.numeric(rownames(corrs[[d]][[j]])[q]),width=12),
                                      format(-9999,width=12),
                                      format(-9999,width=12),
                                      paste0(format(-9999,width=12)),sep=',')
                          } else
                          {
                            row=paste(format(ser_id,width=12),
                                      format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                      format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                      format(ref_ser_id[[j]],width=12),
                                      format(m,width=12),
                                      format(as.numeric(rownames(corrs[[d]][[j]])[q]),width=12),
                                      format(as.numeric(bndrs_ds[[d]][[j]][q,m])*10,width=12),
                                      format((as.numeric(bndrs_ds[[d]][[j]][q+1,m])-0.1)*10,width=12),
                                      paste0(format(-9999,width=12)),sep=',')
                          }
                        }
                      }
                      if(q==10)
                      {
                        if(!is.na(corrs[[d]][[j]][q,m]))
                        {
                          row=paste(format(ser_id,width=12),
                                    format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                    format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                    format(ref_ser_id[[j]],width=12),
                                    format(m,width=12),
                                    format(as.numeric(rownames(corrs[[d]][[j]])[q]),width=12),
                                    format((as.numeric(bndrs_ds[[d]][[j]][q,m])+0.1)*10,width=12),
                                    format((as.numeric(bndrs_ds[[d]][[j]][q+1,m])-0.1)*10,width=12),
                                    paste0(format(corrs[[d]][[j]][q,m]*10,width=12)),sep=',')
                        }else
                        {
                          if(abs(bndrs_ds[[d]][[j]][q,m])>100 & abs(bndrs_ds[[d]][[j]][q+1,m])>100)
                          {
                            row=paste(format(ser_id,width=12),
                                      format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                      format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                      format(ref_ser_id[[j]],width=12),
                                      format(m,width=12),
                                      format(as.numeric(rownames(corrs[[d]][[j]])[q]),width=12),
                                      format(-9999,width=12),
                                      format(-9999,width=12),
                                      paste0(format(-9999,width=12)),sep=',')
                          } else
                          {
                            row=paste(format(ser_id,width=12),
                                      format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                      format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                      format(ref_ser_id[[j]],width=12),
                                      format(m,width=12),
                                      format(as.numeric(rownames(corrs[[d]][[j]])[q]),width=12),
                                      format((as.numeric(bndrs_ds[[d]][[j]][q,m])+0.1)*10,width=12),
                                      format((as.numeric(bndrs_ds[[d]][[j]][q+1,m])-0.1)*10,width=12),
                                      paste0(format(-9999,width=12)),sep=',')
                          }
                        }
                      }
                      if(q>10)
                      {
                        if(!is.na(corrs[[d]][[j]][q,m]))
                        {
                          row=paste(format(ser_id,width=12),
                                    format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                    format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                    format(ref_ser_id[[j]],width=12),
                                    format(m,width=12),
                                    format(as.numeric(rownames(corrs[[d]][[j]])[q]),width=12),
                                    format((as.numeric(bndrs_ds[[d]][[j]][q,m])+0.1)*10,width=12),
                                    format(as.numeric(bndrs_ds[[d]][[j]][q+1,m])*10,width=12),
                                    paste0(format(corrs[[d]][[j]][q,m]*10,width=12)),sep=',')
                        }else
                        {
                          if(abs(bndrs_ds[[d]][[j]][q,m])>100 & abs(bndrs_ds[[d]][[j]][q+1,m])>100)
                          {
                            row=paste(format(ser_id,width=12),
                                      format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                      format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                      format(ref_ser_id[[j]],width=12),
                                      format(m,width=12),
                                      format(as.numeric(rownames(corrs[[d]][[j]])[q]),width=12),
                                      format(-9999,width=12),
                                      format(-9999,width=12),
                                      paste0(format(-9999,width=12)),sep=',')
                          } else
                          {
                            row=paste(format(ser_id,width=12),
                                      format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                      format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                      format(ref_ser_id[[j]],width=12),
                                      format(m,width=12),
                                      format(as.numeric(rownames(corrs[[d]][[j]])[q]),width=12),
                                      format((as.numeric(bndrs_ds[[d]][[j]][q,m])+0.1)*10,width=12),
                                      format(as.numeric(bndrs_ds[[d]][[j]][q+1,m])*10,width=12),
                                      paste0(format(-9999,width=12)),sep=',')
                          }
                        }
                      }
                      write(row,output_a,append=TRUE)
                    }
                  }
                }
                setwd(inpath)
                
              }else #donating series is a subset of the base series
              {
                print(paste('Error in station',sta_id,', series',ser_id,'.Overlapping period between basis and donating!'))
              }
              
              
              
              ### CORRECTION OF THE SERIES
              homdon=donat[c('year','month','day','donat')]
              
              if (error_report==1)
              {
                tryCatch(  {homdon$hom=apply(donat,1,function(x) correct(x,bndrs_ds[[d]],nb,corrs[[d]]))},
                           warning = function(w) {write(c(sprintf('%-4s %7d %7d', nox, blend_sta_id, basis_ser_id),"warning in correction of blended series"), rprt, append=TRUE)},
                           error = function(e) {print('error in correction of blended series');
                             write(c(sprintf('%-4s %7d %7d', nox, blend_sta_id, basis_ser_id),"error in correction of blended series"), rprt, append=TRUE)})
              } else
              {
                homdon$hom=apply(donat,1,function(x) correct(x,bndrs_ds[[d]],nb,corrs[[d]]))
              }
              homdon$diff=homdon$hom-homdon$donat
              colnames(homdon)[5]='basis'
              basis=rbind(homdon[,c('year','month','day','basis')],basis)
            } else #end of if (length(refers)>3)
            {
              print(paste('Corrs calculation of',donat_ser_id,'from',fyd,'to',lyd,'has not enough references'))
              write(paste('Corrs calculation of',donat_ser_id,'from',fyd,'to',lyd,'has not enough references'),rprt,append=TRUE)
              homdon=donat[c('year','month','day','donat')]
              colnames(homdon)[4]='basis'
              basis=rbind(homdon[,c('year','month','day','basis')],basis)
            }
          }else
          {
            if (dim(donat)[1]>0)
            {
              fyd=donat[1,'year']
              lyd=donat[dim(donat)[1],'year']  
              print(paste('Corrs calculation of',donat_ser_id,'from',fyd,'to',lyd,'is impossible due to the short period'))
              write(paste('Corrs calculation of',donat_ser_id,'from',fyd,'to',lyd,'is impossible due to the short period'),rprt,append=TRUE)
              homdon=donat[c('year','month','day','donat')]
              colnames(homdon)[4]='basis'
              basis=rbind(homdon[,c('year','month','day','basis')],basis)
            }
          }
        }#end of the loop on the donating series
        final=basis 
      }else
      {
        print(paste('Corrs calculation of',basis_ser_id,'not possible, since the latest segment is too short'))
        write(paste('Corrs calculation of',basis_ser_id,'not possible, since the latest segment is too short'),rprt,append=TRUE)
        homdon=series
        colnames(homdon)[4]='basis'
        final=homdon
      }
      
      
      
      series_hom=final#basis[,c('year','month','day','basis')]
      colnames(series_hom)=c('year','month','day','hom')
      series_hom=complete(series_hom,fs,ls)
      series_hom$qc=8
      series_hom$qca=8
      series_hom[series_hom$hom<(-90),c('qc','qca')]=9
      setwd(outpath)
      hom4out=series_hom[,1:4]
      #hom4out=data.frame(ser_id=ser_id,ser_date=series_hom$year*10000+series_hom$month*100+series_hom$day,value=series_hom$hom*10,qc=series_hom$qc,
      #                   qca=series_hom$qca,qcm=9)
    }else
    {
      if (length(chp_y_old)>0)
      {
        print(paste('Breaks found in',ser_id,'are out of the range of the series.'))
        write(paste('Breaks found in',ser_id,'are out of the range of the series.'),rprt,append=TRUE)
      }
      else
      {
        print(paste('No breaks have been found in',ser_id))
        write(paste('No breaks have been found in',ser_id),rprt,append=TRUE)
      }
      series$qc=8
      series$qca=8
      series[series$series<(-90),c('qc','qca')]=9
      setwd(outpath)
      hom4out=series[,1:4]
      #hom4out=data.frame(ser_id=ser_id,ser_date=series$year*10000+series$month*100+series$day,value=series$series*10,qc=series$qc,qca=series$qca,qcm=9)
      
     
      
      
      
    }
    setwd(outpath)
    write.fwf(hom4out,output_s,colnames=FALSE,sep=' ',eol='\n',append=FALSE)
    setwd(inpath)
  }
  
}#end of if(length(files.ref)<=3)

