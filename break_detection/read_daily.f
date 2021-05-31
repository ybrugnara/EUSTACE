c ---------------------------------------------------------------------
c     Fortran subroutine RDDAILY version 1.1 (called by Breaks_daily.R)
c     [compile from terminal with 'R CMD SHLIB read_daily.f']
c ---------------------------------------------------------------------
c     Read daily data in KNMI format and calculate yearly means
c     Series with less than 'nmin' years are excluded
c     Series more than 'dmax' km away from the candidate are excluded
c     Series uncorrelated with the candidate (r<rmin) are excluded
c     The correlation is calculated from first differences of annual 
c     means
c
c ---------------------------------------------------------------------
c     INPUT 
c ---------------------------------------------------------------------
c     - files: filename of the list of files to read
c     - nf: number of files
c     - yr1: first year to analyse
c     - yr2: last year to analyse
c     - ddmax: maximum number of missing days in a month
c     - dismax: maximum distance of the reference series (km)
c     - nmin: minimum number of available months
c     - rmin: minimum correlation
c     - nmax: maximum number of reference series
c     - overlap: minimum fraction of overlapping months
c     - ic: index of the candidate series
c     - inpath: path of the input directory
c     - outpath: path of the output directory
c     reads the list of series to analyse from the file <outpath>/list
c
c ---------------------------------------------------------------------
c     OUTPUT
c ---------------------------------------------------------------------
c     - x: matrix of monthly data
c     - rx: vector of correlations
c     - lat: vector of latitudes
c     - lon: vector of longitudes
c     - ele: vector of elevations
c     - sname: name of the candidate station
c     - j: number of available reference series
c     writes the list of available reference series in the file
c     <outpath>/list_ref
c ---------------------------------------------------------------------


      subroutine rddaily(files,nf,x,rx,lat,lon,ele,yr1,yr2,ddmax,
     &dismax,nmin,rmin,nmax,overlap,ic,sname,j,inpath,outpath)

      integer ddmax,nf,yr1,yr2,d,ic,nmin,j,n,k,nmax
      integer year,pyear,month,pmonth,ndays(yr1:yr2,1:12)
      double precision x(1:12,yr1:yr2,1:nf),xmiss,rmin,rx(1:(nf-1))
      double precision fdiff((yr1+1):yr2),fdiffref((yr1+1):yr2)
      double precision xtmp(1:366),tmp,latmp,lotmp,eltmp,overlap
      double precision fy(1:(yr2-yr1+1)),dist,dismax,dmaxk
      double precision ytmp(yr1:yr2,1:12),rdist(1:nf)
      double precision lat(1:nf),lon(1:nf),ele(1:nf)
      character*200 inpath,outpath
      character*50 sname,rname,fname,files,rnames(1:nf)
      data xmiss/-99.99/d,n,k/0,0,0/
      
      do iy=yr1,yr2
        do im=1,12
          if (im.eq.4.or.im.eq.6.or.im.eq.9.or.im.eq.11) then
            ndays(iy,im)=30
          elseif (im.eq.2) then
            if (iy/4.eq.abs(iy/4).and.iy.ne.1900) then
              ndays(iy,im)=29
            else
              ndays(iy,im)=28
            end if
          else
            ndays(iy,im)=31
          end if
        enddo
      enddo

      do iip=1,199
        if(inpath(iip:(iip+1)).eq.'//') goto 51
      enddo
51    continue
      do iop=1,199
        if(outpath(iop:(iop+1)).eq.'//') goto 52
      enddo
52    continue
      open(8,file=outpath(1:iop)//files)


c Read candidate series
      do im=1,12 
        ytmp(yr1:yr2,im)=xmiss
      enddo
     
      do i=1,ic-1
        read(8,*)
      enddo
      read(8,*) fname
      open(9,file=inpath(1:iip)//fname,status='old')
      read(9,*)
      read(9,'(14x,f11.6,2x,f12.6,2x,f9.1)') lat(j),lon(j),ele(j)
      read(9,'(14x,a50)') sname
      read(9,*)
      read(9,*)

      do 100
        read(9,'(i5,i3)',end=299) pyear,pmonth     
        if(pyear.ge.yr1) goto 198        
100   enddo
198   backspace(9)
      yr1=pyear

      do 200

        read(9,'(i5,i3,3x,f10.2)',end=299) year,month,tmp
        if(year.gt.yr2) goto 299

        if(year.ne.pyear.or.month.ne.pmonth) then
          if(d.ge.(ndays(pyear,pmonth)-ddmax)) then
            ytmp(pyear,pmonth)=sum(xtmp(1:d))/d
            n=n+1
          end if
          d=0
          pyear=year
          pmonth=month
        end if
        if(tmp.gt.-90) then
          d=d+1
          xtmp(d)=tmp
        end if

200   enddo

299   if(d.ge.(ndays(pyear,pmonth)-ddmax)) then
        ytmp(pyear,pmonth)=sum(xtmp(1:d))/d
        n=n+1
      end if
      yr2=pyear
      
      close(9)


      write(*,'(a60)') '*****************************************'//
     &'*******************'
      write(*,'(a60)') '*****************************************'//
     &'*******************'
      write(*,'(a60)') '*****************************************'//
     &'*******************'
      write(*,*)
      write(*,'(a60)') '-----------------------------------------'//
     &'-------------------'
      write(*,'(a16)') 'CANDIDATE SERIES'
      write(*,'(20x,5a8)') 'Lat','Long','Elev','Start','End'
      write(*,'(a20,2f8.2,3i8)') sname,lat(j),lon(j),int(ele(j)),yr1
     &,yr2
      write(*,'(a60)') '-----------------------------------------'//
     &'-------------------'


      if(n.lt.nmin) then
        write(*,*) 'Skipped (too short)'       
        goto 999
      end if
      ncand=n

      do iy=yr1,yr2
        x(1:12,iy,j)=ytmp(iy,1:12)
        call yearmean(ytmp(iy,1:12),xmiss,fy(iy-yr1+1))
      enddo
      call frstdiff(fy(1:(yr2-yr1+1)),xmiss,yr2-yr1+1,fdiff)


c Read reference series
      write(*,'(a26)') 'AVAILABLE REFERENCE SERIES'
      write(*,'(20x,5a8)') 'Lat','Long','Elev','Dist','Corr'
      open(10,file=outpath(1:iop)//'list_ref')     
302   i=0
      k=k+1
      dmaxk=k*dismax/10
      if(k.eq.10) ncand=nmin
      rewind(8)

      do 300
        i=i+1
        if(i.eq.ic) then 
          i=i+1
          read(8,*) fname
        end if
        read(8,*,end=399) fname
        open(9,file=inpath(1:iip)//fname,status='old')
        read(9,*)
        read(9,'(14x,f11.6,2x,f12.6,2x,f9.1)') latmp,lotmp,eltmp
        call distance(lotmp,lon(1),latmp,lat(1),dist)

        if(dist.gt.dmaxk) then
          goto 301
        else
          read(9,'(14x,a50)') rname
          read(9,*)
          read(9,*)
          
          do 400
            read(9,'(i5,i3)',end=301) pyear,pmonth
            if(pyear.ge.yr1) goto 498
400       enddo
498       backspace(9)
          d=0
          n=0
          do im=1,12 
            ytmp(yr1:yr2,im)=xmiss
          enddo

          do 500
            read(9,'(i5,i3,3x,f10.2)',end=599) year,month,tmp
            if(year.gt.yr2) goto 599
            if(year.ne.pyear.or.month.ne.pmonth) then
              if(d.ge.(ndays(pyear,pmonth)-ddmax)) then
                ytmp(pyear,pmonth)=sum(xtmp(1:d))/d
                n=n+1
              end if
              d=0
              pyear=year
              pmonth=month
            end if
            if(tmp.gt.-90) then
              d=d+1
              xtmp(d)=tmp
            end if
500       enddo

        end if

599     if(d.ge.(ndays(pyear,pmonth)-ddmax)) then
          ytmp(pyear,pmonth)=sum(xtmp(1:d))/d
          n=n+1
        end if

        if(n.ge.nmin.and.n.ge.int(ncand*overlap)) then

          do iy=yr1,yr2           
            call yearmean(ytmp(iy,1:12),xmiss,fy(iy-yr1+1))
          enddo
          call frstdiff(fy(1:(yr2-yr1+1)),xmiss,yr2-yr1+1,fdiffref)
          call pcorr(fdiff,fdiffref,yr2-yr1,xmiss,10,rx(j))

          if(rx(j).ge.rmin) then
            j=j+1
            do iy=yr1,yr2
              x(1:12,iy,j)=ytmp(iy,1:12)
            enddo
            lat(j)=latmp
            lon(j)=lotmp
            ele(j)=eltmp
            rnames(j-1)=rname
            rdist(j-1)=dist
            write(10,*) fname
          end if

        end if

301   close(9)
      if(j.eq.nf-1) goto 399
         
300   enddo

399   j=j-1
      if(j.lt.nmax.and.k.lt.10) then
        j=1
        rewind(10)
        goto 302
      end if
      if((j.ge.nmax).or.(j.lt.nmax.and.j.ge.3.and.k.eq.10)) then
        do ij=1,j
          write(*,'(a20,2f8.2,i8,f8.1,f8.2)') rnames(ij),lat(ij+1),
     &lon(ij+1),int(ele(ij+1)),rdist(ij),rx(ij)
        enddo
        goto 998
      end if        
      if(j.lt.3.and.k.eq.10) then
        write(*,*) '!!!INSUFFICIENT REFERENCE SERIES AVAILABLE!!!'
      end if
     
998   close(10)

999   write(*,*)
      close(8)
    
   
      end


c ---------------------------------------------------------------------
c Function that calculates the distance between two geographical points
c ---------------------------------------------------------------------
      subroutine distance(xlonA,xlonB,xlatA,xlatB,xdist)

      double precision xlonA,xlonB,xlatA,xlatB,xdist,pi,R,xA,xB,yA,yB
      double precision xfi

      if(xlonA.eq.xlonB.and.xlatA.eq.xlatB) then
        xdist=0
      else
        pi=acos(-1.d0)
        R=6378.d0
        xA=pi*xlonA/180.d0
        xB=pi*xlonB/180.d0
        yA=pi*xlatA/180.d0
        yB=pi*xlatB/180.d0
        xfi=abs(xA-xB)
        xdist=R*acos(sin(yA)*sin(yB)+cos(yA)*cos(yB)*cos(xfi))
      end if

      return
      end


c ---------------------------------------------------------------------
c Function that calculates the Pearson correlation coefficient taking
c into account missing values
c ---------------------------------------------------------------------
      subroutine pcorr(x,y,n,miss,nxymin,r)

      integer n,nxy,nxymin
      double precision x(1:n),y(1:n),r,r1,r2,r3,xavg,yavg,miss
      
      r1=0.0
      r2=0.0
      r3=0.0
      xavg=0.0
      yavg=0.0
      nxy=0
      do i=1,n
        if((x(i).ne.miss).and.(y(i).ne.miss)) then
          xavg=xavg+x(i)
          yavg=yavg+y(i)
          nxy=nxy+1
        end if
      enddo
      if(nxy.lt.nxymin) then
        r=0.0
        goto 98
      end if
      xavg=xavg/float(nxy)
      yavg=yavg/float(nxy)
      
      do i=1,n
        if((x(i).ne.miss).and.(y(i).ne.miss)) then
          r1=r1+(x(i)-xavg)*(y(i)-yavg)
          r2=r2+(x(i)-xavg)**2
          r3=r3+(y(i)-yavg)**2
        end if
      enddo
      r=r1/(sqrt(r2)*sqrt(r3))

98    return
      end


c ---------------------------------------------------------------------
c Function that calculates the first differences
c ---------------------------------------------------------------------
      subroutine frstdiff(x,miss,n,xd)

      integer n
      double precision x(1:n),xd(1:(n-1)),miss

      do i=2,n
        if(x(i).gt.-90.and.x(i-1).gt.-90) then
          xd(i-1)=x(i)-x(i-1)
        else
          xd(i-1)=miss
        end if
      enddo

      return
      end


c ---------------------------------------------------------------------
c Function that calculates annual averages from monthly data
c ---------------------------------------------------------------------
      subroutine yearmean(x,miss,y)

      double precision x(1:12),y,miss

      y=0.0
      do i=1,12
        if(x(i).gt.-90) then
          y=y+x(i)
        else
          y=miss
          goto 99
        end if
      enddo
      y=y/12

99    return
      end
