c     To be compiled from shell with 'R CMD SHLIB myeval.f'  
      subroutine myeval(string,a,b,n1,n2,mylik)
    
      integer n1,n2,string(1:n1),control,indx(1:n2),pos(1:n2,1:n1)
      double precision a(1:n2,1:n2),b(1:n1),p,media(1:n2),sigma(1:n2)
      double precision subpop(1:n2,1:n1),stdobs(1:n1),sigobs(1:n1)
      double precision lolika,lolikb,mylik

      p=a(1,1)
      control=0
      do l=1,n2
        do i=1,n1
          if (string(i).eq.l) control=control+1
        enddo
        if (control.le.1) then
          mylik=1000
          goto 100
        end if
        control=0
      enddo

      indx(1:n2)=0
      do j=1,n2
        do i=1,n1
          if(string(i).eq.j) then
            indx(j)=indx(j)+1
            subpop(j,indx(j))=b(i)
            pos(j,indx(j))=i
          end if
        enddo
        if (indx(j).gt.0) then
          media(j)=sum(subpop(j,1:indx(j)))/indx(j)
          sigma(j)=sqrt(sum((subpop(j,1:indx(j))-media(j))**2)
     &/(indx(j)-1))
        end if
      enddo

      do j=1,n2
        if(indx(j).gt.0) then
          do i=1,indx(j)
            stdobs(pos(j,i))=((subpop(j,i)-media(j))**2)/
     &(2*sigma(j)*sigma(j))
            sigobs(pos(j,i))=log(sigma(j)*sqrt(2*3.141592653589793))
          enddo
        end if
      enddo
      lolika=sum(stdobs)
      lolikb=sum(sigobs)

      mylik=-(-lolika-lolikb+((n2-1)*log(1-p))+(((n1-1)-(n2-1)-
     &(indx(n2)-1))*(log(p))))
      if(mylik.ne.mylik) mylik=1000

100   end
