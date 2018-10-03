      subroutine pcmap(n1,n2,n3,n4,fa,fm)
c  routine to print m,f3,f4 and t,u.
c Written by D. Douglas ca 1982 and modified by Rob Ryne
c and Alex Dragt ca 1986
      use parallel
      use lieaparam, only : monoms
      include 'impli.inc'
      integer colme(6,0:6)
      include 'expon.inc'
      include 'pbkh.inc'
      include 'files.inc'
      dimension fa(monoms),fm(6,6)
      dimension t(monoms),u(monoms),u2(monoms)
      if(idproc.ne.0)return
c
c  test for matrix write
      if(n1.eq.0) goto 20
c
c  procedure for writing out matrix
c  write matrix at terminal
      if(n1.eq.1.or.n1.eq.3)then
        write(jof,13)
   13   format(/1h ,'matrix for map is :'/)
        write(jof,15)((fm(k,i),i=1,6),k=1,6)
   15   format(6(1x,1pe12.5))
        write(jof,*)
        write(jof,*)'nonzero matrix elements in full precision:'
        do k=1,6
        do i=1,6
        if(fm(k,i).ne.0.d0)write(jof,155)k,i,fm(k,i)
  155   format(i1,2x,i1,2x,1pg21.14)
        enddo
        enddo
        write(jof,*)
      endif
c  write matrix on file 12
      if(n1.eq.2.or.n1.eq.3)then
        write(jodf,13)
        write(jodf,15)((fm(k,i),i=1,6),k=1,6)
        write(jodf,*)
        write(jodf,*)'nonzero matrix elements in full precision:'
        do k=1,6
        do i=1,6
        if(fm(k,i).ne.0.d0)write(jodf,155)k,i,fm(k,i)
        enddo
        enddo
        write(jodf,*)
      endif
c
c  test for polynomial write
   20 continue
      if(n2.eq.0)goto 30
c
c  procedure for writing out polynomial
c  write polynomial at terminal
      if(n2.eq.1.or.n2.eq.3)then
        write(jof,22)
   22   format(/1h ,'nonzero elements in generating polynomial are :'/)
        do 25 i=1,monoms
        if(fa(i).eq.0.0d0) goto 25
        write(jof,27)i,(expon(j,i),j=1,6),fa(i)
cryne 6/21/2002   27   format(2x,'f(',i3,')=f( ',3(2i1,1x),')=',d21.14)
   27 format(2x,'f(',i3,')=f( ',3(2i1,1x),')=',1pg21.14)
   25   continue
      endif
c  write polynomial on file 12
      if(n2.eq.2.or.n2.eq.3)then
        write(jodf,22)
        do 26 i=1,monoms
        if(fa(i).eq.0.0d0) goto 26
        write(jodf,27)i,(expon(j,i),j=1,6),fa(i)
   26   continue
      endif
c
c  prepare for higher order matrix write if required
   30 continue
cryne write(6,*)'inside pcmap'
      if(n3.gt.0.or.n4.gt.0)then
cryne   write(6,*)'calling brkts'
        call brkts(fa)
cryne   write(6,*)'returned from brkts'
      endif
c
c  test for t matrix write
      if(n3.eq.0) goto 40
c
c  procedure for writing t matrix
c  write out heading
      if(n3.eq.1.or.n3.eq.3)write(jof,32)
      if(n3.eq.2.or.n3.eq.3)write(jodf,32)
   32   format(/1h ,'nonzero elements in second order matrix are :'/)
c  write out contents
        do 35 i=1,6
        call xform(pbh(1,i),2,fm,i-1,t)
        do 36 n=7,27
        if(t(n).eq.0.0d0) goto 36
        if(n3.eq.1.or.n3.eq.3)
     &  write(jof,38) i,n,i,(expon(j,n),j=1,6),t(n)
        if(n3.eq.2.or.n3.eq.3)
     &  write(jodf,38) i,n,i,(expon(j,n),j=1,6),t(n)
cryne6/21/02 38format(2x,'t',i1,'(',i3,')','=t',i1,'( ',3(2i1,1x),')=',d21.14)
   38 format(2x,'t',i1,'(',i3,')','=t',i1,'( ',3(2i1,1x),')=',1pg21.14)
   36   continue
   35   continue
c
 
c  test for u matrix write
   40 continue
      if(n4.eq.0) goto 50
c
c  procedure for writing u matrix
c  write out heading
      if(n4.eq.1.or.n4.eq.3)write(jof,42)
      if(n4.eq.2.or.n4.eq.3)write(jodf,42)
   42 format(/1h ,'nonzero elements in third order matrix are :'/)
c  write out contents
         do  44 i=1,6
         call xform(pbh(1,i),3,fm,i-1,u)
         call xform(pbh(1,i+6),3,fm,1,u2)
         do  45 n=28,83
         u(n)=u(n)+u2(n)/2.d0
         if(u(n).eq.0.0d0) goto 45
         if(n4.eq.1.or.n4.eq.3)
     &   write(jof,46) i,n,i,(expon(j,n),j=1,6),u(n)
         if(n4.eq.2.or.n4.eq.3)
     &   write(jodf,46) i,n,i,(expon(j,n),j=1,6),u(n)
cryne 6/21/02 46format(2x,'u',i1,'(',i3,')','=u',i1,'( ',3(2i1,1x),')=',d21.14)
   46 format(2x,'u',i1,'(',i3,')','=u',i1,'( ',3(2i1,1x),')=',1pg21.14)
   45    continue
   44    continue
c
c  procedure if all n's are zero or are faulty
   50 continue
c
      return
      end
