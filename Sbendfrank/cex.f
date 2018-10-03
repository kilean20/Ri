      subroutine cex(p,ga,gm,myorder)
c this routine computes t=exp(:f)
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'extalk.inc'
      include 'hmflag.inc'
      include 'buffer.inc'
      character*3 kynd
c
      dimension p(6),ga(monoms),gm(6,6)
c
      dimension ta(monoms)
      dimension tm(6,6),em(6,6),fm(6,6)
cryne 3AM 7/11/2002      dimension y(224)
      dimension y(monoms+15)

c--------------------------------------------
c     write(6,*)'inside routine cex'
c     write(6,*)'input 6-vector p:'
c     do i=1,6
c     write(6,*)p(i)
c     enddo
c     write(6,*)'input matrix:'
c     do i=1,6
c     do j=1,6
c     if(gm(i,j).ne.0.d0)write(69,*)i,j,gm(i,j)
c     enddo
c     enddo
c
c     write(6,*)'input polynomial:'
c     do i=1,209
c     if(ga(i).ne.0.d0)write(69,*)i,ga(i)
c     enddo
c     if(p(1).ne.12345)stop
c--------------------------------------------
c
c set up power and control indices
      power=p(1)
      nmapf=nint(p(2))
      nmapg=nint(p(3))
c
c get map and clear arrays
      if (nmapf.eq.0) call mapmap(ga,gm,fa,fm)
      if (nmapf.ge.1 .and. nmapf.le.5) then
      kynd='gtm'
      call strget(kynd,nmapf,fa,fm)
      endif
      call clear(ta,tm)
c
c perform calculation
c
c set up exponent
      call csmul(power,fa,fa)
c
c compute a scaling factor to bring exponent within range of a
c taylor expansion or GENMAP
      call matify(em,fa)
      call mnorm(em,res)
      kmax=1
      scale=.5d0
   10 continue
      test=res*scale
      if (test.lt..1d0) goto 20
      kmax=kmax+1
      scale=scale/2.d0
      go to 10
   20 continue
c
c select procedure
      itest=1
      do 30 i=28,monoms
      if (fa(i).ne.0.d0) itest=2
      if (itest.ne.1) go to 40
   30 continue
   40 continue
      if (itest.eq.1) go to 50
      if (itest.eq.2) go to 80
c
c procedure using taylor series
   50 continue
      write (12,*) 'exp(:f:) computed using taylor series'
c rescale em
      call smmult(scale,em,em)
c compute taylor series result fm=exp(scale*em)
      call exptay(em,fm)
c raise the result to the 2**kmax (= 1/scale) power
      do 60 i=1,kmax
      call mmult(fm,fm,tm)
      call matmat(tm,fm)
   60 continue
      goto 200
c
c procedure using genmap
   80 continue
      write(12,*) 'exp(:f:) computed using GENMAP'
c rescale fa
      call csmul(scale,fa,fa)
c setup and initialize for GENMAP routines
      iflag=1
      t=0.d0
      ns=50.d0
      h=.02d0
cryne 3AM 7/11/2002      ne=224
cryne 1 August 2004      ne=monoms+15
cryne 1 August 2004 Initialize:
cryne do 90 i=1,ne
      do 90 i=1,monoms+15
   90 y(i)=0.d0
      do 100 i=1,6
      j=7*i
  100 y(j)=1.d0
c call GENMAP routines
cryne there is a better way to do this; fix later.
c  y(7-42) = matrix
c  y(43-98) = f3
c  y(99-224) = f4
c  y(225-476) = f5
c  y(477-938) = f6

      if(myorder.eq.1)then
        ne=42    !36+6
      elseif(myorder.eq.2)then
        ne=98    !83+15
      elseif(myorder.eq.3)then
        ne=224    !209+15
      elseif(myorder.eq.4)then
        ne=476    !461+15
      elseif(myorder.eq.5)then
        ne=938    !923+15
      else
        write(6,*)'ERROR IN CEX: MYORDER IS NOT DEFINED'
        stop
      endif
c
      call adam11(h,ns,'start',t,y,ne)
      call putmap(y,fa,fm)
c
c raise the result to the 2**kmax (= 1/scale) power
      do 110 i=1,kmax
      call concat(fa,fm,fa,fm,ta,tm)
      call mapmap(ta,tm,fa,fm)
  110 continue
      go to 200
c
c decide where to put results
c
  200 continue
      if (nmapg.ge.1 .and. nmapg.le.5) then 
      kynd='stm'
      call strget(kynd,nmapg,ta,tm)      
      endif
c
      if (nmapg.eq.0) call mapmap(ta,tm,ga,gm)
c 
      if (nmapg.eq.-1) call mapmap(ta,tm,buf1a,buf1m)
      if (nmapg.eq.-2) call mapmap(ta,tm,buf2a,buf2m)
      if (nmapg.eq.-3) call mapmap(ta,tm,buf3a,buf3m)
      if (nmapg.eq.-4) call mapmap(ta,tm,buf4a,buf4m)
      if (nmapg.eq.-5) call mapmap(ta,tm,buf5a,buf5m)
c
      return
      end
