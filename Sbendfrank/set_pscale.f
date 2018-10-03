      subroutine set_pscale_mc(h,mh)
c scale a map to use mc as the scale momentum (instead of p0)
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'expon.inc'
      double precision h(monoms),mh(6,6)
c     write(6,*)'inside set_pscale'
c
cryne Jan 30,2003
c this routine does 2 things:
c (1) convert from static to dynamic units (if needed),
c (2) convert variables 5 and 6 to take into account the fact that the
c     scale varables omega and l are independent and do not need to
c     satisfy omega*l/c=1 (as is assumed in the MaryLie library routines)
c
      clite=299792458.d0
c dynamic units:
calsowrong      p2ovrp1=1.d0
      p2ovrp1=p0sc/(gamma*beta*pmass/clite)
c     write(6,*)'initializing p2ovrp1 to ',p2ovrp1
c
      if(lflagdynu)then
        p2ovrp1=1./(beta*gamma)
c       write(6,*)'lflagdynu is true; resetting p2ovrp1 to ',p2ovrp1
      endif
cwrong:      if(lflagdynu)p2ovrp1=(pmass/clite)/p0sc
cryne the reason this was wrong is that the MaryLie library
cryne routines use brho, the p0sc.In other words, everything
cryne is computed relative to the present momentum.
      p1ovrp2=1.d0/p2ovrp1
c omega*l/c:
      scal5=1.d0
      scal6=1.d0
      wlbyc=omegascl*sl/clite
c     write(6,*)'initializing wlbyc to ',wlbyc
      if(wlbyc.ne.1.d0)then
c       write(6,*)'scal5 being set to wlbyc, scal6 set to 1/wlbyc'
        scal5=wlbyc
        scal6=1.d0/wlbyc
      endif
c
      mh(1,2)=mh(1,2)*p2ovrp1
      mh(1,4)=mh(1,4)*p2ovrp1
      mh(1,5)=mh(1,5)/wlbyc
      mh(1,6)=mh(1,6)*p2ovrp1*wlbyc
c
      mh(2,1)=mh(2,1)*p1ovrp2
      mh(2,3)=mh(2,3)*p1ovrp2
      mh(2,5)=mh(2,5)*p1ovrp2/wlbyc
      mh(2,6)=mh(2,6)*wlbyc
c
      mh(3,2)=mh(3,2)*p2ovrp1
      mh(3,4)=mh(3,4)*p2ovrp1
      mh(3,5)=mh(3,5)/wlbyc
      mh(3,6)=mh(3,6)*p2ovrp1*wlbyc
c
      mh(4,1)=mh(4,1)*p1ovrp2
      mh(4,3)=mh(4,3)*p1ovrp2
      mh(4,5)=mh(4,5)*p1ovrp2/wlbyc
      mh(4,6)=mh(4,6)*wlbyc
c
      mh(5,1)=mh(5,1)*wlbyc
      mh(5,2)=mh(5,2)*p2ovrp1*wlbyc
      mh(5,3)=mh(5,3)*wlbyc
      mh(5,4)=mh(5,4)*p2ovrp1*wlbyc
      mh(5,6)=mh(5,6)*p2ovrp1*wlbyc**2
c
      mh(6,1)=mh(6,1)*p1ovrp2/wlbyc
      mh(6,2)=mh(6,2)/wlbyc
      mh(6,3)=mh(6,3)*p1ovrp2/wlbyc
      mh(6,4)=mh(6,4)/wlbyc
      mh(6,5)=mh(6,5)*p1ovrp2/wlbyc**2
c
c now change the units of the nonlinear part of the map:
cryne 12/21/2004 changed "i=28,209" to "i=28,monoms"  !!!!!!!!!!!!!!
c     write(6,*)'modified set_pscale_mc do loop to use monoms=',monoms
      do 200 i=28,monoms
ccc   h(i)=h(i)*(p2ovrp1**nxp246(i))*(x2ovrx1**nxp135(i))
!!!   h(i)=h(i)*(p2ovrp1**nxp24(i))*(scal5**nxp5(i))*(scal6**nxp6(i))
cryne 12/21/2004 note that this is a minus 1 only in nxp13 and nxp24
cryne note that we are assuming that the scale length has been done already
!     h(i)=h(i)*(p2ovrp1**(nxp24(i)+nxp6(i)))*(scal5**(nxp5(i)-nxp6(i)))
      h(i)=h(i)*(p2ovrp1**(nxp24(i)+nxp6(i)))*(scal5**(nxp6(i)-nxp5(i)))
  200 continue
      return
      end
