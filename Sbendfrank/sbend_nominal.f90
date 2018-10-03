c
c====================================================================
c
      subroutine get_sbendmap(pp,xmh,h,jslice,jsltot,slfrac,ihalf)
      use beamdata !contains brho,gamma,gamm1,beta,achg,pmass,bfreq,bcurr,c,
c                  !scaling constants (sl,p0sc,ts,omegascl,freqscl),
c                  !and logicals lflagmagu,lflagdynu
      use lieaparam, only : monoms
      use parallel, only : idproc  !for diagnostic printouts from proc 0
      include 'mpif.h'
      include 'impli.inc'
      include 'map.inc' !contains total transfer map (not used), reftraj, arclen
      include 'previous.inc'  !contains prevlen (needed for strength of sc kick)
      include 'pie.inc'  !contains pi,pi180,twopi
      include 'parset.inc'  !can be used to specify multipole coeffs of sbend
      dimension pp(21)
      dimension xmh(6,6),h(monoms)
c local variables:
      dimension xmht1(6,6),ht1(monoms)
      dimension ptmp(6),coefpst(6),xmsk(6)
c
c Most of the following code from routine lmnt in afro.f
c
c 'sbend   ': MAD sector bend
c  similar to MaryLie cfbd, but with arbitrary entrance/exit angles
c
c  The map for this element is treated as the following:
c  prot*cgfrngg*cgbend*cgfrngg*prot
c where:
c  prot is a pole face rotation
c  cgfrngg is a fringe field (analgous to that used in cfbd)
c  cgbend is a gbend but w/ multipole coefficients
c  note that, in analogy with the fact that gbend=hpf1*nbend*hpf2,
c  in this case we set cgbend=hpf1*(cfbd w/ no fringe fields)*hpf2,
c  i.e.                cgbend=hpf1*cfbend*hpf2
c
      bndang=pp(1)
      rho=brho/pp(2)
      psideg=pp(3)
      phideg=pp(4)
      lfrn=nint(pp(5))
      mtfrn=nint(pp(6))
      gap1=pp(7)
      gap2=pp(8)
      fint1=pp(9)
      fint2=pp(10)
      iopt=nint(pp(11))
      ipset=nint(pp(12))
c
c     if(idproc.eq.0)then
c      write(6,*)'bndang=',bndang
c      write(6,*)'rho=',rho
c      write(6,*)'psideg=',psideg
c      write(6,*)'phideg=',phideg
c      write(6,*)'jslice,nslices=',jslice,jsltot 
c      write(6,*)'jslfrac,ihalf=',slfrac,ihalf
c      write(6,*)'idproc=',idproc
c     endif
c
      do j=1,6
      coefpst(j)=0.d0
      enddo
      if( (ipset.gt.0).and.(ipset.le.maxpst))then
        do j=1,6
        coefpst(j)=pst(j,ipset)
        enddo
      else
        do j=1,6
        if(pp(j+12).ne.0.d0)coefpst(j)=pp(j+12)
        enddo
      endif
c
      tiltang=pp(19)
      if(tiltang.ne.0. .and. idproc.eq.0)then
        write(6,*)'WARNING(sbend): tilt keyword being tested'
      endif
      myorder=nint(pp(20))
      if(idproc.eq.0)then
        if(myorder.ne.5)write(6,*)'(sbend) order = ',myorder
      endif
      azero=0.d0
c----------------
c compute gbend
c     call cgbend(rho,bndang,psideg,phideg,iopt,coefpst,h,xmh)
      aldeg=bndang-psideg-phideg
      ptmp(1)=aldeg+psideg+phideg    ! = bndang
      if(slfrac.ne.1.0)ptmp(1)=ptmp(1)*slfrac
      ptmp(2)=brho/rho
      ptmp(3)=0.
      ptmp(4)=0.
      ptmp(5)=iopt
cryne 1 August 2004 ptmp(6)=0.    !not used
      ptmp(6)=myorder
c
c     ptmp(1)=psideg
c     call cfbend(ptmp,coefpst,ht2,xmht2)
c     ptmp(1)=aldeg
c     call cfbend(ptmp,coefpst,ht3,xmht3)
c     call concat(ht2,xmht2,ht3,xmht3,ht3,xmht3)
c     ptmp(1)=phideg
c     call cfbend(ptmp,coefpst,ht2,xmht2)
c     call concat(ht2,xmht2,ht3,xmht3,ht3,xmht3)
c     if(idproc.eq.0)then
c     write(6,*)'call cfbend:aldeg,slfrac,ptmp(1)=',aldeg,slfrac,ptmp(1)
c     endif
c***
cryne 12/22/2004  mods per Marco Venturini bypass the slow cfbend routines
c if there are no multipole coefficients (note that, at the moment, the
c cfbend routine is only third order)
      scpset=0.d0
      do ii=1,6
        scpset=scpset+coefpst(ii)
      enddo
      if(scpset.eq.0.d0)then    ! use nbend if field is pure dipole
c        write(6,*)'computing pure sbend map (no multipole coefficients)'
        slang=bndang
        if(slfrac.ne.1.0)slang=slang*slfrac
        call nbend(rho,slang,h,xmh)
      else
c        write(6,*)'computing combined fxn sbend (w/ multipole coeffs)'
        call cfbend(ptmp,coefpst,h,xmh)
      endif
c***
c----------------
c
      if(jslice.eq.1.and.(ihalf.eq.0.or.ihalf.eq.1))then
c put on the leading gbend
      call gbend(rho,azero,psideg,azero,ht1,xmht1)
      call concat(ht1,xmht1,h,xmh,h,xmh)
c put on the leading fringe field
      call cgfrngg(1,psideg,rho,lfrn,gap1,fint1,iopt,coefpst,ht1,xmht1)
      call concat(ht1,xmht1,h,xmh,h,xmh)
c put on leading prot
      call prot(psideg,1,ht1,xmht1)
      call concat(ht1,xmht1,h,xmh,h,xmh)
c put on leading arot
      call arot(tiltang,ht1,xmht1)
      call concat(ht1,xmht1,h,xmh,h,xmh)
      endif
c
      if(jslice.eq.jsltot.and.(ihalf.eq.0.or.ihalf.eq.2))then
c put on the trailing gbend
      call gbend(rho,azero,azero,phideg,ht1,xmht1)
      call concat(h,xmh,ht1,xmht1,h,xmh)
c put on trailing fringe field
      call cgfrngg(2,phideg,rho,mtfrn,gap2,fint2,iopt,coefpst,ht1,xmht1)
      call concat(h,xmh,ht1,xmht1,h,xmh)
c put on trailing prot
      call prot(phideg,2,ht1,xmht1)
      call concat(h,xmh,ht1,xmht1,h,xmh)
c put on trailing arot
      tiltang=-tiltang
      call arot(tiltang,ht1,xmht1)
      call concat(h,xmh,ht1,xmht1,h,xmh)
      endif
cryne April 21, 2006:
cryne This is a bit confusing, but is the best I can come up with.
cryne The only other choice is to let myorder equal 2-6, but that
cryne is most confusing because "2" would mean "linear calculation"
cryne So here it is:
cryne If myorder=1, do linear calculation (omit f3-f6)
cryne If myorder=2, do quadratic calculation (omit f4-f6)
cryne If myorder=3, do cubic     calculation (omit f5-f6)
cryne If myorder=4, do quartic   calculation (omit f6-f6)
cryne If myorder=5, do quintic   calculation (do nothing, i.e. keep all)
cryne In other words, myorder is the order of the nonlinear effect,
cryne not the order of the monomial being masked/not masked.
      if(myorder.ge.1 .and. myorder.le.4)then
        xmsk(:)=1.d0
        xmsk(myorder+2:6)=0.d0
        call mask(xmsk,h,xmh)
      endif
      call set_pscale_mc(h,xmh)

      pp1=pp(1)
      if(slfrac.ne.1.0)pp1=pp1*slfrac
      arclen=arclen+rho*pp1*(pi/180.)
      reftraj(5)=reftraj(5)+pp1*rho*(pi/180.)/(beta*c)*(omegascl)
      prevlen=rho*pp1*(pi/180.)

      return
      end
c====================================================================
c
      subroutine myexit
!      use parallel
!      call end_parallel
      stop
      end
c====================================================================
c
      block data misc
c MaryLie and ML/I require that certain variables in common blocks
c be initialized. Here is a subset of "block data misc" from MLI
c
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
c--------------------
      include 'files.inc'
      include 'maxcat.inc'
      include 'lims.inc'
c--------------------
      data ordcat/6/,topcat/monoms/
      data ordlib/6/,toplib/monoms/
c--------------------
      data bottom/0,1, 7,28, 84,210,462, 924,1716,3003,5005, 8008,12376/
      data top   /0,6,27,83,209,461,923,1715,3002,5004,8007,12375,18563/
c--------------------
      data lf,icf,jfcf,mpi,mpo,jif,jof,jodf,ibrief,iquiet/              &
     &     11, 13,  14, 15, 16,  5,  6,  12,     2,     0/
      end
c
