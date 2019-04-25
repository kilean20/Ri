!----------------------------------------------------------------
! (c) Copyright, 2003 by the Regents of the University of California.
! Outputclass: Output class in I/O module of CONTROL layer. 
! Version: 2.0
! Author: Ji Qiang, LBNL, 7/24/03
! Description: This class defines functions to print out the charged
!              particle beam information in the accelerator.
! Comments: We have added 3 more attributes to the particle:
!           x,px,y,py,t,pt,charge/mass,charge weight,id
!----------------------------------------------------------------
      module Outputclass
!        use Timer_class
        use BeamBunchclass
        use Pgrid2dclass
        use Fldmgerclass
        use PhysConstclass
        use Multipoleclass
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Kilean
        private :: get_free_unit,num2str_int,sort
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Kilean
      contains
        ! calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance.
        subroutine diagnostic1_Output(z,this,nchrg,nptlist)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), intent(in) :: this
        integer, intent(in) :: nchrg
        integer*8, dimension(nchrg), intent(in) :: nptlist
        integer :: innp
        integer*8 :: nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: gam,energy,bet,gambet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(27) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax
!        double precision :: alphax,alphay,alphaz
!       for the calculation of bunched beam current.
!        double precision, dimension(3) :: localaa, aa, localbb, bb 
!        double precision :: pi,tmp5,tmp55
!        integer :: ntmp5,nf
        integer :: npctmin,npctmax

        call starttime_Timer(t0)

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        innp = this%Nptlocal
        nptot = this%Npt
        !print*,"pts1: ",this%Pts1(1,1),this%Pts1(2,1),my_rank,innp

        den1 = 1.0/dble(nptot)
        den2 = den1*den1
        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        z0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0

!        localaa = 0.0
!        localbb = 0.0
!        aa = 0.0
!        bb = 0.0
!        pi = 2*asin(1.0)

        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            localmax(i) = abs(this%Pts1(i,1))
          enddo
          lcrmax = this%Pts1(1,1)**2+this%Pts1(3,1)**2
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
          x0lc3 = x0lc3 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)
          x0lc4 = x0lc4 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)*&
                  this%Pts1(1,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*this%Pts1(2,i)
          px0lc = px0lc + this%Pts1(2,i)
          sqsum2local = sqsum2local + this%Pts1(2,i)*this%Pts1(2,i)
          px0lc3 = px0lc3 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)
          px0lc4 = px0lc4 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)*&
                   this%Pts1(2,i)
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
          y0lc3 = y0lc3 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)
          y0lc4 = y0lc4 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)*&
                  this%Pts1(3,i)
          ypylocal = ypylocal + this%Pts1(3,i)*this%Pts1(4,i)
          py0lc = py0lc + this%Pts1(4,i)
          sqsum4local = sqsum4local + this%Pts1(4,i)*this%Pts1(4,i)
          py0lc3 = py0lc3 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)
          py0lc4 = py0lc4 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)*&
                   this%Pts1(4,i)
          z0lc = z0lc + this%Pts1(5,i)
          sqsum5local = sqsum5local + this%Pts1(5,i)*this%Pts1(5,i)
          z0lc3 = z0lc3 + abs(this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i))
          z0lc4 = z0lc4 + this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i)*&
                          this%Pts1(5,i)

          !for the calculation of bunched beam current
!          ntmp5 = this%Pts1(5,i)/pi
!          tmp5 = this%Pts1(5,i) - ntmp5*pi
!          tmp55 = tmp5 - mod(ntmp5,2)*pi
!          do nf = 1, 3
!            localaa(nf) = localaa(nf) + cos(nf*tmp55)
!            localbb(nf) = localbb(nf) + sin(nf*tmp55)
!          enddo

          zpzlocal = zpzlocal + this%Pts1(5,i)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc3 = pz0lc3 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc4 = pz0lc4 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)*&
                            this%Pts1(6,i)
          do j = 1, 6
            if(localmax(j).lt.abs(this%Pts1(j,i))) then
               localmax(j) = abs(this%Pts1(j,i))
            endif
          enddo
          if(lcrmax.lt.(this%Pts1(1,i)**2+this%Pts1(3,i)**2)) then
            lcrmax = this%Pts1(1,i)**2 + this%Pts1(3,i)**2
          endif
        enddo

!        do i = 1, 6
!          localmax(i) = maxval(abs(this%Pts1(i,1:innp)))
!        enddo

        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4
        !print*,"xlc: ",tmplc(1),tmplc(2),sqsum1local,sqsum2local
        
        !for the calculation of bunched beam current.
!        call MPI_REDUCE(localaa,aa,3,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
!        call MPI_REDUCE(localbb,bb,3,MPI_DOUBLE_PRECISION,&
!                        MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call MPI_REDUCE(tmplc,tmpgl,27,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        gam = -this%refptcl(6)
        gambet = sqrt(gam**2-1.0)
        bet = sqrt(gam**2-1.0)/gam
        energy = (gam-1.)*qmc

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
          sqsum5 = sqz - z0*z0
          sqpz = tmpgl(12)*den1
          sqsum6 = sqpz - pz0*pz0
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
          zpz = tmpgl(15)*den1 - z0*pz0
          cubx = tmpgl(16)*den1
          fthx = tmpgl(17)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
                 3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(sqsum1)
          pxrms = sqrt(sqsum2)
          yrms = sqrt(sqsum3)
          pyrms = sqrt(sqsum4)
          zrms = sqrt(sqsum5)
          pzrms = sqrt(sqsum6)
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          !<<<<<<<<<<<<<<<<<<<< free format (kilean) <<<<<<<<<<<<<<<<<<<<<<<<
          write(18,99)z,this%refptcl(5),gam,energy,bet,sqrt(glrmax)*xl
          write(24,100)z,x0*xl,xrms*xl,px0/gambet,pxrms/gambet,-xpx/epx,epx*xl
          write(25,100)z,y0*xl,yrms*xl,py0/gambet,pyrms/gambet,-ypy/epy,epy*xl
          write(26,100)z,z0*xt,zrms*xt,pz0*qmc,pzrms*qmc,-zpz/epz,epz*qmc*xt
          write(27,100)z,glmax(1)*xl,glmax(2)/gambet,glmax(3)*xl,&
                       glmax(4)/gambet,glmax(5)*xt,glmax(6)*qmc
          write(28,101)z,npctmin,npctmax,nptot
          write(29,100)z,x03*xl,px03/gambet,y03*xl,py03/gambet,z03*xt,&
                       pz03*qmc
          write(30,100)z,x04*xl,px04/gambet,y04*xl,py04/gambet,z04*xt,&
                       pz04*qmc
! output bunched beam current.
!          write(31,100)z,aa(1)*den1,aa(2)*den1,aa(3)*den1,bb(1)*den1,&
!                        bb(2)*den1,bb(3)*den1
          write(32,*)z,nptlist(1:nchrg)

!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,ypy*ypyfac,epy*xl
!          write(24,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl/gambet
!          write(25,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl/gambet
!          write(26,100)z,z0*bet*xl,zrms*bet*xl,pz0,pzrms,-zpz/epz,&
!                       epz*xl/gam**3/bet
!          write(26,100)z,z0*bet*xl,zrms*bet*xl,pz0,pzrms,-zpz/epz,epz*qmc*xt
!          write(26,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
!          write(26,100)z,z0*bet*xl,zrms*bet*xl,pz0,pzrms,-zpz/epz,epz*qmc*xt
!          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,-zpz/epz,epz*xl/(gam**3*bet)
!        write(24,99)z,xrms,pxrms,xpx*xpxfac,epx
!        write(25,99)z,yrms,pyrms,ypy*ypyfac,epy
!        write(26,99)z,zrms,pzrms,zpz*zpzfac,epz
!          write(27,100)z,glmax(1)*xl,glmax(2),glmax(3)*xl,glmax(4), &
!                     glmax(5)*xt,glmax(6)
!        write(27,100)z,glmax(1),glmax(2),glmax(3),glmax(4), &
!                     glmax(5),glmax(6)

          call flush(18)
          call flush(24)
          call flush(25)
          call flush(26)
          call flush(27)
          call flush(28)
          call flush(29)
          call flush(30)
          call flush(32)
        endif

99      format(6(1x,e15.8))
100      format(7(1x,e15.8))
101     format(1x,e15.8,3I13)

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic1_Output

        subroutine diagnostic2_Output(this,z,nchrg,nptlist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: z
        integer, intent(in) :: nchrg
        integer*8, dimension(nchrg), intent(in) :: nptlist
        integer, parameter :: nbin=100
        integer, dimension(nbin) :: tmpbin,glbin
        integer :: innp,i,my_rank,ierr,j,npctmin,npctmax,nfreq
        integer :: nii,iitmp
        integer*8 :: nptot
        double precision :: qmc,qmass,xl,xt,tg,pi,t0,gb,bgend,blg
        double precision :: ez1,ezp1,ezpp1,gam,energy,bet,gambet
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                           epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                           xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local,x0lc,x0,px0lc,px0
        double precision:: y0lc,y0,py0lc,py0,z0lc,z0,pz0lc,pz0
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(15) :: tmplc,tmpgl
        double precision :: lcrmax,glrmax,tmp1,tmp2
        double precision :: f90,f95,f99,ex90,ex95,ex99,ey90,ey95,ey99,&
        ez90,ez95,ez99,r90,r95,r99,rrms,ravg,rrmax,rtmp2,rrmaxlc,ravg2,&
        ravglc,ravg2lc
        double precision, dimension(3) :: Ealpha,Ebeta,Egamma
        double precision :: epsmxlc,epsmylc,epsmzlc,xtmp,pxtmp,ytmp,pytmp,&
                   ztmp,pztmp,epsmx,epsmy,epsmz,eps,hxeps,hyeps,hzeps,hreps
        double precision,allocatable,dimension(:) :: epsiontmp

        call starttime_Timer(t0)

        qmc = this%Mass/1.0e6
        qmass =this%Mass
        xl = Scxl
        xt = Rad2deg
        tg = this%refptcl(5)
        pi = 2*asin(1.0)

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        innp = this%Nptlocal
        nptot = this%Npt

        den1 = 1.0/dble(nptot)
        den2 = den1*den1
        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        z0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            localmax(i) = abs(this%Pts1(i,1))
          enddo
          lcrmax = this%Pts1(1,1)**2+this%Pts1(3,1)**2
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
!          tmp1 = this%Pts1(2,i)-(xl*xl*ezp1*this%Pts1(1,i)*&
!            sin(nfreq*(this%Pts1(5,i)+tg)+theta*pi/180)/(2*qmass*nfreq))
          tmp1 = this%Pts1(2,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*tmp1
          px0lc = px0lc + tmp1
          sqsum2local = sqsum2local + tmp1*tmp1
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
!          tmp2 = this%Pts1(4,i)-(xl*xl*ezp1*this%Pts1(3,i)*&
!            sin(nfreq*(this%Pts1(5,i)+tg)+theta*pi/180)/(2*qmass*nfreq))
          tmp2 = this%Pts1(4,i)
          ypylocal = ypylocal + this%Pts1(3,i)*tmp2
          py0lc = py0lc + tmp2
          sqsum4local = sqsum4local + tmp2*tmp2
          z0lc = z0lc + this%Pts1(5,i)
          sqsum5local = sqsum5local + this%Pts1(5,i)*this%Pts1(5,i)
          zpzlocal = zpzlocal + this%Pts1(5,i)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
          do j = 1, 6
            if(localmax(j).lt.abs(this%Pts1(j,i))) then
               localmax(j) = abs(this%Pts1(j,i))
            endif
          enddo
          if(lcrmax.lt.(this%Pts1(1,i)**2+this%Pts1(3,i)**2)) then
            lcrmax = this%Pts1(1,i)**2 + this%Pts1(3,i)**2
          endif
        enddo

        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        
        call MPI_ALLREDUCE(tmplc,tmpgl,15,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        gam = -this%refptcl(6)
        bet = sqrt(gam**2-1.0)/gam
        gambet = sqrt(gam**2-1.0)
        energy = (gam-1.)*qmc

        x0 = tmpgl(1)
        px0 = tmpgl(2)
        y0 = tmpgl(3)
        py0 = tmpgl(4)
        z0 = tmpgl(5)
        pz0 = tmpgl(6)
        sqsum1 = tmpgl(7) - x0*x0*den1
        sqsum2 = tmpgl(8) - px0*px0*den1
        sqsum3 = tmpgl(9) - y0*y0*den1
        sqsum4 = tmpgl(10) - py0*py0*den1
        sqsum5 = tmpgl(11) - z0*z0*den1
        sqsum6 = tmpgl(12) - pz0*pz0*den1
        xpx = tmpgl(13) - x0*px0*den1
        ypy = tmpgl(14) - y0*py0*den1
        zpz = tmpgl(15) - z0*pz0*den1
        epsx2 = (sqsum1*sqsum2-xpx*xpx)*den2
        epsy2 = (sqsum3*sqsum4-ypy*ypy)*den2
        epsz2 = (sqsum5*sqsum6-zpz*zpz)*den2
        epx = sqrt(max(epsx2,0.0d0))
        epy = sqrt(max(epsy2,0.0d0))
        epz = sqrt(max(epsz2,0.0d0))
        xrms = sqrt(sqsum1*den1)
        pxrms = sqrt(sqsum2*den1)
        yrms = sqrt(sqsum3*den1)
        pyrms = sqrt(sqsum4*den1)
        zrms = sqrt(sqsum5*den1)
        pzrms = sqrt(sqsum6*den1)
        xpx = xpx*den1
        ypy = ypy*den1
        zpz = zpz*den1
        xpxfac = 0.0
        ypyfac = 0.0
        zpzfac = 0.0
        if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
        if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
        if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
        x0 = x0*den1
        px0 = px0*den1
        y0 = y0*den1
        py0 = py0*den1
        z0 = z0*den1
        pz0 = pz0*den1

        Ealpha(1) = -xpx/epx
        Ealpha(2) = -ypy/epy
        Ealpha(3) = -zpz/epz
        Ebeta(1) = xrms*xrms*gambet/epx
        Ebeta(2) = yrms*yrms*gambet/epy
        Ebeta(3) = zrms*zrms*bet*bet/(epz/gam**3/bet)
        Egamma(:) = (1.0+Ealpha(:)*Ealpha(:))/Ebeta(:)

        allocate(epsiontmp(innp))
        epsmxlc = -1.0e10
        do i = 1, innp
          xtmp = this%Pts1(1,i) - x0
!          tmp1 = this%Pts1(2,i)-(xl*xl*ezp1*this%Pts1(1,i)*&
!            sin(nfreq*(this%Pts1(5,i)+tg)+theta*pi/180)/(2*qmass*nfreq))
          tmp1 = this%Pts1(2,i)
          pxtmp = (tmp1 - px0)/gambet
          epsiontmp(i)=Egamma(1)*xtmp*xtmp+2*Ealpha(1)*xtmp*pxtmp+&
                     Ebeta(1)*pxtmp*pxtmp
          if(epsmxlc.le.epsiontmp(i)) epsmxlc = epsiontmp(i)
        enddo

        call MPI_ALLREDUCE(epsmxlc,epsmx,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        eps = 1.0e-8
        hxeps = (epsmx+eps)/nbin

        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          iitmp = epsiontmp(i)/hxeps
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1 
        enddo
        call MPI_REDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)

        f90 = 0.9*nptot
        f95 = 0.95*nptot
        f99 = 0.99*nptot
        if(my_rank.eq.0) then
          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              ex90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                    hxeps*(i-1) 
              exit
            endif
          enddo
          !print*,"i1: ",i,nbin,glbin(i-1),glbin(i),f90,f95,f99
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              ex95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                    hxeps*(i-1) 
              exit
            endif
          enddo
          !print*,"i2: ",i,nbin,glbin(i-1),glbin(i)
          do i =1, nbin
            if(glbin(i).gt.f99) then
              ex99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                    hxeps*(i-1) 
              exit
            endif
          enddo
          !print*,"i3: ",i,nbin,glbin(i-1),glbin(i)
          !print*,"pass hxeps: ",glbin(nbin),nptot,epsmx,ex90,ex95,ex99
          ex90 = ex90*gambet
          ex95 = ex95*gambet
          ex99 = ex99*gambet
        endif

        epsmylc = -1.0e10
        do i = 1, innp
          ytmp = this%Pts1(3,i) - y0
!          tmp2 = this%Pts1(4,i)-(xl*xl*ezp1*this%Pts1(3,i)*&
!            sin(nfreq*(this%Pts1(5,i)+tg)+theta*pi/180)/(2*qmass*nfreq))
          tmp2 = this%Pts1(4,i)
          pytmp = (tmp2 - py0)/gambet
          epsiontmp(i)=Egamma(2)*ytmp*ytmp+2*Ealpha(2)*ytmp*pytmp+&
                     Ebeta(2)*pytmp*pytmp
          if(epsmylc.le.epsiontmp(i)) epsmylc = epsiontmp(i)
        enddo
        call MPI_ALLREDUCE(epsmylc,epsmy,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        hyeps = (epsmy+eps)/nbin
        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          iitmp = epsiontmp(i)/hyeps
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1
        enddo
        call MPI_REDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              ey90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                    hyeps*(i-1)
              exit
            endif
          enddo
          !print*,"i4: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              ey95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                    hyeps*(i-1)
              exit
            endif
          enddo
          !print*,"i5: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f99) then
              ey99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                    hyeps*(i-1)
              exit
            endif
          enddo
          !print*,"i6: ",i,nbin,glbin(i-1),glbin(i)
          !print*,"pass hyeps: ",glbin(nbin),nptot,epsmy,ey90,ey95,ey99
          ey90 = ey90*gambet
          ey95 = ey95*gambet
          ey99 = ey99*gambet
        endif

        epsmzlc = -1.0e10
        do i = 1, innp
          ztmp = (this%Pts1(5,i) - z0)*bet
          pztmp = (this%Pts1(6,i) - pz0)/gam**3/bet**2
          epsiontmp(i)=Egamma(3)*ztmp*ztmp+2*Ealpha(3)*ztmp*pztmp+&
                     Ebeta(3)*pztmp*pztmp
          if(epsmzlc.le.epsiontmp(i)) epsmzlc = epsiontmp(i)
        enddo
        call MPI_ALLREDUCE(epsmzlc,epsmz,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        hzeps = (epsmz+eps)/nbin
        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          iitmp = epsiontmp(i)/hzeps
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1
        enddo
        call MPI_REDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              ez90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hzeps+&
                    hzeps*(i-1)
              exit
            endif
          enddo
          !print*,"i7: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              ez95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hzeps+&
                    hzeps*(i-1)
              exit
            endif
          enddo
          !print*,"i8: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f99) then
              ez99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hzeps+&
                    hzeps*(i-1)
              exit
            endif
          enddo
          !print*,"i9: ",i,nbin,glbin(i-1),glbin(i)
          !print*,"pass hzeps: ",glbin(nbin),nptot,epsmz,ez90,ez95,ez99
          ez90 = ez90*gam**3*bet
          ez95 = ez95*gam**3*bet
          ez99 = ez99*gam**3*bet
        endif

        deallocate(epsiontmp)

        ravglc = 0.0
        ravg2lc = 0.0
        rrmaxlc = -1.0e10
        do i = 1, innp
          xtmp = this%Pts1(1,i) - x0
          ytmp = this%Pts1(3,i) - y0
          rtmp2 = xtmp*xtmp + ytmp*ytmp
          ravglc = ravglc + sqrt(rtmp2)
          ravg2lc = ravg2lc + rtmp2
          if(rrmaxlc.le.rtmp2) rrmaxlc = rtmp2
        enddo
        call MPI_REDUCE(ravglc,ravg,1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(ravg2lc,ravg2,1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(rrmaxlc,rrmax,1,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)

        hreps = (sqrt(rrmax)+eps)/nbin
        do i = 1, nbin
          tmpbin(i) = 0
        enddo

        do i = 1, innp
          iitmp = sqrt((this%Pts1(1,i)-x0)**2+(this%Pts1(3,i)-y0)**2)/hreps
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1
        enddo
        call MPI_REDUCE(tmpbin,glbin,nbin,MPI_INTEGER,&
                           MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 2, nbin
            glbin(i) = glbin(i) + glbin(i-1)
          enddo
          do i = 1, nbin
            if(glbin(i).gt.f90) then
              r90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hreps+&
                    hreps*(i-1)
              exit
            endif
          enddo
          !print*,"i10: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f95) then
              r95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hreps+&
                    hreps*(i-1)
              exit
            endif
          enddo
          !print*,"i11: ",i,nbin,glbin(i-1),glbin(i)
          do i = 1, nbin
            if(glbin(i).gt.f99) then
              r99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hreps+&
                    hreps*(i-1)
              exit
            endif
          enddo
          !print*,"i12: ",i,nbin,glbin(i-1),glbin(i)
          !print*,"pass hreps: ",glbin(nbin),nptot,sqrt(rrmax),r90,r95,r99
        endif

        if(my_rank.eq.0) then
          !<<<<<<<<<<<<<<<< free format (kilean) <<<<<<<<<<<<<<<<<<<<<<<<<<<
          write(18,99)z,this%refptcl(5),gam,energy,bet,sqrt(glrmax)*xl
          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl,ex90*xl,&
                       ex95*xl,ex99*xl
          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl,ey90*xl,&
                       ey95*xl,ey99*xl
          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,-zpz/epz,epz*qmc*xt,&
                       ez90*qmc*xt,ez95*qmc*xt,ez99*qmc*xt
          !write(26,100)z,z0*bet*xl,zrms*bet*xl,pz0,pzrms,-zpz/epz,epz*qmc*xt,&
          !             ez90*qmc*xt,ez95*qmc*xt,ez99*qmc*xt
          write(27,102)z,glmax(1)*xl,glmax(2),glmax(3)*xl,glmax(4), &
                     glmax(5)*xt,glmax(6)
          write(28,101)z,npctmin,npctmax,nptot
          ravg = ravg/nptot
          rrms = sqrt(ravg2/nptot - ravg*ravg)
          write(29,102)z,ravg*xl,rrms*xl,r90*xl,r95*xl,r99*xl,sqrt(rrmax)*xl
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          write(32,*)z,nptlist(1:nchrg)

          call flush(18)
          call flush(24)
          call flush(25)
          call flush(26)
          call flush(27)
          call flush(28)
          call flush(29)
          call flush(32)
!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,xpx*xpxfac,epx*xl
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,ypy*ypyfac,epy*xl
!          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,zpz*zpzfac,epz*qmc*xt
!          write(24,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl/gambet
!          write(25,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl/gambet
!          write(26,100)z,z0*bet*xl,zrms*bet*xl,pz0,pzrms,-zpz/epz,&
!                       epz*xl/gam**3/bet
!          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,-zpz/epz,epz*qmc*xt
!          write(26,100)z,z0*xt,zrms*xt,pz0,pzrms,-zpz/epz,epz*xl/(gam**3*bet)
!        write(24,99)z,xrms,pxrms,xpx*xpxfac,epx
!        write(25,99)z,yrms,pyrms,ypy*ypyfac,epy
!        write(26,99)z,zrms,pzrms,zpz*zpzfac,epz
!          write(27,100)z,glmax(1)*xl,glmax(2),glmax(3)*xl,glmax(4), &
!                     glmax(5)*bet*xl,glmax(6)
!        write(27,100)z,glmax(1),glmax(2),glmax(3),glmax(4), &
!                     glmax(5),glmax(6)
        endif

99      format(6(1x,e15.8))
100      format(10(1x,e15.8))
101     format(1x,e15.8,3I13)
102      format(7(1x,e15.8))

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic2_Output

        ! calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance.
        subroutine diagnosticT_Output(z,this)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), intent(in) :: this
        integer :: innp
        integer*8 :: nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: gam,energy,bet,gambet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(27) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax
        integer :: npctmin,npctmax

        call starttime_Timer(t0)

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        innp = this%Nptlocal
        nptot = this%Npt
        !print*,"pts1: ",this%Pts1(1,1),this%Pts1(2,1),my_rank,innp

!        do i = 1, innp
!          write(70,999)this%Pts1(1:9,i)
!        enddo
!999     format(9(1x,e14.6))

        den1 = 1.0/dble(nptot)
        den2 = den1*den1
        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        z0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0

        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            localmax(i) = abs(this%Pts1(i,1))
          enddo
          lcrmax = this%Pts1(1,1)**2+this%Pts1(3,1)**2
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
          x0lc3 = x0lc3 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)
          x0lc4 = x0lc4 + this%Pts1(1,i)*this%Pts1(1,i)*this%Pts1(1,i)*&
                  this%Pts1(1,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*this%Pts1(2,i)
          px0lc = px0lc + this%Pts1(2,i)
          sqsum2local = sqsum2local + this%Pts1(2,i)*this%Pts1(2,i)
          px0lc3 = px0lc3 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)
          px0lc4 = px0lc4 + this%Pts1(2,i)*this%Pts1(2,i)*this%Pts1(2,i)*&
                   this%Pts1(2,i)
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
          y0lc3 = y0lc3 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)
          y0lc4 = y0lc4 + this%Pts1(3,i)*this%Pts1(3,i)*this%Pts1(3,i)*&
                  this%Pts1(3,i)
          ypylocal = ypylocal + this%Pts1(3,i)*this%Pts1(4,i)
          py0lc = py0lc + this%Pts1(4,i)
          sqsum4local = sqsum4local + this%Pts1(4,i)*this%Pts1(4,i)
          py0lc3 = py0lc3 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)
          py0lc4 = py0lc4 + this%Pts1(4,i)*this%Pts1(4,i)*this%Pts1(4,i)*&
                   this%Pts1(4,i)
          z0lc = z0lc + this%Pts1(5,i)
          sqsum5local = sqsum5local + this%Pts1(5,i)*this%Pts1(5,i)
          z0lc3 = z0lc3 + abs(this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i))
          z0lc4 = z0lc4 + this%Pts1(5,i)*this%Pts1(5,i)*this%Pts1(5,i)*&
                          this%Pts1(5,i)

          zpzlocal = zpzlocal + this%Pts1(5,i)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc3 = pz0lc3 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)
          pz0lc4 = pz0lc4 + this%Pts1(6,i)*this%Pts1(6,i)*this%Pts1(6,i)*&
                            this%Pts1(6,i)
          do j = 1, 6
            if(localmax(j).lt.abs(this%Pts1(j,i))) then
               localmax(j) = abs(this%Pts1(j,i))
            endif
          enddo
          if(lcrmax.lt.(this%Pts1(1,i)**2+this%Pts1(3,i)**2)) then
            lcrmax = this%Pts1(1,i)**2 + this%Pts1(3,i)**2
          endif
        enddo

        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4
        
        call MPI_REDUCE(tmplc,tmpgl,27,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        gam = -this%refptcl(6)
        gambet = sqrt(gam**2-1.0)
        bet = sqrt(gam**2-1.0)/gam
        energy = (gam-1.)*qmc

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
          sqsum5 = sqz - z0*z0
          sqpz = tmpgl(12)*den1
          sqsum6 = sqpz - pz0*pz0
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
          zpz = tmpgl(15)*den1 - z0*pz0
          cubx = tmpgl(16)*den1
          fthx = tmpgl(17)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
                 3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(sqsum1)
          pxrms = sqrt(sqsum2)
          yrms = sqrt(sqsum3)
          pyrms = sqrt(sqsum4)
          zrms = sqrt(sqsum5)
          pzrms = sqrt(sqsum6)
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          write(48,100)z,this%refptcl(1)*xl,this%refptcl(2),this%refptcl(3)*xl,&
                      this%refptcl(4),this%refptcl(5)*xl,this%refptcl(6)
          write(34,100)z,x0*xl,xrms*xl,px0,pxrms,-xpx/epx,epx*xl
          write(35,100)z,y0*xl,yrms*xl,py0,pyrms,-ypy/epy,epy*xl
          write(36,100)z,z0*xl,zrms*xl,pz0,pzrms,-zpz/epz,epz*xl

          write(37,100)z,glmax(1)*xl,glmax(2),glmax(3)*xl,&
                       glmax(4),glmax(5)*xl,glmax(6)
          write(38,101)z,npctmin,npctmax,nptot
          write(39,100)z,x03*xl,px03,y03*xl,py03,z03*xl,&
                       pz03
          write(40,100)z,x04*xl,px04,y04*xl,py04,z04*xl,&
                       pz04

          call flush(48)
          call flush(34)
          call flush(35)
          call flush(36)
          call flush(37)
          call flush(38)
          call flush(39)
          call flush(40)
        endif

99      format(6(1x,e15.8))
100      format(7(1x,e15.8))
101     format(1x,e15.8,3I13)

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnosticT_Output

        subroutine ingeom_Output(nfile,z,inb,jstp,nprocrow,nproccol,&
                                geom,nx,ny,nz,myidx,myidy)
        implicit none
        include 'mpif.h'
        double precision, intent(out) :: z
        type(CompDom), intent(out) :: geom
        integer, intent(out) :: inb,jstp
        integer, intent(in) :: nfile,nprocrow,nproccol,nx,ny,nz,myidx,&
                               myidy
        integer, allocatable, dimension(:,:,:) :: Localnum
        double precision, allocatable, dimension(:,:,:) :: Localrange
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range
        integer :: my_rank,ierr
        integer :: i,j,k,l,m,n,ioerr
        character*8 name1
        character*9 name2
        character*10 name3
        character*11 name4
        character*3 name
        integer gethostname

        allocate(Localnum(2,0:nprocrow-1,0:nproccol-1))
        allocate(Localrange(4,0:nprocrow-1,0:nproccol-1))
   
        name1 = 'gm0000sx'
        name2 = 'gm0000sxx'
        name3 = 'gm0000sxxx'
        name4 = 'gm0000sxxxx'
!        ioerr = gethostname(name,3)

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        if((my_rank.ge.0).and.(my_rank.lt.128)) then
          name = 'n14'
        else if((my_rank.ge.128).and.(my_rank.lt.256)) then
          name = 'n15'
        else if((my_rank.ge.256).and.(my_rank.lt.384)) then
          name = 'n13'
        else if((my_rank.ge.384).and.(my_rank.lt.512)) then
          name = 'n16'
        else
        endif
!        else if((my_rank.ge.512).and.(my_rank.lt.640)) then
!          name = 'n14'
!        else if((my_rank.ge.640).and.(my_rank.lt.768)) then
!          name = 'n10'
!        else if((my_rank.ge.768).and.(my_rank.lt.896)) then
!          name = 'n16'
!        else if((my_rank.ge.896).and.(my_rank.lt.1024)) then
!          name = 'n12'
!        else
!        endif
  
        if(nfile < 10) then
          if(my_rank.lt.10) then
            name1(6:6) = char(nfile+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            name2(6:6) = char(nfile+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            name3(6:6) = char(nfile+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            name4(6:6) = char(nfile+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,&
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else if((nfile.ge.10).and.(nfile.lt.100)) then
          if(my_rank.lt.10) then
            i = nfile/10
            j = nfile - 10*i
            name1(5:5) = char(i+48)
            name1(6:6) = char(j+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/10
            j = nfile - 10*i
            name2(5:5) = char(i+48)
            name2(6:6) = char(j+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/10
            j = nfile - 10*i
            name3(5:5) = char(i+48)
            name3(6:6) = char(j+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/10
            j = nfile - 10*i
            name4(5:5) = char(i+48)
            name4(6:6) = char(j+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,&
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else if((nfile.ge.100).and.(nfile.lt.1000)) then
          if(my_rank.lt.10) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name1(4:4) = char(i+48)
            name1(5:5) = char(k+48)
            name1(6:6) = char(l+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name2(4:4) = char(i+48)
            name2(5:5) = char(k+48)
            name2(6:6) = char(l+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name3(4:4) = char(i+48)
            name3(5:5) = char(k+48)
            name3(6:6) = char(l+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name4(4:4) = char(i+48)
            name4(5:5) = char(k+48)
            name4(6:6) = char(l+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,&
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else
          if(my_rank.lt.10) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name1(3:3) = char(i+48)
            name1(4:4) = char(k+48)
            name1(5:5) = char(m+48)
            name1(6:6) = char(n+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name2(3:3) = char(i+48)
            name2(4:4) = char(k+48)
            name2(5:5) = char(m+48)
            name2(6:6) = char(n+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name3(3:3) = char(i+48)
            name3(4:4) = char(k+48)
            name3(5:5) = char(m+48)
            name3(6:6) = char(n+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name4(3:3) = char(i+48)
            name4(4:4) = char(k+48)
            name4(5:5) = char(m+48)
            name4(6:6) = char(n+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,&
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        endif

        read(9)z
        read(9)inb,jstp
        read(9)msize(1:3)
        read(9)range(1:6)
        read(9)Localnum(1:2,0:nprocrow-1,0:nproccol-1)
        read(9)Localrange(1:4,0:nprocrow-1,0:nproccol-1)

        close(9)

        geom%Meshsize = msize
        geom%SpatRange = range
        allocate(geom%LcTabrg(4,0:nprocrow-1,0:nproccol-1))
        allocate(geom%LcTabnm(2,0:nprocrow-1,0:nproccol-1))         
        geom%lcTabnm = Localnum
        geom%LcTabrg = Localrange

        deallocate(Localnum)
        deallocate(Localrange)

        geom%Meshnum(1) = nx
        geom%Meshnum(2) = ny
        geom%Meshnum(3) = nz

        geom%Mshlocal(3) = geom%LcTabnm(1,myidx,myidy)
        geom%Mshlocal(2) = geom%LcTabnm(2,myidx,myidy)
        geom%Mshlocal(1) = geom%Meshnum(1)

        geom%Sptrnglocal(5) = geom%LcTabrg(1,myidx,myidy)
        geom%Sptrnglocal(6) = geom%LcTabrg(2,myidx,myidy)
        geom%Sptrnglocal(3) = geom%LcTabrg(3,myidx,myidy)
        geom%Sptrnglocal(4) = geom%LcTabrg(4,myidx,myidy)
        geom%Sptrnglocal(1) = geom%SpatRange(1)
        geom%Sptrnglocal(2) = geom%SpatRange(2)

        end subroutine ingeom_Output

        subroutine phasein_Output(nfile,this)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(out) :: this
        integer :: my_rank,ierr
        integer :: i,j,k,l,m,n,ioerr
        character*8 name1
        character*9 name2
        character*10 name3
        character*11 name4
        character*3 name
        integer gethostname

        name1 = 'ph0000sx'
        name2 = 'ph0000sxx'
        name3 = 'ph0000sxxx'
        name4 = 'ph0000sxxxx'
!        ioerr = gethostname(name,3)

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
 
!        if((my_rank.ge.0).and.(my_rank.lt.128)) then
!          name = 'n14'
!        else if((my_rank.ge.128).and.(my_rank.lt.256)) then
!          name = 'n15'
!        else if((my_rank.ge.256).and.(my_rank.lt.384)) then
!          name = 'n13'
!        else if((my_rank.ge.384).and.(my_rank.lt.512)) then
!          name = 'n16'
!        else
!        endif
!        else if((my_rank.ge.512).and.(my_rank.lt.640)) then
!          name = 'n14'
!        else if((my_rank.ge.640).and.(my_rank.lt.768)) then
!          name = 'n10'
!        else if((my_rank.ge.768).and.(my_rank.lt.896)) then
!          name = 'n16'
!        else if((my_rank.ge.896).and.(my_rank.lt.1024)) then
!          name = 'n12'
!        else
!        endif

        if(nfile < 10) then
          if(my_rank.lt.10) then
            name1(6:6) = char(nfile+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            name2(6:6) = char(nfile+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            name3(6:6) = char(nfile+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            name4(6:6) = char(nfile+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else if((nfile.ge.10).and.(nfile.lt.100)) then
          if(my_rank.lt.10) then
            i = nfile/10
            j = nfile - 10*i
            name1(5:5) = char(i+48)
            name1(6:6) = char(j+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/10
            j = nfile - 10*i
            name2(5:5) = char(i+48)
            name2(6:6) = char(j+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/10
            j = nfile - 10*i
            name3(5:5) = char(i+48)
            name3(6:6) = char(j+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/10
            j = nfile - 10*i
            name4(5:5) = char(i+48)
            name4(6:6) = char(j+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else if((nfile.ge.100).and.(nfile.lt.1000)) then
          if(my_rank.lt.10) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name1(4:4) = char(i+48)
            name1(5:5) = char(k+48)
            name1(6:6) = char(l+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name2(4:4) = char(i+48)
            name2(5:5) = char(k+48)
            name2(6:6) = char(l+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name3(4:4) = char(i+48)
            name3(5:5) = char(k+48)
            name3(6:6) = char(l+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name4(4:4) = char(i+48)
            name4(5:5) = char(k+48)
            name4(6:6) = char(l+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else
          if(my_rank.lt.10) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name1(3:3) = char(i+48)
            name1(4:4) = char(k+48)
            name1(5:5) = char(m+48)
            name1(6:6) = char(n+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name2(3:3) = char(i+48)
            name2(4:4) = char(k+48)
            name2(5:5) = char(m+48)
            name2(6:6) = char(n+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name3(3:3) = char(i+48)
            name3(4:4) = char(k+48)
            name3(5:5) = char(m+48)
            name3(6:6) = char(n+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name4(3:3) = char(i+48)
            name4(4:4) = char(k+48)
            name4(5:5) = char(m+48)
            name4(6:6) = char(n+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        endif

        read(9)this%refptcl(1:6)
        read(9)this%Nptlocal
        allocate(this%Pts1(9,this%Nptlocal))
        read(9)this%Pts1(1:9,1:this%Nptlocal)

        close(9)

        end subroutine phasein_Output

        subroutine phasein2_Output(nfile,this)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(out) :: this

        open(nfile,status="old",form="unformatted")

        read(nfile)this%Nptlocal
        read(nfile)this%refptcl(1:6)
        allocate(this%Pts1(9,this%Nptlocal))
        read(nfile)this%Pts1(1:9,1:this%Nptlocal)

        close(nfile)

        end subroutine phasein2_Output

        subroutine inpoint_Output(nfile,this,z,inb,jstp,ibal,nstep,nprocrow,nproccol,&
                                  geom,nx,ny,nz,myidx,myidy,nptot)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(out) :: this
        double precision, intent(out) :: z
        type(CompDom), intent(out) :: geom
        integer, intent(out) :: inb,jstp,ibal,nstep
        integer*8, intent(inout) :: nptot
        integer, intent(in) :: nprocrow,nproccol,nx,ny,nz,myidx,myidy
        integer, allocatable, dimension(:,:,:) :: Localnum
        double precision, allocatable, dimension(:,:,:) :: Localrange
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range
        integer :: nplc,ierr
        real*8 :: fnplc,fnptot

        allocate(Localnum(2,0:nprocrow-1,0:nproccol-1))
        allocate(Localrange(4,0:nprocrow-1,0:nproccol-1))

        open(nfile,status="old",form="unformatted")

        read(nfile)z
        read(nfile)inb,jstp,ibal,nstep
        read(nfile)msize(1:3)
        read(nfile)range(1:6)
        read(nfile)Localnum(1:2,0:nprocrow-1,0:nproccol-1)
        read(nfile)Localrange(1:4,0:nprocrow-1,0:nproccol-1)

        read(nfile)this%Nptlocal
        read(nfile)this%refptcl(1:6)
        allocate(this%Pts1(9,this%Nptlocal))
        read(nfile)this%Pts1(1:9,1:this%Nptlocal)

        close(nfile)

        geom%Meshsize = msize
        geom%SpatRange = range
        allocate(geom%LcTabrg(4,0:nprocrow-1,0:nproccol-1))
        allocate(geom%LcTabnm(2,0:nprocrow-1,0:nproccol-1))
        geom%lcTabnm = Localnum
        geom%LcTabrg = Localrange

        deallocate(Localnum)
        deallocate(Localrange)

        geom%Meshnum(1) = nx
        geom%Meshnum(2) = ny
        geom%Meshnum(3) = nz

        geom%Mshlocal(3) = geom%LcTabnm(1,myidx,myidy)
        geom%Mshlocal(2) = geom%LcTabnm(2,myidx,myidy)
        geom%Mshlocal(1) = geom%Meshnum(1)

        geom%Sptrnglocal(5) = geom%LcTabrg(1,myidx,myidy)
        geom%Sptrnglocal(6) = geom%LcTabrg(2,myidx,myidy)
        geom%Sptrnglocal(3) = geom%LcTabrg(3,myidx,myidy)
        geom%Sptrnglocal(4) = geom%LcTabrg(4,myidx,myidy)
        geom%Sptrnglocal(1) = geom%SpatRange(1)
        geom%Sptrnglocal(2) = geom%SpatRange(2)

        nplc = this%Nptlocal
        fnplc = nplc
        !call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,&
        !                   MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(fnplc,fnptot,1,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        this%Npt = fnptot+0.1

        end subroutine inpoint_Output

        subroutine outgeom_Output(nfile,z,inb,jstp,nprocrow,nproccol,&
                                 geom)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type(CompDom), intent(in) :: geom
        integer, intent(in) :: nfile,inb,jstp,nprocrow,nproccol
        integer, allocatable, dimension(:,:,:) :: Localnum
        double precision, allocatable, dimension(:,:,:) :: Localrange
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range
        integer :: my_rank,ierr
        integer :: i,j,k,l,m,n,ioerr
        character*8 name1
        character*9 name2
        character*10 name3
        character*11 name4
        character*3 name
        integer gethostname

        allocate(Localnum(2,0:nprocrow-1,0:nproccol-1))
        call getlctabnm_CompDom(geom,Localnum)

        allocate(Localrange(4,0:nprocrow-1,0:nproccol-1))
        call getlctabrg_CompDom(geom,Localrange)
   
        call getmsize_CompDom(geom,msize)
        call getrange_CompDom(geom,range)

        name1 = 'gm0000sx'
        name2 = 'gm0000sxx'
        name3 = 'gm0000sxxx'
        name4 = 'gm0000sxxxx'
!        ioerr = gethostname(name,3)

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        if(nfile < 10) then
          if(my_rank.lt.10) then
            name1(6:6) = char(nfile+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            name2(6:6) = char(nfile+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            name3(6:6) = char(nfile+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            name4(6:6) = char(nfile+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,&
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else if((nfile.ge.10).and.(nfile.lt.100)) then
          if(my_rank.lt.10) then
            i = nfile/10
            j = nfile - 10*i
            name1(5:5) = char(i+48)
            name1(6:6) = char(j+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/10
            j = nfile - 10*i
            name2(5:5) = char(i+48)
            name2(6:6) = char(j+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/10
            j = nfile - 10*i
            name3(5:5) = char(i+48)
            name3(6:6) = char(j+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/10
            j = nfile - 10*i
            name4(5:5) = char(i+48)
            name4(6:6) = char(j+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,&
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else if((nfile.ge.100).and.(nfile.lt.1000)) then
          if(my_rank.lt.10) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name1(4:4) = char(i+48)
            name1(5:5) = char(k+48)
            name1(6:6) = char(l+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name2(4:4) = char(i+48)
            name2(5:5) = char(k+48)
            name2(6:6) = char(l+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name3(4:4) = char(i+48)
            name3(5:5) = char(k+48)
            name3(6:6) = char(l+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name4(4:4) = char(i+48)
            name4(5:5) = char(k+48)
            name4(6:6) = char(l+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,&
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else
          if(my_rank.lt.10) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name1(3:3) = char(i+48)
            name1(4:4) = char(k+48)
            name1(5:5) = char(m+48)
            name1(6:6) = char(n+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name2(3:3) = char(i+48)
            name2(4:4) = char(k+48)
            name2(5:5) = char(m+48)
            name2(6:6) = char(n+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name3(3:3) = char(i+48)
            name3(4:4) = char(k+48)
            name3(5:5) = char(m+48)
            name3(6:6) = char(n+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name4(3:3) = char(i+48)
            name4(4:4) = char(k+48)
            name4(5:5) = char(m+48)
            name4(6:6) = char(n+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,&
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        endif

        write(9)z
        write(9)inb,jstp
        write(9)msize(1:3)
        write(9)range(1:6)
        write(9)Localnum(1:2,0:nprocrow-1,0:nproccol-1)
        write(9)Localrange(1:4,0:nprocrow-1,0:nproccol-1)

        close(9)

        deallocate(Localnum)
        deallocate(Localrange)

        end subroutine outgeom_Output

        subroutine phaseout_Output(nfile,this)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        integer :: my_rank,ierr
        integer :: i,j,k,l,m,n,ioerr
        character*8 name1
        character*9 name2
        character*10 name3
        character*11 name4
        character*3 name
        integer gethostname
        real*4, dimension(:,:), allocatable :: ptout

        name1 = 'ph0000sx'
        name2 = 'ph0000sxx'
        name3 = 'ph0000sxxx'
        name4 = 'ph0000sxxxx'
!        ioerr = gethostname(name,3)
!        name = 'n01'

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
 
        if(nfile < 10) then
          if(my_rank.lt.10) then
            name1(6:6) = char(nfile+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            name2(6:6) = char(nfile+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            name3(6:6) = char(nfile+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            name4(6:6) = char(nfile+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else if((nfile.ge.10).and.(nfile.lt.100)) then
          if(my_rank.lt.10) then
            i = nfile/10
            j = nfile - 10*i
            name1(5:5) = char(i+48)
            name1(6:6) = char(j+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/10
            j = nfile - 10*i
            name2(5:5) = char(i+48)
            name2(6:6) = char(j+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/10
            j = nfile - 10*i
            name3(5:5) = char(i+48)
            name3(6:6) = char(j+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/10
            j = nfile - 10*i
            name4(5:5) = char(i+48)
            name4(6:6) = char(j+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else if((nfile.ge.100).and.(nfile.lt.1000)) then
          if(my_rank.lt.10) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name1(4:4) = char(i+48)
            name1(5:5) = char(k+48)
            name1(6:6) = char(l+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name2(4:4) = char(i+48)
            name2(5:5) = char(k+48)
            name2(6:6) = char(l+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name3(4:4) = char(i+48)
            name3(5:5) = char(k+48)
            name3(6:6) = char(l+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name4(4:4) = char(i+48)
            name4(5:5) = char(k+48)
            name4(6:6) = char(l+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        else
          if(my_rank.lt.10) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name1(3:3) = char(i+48)
            name1(4:4) = char(k+48)
            name1(5:5) = char(m+48)
            name1(6:6) = char(n+48)
            name1(8:8) = char(my_rank+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name1,&
!             status="unknown",form="unformatted")
            open(9,file=name1,status="unknown",form="unformatted")
          else if((my_rank.ge.10).and.(my_rank.lt.100)) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name2(3:3) = char(i+48)
            name2(4:4) = char(k+48)
            name2(5:5) = char(m+48)
            name2(6:6) = char(n+48)
            i = my_rank/10
            j = my_rank - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name2,&
!             status="unknown",form="unformatted")
            open(9,file=name2,status="unknown",form="unformatted")
          else if((my_rank.ge.100).and.(my_rank.lt.1000)) then
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name3(3:3) = char(i+48)
            name3(4:4) = char(k+48)
            name3(5:5) = char(m+48)
            name3(6:6) = char(n+48)
            i = my_rank/100
            j = my_rank - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name3,&
!             status="unknown",form="unformatted")
            open(9,file=name3,status="unknown",form="unformatted")
          else
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - m*10
            name4(3:3) = char(i+48)
            name4(4:4) = char(k+48)
            name4(5:5) = char(m+48)
            name4(6:6) = char(n+48)
            i = my_rank/1000
            j = my_rank - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(8:8) = char(i+48)
            name4(9:9) = char(k+48)
            name4(10:10) = char(m+48)
            name4(11:11) = char(n+48)
!            open(9,file="/n/"//name//"/scratch/jiqiang/"//name4,& 
!             status="unknown",form="unformatted")
            open(9,file=name4,status="unknown",form="unformatted")
          endif
        endif

!        do i = 1, this%Nptlocal
!          write(9,100)this%Pts1(1:9,i)
!          write(9)this%Pts1(1:6,i)
!        enddo
!         write(9)this%refptcl(1:6)

        allocate(ptout(9,this%Nptlocal))
        do i = 1, this%Nptlocal
           ptout(:,i) = this%Pts1(:,i)
        enddo
        write(9)this%Nptlocal
        write(9)ptout(1:9,1:this%Nptlocal)

        close(9)
100     format(9(1x,e14.7))

        deallocate(ptout)

        end subroutine phaseout_Output

        subroutine outpoint_Output(nfile,this,z,inb,jstp,ibal,nstep,&
                  nprocrow,nproccol,geom)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: z
        type(CompDom), intent(in) :: geom
        integer, intent(in) :: inb,jstp,ibal,nstep,nprocrow,nproccol
        integer, allocatable, dimension(:,:,:) :: Localnum
        double precision, allocatable, dimension(:,:,:) :: Localrange
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range

        allocate(Localnum(2,0:nprocrow-1,0:nproccol-1))
        call getlctabnm_CompDom(geom,Localnum)

        allocate(Localrange(4,0:nprocrow-1,0:nproccol-1))
        call getlctabrg_CompDom(geom,Localrange)
  
        call getmsize_CompDom(geom,msize)
        call getrange_CompDom(geom,range)

        open(nfile,status="unknown",form="unformatted")

        write(nfile)z
        write(nfile)inb,jstp,ibal,nstep
        write(nfile)msize(1:3)
        write(nfile)range(1:6)
        write(nfile)Localnum(1:2,0:nprocrow-1,0:nproccol-1)
        write(nfile)Localrange(1:4,0:nprocrow-1,0:nproccol-1)

        deallocate(Localnum)
        deallocate(Localrange)

        write(nfile)this%Nptlocal
        write(nfile)this%refptcl(1:6)
        write(nfile)this%Pts1(1:9,1:this%Nptlocal)

        close(nfile)

        end subroutine outpoint_Output

        subroutine phaseold_Output(nfile,this)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,sixnpt,mnpt
        integer, allocatable, dimension(:) :: nptlist
        double precision, allocatable,dimension(:,:) :: recvbuf

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        call MPI_ALLREDUCE(this%Nptlocal,mnpt,1,MPI_INTEGER,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)

        allocate(nptlist(0:np-1))
        nptlist = 0
        allocate(recvbuf(9,mnpt))
        sixnpt = 9*this%Nptlocal

        call MPI_GATHER(this%Nptlocal,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        nptlist = 9*nptlist

        if(my_rank.eq.0) then
          open(nfile,status='unknown')
          do i = 1, this%Nptlocal
            write(nfile,100)this%Pts1(1,i),this%Pts1(2,i),this%Pts1(3,i),&
                            this%Pts1(4,i),this%Pts1(5,i),this%Pts1(6,i),&
                            this%Pts1(7,i),this%Pts1(8,i),this%Pts1(9,i)
          enddo
          do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
        
            do j = 1, nptlist(i)/9
              write(nfile,100)recvbuf(1,j),recvbuf(2,j),recvbuf(3,j),&
                              recvbuf(4,j),recvbuf(5,j),recvbuf(6,j),&
                              recvbuf(7,j),recvbuf(8,j),recvbuf(9,j)
            enddo
          enddo
          close(nfile)
        else
          call MPI_SEND(this%Pts1(1,1),sixnpt,MPI_DOUBLE_PRECISION,0,1,&
                        MPI_COMM_WORLD,ierr)
        endif

100     format(9(1x,e14.7))

        deallocate(nptlist)
        deallocate(recvbuf)

        end subroutine phaseold_Output


        subroutine phase_Output(nfile,this,samplePeriod)
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        !   Modified to output binary format when nfile < 0 (Kilean)
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        integer :: samplePeriod
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,sixnpt,mnpt,npt
        integer, allocatable, dimension(:) :: nptlist
        double precision, allocatable,dimension(:,:) :: recvbuf  
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        logical :: flagBinary
        if(nfile < 0 ) then
          flagBinary = .True.
        else
          flagBinary = .False.
        endif
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        if (samplePeriod .eq. 0) then
           samplePeriod = 1
        endif


        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        call MPI_ALLREDUCE(this%Nptlocal,mnpt,1,MPI_INTEGER,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)
        allocate(nptlist(0:np-1))
        nptlist = 0
        allocate(recvbuf(9,mnpt))
        sixnpt = 9*this%Nptlocal

        call MPI_GATHER(this%Nptlocal,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        nptlist = 9*nptlist
        
        if(my_rank.eq.0) then
          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          if(flagBinary) then
            open(-nfile,status='unknown',form='unformatted',&
                        action='write')
            npt = ceiling(real(this%Npt)/real(samplePeriod))
            write(-nfile) npt
            do i = 1, this%Nptlocal,samplePeriod
              write(-nfile) this%Pts1(1:9,i)
            enddo
            do i = 1, np-1
              call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
                            i,1,MPI_COMM_WORLD,status,ierr) 
          
              do j = 1, nptlist(i)/9,samplePeriod
                write(-nfile) recvbuf(1:9,j)
              enddo
            enddo
            close(-nfile)
          else
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            open(nfile,status='unknown')
            do i = 1, this%Nptlocal,samplePeriod
              write(nfile,100)this%Pts1(1:9,i)
            enddo
            do i = 1, np-1
              call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
                            i,1,MPI_COMM_WORLD,status,ierr) 
          
              do j = 1, nptlist(i)/9,samplePeriod
                write(nfile,100)recvbuf(1:9,j)
              enddo
            enddo
            close(nfile)
          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          endif
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
        else
          call MPI_SEND(this%Pts1(1,1),sixnpt,MPI_DOUBLE_PRECISION,0,1,&
                          MPI_COMM_WORLD,ierr)
        endif

!100     format(9(1x,e16.9))
100     format(9(1x,g30.16))

        deallocate(nptlist)
        deallocate(recvbuf)

        end subroutine phase_Output

        subroutine phaseleda_Output(nfile,this)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,sixnpt,mnpt
        integer, allocatable, dimension(:) :: nptlist
        double precision, allocatable,dimension(:,:) :: recvbuf

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        call MPI_ALLREDUCE(this%Nptlocal,mnpt,1,MPI_INTEGER,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)

        allocate(nptlist(0:np-1))
        nptlist = 0
        allocate(recvbuf(9,mnpt))
        sixnpt = 9*this%Nptlocal

        call MPI_GATHER(this%Nptlocal,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        nptlist = 9*nptlist

        if(my_rank.eq.0) then
          !open(nfile,status='unknown')
          open(unit=nfile,file="phase.out",status="unknown",&
               position="append",form="unformatted")
          write(nfile)this%Nptlocal,my_rank
          write(nfile)this%Pts1(1:9,1:this%Nptlocal)
          !do i = 1, this%Nptlocal
          !  write(nfile,100)this%Pts1(1,i),this%Pts1(2,i),this%Pts1(3,i),&
          !                  this%Pts1(4,i),this%Pts1(5,i),this%Pts1(6,i)
          !enddo
          do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
        
            write(nfile)nptlist(i)/9,i
            write(nfile)recvbuf(1:9,1:nptlist(i)/9)
            !do j = 1, nptlist(i)/9
            !  write(nfile,100)recvbuf(1,j),recvbuf(2,j),recvbuf(3,j),&
            !                  recvbuf(4,j),recvbuf(5,j),recvbuf(6,j)
            !enddo
          enddo
          close(nfile)
        else
          call MPI_SEND(this%Pts1(1,1),sixnpt,MPI_DOUBLE_PRECISION,0,1,&
                        MPI_COMM_WORLD,ierr)
        endif

100     format(9(1x,e16.9))

        deallocate(nptlist)
        deallocate(recvbuf)

        end subroutine phaseleda_Output

        !//output the 3D particle number density
        subroutine dens3d_Output(nstep,nfile,this,totnptcls,xmni,xmxi,ymni,&
          ymxi,zmni,zmxi)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nstep
        type (BeamBunch), intent(in) :: this
        integer*8, intent(in) :: totnptcls
        double precision, intent(in) :: xmni,xmxi,ymni,ymxi,zmni,zmxi
        integer, parameter :: nxcell = 100
        integer, parameter :: nycell = 100
        integer, parameter :: nzcell = 100
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,k,ii,jj,kk,icc
        integer :: ix,iy,iz,ic,nsend,req
        integer, dimension(nxcell,nycell,nzcell) :: count
        integer, dimension(4,nxcell*nycell*nzcell) :: sendbuf,recvbuf
        integer, allocatable,dimension(:) :: nptlist
        double precision :: hxx,hyy,hzz,invvol,tmp1,tmp2,tmp3,tmp4,eps
        double precision, dimension(6) :: range
        double precision, dimension(3) :: lcrange1,lcrange2,temp1
        double precision :: xl,xt
        double precision  :: xmnin,xmxin,ymnin,ymxin,zmnin,zmxin,xmin,xmax,&
        ymin,ymax,zmin,zmax

        xl = Scxl
        xt = 90.0/asin(1.0)

        xmnin = xmni/xl
        xmxin = xmxi/xl
        ymnin = ymni/xl
        ymxin = ymxi/xl
        zmnin = zmni/xt
        zmxin = zmxi/xt

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        eps = 1.0e-4
        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange2(i) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i-1) = temp1(i) + eps*temp1(i)
        enddo 
        call MPI_ALLREDUCE(lcrange2,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + eps*temp1(i)
        enddo 

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        zmin = min(zmnin,range(5))
        zmax = max(zmxin,range(6))

        hxx = (xmax-xmin)/nxcell
        hyy = (ymax-ymin)/nycell
        hzz = (zmax-zmin)/nzcell
        count = 0 
        do i = 1, this%Nptlocal
          ix = int((this%Pts1(1,i)-range(1))/hxx) + 1
          iy = int((this%Pts1(3,i)-range(3))/hyy) + 1
          iz = int((this%Pts1(5,i)-range(5))/hzz) + 1
          count(ix,iy,iz) = count(ix,iy,iz) + 1
        enddo

        invvol = 1.0/(hxx*hyy*hzz)
        icc = 0
        if(my_rank.eq.0) then
          do k = 1, nzcell
            do j =1, nycell
              do i =1, nxcell
                ic = (k-1)*nycell*nxcell+(j-1)*nxcell + i
                sendbuf(4,ic) = count(i,j,k)
              enddo
            enddo
          enddo
        else
          do k = 1, nzcell
            do j =1, nycell
              do i =1, nxcell
!                ic = (k-1)*nycell*nxcell+(j-1)*nxcell + i
!                sendbuf(4,ic) = 0
                if(count(i,j,k).gt.0) then
                  icc = icc + 1
                  sendbuf(1,icc) = i
                  sendbuf(2,icc) = j
                  sendbuf(3,icc) = k
                  sendbuf(4,icc) = count(i,j,k)
                endif
              enddo
            enddo
          enddo
        endif
        nsend = 4*icc
        allocate(nptlist(0:np-1))
        nptlist = 0

        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/4
              ii = recvbuf(1,j)
              jj = recvbuf(2,j)
              kk = recvbuf(3,j)
              icc = (kk-1)*nycell*nxcell+(jj-1)*nxcell+ii
              sendbuf(4,icc) = sendbuf(4,icc) + recvbuf(4,j)
            enddo
          enddo

          open(nfile,status='unknown')
          tmp4 = invvol/dble(totnptcls)
          do k = 1, nzcell
            do j = 1, nycell
              do i = 1, nxcell
                tmp1 = (range(1) + i*hxx - 0.5*hxx)*xl
                tmp2 = (range(3) + j*hyy - 0.5*hyy)*xl
                tmp3 = (range(5) + k*hzz - 0.5*hzz)*xt
                ic = (k-1)*nycell*nxcell+(j-1)*nxcell+i
                write(nfile,100)tmp1,tmp2,tmp3,sendbuf(4,ic)*tmp4
              enddo
            enddo
          enddo
          close(nfile)
        else
          call MPI_ISEND(sendbuf(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

100     format(4(1x,e14.7))

        deallocate(nptlist)

        end subroutine dens3d_Output

        !//output 2D particle number density.
        subroutine dens2d_Output(nstep,nfile,this,totnptcls,xmnin,xmxin,&
        pxmnin,pxmxin,ymnin,ymxin,pymnin,pymxin,zmnin,zmxin,pzmnin,pzmxin)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nstep
        type (BeamBunch), intent(in) :: this
        integer*8, intent(in) :: totnptcls
        double precision, intent(in) :: xmnin,xmxin,pxmnin,pxmxin,ymnin,&
        ymxin,pymnin,pymxin,zmnin,zmxin,pzmnin,pzmxin
        integer, parameter :: nxcell = 128
        integer, parameter :: npxcell = 128
        integer, parameter :: nycell = 128
        integer, parameter :: npycell = 128
        integer, parameter :: nzcell = 128
        integer, parameter :: npzcell = 128
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,k,ii,req,jj,kk,icc
        integer :: ix,ipx,iy,ipy,iz,ipz,ic,nsend,nfile1,nfile2
        integer :: nfile3,nfile4,nfile5
        integer, dimension(nxcell,npxcell) :: count1
        integer, dimension(nxcell,npxcell) :: count2
        integer, dimension(nxcell,npxcell) :: count3
        integer, dimension(nxcell,nycell) :: count4
        integer, dimension(nxcell,nzcell) :: count5
        integer, dimension(nycell,nzcell) :: count6
        integer, dimension(3,nxcell*npxcell) :: sendbuf1,recvbuf1
        integer, dimension(3,nycell*npycell) :: sendbuf2,recvbuf2
        integer, dimension(3,nzcell*npzcell) :: sendbuf3,recvbuf3
        integer, dimension(3,nxcell*nycell) :: sendbuf4,recvbuf4
        integer, dimension(3,nxcell*nzcell) :: sendbuf5,recvbuf5
        integer, dimension(3,nycell*nzcell) :: sendbuf6,recvbuf6
        integer, allocatable, dimension(:) :: nptlist
        double precision, dimension(6) :: range,prange
        double precision, dimension(6) :: lcrange1,lcrange2,temp1
        double precision :: hxx,hyy,hzz,hpxx,hpyy,hpzz,invvol
        double precision ::eps,tmp1,tmp2,tmp3,xmax,xmin,pxmax,pxmin,&
        ymax,ymin,pymax,pymin,zmin,zmax,pzmin,pzmax
        character*3 name
        integer :: ioerr
        integer gethostname

!        ioerr = gethostname(name,3)
        eps = 1.0e-4
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange1(i+3) = this%Pts1(2*i,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
            lcrange2(i+3) = this%Pts1(2*i,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange1(i+3) = 5000.0
            lcrange2(i) = -5000.0
            lcrange2(i+3) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
            ii = i + 3
            if(lcrange1(ii).gt.this%Pts1(2*i,j)) then
              lcrange1(ii) = this%Pts1(2*i,j)
            endif
            if(lcrange2(ii).lt.this%Pts1(2*i,j)) then
              lcrange2(ii) = this%Pts1(2*i,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,6,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)

        do i = 1, 3
          range(2*i-1) = temp1(i) + temp1(i)*eps
        enddo 
        do i = 1, 3
          prange(2*i-1) = temp1(i+3) + temp1(i+3)*eps
        enddo

        call MPI_ALLREDUCE(lcrange2,temp1,6,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + temp1(i)*eps
        enddo 
        do i = 1, 3
          prange(2*i) = temp1(i+3) + temp1(i+3)*eps
        enddo 

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        zmin = min(zmnin,range(5))
        zmax = max(zmxin,range(6))
        pxmin = min(pxmnin,prange(1))
        pxmax = max(pxmxin,prange(2))
        pymin = min(pymnin,prange(3))
        pymax = max(pymxin,prange(4))
        pzmin = min(pzmnin,prange(5))
        pzmax = max(pzmxin,prange(6))

        hxx = (xmax-xmin)/nxcell
        hpxx = (pxmax-pxmin)/npxcell
        hyy = (ymax-ymin)/nycell
        hpyy = (pymax-pymin)/npycell
        hzz = (zmax-zmin)/nzcell
        hpzz = (pzmax-pzmin)/npzcell

        count1 = 0 
        count2 = 0 
        count3 = 0 
        count4 = 0
        count5 = 0
        count6 = 0
        do i = 1, this%Nptlocal
          ix = int((this%Pts1(1,i)-range(1))/hxx) + 1
          ipx = int((this%Pts1(2,i)-prange(1))/hpxx) + 1
          count1(ix,ipx) = count1(ix,ipx) + 1
          iy = int((this%Pts1(3,i)-range(3))/hyy) + 1
          ipy = int((this%Pts1(4,i)-prange(3))/hpyy) + 1
          count2(iy,ipy) = count2(iy,ipy) + 1
          iz = int((this%Pts1(5,i)-range(5))/hzz) + 1
          ipz = int((this%Pts1(6,i)-prange(5))/hpzz) + 1
          count3(iz,ipz) = count3(iz,ipz) + 1
          count4(ix,iy) = count4(ix,iy) + 1
          count5(ix,iz) = count5(ix,iz) + 1
          count6(iy,iz) = count6(iy,iz) + 1
        enddo

        invvol = 1.0/(hxx*hpxx)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npxcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf1(3,ic) = count1(i,j)
            enddo
          enddo
        else
          do j = 1, npxcell
            do i = 1, nxcell
              if(count1(i,j).gt.0) then
                icc = icc + 1
                sendbuf1(1,icc) = i
                sendbuf1(2,icc) = j
                sendbuf1(3,icc) = count1(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        allocate(nptlist(0:np-1))
        nptlist = 0

        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf1(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf1(1,j)
              jj = recvbuf1(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf1(3,icc) = sendbuf1(3,icc) + recvbuf1(3,j)
            enddo
          enddo

!          open(unit=nfile,file="/n/"//name//&
!               "/scratch/jiqiang/XPx.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile,file='XPx.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,prange(1)+0.5*hpxx,prange(1)+npxcell*hpxx-&
           0.5*hpxx
          do j = 1, npxcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = prange(1) + j*hpxx - 0.5*hpxx
!              write(nfile)sendbuf1(3,ic)
!              write(nfile)tmp1,tmp2,sendbuf1(3,ic)*tmp3
              write(nfile,100)tmp1,tmp2,sendbuf1(3,ic)*tmp3
            enddo
          enddo
          close(nfile)
        else
          call MPI_ISEND(sendbuf1(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hyy*hpyy)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npycell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              sendbuf2(3,ic) = count2(i,j) 
            enddo
          enddo
        else
          do j = 1, npycell
            do i = 1, nycell
              if(count2(i,j).gt.0) then
                icc = icc + 1
                sendbuf2(1,icc) = i
                sendbuf2(2,icc) = j
                sendbuf2(3,icc) = count2(i,j) 
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf2(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf2(1,j)
              jj = recvbuf2(2,j)
              icc = (jj-1)*nycell+ii
              sendbuf2(3,icc) = sendbuf2(3,icc) + recvbuf2(3,j)
            enddo
          enddo

          nfile1 = nfile + 1
!          open(unit=nfile1,file="/n/"//name//&
!               "/scratch/jiqiang/YPy.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile1,file='YPy.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile1,*)nstep,range(3)+0.5*hyy,range(3)+nycell*hyy &
           -0.5*hyy,prange(3)+0.5*hpyy,prange(3)+npycell*hpyy-&
           0.5*hpyy
          do j = 1, npycell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              tmp1 = range(3) + i*hyy - 0.5*hyy
              tmp2 = prange(3) + j*hpyy - 0.5*hpyy
!              write(nfile1)sendbuf2(3,ic)
!              write(nfile1)tmp1,tmp2,sendbuf2(3,ic)*tmp3
              write(nfile1,100)tmp1,tmp2,sendbuf2(3,ic)*tmp3
            enddo
          enddo
          close(nfile1)
        else
          call MPI_ISEND(sendbuf2(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hzz*hpzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, npzcell
            do i = 1, nzcell
              ic = (j-1)*nzcell + i
              sendbuf3(3,ic) = count3(i,j) 
            enddo
          enddo
        else
          do j = 1, npzcell
            do i = 1, nzcell
              if(count3(i,j).gt.0) then
                icc = icc + 1
                sendbuf3(1,icc) = i
                sendbuf3(2,icc) = j
                sendbuf3(3,icc) = count3(i,j) 
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf3(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf3(1,j)
              jj = recvbuf3(2,j)
              icc = (jj-1)*nzcell + ii
              sendbuf3(3,icc) = sendbuf3(3,icc) + recvbuf3(3,j)
            enddo
          enddo

          nfile2 = nfile + 2
!          open(unit=nfile2,file="/n/"//name//&
!               "/scratch/jiqiang/ZPz.data",status="unknown",&
!               position="append",form="unformatted")
         open(unit=nfile2,file='ZPz.data',status='unknown',&
              position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile2,*)nstep,range(5)+0.5*hzz,range(5)+nzcell*hzz &
           -0.5*hzz,prange(5)+0.5*hpzz,prange(5)+npzcell*hpzz-&
           0.5*hpzz
          do j = 1, npzcell
            do i = 1, nzcell
              ic = (j-1)*nzcell + i
              tmp1 = range(5) + i*hzz - 0.5*hzz
              tmp2 = prange(5) + j*hpzz - 0.5*hpzz
!              write(nfile2)sendbuf3(3,ic)
!              write(nfile2)tmp1,tmp2,sendbuf3(3,ic)*tmp3
              write(nfile2,100)tmp1,tmp2,sendbuf3(3,ic)*tmp3
            enddo
          enddo
          close(nfile2)
        else
          call MPI_ISEND(sendbuf3(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hxx*hyy)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf4(3,ic) = count4(i,j)
            enddo
          enddo
        else
          do j = 1, nycell
            do i = 1, nxcell
              if(count4(i,j).gt.0) then
                icc = icc + 1
                sendbuf4(1,icc) = i
                sendbuf4(2,icc) = j
                sendbuf4(3,icc) = count4(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf4(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf4(1,j)
              jj = recvbuf4(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf4(3,icc) = sendbuf4(3,icc) + recvbuf4(3,j)
            enddo
          enddo

          nfile3 = nfile + 3
!          open(unit=nfile3,file="/n/"//name//&
!                "/scratch/jiqiang/XY.data",status="unknown",&
!                position="append",form="unformatted")
          open(unit=nfile3,file='XY.data',status='unknown',&
                position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile3,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,range(3)+0.5*hyy,range(3)+nycell*hyy-&
           0.5*hyy
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = range(3) + j*hyy - 0.5*hyy
!              write(nfile3)sendbuf4(3,ic)
!              write(nfile3)tmp1,tmp2,sendbuf4(3,ic)*tmp3
              write(nfile3,100)tmp1,tmp2,sendbuf4(3,ic)*tmp3
            enddo
          enddo
          close(nfile3)
        else
          call MPI_ISEND(sendbuf4(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hxx*hzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nzcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf5(3,ic) = count5(i,j)
            enddo
          enddo
        else
          do j = 1, nzcell
            do i = 1, nxcell
              if(count5(i,j).gt.0) then
                icc = icc + 1
                sendbuf5(1,icc) = i
                sendbuf5(2,icc) = j
                sendbuf5(3,icc) = count5(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf5(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf5(1,j)
              jj = recvbuf5(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf5(3,icc) = sendbuf5(3,icc) + recvbuf5(3,j)
            enddo
          enddo

          nfile4 = nfile + 4
!          open(unit=nfile4,file="/n/"//name//&
!               "/scratch/jiqiang/XZ.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile4,file='XZ.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile4,*)nstep,range(1)+0.5*hxx,range(1)+nxcell*hxx &
           -0.5*hxx,range(5)+0.5*hzz,range(5)+nzcell*hzz-&
           0.5*hzz
          do j = 1, nzcell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              tmp1 = range(1) + i*hxx - 0.5*hxx
              tmp2 = range(5) + j*hzz - 0.5*hzz
!              write(nfile4)sendbuf5(3,ic)
!              write(nfile4)tmp1,tmp2,sendbuf5(3,ic)*tmp3
              write(nfile4,100)tmp1,tmp2,sendbuf5(3,ic)*tmp3
            enddo
          enddo
          close(nfile4)
        else
          call MPI_ISEND(sendbuf5(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        invvol = 1.0/(hyy*hzz)
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nzcell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              sendbuf6(3,ic) = count6(i,j)
            enddo
          enddo
        else
          do j = 1, nzcell
            do i = 1, nycell
              if(count6(i,j).gt.0) then
                icc = icc + 1
                sendbuf6(1,icc) = i
                sendbuf6(2,icc) = j
                sendbuf6(3,icc) = count6(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf6(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
            do j = 1, nptlist(i)/3
              ii = recvbuf6(1,j)
              jj = recvbuf6(2,j)
              icc = (jj-1)*nycell + ii
              sendbuf6(3,icc) = sendbuf6(3,icc) + recvbuf6(3,j)
            enddo
          enddo

          nfile5 = nfile + 5
!          open(unit=nfile5,file="/n/"//name//&
!               "/scratch/jiqiang/YZ.data",status="unknown",&
!               position="append",form="unformatted")
          open(unit=nfile5,file='YZ.data',status='unknown',&
               position='append')
          tmp3 = invvol/dble(totnptcls)
          write(nfile5,*)nstep,range(3)+0.5*hyy,range(3)+nycell*hyy &
           -0.5*hyy,range(5)+0.5*hzz,range(5)+nzcell*hzz-&
           0.5*hzz
          do j = 1, nzcell
            do i = 1, nycell
              ic = (j-1)*nycell + i
              tmp1 = range(3) + i*hyy - 0.5*hyy
              tmp2 = range(5) + j*hzz - 0.5*hzz
!              write(nfile5)sendbuf6(3,ic)
!              write(nfile5)tmp1,tmp2,sendbuf6(3,ic)*tmp3
              write(nfile5,100)tmp1,tmp2,sendbuf6(3,ic)*tmp3
            enddo
          enddo
          close(nfile5)
        else
          call MPI_ISEND(sendbuf6(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

        deallocate(nptlist)

100     format(3(1x,e14.7))

        end subroutine dens2d_Output

        !//output the 1D accumulated particle density
        subroutine accdens1d_Output(nstep,nfile,this,nptot,rmxi,xmni,xmxi,&
                                 ymni,ymxi)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nstep
        integer*8, intent(in) :: nptot
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: rmxi,xmni,xmxi,ymni,ymxi
        integer, parameter :: nrcell = 128
        integer, parameter :: nxcell = 128
        integer, parameter :: nycell = 128
        integer, dimension(nxcell,nycell) :: count4
        integer, dimension(3,nxcell*nycell) :: sendbuf4,recvbuf4
        integer, allocatable, dimension(:) :: nptlist
        integer, dimension(nxcell) :: densx,acdensx
        integer, dimension(nycell) :: densy,acdensy
        integer :: np,my_rank,ierr,i,ir,ic,icc,sumx,sumy,ix,iy,nsend,&
                   nfile3,ii,jj,req,j,testsum
        integer status(MPI_STATUS_SIZE)
        integer, dimension(nrcell) :: count,glcount
        double precision :: hrr,rmax,rad,xl,hxx,hyy,xmin,xmax,ymin,ymax
        double precision, dimension(6) :: range
        double precision, dimension(3) :: lcrange1,lcrange2,temp1
        double precision :: rmxin,xmnin,xmxin,ymnin,ymxin,eps

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        xl = Scxl

        rmxin = rmxi/xl
        xmnin = xmni/xl
        xmxin = xmxi/xl
        ymnin = ymni/xl
        ymxin = ymxi/xl

        eps = 0.001
        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange2(i) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i-1) = temp1(i) + temp1(i)*eps
        enddo

        call MPI_ALLREDUCE(lcrange2,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + temp1(i)*eps
        enddo

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        rmax = max(rmxin,sqrt(xmax*xmax+ymax*ymax))
        call MPI_BARRIER(MPI_COMM_WORLD,ierr);
        !print*,"xmin: ",xmin,xmax,ymin,ymax,rmax

        hxx = (xmax-xmin)/nxcell
        hyy = (ymax-ymin)/nycell
        hrr = rmax/nrcell
        count = 0 
        count4 = 0
        do i = 1, this%Nptlocal
          rad = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+this%Pts1(3,i)*&
                     this%Pts1(3,i))
          ir = int(rad/hrr) + 1
          count(ir) = count(ir) + 1
          ix = int((this%Pts1(1,i)-xmin)/hxx) + 1
          iy = int((this%Pts1(3,i)-ymin)/hyy) + 1
          count4(ix,iy) = count4(ix,iy) + 1
        enddo

        call MPI_REDUCE(count,glcount,nrcell,MPI_INTEGER,MPI_SUM,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 2, nrcell
            glcount(i) = glcount(i) + glcount(i-1)
          enddo
          !print*,"glcount: ",glcount(nrcell)

!          open(unit=nfile,file='RadDens.data',status='unknown',&
!               position='append',form='unformatted')
          open(unit=nfile,file='RadDens.data',status='unknown',&
               position='append')

          write(nfile,*)nstep
          do i = 1, nrcell
            !write(nfile,100)dble(glcount(i))/glcount(nrcell)
            write(nfile,100)(i*hrr-0.5*hrr)*xl,dble(glcount(i))/nptot
          enddo
          close(nfile)
        endif

        allocate(nptlist(0:np-1))
        nptlist = 0
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf4(3,ic) = count4(i,j)
            enddo
          enddo
        else
          do j = 1, nycell
            do i = 1, nxcell
              if(count4(i,j).gt.0) then
                icc = icc + 1
                sendbuf4(1,icc) = i
                sendbuf4(2,icc) = j
                sendbuf4(3,icc) = count4(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf4(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf4(1,j)
              jj = recvbuf4(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf4(3,icc) = sendbuf4(3,icc) + recvbuf4(3,j)
            enddo
          enddo

          nfile3 = nfile + 3
          open(unit=nfile3,file='Yprof.data',status='unknown',&
                position='append')
          !testsum = 0
          do j = 1, nycell
            sumy = 0
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sumy = sumy + sendbuf4(3,ic)
            enddo
            !write(nfile3,100)dble(sumy)/nptot
            !testsum = testsum + sumy
            densy(j) = sumy
          enddo
          acdensy = 0
          do j = nycell/2,1,-1
            acdensy(j) = densy(j) + densy(nycell-j+1) + acdensy(j+1)
          enddo
          do j = 1,nycell/2
            acdensy(j+nycell/2) = acdensy(nycell/2-j+1)
          enddo
          write(nfile3,*)nstep
          do j = 1, nycell
            write(nfile3,100)((j*hyy-0.5*hyy)+ymin)*xl,dble(acdensy(j))/nptot
          enddo

          close(nfile3)
          !print*,"testsumy: ",testsum
        
          open(unit=nfile3,file='Xprof.data',status='unknown',&
                position='append')
          !testsum = 0
          do i = 1, nxcell
            sumx = 0
            do j = 1, nycell
              ic = (j-1)*nxcell + i
              sumx = sumx + sendbuf4(3,ic)
            enddo
            !write(nfile3,100)dble(sumx)/nptot
            !testsum = testsum + sumx
            densx(i) = sumx
          enddo
          acdensx = 0
          do j = nxcell/2,1,-1
            acdensx(j) = densx(j) + densx(nxcell-j+1) + acdensx(j+1)
          enddo
          do j = 1,nxcell/2
            acdensx(j+nxcell/2) = acdensx(nxcell/2-j+1)
          enddo
          write(nfile3,*)nstep
          do j = 1, nxcell
            write(nfile3,100)((i*hxx-0.5*hxx)+xmin)*xl,dble(acdensx(j))/nptot
          enddo
          close(nfile3)
          !print*,"testsumx: ",testsum
        else
          call MPI_ISEND(sendbuf4(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

100     format(2(1x,e14.7))

        deallocate(nptlist)

        end subroutine accdens1d_Output

        !//output 1d particle density
        subroutine dens1d_Output(nstep,nfile,this,nptot,rmxi,xmni,xmxi,&
                                  ymni,ymxi)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,nstep
        integer*8, intent(in) :: nptot
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: rmxi,xmni,xmxi,ymni,ymxi
        integer, parameter :: nrcell = 64
        integer, parameter :: nxcell = 64
        integer, parameter :: nycell = 64
        integer, dimension(nxcell,nycell) :: count4
        integer, dimension(3,nxcell*nycell) :: sendbuf4,recvbuf4
        integer, allocatable, dimension(:) :: nptlist
        integer :: np,my_rank,ierr,i,ir,ic,icc,sumx,sumy,ix,iy,nsend,&
                   nfile3,ii,jj,req,j,testsum
        integer status(MPI_STATUS_SIZE)
        integer, dimension(nrcell) :: count,glcount
        double precision :: hrr,rmax,rad,xl,hxx,hyy,xmin,xmax,ymin,ymax,eps
        double precision, dimension(6) :: range
        double precision, dimension(3) :: lcrange1,lcrange2,temp1
        double precision :: rmxin,xmnin,xmxin,ymnin,ymxin

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        xl = Scxl
        rmxin = rmxi/xl
        xmnin = xmni/xl
        xmxin = xmxi/xl
        ymnin = ymni/xl
        ymxin = ymxi/xl

        eps = 0.0001
        if(this%Nptlocal.gt.0) then
          do i = 1, 3
            lcrange1(i) = this%Pts1(2*i-1,1)
            lcrange2(i) = this%Pts1(2*i-1,1)
          enddo
        else
          do i = 1, 3
            lcrange1(i) = 5000.0
            lcrange2(i) = -5000.0
          enddo
        endif

        do j = 1, this%Nptlocal
          do i = 1, 3
            if(lcrange1(i).gt.this%Pts1(2*i-1,j)) then
              lcrange1(i) = this%Pts1(2*i-1,j)
            endif
            if(lcrange2(i).lt.this%Pts1(2*i-1,j)) then
              lcrange2(i) = this%Pts1(2*i-1,j)
            endif
          enddo
        enddo

        call MPI_ALLREDUCE(lcrange1,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i-1) = temp1(i) + temp1(i)*eps
        enddo

        call MPI_ALLREDUCE(lcrange2,temp1,3,MPI_DOUBLE_PRECISION,&
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        do i = 1, 3
          range(2*i) = temp1(i) + temp1(i)*eps
        enddo

        xmin = min(xmnin,range(1))
        xmax = max(xmxin,range(2))
        ymin = min(ymnin,range(3))
        ymax = max(ymxin,range(4))
        rmax = max(rmxin,sqrt(xmax*xmax+ymax*ymax))

        hxx = (xmax-xmin)/nxcell
        hyy = (ymax-ymin)/nycell
        hrr = rmax/nrcell
        count = 0 
        count4 = 0
        do i = 1, this%Nptlocal
          rad = sqrt(this%Pts1(1,i)*this%Pts1(1,i)+this%Pts1(3,i)*&
                     this%Pts1(3,i))
          ir = int(rad/hrr) + 1
          count(ir) = count(ir) + 1
          ix = int((this%Pts1(1,i)-xmin)/hxx) + 1
          iy = int((this%Pts1(3,i)-ymin)/hyy) + 1
          count4(ix,iy) = count4(ix,iy) + 1
        enddo

        call MPI_REDUCE(count,glcount,nrcell,MPI_INTEGER,MPI_SUM,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
!          do i = 2, nrcell
!            glcount(i) = glcount(i) + glcount(i-1)
!          enddo
!          print*,"glcount: ",glcount(nrcell)

!          open(unit=nfile,file='RadDens2.data',status='unknown',&
!               position='append',form='unformatted')
          open(unit=nfile,file='RadDens2.data',status='unknown',&
               position='append')
          write(nfile,*)nstep
          do i = 1, nrcell
            !write(nfile,100)dble(glcount(i))/glcount(nrcell)
            write(nfile,100)(i*hrr-0.5*hrr)*xl,dble(glcount(i))/nptot
          enddo
          close(nfile)
        endif

        allocate(nptlist(0:np-1))
        nptlist = 0
        icc = 0
        if(my_rank.eq.0) then
          do j = 1, nycell
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sendbuf4(3,ic) = count4(i,j)
            enddo
          enddo
        else
          do j = 1, nycell
            do i = 1, nxcell
              if(count4(i,j).gt.0) then
                icc = icc + 1
                sendbuf4(1,icc) = i
                sendbuf4(2,icc) = j
                sendbuf4(3,icc) = count4(i,j)
              endif
            enddo
          enddo
        endif
        nsend = 3*icc
        call MPI_GATHER(nsend,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          do i = 1, np-1
            call MPI_RECV(recvbuf4(1,1),nptlist(i),MPI_INTEGER,&
                          i,1,MPI_COMM_WORLD,status,ierr)
            do j = 1, nptlist(i)/3
              ii = recvbuf4(1,j)
              jj = recvbuf4(2,j)
              icc = (jj-1)*nxcell + ii
              sendbuf4(3,icc) = sendbuf4(3,icc) + recvbuf4(3,j)
            enddo
          enddo

          nfile3 = nfile + 3
          open(unit=nfile3,file='Yprof2.data',status='unknown',&
                position='append')
          testsum = 0
          write(nfile3,*)nstep
          do j = 1, nycell
            sumy = 0
            do i = 1, nxcell
              ic = (j-1)*nxcell + i
              sumy = sumy + sendbuf4(3,ic)
            enddo
            write(nfile3,100)((j*hyy-0.5*hyy)+ymin)*xl,dble(sumy)/nptot
            testsum = testsum + sumy
          enddo
          close(nfile3)
!          print*,"testsumy: ",testsum
        
          open(unit=nfile3,file='Xprof2.data',status='unknown',&
                position='append')
          testsum = 0
          write(nfile3,*)nstep
          do i = 1, nxcell
            sumx = 0
            do j = 1, nycell
              ic = (j-1)*nxcell + i
              sumx = sumx + sendbuf4(3,ic)
            enddo
            write(nfile3,100)((i*hxx-0.5*hxx)+xmin)*xl,dble(sumx)/nptot
            testsum = testsum + sumx
          enddo
          close(nfile3)
!          print*,"testsumx: ",testsum
        else
          call MPI_ISEND(sendbuf4(1,1),nsend,MPI_INTEGER,0,1,&
                        MPI_COMM_WORLD,req,ierr)
          call MPI_WAIT(req,status,ierr)
        endif

100     format(2(1x,e14.7))

        deallocate(nptlist)

        end subroutine dens1d_Output

        ! deposit particles onto grid and interpolate onto particles.
        subroutine depint_Output(innp,innx,inny,innz,rays,rhopts,&
                                ptsgeom,npx,npy,myidx,myidy,grid)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,innx,inny,innz,npx,npy,myidx,myidy
        double precision, intent (in), dimension (9, innp) :: rays
        type (CompDom), intent(in) :: ptsgeom
        type (Pgrid2d), intent(in) :: grid
        double precision, intent (out), dimension (innp) :: rhopts
        double precision, dimension (innx,inny,innz) :: rho
        integer :: indx,jndx,kndx
        integer, dimension(2,0:npx-1,0:npy-1) :: table
        integer, dimension(0:npx-1) :: xtable
        integer, dimension(0:npy-1) :: ytable
        double precision :: ab,cd,ef
        double precision :: vol
        integer :: i, j, k
        double precision :: hx,hxi,hy,hyi,hz,hzi,xmin,ymin,zmin
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range
        integer :: ix,ix1,jx,jx1,kx,kx1,n
        double precision :: sumlocal, sumtot
        integer :: my_rank,ierr,kadd,jadd 

        call getmsize_CompDom(ptsgeom,msize)
        hxi = 1.0/msize(1)
        hyi = 1.0/msize(2)
        hzi = 1.0/msize(3)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        call getrange_CompDom(ptsgeom,range)
        xmin = range(1)
        ymin = range(3)
        zmin = range(5)

        call getlctabnm_CompDom(ptsgeom,table)
        xtable(0) = 0
        do i = 1, npx - 1
          xtable(i) = xtable(i-1) + table(1,i-1,0)
        enddo

        ytable(0) = 0
        do i = 1, npy - 1
          ytable(i) = ytable(i-1) + table(2,0,i-1)
        enddo

        rho = 0.0

        if(npx.gt.1) then
          kadd = -xtable(myidx) + 1
        else
          kadd = 0
        endif
        if(npy.gt.1) then
          jadd = -ytable(myidy) + 1
        else
          jadd = 0
        endif
        do i = 1, innp   
          kx=(rays(5,i)-zmin)*hzi + 1 + kadd
          jx=(rays(3,i)-ymin)*hyi + 1 + jadd
          ix=(rays(1,i)-xmin)*hxi + 1
          ix1=ix+1
          jx1=jx+1
          kx1=kx+1
          ef=((zmin-rays(5,i))+(kx-kadd)*hz)*hzi
          cd=((ymin-rays(3,i))+(jx-jadd)*hy)*hyi
          ab=((xmin-rays(1,i))+ix*hx)*hxi
          ! (i,j,k):
          rho(ix,jx,kx) = rho(ix,jx,kx) + ab*cd*ef
          ! (i,j+1,k):
          rho(ix,jx1,kx) = rho(ix,jx1,kx) + ab*(1.0-cd)*ef
          ! (i,j+1,k+1):
          rho(ix,jx1,kx1) = rho(ix,jx1,kx1)+ab*(1.0-cd)*(1.0-ef)
          ! (i,j,k+1):
          rho(ix,jx,kx1) = rho(ix,jx,kx1)+ab*cd*(1.0-ef)
          ! (i+1,j,k+1):
          rho(ix1,jx,kx1) = rho(ix1,jx,kx1)+(1.0-ab)*cd*(1.0-ef)
          ! (i+1,j+1,k+1):
          rho(ix1,jx1,kx1) = rho(ix1,jx1,kx1)+(1.0-ab)*(1.0-cd)*(1.0-ef)
          ! (i+1,j+1,k):
          rho(ix1,jx1,kx) = rho(ix1,jx1,kx)+(1.0-ab)*(1.0-cd)*ef
          ! (i+1,j,k):
          rho(ix1,jx,kx) = rho(ix1,jx,kx)+(1.0-ab)*cd*ef
        enddo
          
        sumlocal = sum(rho)
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_REDUCE(sumlocal,sumtot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
                        MPI_COMM_WORLD,ierr)
        if(my_rank.eq.0) then
!          print*,"sumrho: ",sumtot
        endif

        call guardsum2_Fldmger(rho,innx,inny,innz,grid)

!        sumlocal = sum(rho)
!        call MPI_REDUCE(sumlocal,sumtot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
!                        MPI_COMM_WORLD,ierr)
!        if(my_rank.eq.0) then
!          print*,"sumrho1: ",sumtot
!        endif

        call boundint_Fldmger(rho,innx,inny,innz,grid)

!        sumlocal = sum(rho)
!        call MPI_REDUCE(sumlocal,sumtot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
!                        MPI_COMM_WORLD,ierr)
!        if(my_rank.eq.0) then
!          print*,"sumrho2: ",sumtot
!        endif

        do n = 1, innp
          kx=(rays(5,n)-zmin)*hzi + 1 + kadd
          jx=(rays(3,n)-ymin)*hyi + 1 + jadd
          ix=(rays(1,n)-xmin)*hxi + 1
          ix1 = ix + 1
          jx1 = jx + 1
          kx1 = kx + 1
          ef=((zmin-rays(5,n))+(kx-kadd)*hz)*hzi
          cd=((ymin-rays(3,n))+(jx-jadd)*hy)*hyi
          ab=((xmin-rays(1,n))+ix*hx)*hxi
          rhopts(n) = rho(ix,jx,kx)*ab*cd*ef  &
                  +rho(ix,jx1,kx)*ab*(1.-cd)*ef &
                  +rho(ix,jx1,kx1)*ab*(1.-cd)*(1.0-ef) &
                  +rho(ix,jx,kx1)*ab*cd*(1.0-ef) &
                  +rho(ix1,jx,kx1)*(1.0-ab)*cd*(1.0-ef) &
                  +rho(ix1,jx1,kx1)*(1.0-ab)*(1.0-cd)*(1.0-ef)&
                  +rho(ix1,jx1,kx)*(1.0-ab)*(1.0-cd)*ef &
                  +rho(ix1,jx,kx)*(1.0-ab)*cd*ef
        enddo

!        sumlocal = sum(rhopts)
!        call MPI_REDUCE(sumlocal,sumtot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
!                        MPI_COMM_WORLD,ierr)
!        if(my_rank.eq.0) then
!          print*,"sumrho3: ",sumtot
!        endif

        end subroutine depint_Output

        subroutine phasecut_Output(nfile,ray,npout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile,npout
        double precision, dimension(10,npout), intent(in) :: ray
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,sixnpt,mnpt
        integer, allocatable, dimension(:) :: nptlist
        double precision, allocatable,dimension(:,:) :: recvbuf

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        call MPI_ALLREDUCE(npout,mnpt,1,MPI_INTEGER,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)

        allocate(nptlist(0:np-1))
        nptlist = 0
        allocate(recvbuf(10,mnpt))
        sixnpt = 10*npout

        call MPI_GATHER(npout,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        nptlist = 10*nptlist

        if(my_rank.eq.0) then
          open(nfile,status='unknown')
          do i = 1, npout
            write(nfile,100)ray(1,i),ray(2,i),ray(3,i),&
                            ray(4,i),ray(5,i),ray(6,i),ray(7,i)
          enddo
          do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
                          i,1,MPI_COMM_WORLD,status,ierr)

            do j = 1, nptlist(i)/10
              write(nfile,100)recvbuf(1,j),recvbuf(2,j),recvbuf(3,j),&
                  recvbuf(4,j),recvbuf(5,j),recvbuf(6,j),recvbuf(7,j)
            enddo
          enddo
          close(nfile)
        else
          call MPI_SEND(ray(1,1),sixnpt,MPI_DOUBLE_PRECISION,0,1,&
                        MPI_COMM_WORLD,ierr)
        endif

100     format(7(1x,e11.4))

        deallocate(nptlist)
        deallocate(recvbuf)

        end subroutine phasecut_Output

        !//output 3D particle density
        subroutine density_Output(nfile,this,myidx,myidy,z,innx,inny,&
                                 innz)
        implicit none
        include 'mpif.h'
        double precision, dimension(:,:,:)  :: this
!        type (FieldQuant), intent(in) :: this
        integer,intent(in) :: nfile,myidx,myidy,innx,inny,innz
        double precision, intent(in) :: z
        double precision, allocatable, dimension(:,:,:) :: charge,recvbuf
        double precision :: tmp
        integer :: i,j,k,my_rank,np,ierr
        integer, allocatable, dimension(:) :: ngxlist,ngylist,idx,idy
        integer :: nlcgrid,numgrid,ii
        integer status(MPI_STATUS_SIZE)

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

!        call getlcgrid_FieldQuant(this,innx,inny,innz)
!        innx = size(this,1)
!        inny = size(this,2)
!        innz = size(this,3)
        allocate(charge(innx,inny-2,innz-2))

        nlcgrid = innx*(inny-2)*(innz-2)
     
        do k = 2, innz-1
          do j = 2, inny-1
            do i = 1, innx
              tmp = this(i,j,k)
!              tmp = get_FieldQuant(this,i,j-1,k-1)
              charge(i,j-1,k-1) = tmp
            enddo
          enddo
        enddo

        allocate(ngxlist(0:np-1))
        ngxlist = 0
        call MPI_GATHER(innz,1,MPI_INTEGER,ngxlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        allocate(ngylist(0:np-1))
        ngylist = 0
        call MPI_GATHER(inny,1,MPI_INTEGER,ngylist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        allocate(idx(0:np-1))
        idx = 0
        call MPI_GATHER(myidx,1,MPI_INTEGER,idx,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        allocate(idy(0:np-1))
        idy = 0
        call MPI_GATHER(myidy,1,MPI_INTEGER,idy,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          open(nfile,status='unknown')
          write(nfile,*)z
          write(nfile,*)myidx, myidy,innx,inny-2,innz-2
          do k = 1, innz-2
            write(nfile,*)k
            do j = 1, inny-2
              do i = 1, innx
                write(nfile,100)charge(i,j,k)
              enddo
            enddo
          enddo
          do i = 1, np-1
            numgrid = (ngxlist(i)-2)*(ngylist(i)-2)*innx
            allocate(recvbuf(innx,ngylist(i)-2,ngxlist(i)-2))
            call MPI_RECV(recvbuf(1,1,1),numgrid,MPI_DOUBLE_PRECISION,i,1,&
                          MPI_COMM_WORLD,status,ierr) 
            write(nfile,*)idx(i),idy(i),innx,ngylist(i)-2,ngxlist(i)-2
            do k = 1, ngxlist(i)-2
              write(nfile,*)k
              do j = 1, ngylist(i)-2
                do ii = 1, innx
                  write(nfile,100)recvbuf(ii,j,k)
                enddo
              enddo
            enddo
            deallocate(recvbuf)
          enddo
          close(nfile)
        else
          call MPI_SEND(charge(1,1,1),nlcgrid,MPI_DOUBLE_PRECISION,0,1,&
                        MPI_COMM_WORLD,ierr)
        endif

100     format(2(1x,e14.7))

        deallocate(charge)
        deallocate(idx)
        deallocate(idy)
        deallocate(ngxlist)
        deallocate(ngylist)
 
        end subroutine density_Output
        
        ! Terminate MPI
        subroutine end_Output(time)
        implicit none
        include 'mpif.h'
        double precision, intent(inout) :: time
        double precision :: endtime, mtime
        integer :: my_rank,ierr

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        endtime = MPI_WTIME()
        time = endtime - time
        call MPI_REDUCE(time,mtime,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          close(24)
          close(25)
          print*,"time: ",mtime
        endif

        !for measurement of memory
        !call system_stats()

        call MPI_Finalize(ierr)

        end subroutine end_Output

!<<<<<<<<<<<<<<<<<<<<<<<<< TBT output(Kilean) <<<<<<<<<<<<<<<<<<<<<<<<<<
        integer function get_free_unit(init)
        ! get a unit numer not occupied by others
        implicit none 
        integer,optional,intent(in) :: init
        integer :: lun 
        logical :: file_open 
        if(present(init)) then
          lun = init
        else
          lun = 4692 
        endif
        file_open = .true. 
        do while ( file_open ) 
          lun = lun + 1 
          inquire( lun, opened = file_open ) 
        end do 
        get_free_unit = lun 
        return 
        end function get_free_unit
        
        character(len=6) function num2str_int(num)
        ! type cast : from integer to string
        implicit none
        integer, intent(in) :: num
        character(len=8), parameter  :: fmt_ = "(I0)"
        write(num2str_int,fmt_) num
        return
        end function
        
        recursive subroutine sort(data2sort, iRow4Sort, nRow, nCol, iColStart, iColEnd)
        !==============================================================================
        ! quicksort algorithm of a maxtrix sorting column order
        ! input : 
        !   iRow4Sort : id of the row to be used for sorting criteria
        !   nRow      : # of rows
        !   nCol      : # of cols
        !   data2sort : matrix of dimension nRow-by-nCol
        !   iColStart/iColEnd : start/end col index of an interval to be sorted.
        !------------------------------------------------------------------------------
          implicit none
          real*8 data2sort(nRow, nCol), mid, temp(nRow)
          integer iRow4Sort, nCol, nRow, iColStart, iColEnd
          integer i, j
         
          mid = data2sort(iRow4Sort,(iColStart+iColEnd)/2)
          i = iColStart
          j = iColEnd
          do
            do while (data2sort(iRow4Sort,i) < mid)
              i=i+1
            end do
            do while (mid < data2sort(iRow4Sort,j))
              j=j-1
            end do
            if (i >= j) exit
            temp = data2sort(:,i)
            data2sort(:,i) = data2sort(:,j)
            data2sort(:,j) = temp
            i=i+1
            j=j-1
          end do
          if (iColStart < i-1) call sort(data2sort, iRow4Sort, nRow, nCol, iColStart, i-1)
          if (j+1 < iColEnd)   call sort(data2sort, iRow4Sort, nRow, nCol, j+1, iColEnd)
        end subroutine sort


        subroutine turn_by_turn_phasespace(BB,fileID)
        ! write phase-space of test particles (q=0) turn-by-turn
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: BB
        integer, intent(in) :: fileID
        integer, save :: unitfID(2,1000)
        logical, save :: isOn(1000)=.false.
        integer :: i,j,ifail,np,my_rank,ierr,tpt,mtpt,iUnit
        integer :: status(MPI_STATUS_SIZE)
        logical :: isTest(BB%Nptlocal)
        integer, allocatable, dimension(:) :: nptlist,nptdisp
        double precision, allocatable,dimension(:,:) :: recvbuf, sendbuf
             
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        
        isTest = BB%Pts1(8,1:BB%Nptlocal) < tiny(0.0)  ! inteded type cast. ignore compiler warining.
        tpt = count(isTest)
        allocate(sendbuf(7,tpt))
        tpt = 0
        do i=1,BB%Nptlocal
          if(isTest(i)) then
            tpt = tpt+1
            sendbuf(1:6,tpt) = BB%Pts1(1:6,i)
            sendbuf(7,tpt) = BB%Pts1(9,i)
          endif
        enddo
        
        allocate(nptlist(0:np-1))
        allocate(nptdisp(0:np-1))
        nptlist = 0
        nptdisp = 0
        call MPI_GATHER(tpt,1,MPI_INTEGER,nptlist,1,&
                           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        mtpt = sum(nptlist)
        nptlist = nptlist*7
        do i=0,np-2
          nptdisp(i+1) = nptlist(i)+nptdisp(i)
        enddo
        allocate(recvbuf(7,mtpt))
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(my_rank==0) print*, '[TBTphase]MPI_GATHERV, mtpt=',mtpt
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        call MPI_GATHERV(sendbuf,tpt*7,MPI_DOUBLE_PRECISION,&
                         recvbuf,nptlist,nptdisp,MPI_DOUBLE_PRECISION,&
                         0,MPI_COMM_WORLD,ierr)
        if(my_rank.eq.0) then
          call sort(recvbuf, 7, 7, mtpt, 1, mtpt)
          do i=1,1000
            if(isOn(i)) then
              if(UnitfID(2,i)==fileID) then
                iUnit = UnitfID(1,i)
                exit
              endif
            else
              isOn(i) = .true.
              UnitfID(2,i)=fileID
              UnitfID(1,i)=get_free_unit(3582)
              iUnit = UnitfID(1,i)
              open(iUnit,file='TBT.'//trim(num2str_int(fileID)),form='unformatted',&
                   action='write', iostat=ifail)
              if(ifail /= 0)  STOP '--- Error in opening TBT file ---'
              exit
            endif
          enddo
          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          print*, '[TBTphase]i,UnitfID(:,i)=',i,UnitfID(:,i)
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          if(i==1000) then
            STOP 'Error : maximum number of TBT file reached'
          endif
          
          write(iUnit) mtpt
          write(iUnit) int(recvbuf(7,:))
          write(iUnit) recvbuf(1:6,:)
          flush(iUnit)
!          do i = 1, mtpt
!            if(i>mtpt-5) print*, i,'-th recvbuf',recvbuf(:,i)
!            write(iUnit,*) recvbuf(:,i)
!          enddo
        endif
        end subroutine turn_by_turn_phasespace
        
        subroutine turn_by_turn_integral(BB,fileID,beta,alfa,tn,cn)
        ! write phase-space of test particles (q=0) turn-by-turn
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: BB
        integer, intent(in) :: fileID
        double precision, intent(in) :: beta,alfa,cn,tn
        integer, save :: unitfID(2,1000)
        logical, save :: isOn(1000)=.false.
        integer :: i,j,ifail,np,my_rank,ierr,tpt,mtpt,iUnit
        integer :: status(MPI_STATUS_SIZE)
        logical :: isTest(BB%Nptlocal)
        integer, allocatable, dimension(:) :: nptlist,nptdisp
        double precision :: xn,pxn,yn,pyn,Hinv,Iinv,gambet0,p_p0
        double precision, allocatable,dimension(:,:) :: recvbuf, sendbuf
        
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        gambet0 = sqrt(BB%refptcl(6)**2-1.0d0) 
        isTest = BB%Pts1(8,1:BB%Nptlocal) < tiny(0.0)  ! inteded type cast. ignore compiler warining.
        tpt = count(isTest)
        allocate(sendbuf(3,tpt))
        tpt = 0
        do i=1,BB%Nptlocal
          if(isTest(i)) then
            tpt = tpt+1
            xn  = BB%Pts1(1,i)*Scxl/cn/sqrt(beta)
            yn  = BB%Pts1(3,i)*Scxl/cn/sqrt(beta)
            p_p0 = sqrt((-BB%Pts1(6,i)-BB%refptcl(6))**2-1.0d0)/gambet0
            pxn = BB%Pts1(2,i)/(p_p0*gambet0)*sqrt(beta)/cn + alfa*xn
            pyn = BB%Pts1(4,i)/(p_p0*gambet0)*sqrt(beta)/cn + alfa*yn
            call InvariantPotentials(xn,yn,Hinv,Iinv)
            sendbuf(1,tpt) = (xn**2+yn**2+pxn**2+pyn**2)/2.d0+tn*Hinv/p_p0
            sendbuf(2,tpt) = sqrt((xn*pyn-yn*pxn)**2+pxn**2+xn**2+tn*Iinv/p_p0)
            sendbuf(3,tpt) = BB%Pts1(9,i)
          endif
        enddo
        
        allocate(nptlist(0:np-1))
        allocate(nptdisp(0:np-1))
        nptlist = 0
        nptdisp = 0
        call MPI_GATHER(tpt,1,MPI_INTEGER,nptlist,1,&
                           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        mtpt = sum(nptlist)
        nptlist = nptlist*3
        do i=0,np-2
          nptdisp(i+1) = nptlist(i)+nptdisp(i)
        enddo
        allocate(recvbuf(3,mtpt))
        call MPI_GATHERV(sendbuf,tpt*3,MPI_DOUBLE_PRECISION,&
                         recvbuf,nptlist,nptdisp,MPI_DOUBLE_PRECISION,&
                         0,MPI_COMM_WORLD,ierr)
        if(my_rank.eq.np-1) then
          do i=1,1000
            if(isOn(i)) then
              if(UnitfID(2,i)==fileID) then
                iUnit = UnitfID(1,i)
                exit
              endif
            else
              isOn(i) = .true.
              UnitfID(2,i)=fileID
              UnitfID(1,i)=get_free_unit(7251)
              iUnit = UnitfID(1,i)
              open(iUnit,file='TBT.integral.'//trim(num2str_int(fileID)),form='unformatted',&
                   action='write', iostat=ifail)
              if(ifail /= 0)  STOP '--- Error in opening TBT.integral file ---'
              exit
            endif
          enddo
          if(i==1000) then
            STOP 'Error : maximum number of TBT.integral file reached'
          endif
        
          call sort(recvbuf, 3, 3, mtpt, 1, mtpt)
          
          write(iUnit) mtpt
          write(iUnit) int(recvbuf(3,:))
          write(iUnit) recvbuf(1:2,:)
          flush(iUnit)
        endif
        end subroutine turn_by_turn_integral
        
        subroutine turn_by_turn_integral_on_momentum(BB,fileID,beta,alfa,tn,cn)
        ! write phase-space of test particles (q=0) turn-by-turn
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: BB
        integer, intent(in) :: fileID
        double precision, intent(in) :: beta,alfa,cn,tn
        integer, save :: unitfID(2,1000)
        logical, save :: isOn(1000)=.false.
        integer :: i,j,ifail,np,my_rank,ierr,tpt,mtpt,iUnit
        integer :: status(MPI_STATUS_SIZE)
        logical :: isTest(BB%Nptlocal)
        integer, allocatable, dimension(:) :: nptlist,nptdisp
        double precision :: xn,pxn,yn,pyn,Hinv,Iinv,gambet0
        double precision, allocatable,dimension(:,:) :: recvbuf, sendbuf
        
        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        gambet0 = sqrt(BB%refptcl(6)**2-1.0d0) 
        isTest = BB%Pts1(8,1:BB%Nptlocal) < tiny(0.0)  ! inteded type cast. ignore compiler warining.
        tpt = count(isTest)
        allocate(sendbuf(3,tpt))
        tpt = 0
        do i=1,BB%Nptlocal
          if(isTest(i)) then
            tpt = tpt+1
            xn  = BB%Pts1(1,i)*Scxl/cn/sqrt(beta)
            yn  = BB%Pts1(3,i)*Scxl/cn/sqrt(beta)
            pxn = BB%Pts1(2,i)/gambet0*sqrt(beta)/cn + alfa*xn
            pyn = BB%Pts1(4,i)/gambet0*sqrt(beta)/cn + alfa*yn
            call InvariantPotentials(xn,yn,Hinv,Iinv)
            sendbuf(1,tpt) = (xn**2+yn**2+pxn**2+pyn**2)/2.d0+tn*Hinv
            sendbuf(2,tpt) = sqrt((xn*pyn-yn*pxn)**2+pxn**2+xn**2+tn*Iinv)
            sendbuf(3,tpt) = BB%Pts1(9,i)
          endif
        enddo
        
        allocate(nptlist(0:np-1))
        allocate(nptdisp(0:np-1))
        nptlist = 0
        nptdisp = 0
        call MPI_GATHER(tpt,1,MPI_INTEGER,nptlist,1,&
                           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        mtpt = sum(nptlist)
        nptlist = nptlist*3
        do i=0,np-2
          nptdisp(i+1) = nptlist(i)+nptdisp(i)
        enddo
        allocate(recvbuf(3,mtpt))
        call MPI_GATHERV(sendbuf,tpt*3,MPI_DOUBLE_PRECISION,&
                         recvbuf,nptlist,nptdisp,MPI_DOUBLE_PRECISION,&
                         0,MPI_COMM_WORLD,ierr)
        if(my_rank.eq.np-1) then
          do i=1,1000
            if(isOn(i)) then
              if(UnitfID(2,i)==fileID) then
                iUnit = UnitfID(1,i)
                exit
              endif
            else
              isOn(i) = .true.
              UnitfID(2,i)=fileID
              UnitfID(1,i)=get_free_unit(7251)
              iUnit = UnitfID(1,i)
              open(iUnit,file='TBT.integral.onMomentum.'//trim(num2str_int(fileID)),form='unformatted',&
                   action='write', iostat=ifail)
              if(ifail /= 0)  STOP '--- Error in opening TBT.integral file ---'
              exit
            endif
          enddo
          if(i==1000) then
            STOP 'Error : maximum number of TBT.integral file reached'
          endif
        
          call sort(recvbuf, 3, 3, mtpt, 1, mtpt)
          
          write(iUnit) mtpt
          write(iUnit) int(recvbuf(3,:))
          write(iUnit) recvbuf(1:2,:)
          flush(iUnit)
        endif
        end subroutine turn_by_turn_integral_on_momentum
!>>>>>>>>>>>>>>>>>>>> end of TBToutput(Kilean) >>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<< write lost particle information <<<<<<<<<<<<<<<<<<<
        subroutine write_lost_pData(lost_pdata,lost_pID,nlost,elemType,elemID)
        implicit none
        integer, intent(in) :: elemType,elemID
        integer, intent(inout) :: nlost
        integer*8, allocatable, intent(inout) :: lost_pID(:)
        real*8, allocatable, intent(inout) :: lost_pdata(:,:)
        
        logical,save :: fileCreated = .false.
        integer,save :: iUnit
        integer :: i,myrank,ierr,np,nlost_tot
        integer,  allocatable :: nlost_list,nlost_disp
        integer*8,allocatable :: lost_pID_tot
        double precision, allocatable :: lost_pdata_tot

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        allocate(nlost_list(0:np-1))
        allocate(nlost_disp(0:np-1))
        nlost_list = 0

        call MPI_GATHER(nlost,1,MPI_INTEGER,nlost_list,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        nlost_tot = sum(nptlist)
        if (nlost_tot=0) return
        nlost_disp = 0
        do i=0,np-2
          nlost_disp(i+1) = nlost_list(i)+nlost_disp(i)
        enddo
        allocate(lost_pdata_tot(3,nlost_tot))
        allocate(lost_pID_tot(nlost_tot))
        call MPI_GATHERV(lost_pdata,nlost*3,MPI_DOUBLE_PRECISION,&
                         lost_pdata_tot,nlost_list*3,nlost_disp*3,MPI_DOUBLE_PRECISION,&
                         0,MPI_COMM_WORLD,ierr)
        call MPI_GATHERV(lost_pID,nlost*2,MPI_LONG,&
                         lost_pID_tot,nlost_list,nlost_disp,MPI_LONG,&
                         0,MPI_COMM_WORLD,ierr)
        if(my_rank.eq.np-1) then
          if(.not. fileCreated) then
            fileCreated = .true.
            iUnit=get_free_unit(4351)
            open(iUnit,file='lost_partcl.data',action='write', iostat=ifail)
            if(ifail /= 0)  STOP '--- Error in opening lost_partcl.data file ---'
            write(iUnit,*) 'z','x','y','partcl ID','element type ID','element sequence number'
          endif
          do i=1,nlost
            write(iUnit,*) lost_pdata_tot(1:3,i),lost_pID_tot(i),elemType,elemID
          enddo
          flush(iUnit)
        endif
        
        nlost = 0
        deallocate(lost_pdata,lost_pID)
        
        end subroutine
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      end module Outputclass
