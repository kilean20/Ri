!----------------------------------------------------------------
! (c) Copyright, 2006 by the Regents of the University of California.
! BPMclass: Beam position monitor class in Lattice module of APPLICATION 
!           layer.
! Version: 2.0
! Author: Ji Qiang, LBNL, 1/17/06
! Description: This class defines the different beam diagnostics at given
!              beam position.
! Comments:
!  1) Itype = -1, shift the transverse centroid position to 0.
!  2) Itype = -2, shift the transverse centroid position and angle to 0.
!                 (this one not work yet due to conflict of definition)
!  *) Itype = -9,  pipe_override   !<<<<<<< Kilean <<<<<<<<
!                  Param(2) : pipe_id(1=rectangular, 2=ellipse)
!                  Param(3:4) : rad_x, rad_y  ! >>>>>>>>>>>>
!  3) Itype = -10, mismatch the beam distribution by the amount given in
!                  Param(3) - Param(8).  
!  4) Itype = -14, lumped space-charge and wake field kick
!                  Param(3): lumped kick length, Param(4): aawk, Param(5): ggwk,
!                  Param(6): lengwk. If lengwk < 0, no wake field kick
!  5) Itype = -15, switch for lumped space-charge and wake field kick
!                  Param(3) > 0, lumped kick, otherwise, separated one
!  6) Itype = -21, shift the beam centroid in 6D phase space by the amount
!                  given in Param(3) - Param(8).
!  7) Itype = -40:
        !kick the beam longitudinally by the rf nonlinearity (the linear
        !part has been included in the map integrator and substracted.)
        !Param(3); vmax (V)
        !Param(4); phi0 (degree)
!  8) Itype = -45, apply a zero-length thin-lens focusing kick
!                  Param(3): horizontal focusing strength kx
!                  Param(4): vertical focusing strength ky
!  9) Itype = -46, apply a zero-length linear transformation
!                  corresponding to a phase advance rotation
!                  inside the IOTA toy lattice
!                  Param(3): length of the NLI (in m)
!                  Param(4): phase advance/2pi across the NLI
!                  Param(5): phase advance/2pi across the arc
! 10) Itype = -47, apply a zero-length linear transformation
!                  inside the IOTA toy lattice corresponding
!                  to a small phase times the generator of
!                  the linear D&N map
!                  Param(3): length of the NLI (in m)
!                  Param(4): phase advance/2pi across the NLI
!                  Param(5): phase advance/2pi across the arc
!                  Param(6): dimensionless strength of the NLI
!----------------------------------------------------------------
      module BPMclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 8
        type BPM
          !Itype < 0
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : radius
          !      (3) : xmax
          !      (4) : pxmax
          !      (5) : ymax
          !      (6) : pymax
          !      (7) : zmax
          !      (8) : pzmax
        end type BPM
        interface getparam_BPM
          module procedure getparam1_BPM,  &
                          getparam2_BPM,   &
                          getparam3_BPM
        end interface
        interface setparam_BPM
          module procedure setparam1_BPM,  &
                           setparam2_BPM, setparam3_BPM
        end interface
      contains
        subroutine construct_BPM(this,numseg,nmpstp,type,blength)
        implicit none
        type (BPM), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_BPM
   
        subroutine setparam1_BPM(this,i,value)
        implicit none
        type (BPM), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_BPM

        subroutine setparam2_BPM(this,values)
        implicit none
        type (BPM), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_BPM

        subroutine setparam3_BPM(this,numseg,nmpstp,type,blength)
        implicit none
        type (BPM), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_BPM
   
        subroutine getparam1_BPM(this,i,blparam) 
        implicit none 
        type (BPM), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_BPM
  
        subroutine getparam2_BPM(this,blparams)
        implicit none
        type (BPM), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams
        !<<<<<<<< Kilean <<<<<<<<<<
        !blparams = this%Param
        blparams(1:Nparam) = this%Param   ! Nparam = 8 in ver. 11/6/2018
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>
        end subroutine getparam2_BPM

        subroutine getparam3_BPM(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (BPM), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_BPM

        subroutine shift_BPM(Pts1,itype,innp,nptot)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: itype,innp
        integer*8, intent(in) :: nptot
        double precision, pointer, dimension(:,:) :: Pts1
        double precision:: x0lc,px0lc,y0lc,py0lc
        double precision, dimension(4) :: tmplc,tmpgl
        integer :: i,j,ierr

        tmplc = 0.0
        tmpgl = 0.0
        if(itype.eq.(-2)) then
          x0lc = 0.0
          px0lc = 0.0
          y0lc = 0.0
          py0lc = 0.0
          do i = 1, innp
            x0lc = x0lc + Pts1(1,i)
            px0lc = px0lc + Pts1(2,i)
            y0lc = y0lc + Pts1(3,i)
            py0lc = py0lc + Pts1(4,i)
          enddo

          tmplc(1) = x0lc
          tmplc(2) = px0lc
          tmplc(3) = y0lc
          tmplc(4) = py0lc
        
          call MPI_ALLREDUCE(tmplc,tmpgl,4,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)

          tmpgl(1) = tmpgl(1)/nptot
          tmpgl(2) = tmpgl(2)/nptot
          tmpgl(3) = tmpgl(3)/nptot
          tmpgl(4) = tmpgl(4)/nptot

          do i = 1, innp
            Pts1(1,i) = Pts1(1,i) - tmpgl(1)
            Pts1(2,i) = Pts1(2,i) - tmpgl(2)
            Pts1(3,i) = Pts1(3,i) - tmpgl(3)
            Pts1(4,i) = Pts1(4,i) - tmpgl(4)
          enddo
        else if(itype.eq.(-1)) then
          x0lc = 0.0
          y0lc = 0.0
          do i = 1, innp
            x0lc = x0lc + Pts1(1,i)
            y0lc = y0lc + Pts1(3,i)
          enddo

          tmplc(1) = x0lc
          tmplc(3) = y0lc
        
          call MPI_ALLREDUCE(tmplc,tmpgl,4,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)

          tmpgl(1) = tmpgl(1)/nptot
          tmpgl(3) = tmpgl(3)/nptot

          do i = 1, innp
            Pts1(1,i) = Pts1(1,i) - tmpgl(1)
            Pts1(3,i) = Pts1(3,i) - tmpgl(3)
          enddo
        else
        endif

        end subroutine shift_BPM

        !mismatch the beam at given location.
        !Here, the storage Param(3:8) is used to store the mismatch factors
        subroutine scale_BPM(Pts1,innp,xmis,pxmis,ymis,pymis,zmis,pzmis)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: xmis,pxmis,ymis,pymis,zmis,pzmis
        integer :: i
 
        do i = 1, innp
            Pts1(1,i) = xmis*Pts1(1,i)
            Pts1(2,i) = pxmis*Pts1(2,i)
            Pts1(3,i) = ymis*Pts1(3,i)
            Pts1(4,i) = pymis*Pts1(4,i)
            Pts1(5,i) = zmis*Pts1(5,i)
            Pts1(6,i) = pzmis*Pts1(6,i)
        enddo
 
        end subroutine scale_BPM

        !added by M. I. of KEK
        !shift the beam centroid in the 6D phase space.
        !This element can be used to model steering magnet etc.
        !Here, the storage Param(3:8) is used to store the amount of shift.
        !drange(3); shift in x (m)
        !drange(4); shift in Px (rad)
        !drange(5); shift in y (m)
        !drange(6); shift in Py (rad)
        !drange(7); shift in z (deg)
        !drange(8); shift in Pz (MeV)
        !gam; relativistic gamma factor for design particle
        !mass; mass in eV/c^2
        subroutine kick_BPM(Pts1,innp,xshift,pxshift,yshift,pyshift,zshift,&
                   pzshift,gam,mass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: xshift,pxshift,yshift,pyshift,zshift,pzshift
        double precision, intent(in) :: gam,mass
        integer :: i
        double precision :: gambetz
 
        do i = 1, innp
            gambetz = sqrt((gam-Pts1(6,i))**2-Pts1(2,i)**2-Pts1(4,i)**2-1.0)
            Pts1(1,i) = Pts1(1,i)+xshift/Scxl
            Pts1(2,i) = Pts1(2,i)+pxshift*gambetz
            Pts1(3,i) = Pts1(3,i)+yshift/Scxl
            Pts1(4,i) = Pts1(4,i)+pyshift*gambetz
            Pts1(5,i) = Pts1(5,i)+zshift/Rad2deg
            Pts1(6,i) = Pts1(6,i)+pzshift*1.0e6/mass
        enddo
 
        end subroutine kick_BPM

        !kick the beam longitudinally by the rf nonlinearity (the linear
        !part has been included in the map integrator and substracted.)
        !drange(3); vmax (V)
        !drange(4); phi0 (degree)
        !drange(5); horm number of rf
        subroutine kickRF_BPM(Pts1,innp,vmax,phi0,horm,mass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: vmax,phi0,mass,horm
        integer :: i
        real*8 :: vtmp,phi0lc,sinphi0,cosphi0  
 
        vtmp = vmax/mass 
        phi0lc = phi0*asin(1.0)/90
        sinphi0 = sin(phi0lc)
        cosphi0 = cos(phi0lc)
        do i = 1, innp
            Pts1(6,i) = Pts1(6,i)+ vtmp*(sinphi0-sin(Pts1(5,i)/horm+phi0lc)+ &
                                         cosphi0*Pts1(5,i)/horm)
            !Pts1(6,i) = Pts1(6,i)+ vtmp*(sinphi0-sin(Pts1(5,i)/horm+phi0lc))
        enddo
 
        end subroutine kickRF_BPM

        !kick the beam longitudinally by the rf nonlinearity (the linear
        !part has been included in the map integrator and substracted.)
        !drange(3); vmax (V)
        !drange(4); phi0 (degree)
        !drange(5); horm number of rf
        subroutine kickRF2_BPM(Pts1,innp,vmax,phi0,horm,mass,rpt)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: vmax,phi0,mass,horm
        real*8, dimension(6) :: rpt 
        integer :: i
        real*8 :: vtmp,phi0lc,sinphi0,cosphi0  
 
        vtmp = vmax/mass 
        phi0lc = phi0*asin(1.0d0)/90
        sinphi0 = sin(phi0lc)
        cosphi0 = cos(phi0lc)
        do i = 1, innp
            !Pts1(6,i) = Pts1(6,i)+ vtmp*(sinphi0-sin(Pts1(5,i)/horm+phi0lc))
            !linear approximation
            Pts1(6,i) = Pts1(6,i)- vtmp*Pts1(5,i)
        enddo
        rpt(6) = rpt(6) - vtmp*sinphi0

!        print*,"vtmp: ",vmax,phi0,horm,mass
 
        end subroutine kickRF2_BPM


        !Apply a zero-length thin lens kick in the 4D phase space.
        !drange(3); focusing strength in x (1/m)
        !drange(4); focusing strength in y (1/m)
        !gam; relativistic gamma factor for design particle
        !mass; mass in eV/c^2
        subroutine kick_focusing(Pts1,innp,kxfoc,kyfoc,gam,mass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: kxfoc,kyfoc
        double precision, intent(in) :: gam,mass
        integer :: i
        double precision :: gambet0

        gambet0 = sqrt(gam**2-1.0d0)

        !print*, 'gambet0,kx,ky:',gambet0,kxfoc,kyfoc
        !print*, 'Scxl = ',Scxl
        !print*, 'x, y = ',Pts1(1,1)*Scxl,Pts1(3,1)*Scxl
        do i = 1, innp
            Pts1(2,i) = Pts1(2,i)-kxfoc*gambet0*Pts1(1,i)*Scxl
            Pts1(4,i) = Pts1(4,i)-kyfoc*gambet0*Pts1(3,i)*Scxl
        enddo

        end subroutine kick_focusing


        !Apply a zero-length phase advance in the 4D phase space
        !for the arc in the IOTA toy lattice
        !drange(3); length of the nonlinear insert (m)
        !drange(4); phase advance across the nonlinear insert (2pi)
        !drange(5); phase advance across the arc (2pi)
        !gam; relativistic gamma factor for design particle
        !mass; mass in eV/c^2
        subroutine kick_phaseadvance(Pts1,innp,L,muNLI,muArc,gam,mass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision:: L,muNLI,muArc
!        double precision, intent(in) :: L,muNLI,muArc
        double precision, intent(in) :: gam,mass
        integer :: i
        double precision :: gambet0,m11,m12,m21,m22,psi,u
        double precision:: x0,px0,y0,py0,p_p0

        gambet0 = sqrt(gam**2-1.0d0)
        psi = 2.0d0*pi*muNLI
        u = 2.0d0*pi*muArc
!  Test only
!        gambet0 = 7.304823257567683d-2
!        Scxl = 1.59044838641231d0
!        u = 0.0d0
!        L = 1.8d0
!        muNLI = 0.3034496449165134d0
!        psi = 1.9066303504082995d0
!        pi = 4.0d0*datan(1.0d0)

        if(u.eq.0.0d0) then
          m11 = 1.0d0
          m12 = 0.0d0
          m21 = -4.0d0/L*(dsin(pi*muNLI))**2
          m22 = 1.0d0
        else
          m11 = dcos(u+psi/2.0d0)/dcos(psi/2.0d0)
          m12 = L*dsin(u)/dsin(psi)
          m21 = -2.0d0/L*dsin(u+psi)*dtan(psi/2.0d0)
          m22 = dcos(u+psi/2.0d0)/dcos(psi/2.0d0)
        endif

        !print*, 'gambet0,kx,ky:',gambet0,kxfoc,kyfoc
        !print*, 'Scxl = ',Scxl
        !print*, 'x, y = ',Pts1(1,1)*Scxl,Pts1(3,1)*Scxl
        !print*, 'pi = ',pi
        !print*, 'L,muNLI,muArc = ',L,muNLI,muArc
        !print*, 'gambet = ',gambet0
        !print*, 'm11,m12,m21,m22=',m11,m12,m21,m22
        do i = 1, innp
            x0 = Pts1(1,i)
            px0 = Pts1(2,i)
            y0 = Pts1(3,i)
            py0 = Pts1(4,i)
            !<<<<<<<<<<<<<<< p_p0 factor (Kilean) <<<<<<<<<<<<<<<<
            p_p0 = sqrt((gam-Pts1(6,i))**2-1.0d0)/gambet0
            Pts1(1,i) = m11*x0 + m12/p_p0*px0/(gambet0*Scxl)
            Pts1(2,i) = m21*x0*p_p0*gambet0*Scxl + m22*px0
            Pts1(3,i) = m11*y0 + m12/p_p0*py0/(gambet0*Scxl)
            Pts1(4,i) = m21*y0*p_p0*gambet0*Scxl + m22*py0
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        enddo

        end subroutine kick_phaseadvance

        !Apply a zero-length linear advance in the 4D phase space
        !for the arc in the IOTA toy lattice (using linear DN map)
        !drange(3); length of the nonlinear insert (m)
        !drange(4); phase advance across the nonlinear insert (2pi)
        !drange(5); phase advance across the arc (2pi)
        !drange(6); dimensionless strength of NLI
        !gam; relativistic gamma factor for design particle
        !mass; mass in eV/c^2
        subroutine kick_linDNadvance(Pts1,innp,L,muNLI,muArc,tau,gam,mass)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: L,muNLI,muArc,tau
        double precision, intent(in) :: gam,mass
        integer :: i
        double precision :: gambet0,m11,m12,m21,m22,m44,psi,u
        double precision:: x0,px0,y0,py0,kx,ky,m33,m34,m43,t

        gambet0 = sqrt(gam**2-1.0d0)
        psi = 2.0*pi*muNLI
        u = 2.0*pi*muArc
        t = -tau
        kx = dsqrt(1.0d0-2.0d0*t)
        ky = dsqrt(1.0d0+2.0d0*t)

        if(u.eq.0.0d0) then
          m11 = 1.0d0
          m12 = 0.0d0
          m21 = -4.0d0/L*(dsin(pi*muNLI))**2
          m22 = 1.0d0
        else
          m11 = dcos(kx*u)-dsin(kx*u)*dtan(psi/2.0d0)/kx
          m12 = L*dsin(kx*u)/dsin(psi)/kx
          m21 = -2.0d0*kx*dcos(kx*u)*dtan(psi/2.0d0)
          m21 = m21 + dsin(kx*u)*(-kx**2+dtan(psi/2.0d0)**2)
          m21 = m21*dsin(psi)/(kx*L)
          m22 = m11
          m33 = dcos(ky*u)-dsin(ky*u)*dtan(psi/2.0d0)/ky
          m34 = L*dsin(ky*u)/dsin(psi)/ky
          m43 = -2.0d0*ky*dcos(ky*u)*dtan(psi/2.0d0)
          m43 = m43 + dsin(ky*u)*(-ky**2+dtan(psi/2.0d0)**2)
          m43 = m43*dsin(psi)/(ky*L)
          m44 = m33
        endif

        !print*, 'gambet0,kx,ky:',gambet0,kxfoc,kyfoc
        !print*, 'Scxl = ',Scxl
        !print*, 'x, y = ',Pts1(1,1)*Scxl,Pts1(3,1)*Scxl
        !print*, 'pi = ',pi
        !print*, 'L,muNLI,muArc = ',L,muNLI,muArc
        !print*, 'gambet = ',gambet0
        !print*, 'm11,m12,m21,m22=',m11,m12,m21,m22
        do i = 1, innp
            x0 = Pts1(1,i)
            px0 = Pts1(2,i)
            y0 = Pts1(3,i)
            py0 = Pts1(4,i)
            Pts1(1,i) = m11*x0 + m12*px0/(gambet0*Scxl)
            Pts1(2,i) = m21*x0*gambet0*Scxl + m22*px0
            Pts1(3,i) = m33*y0 + m34*py0/(gambet0*Scxl)
            Pts1(4,i) = m43*y0*gambet0*Scxl + m44*py0
        enddo

        end subroutine kick_linDNadvance


      end module BPMclass
