!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! AccSimulatorclass: Linear accelerator simulator class in CONTROL layer.
! Version: beta
! Author: Ji Qiang, LBNL
! Description: This class defines functions to set up the initial beam 
!              particle distribution, field information, computational
!              domain, beam line element lattice and run the dynamics
!              simulation through the system.
! Comments:
!----------------------------------------------------------------
      module AccSimulatorclass
!from MaryLie
        use beamdata   !contains beta,gamma,brho, also contains scaling parameters
        use lieaparam, only : monoms
        use parallel
        use Pgrid2dclass
        use CompDomclass
        use FieldQuantclass
        use BeamLineElemclass
        use Ptclmgerclass
        use BeamBunchclass
        use Timerclass
        use Inputclass
        use Outputclass
        use Dataclass
        use PhysConstclass
        use NumConstclass
        use Distributionclass
        use Besselclass
        use Filterclass
        use SpaceChargeSF
        use PipeInfoClass

!        implicit none
        !# of phase dim., num. total and local particles, int. dist. 
        !and restart switch, error study switch, substep for space-charge
        !switch
        integer, private :: Dim, Nplocal,Flagdist,Rstartflg,Flagerr,&
                            Flagsubstep 
        integer*8, private :: Np

        !# of num. total x, total and local y mesh pts., type of BC, 
        !# of beam elems, type of integrator.
        integer, private :: Nx,Ny,Nz,Nxlocal,Nylocal,Nzlocal,Flagbc,&
                            Nblem,Flagmap,Flagdiag

        !# of processors in column and row direction.
        integer, private :: npcol, nprow

        !beam current, kin. energy, part. mass, and charge.
        double precision, private :: BcurrImp,Bkenergy,Bmass,Bcharge,BfreqImp,&
                                     Perdlen,xrad,yrad

        !conts. in init. dist.
        double precision, private, dimension(21) :: distparam

        !1d logical processor array.
        type (Pgrid2d), private :: grid2d

        !beam particle object and array.
        type (BeamBunch), private :: Bpts

        !beam charge density and field potential arrays.
        type (FieldQuant), private :: Potential

        !geometry object.
        type (CompDom), private :: Ageom

        !the following variables are used for restart purpose
        integer :: iend,jend,ibalend,nstepend
        double precision :: zend

        !beam line element array.
        type (BPM),target,dimension(Nbpmmax) :: beamln0
        type (DriftTube),target,dimension(Ndriftmax) :: beamln1
        type (Quadrupole),target,dimension(Nquadmax) :: beamln2
        type (DTL),target,dimension(Ndtlmax) :: beamln3
        type (CCDTL),target,dimension(Nccdtlmax) :: beamln4
        type (CCL),target,dimension(Ncclmax) :: beamln5
        type (SC),target,dimension(Nscmax) :: beamln6
        type (ConstFoc),target,dimension(Ncfmax) :: beamln7
        type (SolRF),target,dimension(Nslrfmax) :: beamln8
        type (Sol),target,dimension(Nslmax) :: beamln9
        type (Dipole),target,dimension(Ndipolemax) :: beamln10
        type (EMfld),target,dimension(Ncclmax) :: beamln11
        type (Multipole),target,dimension(Nquadmax) :: beamln12
        type (TWS),target,dimension(Nscmax) :: beamln13
        type (NonlinearLens),target,dimension(Nnllmax) :: beamln14
        type (BeamLineElem),private,dimension(Nblemtmax)::Blnelem
        !beam line element period.
        interface construct_AccSimulator
          module procedure init_AccSimulator
        end interface

        !//total # of charge state
        integer :: nchrg
        !//current list of charge state.
        !//charge/mass list of charge state.
        double precision, dimension(100) :: currlist,qmcclist
        !//number of particles of charge state.
        integer*8, dimension(100) :: nptlist
        integer*8, allocatable, dimension(:) :: nptlist0
        double precision, allocatable, dimension(:) :: currlist0,qmcclist0
      contains
        !set up objects and parameters.
        subroutine init_AccSimulator(time)
        implicit none
        include 'mpif.h'
        integer :: i,test1,test2,j
        integer :: myid,myidx,myidy,ierr,inb,jstp
        integer, allocatable, dimension(:) :: bnseg,bmpstp,bitype
        double precision, allocatable, dimension(:) :: blength,val1,&
        val2, val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,&
        val15,val16,val17,val18,val19,val20,val21,val22,val23,val24
        double precision :: time
        double precision :: t0
        double precision :: z,phsini
        double precision, dimension(2) :: tmpdr 
        double precision, dimension(5) :: tmpcf 
        double precision, dimension(8) :: tmpbpm 
        double precision, dimension(9) :: tmpquad
        double precision, dimension(10) :: tmpdipole 
        double precision, dimension(10) :: tmpnll
        double precision, dimension(11) :: tmprf
        double precision, dimension(15) :: tmpslrf
        double precision, dimension(14) :: tmp13
        double precision, dimension(25) :: tmpdtl
        integer :: iqr,idr,ibpm,iccl,iccdtl,idtl,isc,icf,islrf,isl,idipole,&
                   iemfld,myrank,imultpole,itws,inll

        !start up MPI.
        call init_Input(time)

        ! initialize Timer.
        call construct_Timer(0.0d0)

        call starttime_Timer(t0)

        !Flagmap = 0

        nptlist = 0
        currlist = 0.0
        qmcclist = 0.0
!-------------------------------------------------------------------
! get all global input parameters.
        call in_Input(Dim,Np,Nx,Ny,Nz,Flagbc,Flagdist,Rstartflg,&
              Flagmap,distparam,21,BcurrImp,Bkenergy,Bmass,Bcharge,&
        BfreqImp,xrad,yrad,Perdlen,Nblem,npcol,nprow,Flagerr,Flagdiag,&
        Flagsubstep,phsini,nchrg,nptlist,currlist,qmcclist)
        
        !<<<<<<<<<<< SmoothFocusing init (Chad) <<<<<<<<<<<<
        if(FlagBc.eq.9) call initializeSpaceChargeSF
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
 
        allocate(nptlist0(nchrg))
        allocate(currlist0(nchrg))
        allocate(qmcclist0(nchrg))
        do i = 1, nchrg
          nptlist0(i) =  nptlist(i)
          currlist0(i) =  currlist(i)
          qmcclist0(i) =  qmcclist(i)
        enddo
!-------------------------------------------------------------------
! construct 2D logical processor Cartesian coordinate
        call construct_Pgrid2d(grid2d,MPI_COMM_WORLD,nprow,npcol)
        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        if(myid.eq.0) then
          print*,"Start ImpactZ simulation:"
        endif

        !construct Constants class.
        call construct_PhysConst(BfreqImp)

!-------------------------------------------------------------------
! construct computational domain CompDom class and get local geometry 
! information on each processor.
        !if(Rstartflg.eq.1) then
        !  call ingeom_Output(1500,z,inb,jstp,nprow,npcol,Ageom,Nx,Ny,Nz,&
        !                    myidx,myidy)
        !  if(myid.eq.0) print*,"rstart at: ",z,inb,jstp
        !  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !else
          !xrad = 0.1363243029*0.2
          call construct_CompDom(Ageom,distparam,21,Flagdist,&
               Nx,Ny,Nz,grid2d,nprow,npcol,Flagbc,xrad,yrad,Perdlen)
        !endif

!-------------------------------------------------------------------
! initialize Data class.
        call init_Data()
!-------------------------------------------------------------------
! construct BeamBunch class.
        call construct_BeamBunch(Bpts,BcurrImp,Bkenergy,Bmass,Bcharge,&
                            Np,phsini)

!-------------------------------------------------------------------
! sample initial particle distribution.
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
        iend = 0
        jend = 0
        ibalend = 0
        nstepend = 0
        zend = 0.0
        if(Rstartflg.eq.1) then
          !call phasein_Output(1500,Bpts)
          !call phasein2_Output(myrank+31,Bpts)
          !call inpoint_Output(myid+31,Bpts,z,inb,jstp,nprow,npcol,&
          !     Ageom,Nx,Ny,Nz,myidx,myidy,Np)
          !if(myid.eq.0) print*,"rstart at: ",z,inb,jstp
          call inpoint_Output(myid+51,Bpts,zend,iend,jend,ibalend,nstepend,&
               nprow,npcol,Ageom,Nx,Ny,Nz,myidx,myidy,Np)
          if(myid.eq.0)print*,"rstart at: ",zend,iend,jend,ibalend,nstepend
        else
          call sample_Dist(Bpts,distparam,21,Flagdist,Ageom,grid2d,Flagbc,&
                           nchrg,nptlist0,qmcclist0,currlist0)
          !<<<<<<<<<<<<<<<<<<<<
          nptlist0(1) = Bpts%Npt
          !>>>>>>>>>>>>>>>>>>>>
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid.eq.0) print*,"pass generating initial distribution..."

        !get local particle number and mesh number on each processor.
        call getnpt_BeamBunch(Bpts,Nplocal)

!-------------------------------------------------------------------
! construct FieldQuant class objects.
        call construct_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d)

!-------------------------------------------------------------------
! construct beam line elements.
        allocate(blength(Nblem),bnseg(Nblem),bmpstp(Nblem))
        allocate(bitype(Nblem))
        allocate(val1(Nblem),val2(Nblem),val3(Nblem),val4(Nblem))
        allocate(val5(Nblem),val6(Nblem),val7(Nblem),val8(Nblem))
        allocate(val9(Nblem),val10(Nblem),val11(Nblem),val12(Nblem))
        allocate(val13(Nblem),val14(Nblem),val15(Nblem),val16(Nblem))
        allocate(val17(Nblem),val18(Nblem),val19(Nblem),val20(Nblem))
        allocate(val21(Nblem),val22(Nblem),val23(Nblem),val24(Nblem))

        call in_Input(Nblem,blength,bnseg,bmpstp,bitype,val1,val2,val3,&
        val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,&
        val15,val16,val17,val18,val19,val20,val21,val22,val23,val24)

        iccl = 0
        iccdtl = 0
        idtl = 0
        isc = 0
        idr = 0
        iqr = 0
        ibpm = 0
        icf = 0
        islrf = 0
        isl = 0
        idipole = 0
        iemfld = 0
        imultpole = 0
        itws = 0
        inll = 0
        do i = 1, Nblem
          if(bitype(i).lt.0) then
            ibpm = ibpm + 1
            call construct_BPM(beamln0(ibpm),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpbpm(1) = 0.0
            tmpbpm(2) = val1(i)
            tmpbpm(3) = val2(i)
            tmpbpm(4) = val3(i)
            tmpbpm(5) = val4(i)
            tmpbpm(6) = val5(i)
            tmpbpm(7) = val6(i)
            tmpbpm(8) = val7(i)
            call setparam_BPM(beamln0(ibpm),tmpbpm)
            Blnelem(i) = assign_BeamLineElem(beamln0(ibpm))
          else if(bitype(i).eq.0) then
            idr = idr + 1
            call construct_DriftTube(beamln1(idr),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpdr(1) = 0.0
            tmpdr(2) = val1(i)
            call setparam_DriftTube(beamln1(idr),tmpdr)
            Blnelem(i) = assign_BeamLineElem(beamln1(idr))
          else if(bitype(i).eq.1) then
            iqr = iqr + 1
            call construct_Quadrupole(beamln2(iqr),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpquad(1) = 0.0
            tmpquad(2) = val1(i)
            tmpquad(3) = val2(i)
            tmpquad(4) = val3(i)
            tmpquad(5) = val4(i)
            tmpquad(6) = val5(i)
            tmpquad(7) = val6(i)
            tmpquad(8) = val7(i)
            tmpquad(9) = val8(i)
            call setparam_Quadrupole(beamln2(iqr),tmpquad)
            Blnelem(i) = assign_BeamLineElem(beamln2(iqr))
          else if(bitype(i).eq.2) then
            icf = icf + 1
            call construct_ConstFoc(beamln7(icf),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpcf(1) = 0.0
            tmpcf(2) = val1(i)
            tmpcf(3) = val2(i)
            tmpcf(4) = val3(i)
            tmpcf(5) = val4(i)
            call setparam_ConstFoc(beamln7(icf),tmpcf)
            Blnelem(i) = assign_BeamLineElem(beamln7(icf))
          else if(bitype(i).eq.3) then
            isl = isl + 1
            call construct_Sol(beamln9(isl),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpquad(1) = 0.0
            tmpquad(2) = val1(i)
            tmpquad(3) = val2(i)
            tmpquad(4) = val3(i)
            tmpquad(5) = val4(i)
            tmpquad(6) = val5(i)
            tmpquad(7) = val6(i)
            tmpquad(8) = val7(i)
            tmpquad(9) = val8(i)
            call setparam_Sol(beamln9(isl),tmpquad)
            Blnelem(i) = assign_BeamLineElem(beamln9(isl))
          else if(bitype(i).eq.4) then
            idipole = idipole + 1
            call construct_Dipole(beamln10(idipole),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpdipole(1) = 0.0
            tmpdipole(2) = val1(i)
            tmpdipole(3) = val2(i)
            tmpdipole(4) = val3(i)
            tmpdipole(5) = val4(i)
            tmpdipole(6) = val5(i)
            tmpdipole(7) = val6(i)
            tmpdipole(8) = val7(i)
            tmpdipole(9) = val8(i)
            tmpdipole(10) = val9(i)
            call setparam_Dipole(beamln10(idipole),tmpdipole)
            Blnelem(i) = assign_BeamLineElem(beamln10(idipole))
          else if(bitype(i).eq.5) then
            imultpole = imultpole + 1
            call construct_Multipole(beamln12(imultpole),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpdipole(1) = 0.0
            tmpdipole(2) = val1(i)
            tmpdipole(3) = val2(i)
            tmpdipole(4) = val3(i)
            tmpdipole(5) = val4(i)
            tmpdipole(6) = val5(i)
            tmpdipole(7) = val6(i)
            tmpdipole(8) = val7(i)
            tmpdipole(9) = val8(i)
            tmpdipole(10) = val9(i)
            call setparam_Multipole(beamln12(imultpole),tmpdipole)
            Blnelem(i) = assign_BeamLineElem(beamln12(imultpole))
          else if(bitype(i).eq.6 .or. bitype(i).eq.9) then
            inll = inll + 1
            call construct_NonlinearLens(beamln14(inll),bnseg(i),bmpstp(i),&
            bitype(i),blength(i))
            tmpnll(1) = 0.0
            tmpnll(2) = val1(i)
            tmpnll(3) = val2(i)
            tmpnll(4) = val3(i)
            tmpnll(5) = val4(i)
            tmpnll(6) = val5(i)
            tmpnll(7) = val6(i)
            tmpnll(8) = val7(i)
            tmpnll(9) = val8(i)
            tmpnll(10) = val9(i)
            call setparam_NonlinearLens(beamln14(inll),tmpnll)
!   This line added to treat the analytical smooth focusing SC potential case
            Blnelem(i) = assign_BeamLineElem(beamln14(inll))
          else if(bitype(i).eq.101) then
            idtl = idtl + 1
            call construct_DTL(beamln3(idtl),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpdtl(1) = 0.0
            tmpdtl(2) = val1(i) 
            tmpdtl(3) = val2(i) 
            tmpdtl(4) = val3(i) 
            tmpdtl(5) = val4(i) 
            tmpdtl(6) = val5(i) 
            tmpdtl(7) = val6(i) 
            tmpdtl(8) = val7(i) 
            tmpdtl(9) = val8(i) 
            tmpdtl(10) = val9(i) 
            tmpdtl(11) = val10(i) 
            tmpdtl(12) = val11(i) 
            tmpdtl(13) = val12(i) 
            tmpdtl(14) = val13(i) 
            tmpdtl(15) = val14(i) 
            tmpdtl(16) = val15(i) 
            tmpdtl(17) = val16(i) 
            tmpdtl(18) = val17(i) 
            tmpdtl(19) = val18(i) 
            tmpdtl(20) = val19(i) 
            tmpdtl(21) = val20(i) 
            tmpdtl(22) = val21(i) 
            tmpdtl(23) = val22(i) 
            tmpdtl(24) = val23(i) 
            tmpdtl(25) = val24(i) 
            call setparam_DTL(beamln3(idtl),tmpdtl)
            Blnelem(i) = assign_BeamLineElem(beamln3(idtl))
          else if(bitype(i).eq.102) then
            iccdtl = iccdtl + 1
            call construct_CCDTL(beamln4(iccdtl),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmprf(1) = 0.0
            tmprf(2) = val1(i) 
            tmprf(3) = val2(i) 
            tmprf(4) = val3(i) 
            tmprf(5) = val4(i) 
            tmprf(6) = val5(i) 
            tmprf(7) = val6(i)
            tmprf(8) = val7(i) 
            tmprf(9) = val8(i) 
            tmprf(10) = val9(i) 
            tmprf(11) = val10(i) 
            call setparam_CCDTL(beamln4(iccdtl),tmprf)
            Blnelem(i) = assign_BeamLineElem(beamln4(iccdtl))
          else if(bitype(i).eq.103) then
            iccl = iccl + 1
            call construct_CCL(beamln5(iccl),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmprf(1) = 0.0
            tmprf(2) = val1(i) 
            tmprf(3) = val2(i) 
            tmprf(4) = val3(i) 
            tmprf(5) = val4(i) 
            tmprf(6) = val5(i) 
            tmprf(7) = val6(i) 
            tmprf(8) = val7(i) 
            tmprf(9) = val8(i) 
            tmprf(10) = val9(i) 
            tmprf(11) = val10(i) 
            call setparam_CCL(beamln5(iccl),tmprf)
            Blnelem(i) = assign_BeamLineElem(beamln5(iccl))
          else if(bitype(i).eq.104) then
            isc = isc + 1
            call construct_SC(beamln6(isc),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmprf(1) = 0.0
            tmprf(2) = val1(i) 
            tmprf(3) = val2(i) 
            tmprf(4) = val3(i) 
            tmprf(5) = val4(i) 
            tmprf(6) = val5(i) 
            tmprf(7) = val6(i)
            tmprf(8) = val7(i) 
            tmprf(9) = val8(i) 
            tmprf(10) = val9(i) 
            tmprf(11) = val10(i) 
            call setparam_SC(beamln6(isc),tmprf)
            Blnelem(i) = assign_BeamLineElem(beamln6(isc))
          else if(bitype(i).eq.105) then
            islrf = islrf + 1
            call construct_SolRF(beamln8(islrf),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpslrf(1) = 0.0
            tmpslrf(2) = val1(i) 
            tmpslrf(3) = val2(i) 
            tmpslrf(4) = val3(i) 
            tmpslrf(5) = val4(i) 
            tmpslrf(6) = val5(i) 
            tmpslrf(7) = val6(i) 
            tmpslrf(8) = val7(i) 
            tmpslrf(9) = val8(i) 
            tmpslrf(10) = val9(i) 
            tmpslrf(11) = val10(i) 
            tmpslrf(12) = val11(i) 
            tmpslrf(13) = val12(i) 
            tmpslrf(14) = val13(i) 
            tmpslrf(15) = val14(i) 
            call setparam_SolRF(beamln8(islrf),tmpslrf)
            Blnelem(i) = assign_BeamLineElem(beamln8(islrf))
          else if(bitype(i).eq.106) then
            itws = itws + 1
            call construct_TWS(beamln13(itws),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpslrf(1) = 0.0
            tmpslrf(2) = val1(i) 
            tmpslrf(3) = val2(i) 
            tmpslrf(4) = val3(i) 
            tmpslrf(5) = val4(i) 
            tmpslrf(6) = val5(i) 
            tmpslrf(7) = val6(i) 
            tmpslrf(8) = val7(i) 
            tmpslrf(9) = val8(i) 
            tmpslrf(10) = val9(i) 
            tmpslrf(11) = val10(i) 
            tmpslrf(12) = val11(i) 
            tmpslrf(13) = val12(i) 
            tmpslrf(14) = val13(i) 
            tmpslrf(15) = val14(i) 
            call setparam_TWS(beamln13(itws),tmpslrf)
            Blnelem(i) = assign_BeamLineElem(beamln13(itws))
          else if(bitype(i).eq.110) then
            iemfld = iemfld + 1
            call construct_EMfld(beamln11(iemfld),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmp13(1) = 0.0
            tmp13(2) = val1(i) 
            tmp13(3) = val2(i) 
            tmp13(4) = val3(i) 
            tmp13(5) = val4(i) 
            tmp13(6) = val5(i) 
            tmp13(7) = val6(i) 
            tmp13(8) = val7(i) 
            tmp13(9) = val8(i) 
            tmp13(10) = val9(i) 
            tmp13(11) = val10(i) 
            tmp13(12) = val11(i) 
            tmp13(13) = val12(i) 
            tmp13(14) = val13(i) 
            call setparam_EMfld(beamln11(iemfld),tmp13)
            Blnelem(i) = assign_BeamLineElem(beamln11(iemfld))
          else
          endif 
        enddo
!-------------------------------------------------------------------
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid.eq.0) print*,"pass setting up lattice..."

        deallocate(blength,bnseg,bmpstp,bitype)
        deallocate(val1,val2,val3,val4,val5,val6,val7,val8,val9)
        deallocate(val10,val11,val12,val13,val14,val15,val16)
        deallocate(val17,val18,val19,val20,val21,val22,val23)
        deallocate(val24)

        t_init = t_init + elapsedtime_Timer(t0)

        end subroutine init_AccSimulator

        !Run beam dynamics simulation through accelerator.
        subroutine run_AccSimulator()
        !implicit none
        include 'impli.inc'
        include 'mpif.h'
        !from MaryLie
        include 'map.inc' !contains total transfer map (not used), reftraj, arcle
        !....
        integer :: i,j,bnseg,bmpstp,bitype
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,npx,npy,&
                   totnp,ierr,nylcr,nzlcr
        integer :: nbal,ibal,nstep,ifile
        integer :: tmpfile,nfile
        integer, dimension(3) :: lcgrid
        integer, allocatable, dimension(:) :: lctabnmx,lctabnmy,ytable,&
                                              ydisp,ztable,zdisp
        integer, allocatable, dimension(:,:,:) :: temptab
        double precision :: z,tau1,tau2,blength,t0
        double precision, allocatable, dimension(:,:) :: lctabrgx, lctabrgy
        double precision, dimension(6) :: lcrange, range, ptrange,ptref
        double precision, dimension(3) :: msize
        double precision :: hy,hz,ymin,zmin,piperad,zedge,circumference
        double precision :: tmp1,tmp2,tmp3,tmp4,rfile
        double precision, allocatable, dimension(:,:,:) :: chgdens
        double precision, allocatable, dimension(:,:,:) :: besscoef
        double precision, allocatable, dimension(:,:) :: bessnorm,gml
        integer, allocatable, dimension(:) :: modth,pydisp
        integer :: nmod,k,ii,jj
        !double precision :: sumtest, sumtest2, sumtest3
        double precision, dimension(8) :: drange
        double precision, dimension(3) :: al0,ga0,epson0
        double precision :: tg,tv,gam,piperad2!,piperad_x, piperad_y
        integer :: nsubstep,Flagbctmp
        double precision :: zz,vref
        !parameters for stripper modeling
        double precision :: beta0,gamma0,gambetz
        double precision, allocatable, dimension(:,:) :: tmpptcs
        double precision :: rwkinq,rwkinf,avgw
        integer :: itmp,jtmp
        !for WAKEFIELD purpose -----------------------------------------
        !double precision, allocatable, dimension(:,:,:) :: exg,eyg,ezg
        !double precision, allocatable, dimension(:) :: ztable,zdisp,&
        double precision, allocatable, dimension(:) :: &
            denszlc,densz,exwake,eywake,ezwake,xwakelc,xwakez,ywakelc,&
            ywakez,denszp,denszpp,csg,ans
        double precision, allocatable, dimension(:,:) :: sendensz,&
                                                         recvdensz
        double precision :: xx,yy,t3dstart,rr,tmger,tmpwk
        double precision  :: aawk,ggwk,lengwk,hzwake,ab,tbtwstart
        integer :: flagwake,iizz,iizz1,kz,kadd,ipt,flagcsr,flagbtw
        !for bending magnet Transport transfer matrix implementation
        double precision :: hd0,hd1,dstr1,dstr2,angF,tanphiF,tanphiFb,&
            angB,tanphiB,tanphiBb,hF,hB,qm0,qmi,psi1,psi2,r0,gamn,gambet,&
            angz,dpi,ang0,hgap,betai,rcpgammai,ssll
        double precision, dimension(6) :: ptarry
        double precision, dimension(10) :: dparam
        real*8 :: tmpgamlc,tmpgamgl
        real*8, dimension(2) :: tmp56lc,tmp56gl
        real*8 :: tmppx,tmppy,tmppt,tmph
!for MaryLie S-bend ---------------------
! pp is an array that contains the sbend parameters:
        dimension pp(21)
! xmh is the 6x6 map, h contains the nonlinear map (Lie polynomial coefficients)
        dimension xmh(6,6),h(monoms)
        real*8, dimension(6) :: tmp6
        integer :: ntrace,ntaysym,norder
        !<<<<<<<<<<<<<<<<<<<<<<<<< kilean <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        logical :: pipe_override, flagExternalPipe
        integer :: ihalf, jslice, nslices, jadd, jend1, kend1, pipeID, nlost
        double precision, allocatable :: lost_pdata(:,:)
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !integer :: ihalf,jslice,nslices
        !real*8 :: slfrac,pp,reftraj,tmh,th,angle,slen,arclen
        double precision, allocatable, dimension(:,:):: tmpPts
!for csr wake ------------------
        integer :: bitypeold 
        real*8 :: bendlen,zz1,zz2,zwkmin,zfmin,zfmax,blengthold,zbleng
        integer :: ldsg,nlsg,nrsg,npsg,msg,ncoefreal
        integer :: flagcsrTr
!for lumped space charge --------------------------------
        !flagtmp is set by "-15" to turn on(-)/off(+) the SC for ordinary elements
        !other than "-14".
        !flagcoll is set by "-14" for lumped SC or Wake or both
        !depending on the input parameters
        integer :: flagcoll,flagtmp
        real*8 :: tmplump,scwk
! ----------------------------------------------
! for rescaling the uncorrelated energy spread
!        integer :: Nplagrange
!        real*8 ::fract1,fract2,dEscale
!        real*8, allocatable, dimension(:) :: tphase,Ptenergy
!----------------------------------------------
!for global decomposition
        real*8, allocatable, dimension(:) :: tmppot,recvtmppot,gltmppot
        real*8, allocatable, dimension(:,:,:) :: glpot
        integer :: nsendy,nsendz,nsendz2,kk,ijk,flagdecomp
        integer :: k0,j0,ijk0,ijk00,nxylc
!------------------------------------------------
!for laser heater
        real*8 :: rklaser,a0laser,csiglaser,hzl,hrl,ezlaser,ss,cd,betz,&
                  ezamp
        integer :: ir,ir1,iz,iz1
!------------------------------------------------
!for longitudinal mesh adjustment
        real*8 :: zrandi,zadjust,zadjmax
!for longitudinal energy shift before dipole
        real*8 :: y1,t1,c1,pz0lc,pz0 
!--------------------------------
!for ring model.
!whenever we use wake field, we need to insert "-41" element before
!for read in discrete wake function.
!flagbtw = 1 - both transverse and longitude wake, 2 or 3 - only longitude
!          4 - only transverse  
!flagcoll = 1 whenever sc or wake is used
!-15 is used to set flagcol for none "-14" lump element
        integer :: Nturn,iturn,flagsc
        real*8 :: qmass 
        integer :: iend0,jend0,kend0
        real*8 :: t11,t2,t3,t4,t5,t6

        qmass = Bpts%Charge/Bpts%Mass

        !zadjmax = 0.15 !30% increase of z domain
        zadjmax = 0.0d0 !0% increase of z domain

!-------------------------------------------
! initialization: MaryLie S-bend
! note: since this code does not involve particles, the only need for
! module parallel is to specify idproc (needed to print from just proc 0)
!      call init_parallel
! setup intializes Lie algebraic routines for map concatenation, etc.
       call setup
!--------------------------------

!-------------------------------------------------------------------
! prepare initial parameters, allocate temporary array.
        flagbtw = 0 !initial no BTW
        !iend = 0
        !ibalend = 0
        !nstepend = 0
        !zend = 0.0

        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        idproc = myid
        call getcomm_Pgrid2d(grid2d,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid2d,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        allocate(lctabnmx(0:npx-1))
        allocate(ztable(0:npx-1))
        allocate(zdisp(0:npx-1))
        allocate(lctabnmy(0:npy-1))
        allocate(ytable(0:npy-1))
        allocate(ydisp(0:npy-1))
        allocate(lctabrgx(2,0:npx-1))
        allocate(lctabrgy(2,0:npy-1))
        allocate(temptab(2,0:npx-1,0:npy-1))
        allocate(tmpptcs(1,1))
        allocate(tmpPts(6,Nplocal))

        nbal = 5
        ibal = ibalend
        nstep = nstepend
        z = zend

        if(Rstartflg.eq.1) then
          call restart_AccSimulator()
        else
          if(Flagdiag.eq.1) then
            call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
          else
            call diagnostic2_Output(Bpts,z,nchrg,nptlist0)
          endif
        endif

        flagdecomp = 1

        iend0 = 0
        jend0 = 0
        kend0 = 0
        jend1 = 0
        kend1 = 0
        if((Flagbc.eq.4) .or. (Flagbc.eq.6)) then !(odd,odd,odd)
               iend0 = 1
               if((myidy.eq.npy-1) .and. (myidx.eq.npx-1)) then
                 jend0 = 1
                 kend0 = 1
               else if(myidy.eq.npy-1) then
                 jend0 = 1
               else if (myidx.eq.npx-1) then
                 kend0 = 1
               else
               endif
               jend1 = 1
               kend1 = 1
        else if((Flagbc.eq.3) .or. (Flagbc.eq.5)) then !(odd,odd,even)
               iend0 = 1  
               if(myidy.eq.npy-1) then
                   jend0 = 1
               else
               endif
               jend1 = 1
        else if(Flagbc.eq.2) then !(even,even,odd)
               if(myidx.eq.npx-1) then
                   kend0 = 1
               else
               endif
               kend1 = 1
        else !(even,even,even)
        endif

        call getlcmnum_CompDom(Ageom,lcgrid)
        Nxlocal = lcgrid(1)
        if(npy.gt.1) then
          Nylocal = lcgrid(2) + 2
        else
          Nylocal = lcgrid(2)
        endif
        if(npx.gt.1) then
          Nzlocal = lcgrid(3) + 2
        else
          Nzlocal = lcgrid(3)
        endif
        if(flagdecomp.eq.1) then
          allocate(chgdens(1,1,1))
          allocate(glpot(1,1,1))
        else
          allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
          allocate(glpot(Nx,Ny,Nz))
        endif

        if(flagdecomp.eq.1) then
          allocate(tmppot(1))
          allocate(recvtmppot(1))
          allocate(gltmppot(1))
        else
          nsendy = (lcgrid(1)-iend0)*(lcgrid(2)-jend0)*(lcgrid(3)-kend0)
          allocate(tmppot(nsendy))
          nsendz = (lcgrid(1)-iend0)*(Ny-jend1)*(lcgrid(3)-kend0)
          allocate(recvtmppot(nsendz))
          nsendz2 = (lcgrid(1)-iend0)*(Ny-jend1)*(Nz-kend1)
          allocate(gltmppot(nsendz2))
        endif
!-------------------------------------------------------------------
! prepare for round pipe, open longitudinal
        if(Flagbc.eq.3) then
          allocate(pydisp(0:npy-1))
          call getlctabnm_CompDom(Ageom,temptab)
          pydisp(0) = 0
          do i = 1, npy-1
            pydisp(i) = pydisp(i-1) + temptab(2,0,i-1)
          enddo
          call getlcmnum_CompDom(Ageom,lcgrid)
          Nxlocal = lcgrid(1)
          if(npcol.gt.1) then
            Nylocal = lcgrid(2) + 2
          else
            Nylocal = lcgrid(2)
          endif
          if(nprow.gt.1) then
            Nzlocal = lcgrid(3) + 2
          else
            Nzlocal = lcgrid(3)
          endif

          if(nprow.gt.1) then
            nzlcr = Nzlocal-2
          else
            nzlcr = Nzlocal
          endif
          if(npcol.gt.1) then
            nylcr = Nylocal-2
          else
            nylcr = Nylocal
          endif
          if(myidy.eq.(npcol-1)) nylcr = nylcr - 1
          allocate(modth(nylcr))
          do i = 1, nylcr
            modth(i) = (pydisp(myidy)+i-2+1)/2
          enddo
          if(myidy.eq.0) then
            modth(1) = 0
            modth(2) = (ny-1)/2
            nmod = modth(nylcr) - modth(1) + 2
          else
          nmod = modth(nylcr) - modth(1) + 1
          endif
          allocate(besscoef(lcgrid(1),lcgrid(1),nmod))
          allocate(bessnorm(lcgrid(1),nmod))
          allocate(gml(lcgrid(1),nmod))

          call getmsize_CompDom(Ageom,msize)
          call Bessprep_Bessel(msize(1),lcgrid(1),nylcr,nmod,modth,&
                        besscoef,bessnorm,gml)
        endif

        allocate(denszlc(Nz))
        allocate(densz(Nz))
        allocate(denszp(Nz))
        allocate(denszpp(Nz))
        allocate(csg(Nz))
        allocate(ans(Nz))
        allocate(xwakelc(Nz))
        allocate(xwakez(Nz))
        allocate(ywakelc(Nz))
        allocate(ywakez(Nz))
        allocate(recvdensz(Nz,2))
        allocate(sendensz(Nz,2))
        allocate(exwake(Nz))
        allocate(eywake(Nz))
        allocate(ezwake(Nz))
        exwake = 0.0
        eywake = 0.0
        ezwake = 0.0
        flagwake = 0
        flagcsr = 0
        flagcsrTr = 0
        aawk = 0.05
        ggwk = 0.05
        lengwk = 0.1
        scwk = 1.0
        flagsc = 1

        densz = 0.0
        denslc = 0.0
        exwake = 0.0
        eywake = 0.0
        ezwake = 0.0
        xwakelc = 0.0
        xwakez = 0.0
        ywakelc = 0.0
        ywakez = 0.0

        Nturn = 100000
        !set initial flagcoll
        flagcoll = 1
        flagtmp = -1
        do i = iend+1, Nblem
          call getparam_BeamLineElem(Blnelem(i),blength,bnseg,bmpstp,&
                                     bitype)
        !The flag to set the ring simulation parameters
          if(bitype.eq.-16) then !ring simulation
            call getparam_BeamLineElem(Blnelem(i),drange)
            Nturn = drange(3)+0.1
          endif
          if(bitype.eq.-17) then !ring simulation
            call readPipeInfo()
            flagExternalPipe = .True.
          endif

        enddo
      
      !<<<<<<<<<<<<<<<<< prepare particle lost (Kilean) <<<<<<<<<<<<<<<<
      nlost = 0
      pipe_override = .false.
      pipeID = 1
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !-----------------------------------------------------------------
      !start looping through 'Nturn'
      do iturn = 1, Nturn  
        if(iturn == 2) then
          circumference = z
        endif
        tmpfile = 0
        bitypeold = 0
        blengthold = 0.0d0
        zbleng = 0.0d0
        !---------------------------------------------------------------
        ! start looping through 'Nblem' beam line elements.
        do i = iend+1, Nblem
          call starttime_Timer(t11)
          call getparam_BeamLineElem(Blnelem(i),blength,bnseg,bmpstp,&
                                     bitype)
          if(.not. pipe_override) call getradius_BeamLineElem(Blnelem(i),piperad,piperad2)
          nsubstep = bmpstp
          if(myid.eq.0) print*,"enter elment: ",i,bitype
          nfile = 0
          tau1 = 0.0
          if(bitype.ge.0) tau1 = 0.5*blength/bnseg
          if(bitype.ne.-14) then
            tau2 = 2.0*tau1
            flagwake = 0
            if(flagtmp.gt.0) then
              flagcoll = 0
            else
              flagcoll = 1
            endif
            scwk = 1.0d0
          else !instanteous lumped space-charge and wake field kick
            call getparam_BeamLineElem(Blnelem(i),drange)
            tau2 = drange(3)
            aawk = drange(4)
            ggwk = drange(5)
            lengwk = drange(6)
            scwk = drange(7)
             
            flagcoll = 1

            if(lengwk.gt.0.0) then !flag for lumped wake
              flagwake = 1
              if(aawk.gt.100) then
                flagbtw = 1
              else if(aawk.lt.0.0) then
                if(aawk.gt.-10.0) then
                  flagbtw = 2
                else
                  flagbtw = 3
                endif
              else
                flagbtw = 4
              endif
            else
              flagwake = 0
            endif

            if(ggwk.gt.0.0) then !flag for sc
              flagsc = 1
            else
              flagsc = 0
            endif

          endif
          if(bitype.eq.-15) then
            call getparam_BeamLineElem(Blnelem(i),3,tmplump)
            if(tmplump.gt.0.0) then
              flagtmp = 1
            else
              flagtmp = -1
            endif
          endif

!-------------------------------------------------------------------
! read in the on axis E field for rf cavities.
          if(bitype.gt.100) then
            call getparam_BeamLineElem(Blnelem(i),5,rfile)
            nfile = int(rfile + 0.1)
            ifile = nfile
            if(ifile.ne.tmpfile)then
              !for linear map integrator
              if(Flagmap.eq.1) then
                !read in Ez, Ez', Ez'' on the axis
                !for SolRF and TWS, use Fourier coefficients
                !otherwise, use discrete data on axis.
                if(bitype.eq.105 .or. bitype.eq.106) then
                  call read3_Data(ifile)
                else
                  call read1_Data(ifile)
                endif
              else
                !read in Er, Ez, H_theta on r-z grid 
                !call read2_Data(ifile)
                !read in Fourier coefficients of Ez on the axis
                !call read3_Data(ifile)
                if(bitype.eq.110) then
                  call read4_Data(ifile)
                  !call read3_Data(ifile)
                else
                  call read3_Data(ifile)
                endif
              endif
              tmpfile=ifile
            endif
          endif

          if(bitype.eq.1) then
            call getparam_BeamLineElem(Blnelem(i),3,rfile)
            if(rfile.gt.1.0e-5) then
               nfile = int(rfile+0.1)
               ifile = nfile
               if(ifile.ne.tmpfile) then
                  call read1_Data(ifile)
                  tmpfile = ifile
               endif
            endif
          endif

          if(bitype.eq.105 .or. bitype.eq.106) then
            call getparam_BeamLineElem(Blnelem(i),13,aawk)
            call getparam_BeamLineElem(Blnelem(i),14,ggwk)
            call getparam_BeamLineElem(Blnelem(i),15,lengwk)
            if(lengwk.gt.0.0 .and. flagcoll.eq.1) then
              flagwake = 1
              if(aawk.gt.100) then
                flagbtw = 1
              else if(aawk.lt.0.0) then
                if(aawk.gt.-10) then
                  flagbtw = 2
                else
                  flagbtw = 3
                endif
              endif
            else
              flagwake = 0
            endif
          else
            if(bitype.ne.-14) then
              flagwake = 0
            endif
          endif

!-------------------------------------------------------------------
! print out beam information using BPM 
          if(bitype.eq.-1) then
            call shift_BPM(Bpts%Pts1,bitype,Nplocal,Np)
          endif
!<<<<<<<<<<<<<<<<<<<<<< TBToutput(Kilean) <<<<<<<<<<<<<<<<<<<<<<<<<<<
          if(bitype.eq.-90) then
            call getparam_BeamLineElem(Blnelem(i),dparam)
            call turn_by_turn_phasespace_split(Bpts,bmpstp,int(dparam(2),8),int(dparam(3),8),int(dparam(4)))
          endif
          if(bitype.eq.-89) then
            call getparam_BeamLineElem(Blnelem(i),dparam)
            call turn_by_turn_phasespace(Bpts,bmpstp,int(dparam(2),8),int(dparam(3),8))
          endif
          if(bitype.eq.-88) then
            call getparam_BeamLineElem(Blnelem(i),dparam)
            call turn_by_turn_integral(Bpts,bmpstp,dparam(2),dparam(3),dparam(4)&
                                      ,dparam(5),int(dparam(6),8),int(dparam(7),8))
          endif
          if(bitype.eq.-87) then
            call getparam_BeamLineElem(Blnelem(i),dparam)
            call turn_by_turn_integral_on_momentum(Bpts,bmpstp,dparam(2),dparam(3),dparam(4)&
                                      ,dparam(5),int(dparam(6),8),int(dparam(7),8))
          endif
!>>>>>>>>>>>>>>>>>>>>>>> end of TBToutput >>>>>>>>>>>>>>>>>>>>>>>>>>>
          if(bitype.eq.-2) then
            !call phase_Output(bmpstp,Bpts)
            !call phaseleda_Output(bmpstp,Bpts)
            call getparam_BeamLineElem(Blnelem(i),dparam)
            if(int(dparam(3))==iturn) then
              call phase_Output(Bpts,int(dparam(2)),bmpstp,iturn,int(dparam(4)))
            endif
          else if(bitype.eq.-3) then
            call getparam_BeamLineElem(Blnelem(i),drange)
            call accdens1d_Output(nstep,8,Bpts,Np,drange(2),-drange(3),&
            drange(3),-drange(5),drange(5))
          else if(bitype.eq.-4) then
            call getparam_BeamLineElem(Blnelem(i),drange)
            call dens1d_Output(nstep,8,Bpts,Np,drange(2),-drange(3),&
            drange(3),-drange(5),drange(5))
          else if(bitype.eq.-5) then
            call getparam_BeamLineElem(Blnelem(i),drange)
            call dens2d_Output(nstep,8,Bpts,Np,-drange(3),drange(3),&
            -drange(4),drange(4),-drange(5),drange(5),-drange(6),drange(6),&
            -drange(7),drange(7),-drange(8),drange(8))
          else if(bitype.eq.-6) then
            call dens3d_Output(nstep,8,Bpts,Np,-drange(3),drange(3),&
              -drange(5),drange(5),-drange(7),drange(7))
          else if(bitype.eq.-7) then
            !output all particles in 6d phase space in file "xxnstep"
            !call phaseout_Output(nstep,Bpts)
            !output all geomtry information in file "xxnstep".
            !call outgeom_Output(nstep,z,i,j,npx,npy,Ageom)
            !call outpoint_Output(myid+31,Bpts,z,i,j,npx,npy,Ageom)
            call outpoint_Output(myid+51,Bpts,z,i,j,ibal,nstep,npx,npy,Ageom)
          else if(bitype.eq.-9) then
            call getparam_BeamLineElem(Blnelem(i),dparam)
            pipe_override = .true.
            pipeID  = int(dparam(2))
            piperad = dparam(3)
            piperad2 = dparam(4)
            cycle
          else if(bitype.eq.-10) then
            !mismatch the beam at given location.
            !here, drange(3:8) stores the mismatch factor.
            call getparam_BeamLineElem(Blnelem(i),drange)
            call scale_BPM(Bpts%Pts1,Nplocal,drange(3),&
            drange(4),drange(5),drange(6),drange(7),drange(8))
          else if(bitype.eq.-11) then
            deallocate(tmpptcs)
            allocate(tmpptcs(8,Nplocal))
            gam = -Bpts%refptcl(6)
            nfile = bmpstp
            beta0 = sqrt(1.0-(1.0/gam)**2)
            do ii = 1, Nplocal
              tmpptcs(1,ii) = Bpts%Pts1(9,ii)
              tmpptcs(2,ii) = Bpts%Pts1(1,ii)*Scxl*100
              tmpptcs(3,ii) = Bpts%Pts1(3,ii)*Scxl*100
              tmpptcs(4,ii) = Bpts%Pts1(5,ii)*90/asin(1.0)
              gambetz = sqrt((gam-Bpts%Pts1(6,ii))**2-Bpts%Pts1(2,ii)**2-&
                             Bpts%Pts1(4,ii)**2-1.0)
              tmpptcs(5,ii) = Bpts%Pts1(2,ii)/gambetz * 1000
              tmpptcs(6,ii) = Bpts%Pts1(4,ii)/gambetz * 1000
              tmpptcs(7,ii) = -Bpts%Pts1(6,ii)/(gam-1)*100
              !tmpptcs(8,ii) = Bpts%Pts1(7,ii)*Bpts%mass
              tmpptcs(8,ii) = Bpts%Pts1(7,ii)*931.494326d6
            enddo
            deallocate(nptlist0)
            deallocate(currlist0)
            deallocate(qmcclist0)
            allocate(nptlist0(nchrg))
            allocate(currlist0(nchrg))
            allocate(qmcclist0(nchrg))
            do ii = 1, nchrg
              nptlist0(ii) =  nptlist(ii)
              currlist0(ii) =  currlist(ii)
!              qmcclist0(ii) =  qmcclist(ii)/Bpts%mass
              qmcclist0(ii) =  qmcclist(ii)/931.49326d6
            enddo

            !gam = (1.0+rwkinq/Bmass)
            gam = (1.0+rwkinq/931.49326d6)
            gambetz = sqrt(gam**2-1.0)
            avgw = 0.0
            do ii = 1, Nplocal
               Bpts%Pts1(6,ii) = -(gam-1)*tmpptcs(7,ii)/100
               gambetz = sqrt((gam-Bpts%Pts1(6,ii))**2-1.0) /  &
               sqrt(1.0+(tmpptcs(5,ii)/1000)**2+(tmpptcs(6,ii)/1000)**2)
               Bpts%Pts1(2,ii) = tmpptcs(5,ii)/1000*gambetz
               Bpts%Pts1(4,ii) = tmpptcs(6,ii)/1000*gambetz
!               Bpts%Pts1(7,ii) = tmpptcs(8,ii)/Bpts%mass
               Bpts%Pts1(7,ii) = tmpptcs(8,ii)/931.49326d6
               avgw = avgw + tmpptcs(7,ii)
            enddo
            avgw = avgw/Nplocal
            call chgupdate_BeamBunch(Bpts,nchrg,nptlist0,qmcclist0)
            !update reference particle information
            !Bcharge = Bmass*qmcclist0(1)
            Bcharge = 931.49326d6*qmcclist0(1)
            Bpts%Charge = Bcharge
            !Bpts%refptcl(6) = -(1.0+rwkinq/Bmass)
            Bpts%refptcl(6) = -(1.0+rwkinq/931.49326d6)
          else if(bitype.eq.-21)then
            !shift the beam centroid in the 6D phase space.
            !This element can be used to model steering magnets etc.
            !here, drange(3:8) stores the amount of shift.
            !drange(3); shift in x (m)
            !drange(4); shift in Px (rad)
            !drange(5); shift in y (m)
            !drange(6); shift in Py (rad)
            !drange(7); shift in z (deg)
            !drange(8); shift in Pz (MeV)
            call getparam_BeamLineElem(Blnelem(i),drange)
            call kick_BPM(Bpts%Pts1,Nplocal,drange(3),drange(4),drange(5),&
            drange(6),drange(7),drange(8),-Bpts%refptcl(6),Bpts%Mass)
          else if(bitype.eq.-30) then
            !rescale the uncorrelated energy spread by a factor of dEscale
            !number of section points for Lagrange fitting,typically between 8 and 10
!            Nplagrange = drange(3) + 0.01 
!            !fraction of longitudinal length used for fitting
!            fract1 = drange(4) 
!            !fraction of longitudinal length used for computing rms information
!            fract2 = drange(5)
!            !scaling factor
!            dEscale = drange(6)
!            allocate(tphase(Nplocal))
!            tphase(:) = Bpts%Pts1(:,5)
!            allocate(Ptenergy(Nplocal))
!            Ptenergy(:) = Bpts%Pts1(:,6)
!            call rescaler(tphase,Ptenergy,dble(Nplagrange),fract1,fract2,dEscale)
!            deallocate(tphase)
!            deallocate(Ptenergy)
          else if(bitype.eq.-40)then
            call getparam_BeamLineElem(Blnelem(i),drange)
            call kickRF_BPM(Bpts%Pts1,Nplocal,drange(3),drange(4),drange(5),&
                            Bpts%Mass)
          !read in discrete wakefield
          else if(bitype.eq.-41)then
            call getparam_BeamLineElem(Blnelem(i),3,rfile)
            ifile = int(rfile + 0.1)
            call read1wk_Data(ifile)
          else if(bitype.eq.-42)then !instant energy kick in ring
            call getparam_BeamLineElem(Blnelem(i),drange)
            call kickRF2_BPM(Bpts%Pts1,Nplocal,drange(3),drange(4),drange(5),&
                            Bpts%Mass,Bpts%refptcl)
          !laser heater
          else if(bitype.eq.-45)then
            !Apply zero-length thin lens focusing kick in x,y.
            !here, drange(3:8) stores the amount of shift.
            !drange(3); focusing strength in x (1/m)
            !drange(4); focusing strength in y (1/m)
            call getparam_BeamLineElem(Blnelem(i),drange)
            call kick_focusing(Bpts%Pts1,Nplocal,drange(3),drange(4),&
                         -Bpts%refptcl(6),Bpts%Mass)
          else if(bitype.eq.-46)then
            !Apply zero-length 4D phase advance rotation 
            !here, drange(3:8) stores the amount of shift.
            !drange(3); length of nonlinear insert (m)
            !drange(4); phase advance/2pi across NLI 
            !drange(5); phase advance/2pi across arc
            call getparam_BeamLineElem(Blnelem(i),drange)
            call kick_phaseadvance(Bpts%Pts1,Nplocal,drange(3),drange(4),&
                         drange(5),-Bpts%refptcl(6),Bpts%Mass)
          else if(bitype.eq.-47)then
            !Apply zero-length 4D phase advance rotation
            !here, drange(3:8) stores the amount of shift.
            !drange(3); length of nonlinear insert (m)
            !drange(4); phase advance/2pi across NLI
            !drange(5); phase advance/2pi across arc
            !drange(6); dimensionless strength of nonlinear insert
            call getparam_BeamLineElem(Blnelem(i),drange)
            call kick_linDNadvance(Bpts%Pts1,Nplocal,drange(3),drange(4),&
                         drange(5),drange(6),-Bpts%refptcl(6),Bpts%Mass)
          else if(bitype.eq.-50)then
            call getparam_BeamLineElem(Blnelem(i),3,rfile)
            call getparam_BeamLineElem(Blnelem(i),4,a0laser)
            call getparam_BeamLineElem(Blnelem(i),5,csiglaser)
            call getparam_BeamLineElem(Blnelem(i),6,rklaser)
            call getparam_BeamLineElem(Blnelem(i),7,ezamp)
            ifile = int(rfile + 0.1)
            call read5_Data(ifile)
            hrl = (RmaxRf - RminRf)/(NrIntvRf-1) 
            hzl = (ZmaxRf - ZminRf)/(NzIntvRf-1) 
            do ii = 1, Nplocal
              rr = sqrt(Bpts%Pts1(1,ii)**2+Bpts%Pts1(3,ii)**2)*Scxl/a0laser
              gam = -Bpts%refptcl(6)-Bpts%Pts1(6,ii)
              betz = sqrt(gam**2-1.0d0-Bpts%Pts1(2,ii)**2-Bpts%Pts1(4,ii)**2)/gam
              ss = -Bpts%Pts1(5,ii)*betz*Scxl/csiglaser
              ir = (rr-RminRf)/hrl + 1
              ir1 = ir + 1 
              if(ir1.le.NrIntvRf) then
                iz = (ss-ZminRf)/hzl + 1
                iz1 = iz + 1
                ab=(ir*ir*hrl*hrl-(rr-RminRf)*(rr-RminRf))/ &
                 (ir*ir*hrl*hrl-(ir-1)*(ir-1)*hrl*hrl)
                cd = (iz*hzl-(ss-ZminRf))/hzl
                ezlaser = ezdata(ir,iz)*ab*cd+ezdata(ir1,iz)*(1.0d0-ab)*cd+&
                        ezdata(ir,iz1)*ab*(1.0d0-cd)+&
                        ezdata(ir1,iz1)*(1.0d0-ab)*(1.0d0-cd)
                Bpts%Pts1(6,ii) = Bpts%Pts1(6,ii) - &
                                ezamp*ezlaser*cos(rklaser*ss*csiglaser)
              endif
            enddo
          else if(bitype.eq.5)then
           qmass = Bpts%Charge/Bpts%Mass
           call kickmultthinK(Blnelem(i)%pmult,Bpts%refptcl,&
                Bpts%Pts1,Nplocal,qmass)
          elseif(bitype.eq.-99) then
            exit
          endif


!-------------------------------------------------------------------
! loop through 'bnseg' numerical segments in each beam element
! using 2 step symplectic integeration (ie. leap frog).
          zedge = z
          call setparam_BeamLineElem(Blnelem(i),1,zedge)
          if(bitype.eq.4) then
            call getparam_BeamLineElem(Blnelem(i),dparam)
            !add the tranisent drift effects
            if((BcurrImp.gt.0.0) .and. (dparam(4).gt.500)) then 
              flagcsrTr = 1
            else
              flagcsrTr = 0
            endif
          endif
          !check whether the preceding beam line elment is bend or not
!comment out for test purpose
          if((bitype.eq.0) .and. (bitypeold.eq.4) .and. &
             (flagcsrTr.eq.1)) then
            flagcsr = 1
          endif

          !use Z as independent variable for no bend magnet or bend using transfer map
          !//no bend or bend using Transport transfer map
          !-----------------------------------------------------------------
          if((bitype.ne.4).or.(bitype.eq.4 .and. dparam(4).gt.100) ) then  

            !use TRANSPORT bend transfer map
            if(bitype.eq.4 .and. dparam(4).le.300) then
              !transfter to the Transport coordinate (except the 5th coordinate 
              !is -v dt instead of v dt) and apply front edge transfer map
              dpi = 2*asin(1.0)

              gamma0 = -Bpts%refptcl(6)
              
              beta0 = sqrt(1.0-1.0/gamma0/gamma0)
              Bpts%refptcl(5) = Bpts%refptcl(5)+blength/(Scxl*beta0)
              gambet = beta0*gamma0
              hgap = 2*dparam(5)
              ang0 = dparam(2)
              hd1 = dparam(3) !k1
              angF = dparam(6) !e1
              angB = dparam(7) !e2
              hF = dparam(8) !1/r1
              hB = dparam(9) !/1/r2
              dstr1 = dparam(10) !fringe field K of entrance
              dstr2 = dstr1 !fringe field K of exit. here, assume Kb = Kf

              hd0 = ang0/blength !k0
              tanphiF = tan(angF)
              psi1 = hd0*hgap*dstr1*(1.0+sin(angF)*sin(angF))/cos(angF)
              tanphiFb = tan(angF-psi1)
              tanphiB = tan(angB)
              psi2 = hd0*hgap*dstr2*(1.0+sin(angB)*sin(angB))/cos(angB)
              tanphiBb = tan(angB-psi2)
              qm0 = Bpts%Charge/Bpts%Mass
              r0  = abs(1.0/hd0)

              do ipt = 1, Nplocal
                ptarry(1) = Bpts%Pts1(1,ipt)*Scxl
                gamn = gamma0 - Bpts%Pts1(6,ipt)
                gambetz = sqrt(gamn**2-1.0-Bpts%Pts1(2,ipt)**2-&
                               Bpts%Pts1(4,ipt)**2)
                ptarry(2) = Bpts%Pts1(2,ipt)/gambetz
                ptarry(3) = Bpts%Pts1(3,ipt)*Scxl
                ptarry(4) = Bpts%Pts1(4,ipt)/gambetz
                ptarry(5) = -Bpts%Pts1(5,ipt)*beta0*Scxl
                ptarry(6) = -Bpts%Pts1(6,ipt)/beta0/gambet - &
                            (Bpts%Pts1(7,ipt)-qm0)/qm0
                Bpts%Pts1(1,ipt) = ptarry(1)
                Bpts%Pts1(2,ipt) = ptarry(2)
                Bpts%Pts1(3,ipt) = ptarry(3)
                Bpts%Pts1(4,ipt) = ptarry(4)
                Bpts%Pts1(5,ipt) = ptarry(5)
                Bpts%Pts1(6,ipt) = ptarry(6)
              enddo
              call Fpol_Dipole(hd0,hF,tanphiF,tanphiFb,hd1,&
                                 psi1,Bpts%Pts1,angF,Nplocal)
              angz = 0.0
            endif
! The following uses S-bend from MaryLie
!      
            if(bitype.eq.4 .and. dparam(4).gt.300) then

              gamma0 = -Bpts%refptcl(6)
              beta0 = sqrt(1.0-1.0/gamma0/gamma0)
              Bpts%refptcl(5) = Bpts%refptcl(5)+blength/(Scxl*beta0)

              bcurr=BcurrImp !not needed here, needed later (but this is in module beamdata)
              achg=1.d0     !not used
              pmass=Bmass !in eV/c^2
              bfreq=BfreqImp
              gamma = -Bpts%refptcl(6)
              beta=sqrt((gamma+1.d0)*(gamma-1.d0))/gamma
              brho=gamma*beta/c*pmass
!     
! scaling parameters:
              p0sc=pmass/c    !"dynamic" units (scale momenta by mc)
              freqscl=BfreqImp
              omegascl=4.d0*asin(1.d0)*freqscl
              ts=1.d0/omegascl  !not used
              sl=c/omegascl     !IMPACT units
!
! initialize the reference trajectory (not needed to compute bend map, but
! included here for completeness):
              reftraj(1:5)=0.d0
              reftraj(6)=-gamma*pmass/(omegascl*sl*p0sc)
              arclen=0.d0
!
! example bend parameters:
! B_BC1.1: SBEND,L=0.5,ANGLE=0.06774,E2=0.06774  !these angles are in radians
! B_BC1.2: SBEND,L=0.5,ANGLE=-0.06774,E1=-0.06774  !these angles are in radians
! B_BC1.3: SBEND,L=0.5,ANGLE=-0.06774,E2=-0.06774  !these angles are in radians
! B_BC1.4: SBEND,L=0.5,ANGLE=0.06774,E1=0.06774  !these angles are in radians
!
              slen=blength
              angle=dparam(2)  !angle in radians     [but pp(2) is in degrees]
              r0 = abs(slen/angle)
              pp(1)=angle*90.d0/asin(1.d0)  !bend angle (deg)
              pp(2)=brho*angle/slen  !B field (Tesla)
              pp(3)=dparam(6)*90.d0/asin(1.d0)      !e1 entry pole face rotation angle (deg)
              pp(4)=dparam(7)*90.d0/asin(1.d0) !e2 exit  pole face rotation angle (deg)
              pp(5)=3.d0  !entry fringe flag (0=none,1=quad fringe,2=dipole frng,3=both)
              pp(6)=3.d0  !exit  fringe flag
              pp(7)=dparam(5)  !gap size for entry fringe field calculation
              pp(8)=dparam(5)  !gap size for exit  fringe field calculation
              pp(9)=dparam(10)  !normalized field integral for entry fringe field
              pp(10)=dparam(10) !normalized field integral for exit  fringe field
              pp(11)=3.d0 !iopt (controls handling of multipole coeffs;3=same as MAD)
              pp(12)=0.d0 !ipset (legacy MaryLie option to specificy coeffs in a pset)
              !pp(13)=0.d0 ! (13-18)=
              pp(13)=dparam(3) ! (13-18)=
              pp(14)=0.d0 ! BQD,AQD,BSEX,ASEX,BOCT,AOCT  (if iopt=1)
              pp(15)=0.d0 ! Tay1,AQD,Tay2,ASEX,Tay3,AOCT (if iopt=2)
              pp(16)=0.d0 ! Tay1/brho,AQD/brho,Tay2/brho,ASEX/brho,Tay3/brho,AOCT/brho (=3)
              pp(17)=0.d0
              pp(18)=0.d0
              pp(19)=0.d0 ! axial rotation angle ("TILT" in MAD) [this feature untested]
              pp(20)=1. ! linear order
!              pp(20)=3.  ! 3rd order 
              pp(21)=1. ! slices
!
! jslice=present slice number; nslices=total number of slices,
! slfrac=fraction by which the element length is multiplied for slicing
! ihalf = 1 (first half) or 2 (second half) if a slice is split in half,
!         as when performing space-charge kicks
! ihalf = 0 if slices are not split in half.
! here is the logic to do multiple slices with a space-charge kick in
! the middle of each slice:
              nslices=bnseg
              slfrac=1.d0/nslices
              slfrac=slfrac*0.5d0  !cut each slice in half, do sc kick in the middle
! Track particles through the map:
              ntaysym=1 !ntaysym=1 for taylor, =2 for symplectic
              norder=1  !order of tracking (1=linear, ..., 5=5th order nonlinear)
!              norder=3  !3rd order !<<kilean>>
              ntrace=1  !number of times to apply the map
            endif

            t_transp = t_transp + elapsedtime_Timer(t11)

          do j = 1, bnseg
            call starttime_Timer(t2)

            if((Flagerr.eq.1).and.(Flagmap.eq.1)) then
              call geomerrL_BeamBunch(Bpts,Blnelem(i)) 
            end if
            
!-------------------------------------------------------------------
! use linear map or nonlinear Lorentz integrator to advance particles.
            if(bitype.ne.4) then
              if(bitype.ne.-14) then
              ! spatial drift.
              !linear map integrator
                if(Flagmap.eq.1) then
                  call map1_BeamBunch(Bpts,Blnelem(i),z,tau1,bitype)
                else
                  call map1_BeamBunch(Bpts,z,tau1)
                endif
              endif
            else
              if(dparam(4).gt.100 .and. dparam(4).le.300) then

              call Sector_Dipole(tau1,beta0,hd0,hd1,Bpts%Pts1,&
                                   Nplocal,qm0)
              z = z + tau1
              angz = angz + tau1*hd0
              !convert the coordinate (2), (4), (6) from the Transport back 
              !to the Impact coordinate, and also add gamma to (5) to go to beam frame.
              do ipt = 1, Nplocal
                ptarry(1) = Bpts%Pts1(1,ipt)
                ptarry(2) = Bpts%Pts1(2,ipt)
                ptarry(3) = Bpts%Pts1(3,ipt)
                ptarry(4) = Bpts%Pts1(4,ipt)
                ptarry(5) = Bpts%Pts1(5,ipt)
                ptarry(6) = Bpts%Pts1(6,ipt)
                Bpts%Pts1(6,ipt) = -(ptarry(6)+(Bpts%Pts1(7,ipt)-qm0)/qm0)*beta0*gambet
                gamn = gamma0 - Bpts%Pts1(6,ipt)
                Bpts%Pts1(5,ipt) = ptarry(5)*gamma0
                gambetz = sqrt((gamn**2-1)/(1+ptarry(2)**2+ptarry(4)**2))
                Bpts%Pts1(2,ipt) = ptarry(2)*gambetz
                Bpts%Pts1(4,ipt) = ptarry(4)*gambetz
              enddo
              else if(dparam(4).gt.300) then
                jslice = j
                ihalf=1
                call get_sbendmap(pp,xmh,h,jslice,nslices,slfrac,ihalf)
                if(ntaysym.eq.1)then
                  call brkts(h)
                  call eval(Bpts%Pts1,Nplocal,xmh,norder,ntrace)
                else
                  call canx(xmh,h,norder)
                  call rearr
!                    call evalsr(Bpts%Pts1,Nplocal,xmh,norder,ntrace)
                    tmpPts(1:6,1:Nplocal)=Bpts%Pts1(1:6,1:Nplocal)
                    call evalsr(tmpPts,Nplocal,xmh,norder,ntrace)
                    Bpts%Pts1(1:6,1:Nplocal)=tmpPts(1:6,1:Nplocal)
                endif
                z = z + tau1
              endif
            endif

            t_enlarge = t_enlarge + elapsedtime_Timer(t2)

            call starttime_Timer(t3)
            
!-------------------------------------------------------------------
! escape the space charge calculation for 0 current case
            if(flagExternalPipe) then
              call getPipeInfo(modulo(z,circumference),piperad,piperad2)
              pipeID = 2
            endif
            if(BcurrImp.lt.1.0e-30)  then !no space-charge
            !<<<<<<<<<<<<<< check particle loss (Kilean) <<<<<<<<<<<<<<<
              call lostcount_BeamBunch(Bpts,Nplocal,Np,&
                                       pipeID,piperad,piperad2,&
                                       lost_pdata,z,nlost)
            else if(Flagbc.eq.7 .or. Flagbc.eq.8 .or. Flagbc.eq.9) then
              call lostcount_BeamBunch(Bpts,Nplocal,Np,&
                                       pipeID,piperad,piperad2,&
                                       lost_pdata,z,nlost)
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              goto 200
            else if(flagcoll.eq.1) then !calculate space charge forces
            !start sc calculation-------------------
              if(bitype.ne.4 .or. (bitype.eq.4 .and. dparam(4).gt.300) ) then
                call conv0th_BeamBunch(Bpts,tau2,Nplocal,Np,ptrange,&
                                   Flagbc,Perdlen,piperad,piperad2)
                call chgupdate_BeamBunch(Bpts,nchrg,nptlist0,qmcclist0)
!comment out for test purpose
                if((bitype.eq.4 .or. (bitypeold.eq.4 .and. flagcsrTr.eq.1))&
                    .and. dparam(4).gt.400) then
                  flagcsr = 1
                else
                  flagcsr = 0
                endif
              else
                if(dparam(4).gt.200) then
                  flagcsr = 1
                else
                  flagcsr = 0
                endif
                ptrange(1) = 1.0e10
                ptrange(2) = -1.0e10
                ptrange(3) = 1.0e10
                ptrange(4) = -1.0e10
                ptrange(5) = 1.0e10
                ptrange(6) = -1.0e10
                do ipt = 1, Nplocal
                  if(ptrange(1).ge.Bpts%Pts1(1,ipt)) then
                    ptrange(1) = Bpts%Pts1(1,ipt)
                  endif
                  if(ptrange(2).le.Bpts%Pts1(1,ipt)) then
                    ptrange(2) = Bpts%Pts1(1,ipt)
                  endif
                  if(ptrange(3).ge.Bpts%Pts1(3,ipt)) then
                    ptrange(3) = Bpts%Pts1(3,ipt)
                  endif
                  if(ptrange(4).le.Bpts%Pts1(3,ipt)) then
                    ptrange(4) = Bpts%Pts1(3,ipt)
                  endif
                  if(ptrange(5).ge.Bpts%Pts1(5,ipt)) then
                    ptrange(5) = Bpts%Pts1(5,ipt)
                  endif
                  if(ptrange(6).le.Bpts%Pts1(5,ipt)) then
                    ptrange(6) = Bpts%Pts1(5,ipt)
                  endif
                enddo
              endif
              !fix the global range for sub-cycle of space charge potential.
              if(Flagsubstep.eq.1) then
                if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                  ptrange(1) = 0.0
                  ptrange(2) = piperad
                  ptrange(3) = 0.0
                  ptrange(4) = 4*asin(1.0)
                else
                  ptrange(1) = -piperad
                  ptrange(2) = piperad
                  ptrange(3) = -piperad2
                  ptrange(4) = piperad2
                endif
                ptrange(5) = -Perdlen/2
                ptrange(6) = Perdlen/2
              endif

            call random_number(zrandi)
            zadjust = zadjmax*zrandi
            ptrange(5) = ptrange(5)-(ptrange(6)-ptrange(5))*zadjust
            ptrange(6) = ptrange(6)+(ptrange(6)-ptrange(5))*zadjust
            ! get new boundary from the range of beam particles.
            if(Flagbc.eq.4 .or. Flagbc.eq.6) then
            !if(Flagbc.eq.4) then
            else
              if(flagdecomp.eq.1) then
                call update_CompDom(Ageom,ptrange,grid2d,Flagbc)
              else
                call updateglb_CompDom(Ageom,ptrange,grid2d,Flagbc)
              endif
            endif
            call getlcmnum_CompDom(Ageom,lcgrid)
            Nxlocal = lcgrid(1) 
            if(npy.gt.1) then
              Nylocal = lcgrid(2) + 2
            else
              Nylocal = lcgrid(2) 
            endif
            if(npx.gt.1) then
              Nzlocal = lcgrid(3) + 2
            else
              Nzlocal = lcgrid(3)
            endif
            call getlcrange_CompDom(Ageom,lcrange)
            call getrange_CompDom(Ageom,range)

            if((totnp.gt.1) .and. (flagdecomp.eq.1)) then
              ! pass particles to local space domain on new processor.
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                      Nplcmax,lcrange)
              else
                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                Nplcmax,lcrange)
              endif
            endif
            ! assign new 'Nplocal' local particles on each processor.
            call setnpt_BeamBunch(Bpts,Nplocal)

            if((mod(j-1,nsubstep).eq.0).or.(Flagsubstep.ne.1)) then
            !start sc for substep

            if(flagdecomp.eq.1) then
              call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,&
                                npy)
              deallocate(chgdens)
              allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
              ! deposit particles onto grid to obtain charge density.
              call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,Ageom,&
                                  grid2d,chgdens,Flagbc,Perdlen)
            else
            endif

!-------------------------------------------------------------------
! start load balance. (at the location of new space charge calculation)
!                if((mod(ibal,nbal).eq.0).and.(totnp.gt.1).and.(bitype.ne.4)&
!                    .and.(bitypeold.ne.4)) then

            if((mod(ibal,nbal).eq.0).and.(totnp.gt.1).and.&
               (flagdecomp.eq.1)) then
              call MPI_BARRIER(comm2d,ierr)
              call getlctabnm_CompDom(Ageom,temptab)
              lctabnmx(0:npx-1) = temptab(1,0:npx-1,0)
              lctabnmy(0:npy-1) = temptab(2,0,0:npy-1)
              call getmsize_CompDom(Ageom,msize)
              hy = msize(2)
              hz = msize(3) 
              call getrange_CompDom(Ageom,range)
              ymin = range(3)
              zmin = range(5)
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                call balance_CompDom(chgdens,lctabnmx,&
                lctabrgx,npx,npy,commrow,commcol, &
                lcgrid(1),lcgrid(2),lcgrid(3),Ny,Nz,hz,zmin)
                call setlctab_CompDom(Ageom,lctabnmx,lctabrgx,&
                                     npx,npy,myidx,myidy)
              else
                call balance_CompDom(chgdens,lctabnmx,&
                lctabnmy,lctabrgx,lctabrgy,npx,npy,commrow,commcol, &
                lcgrid(1),lcgrid(2),lcgrid(3),Ny,Nz,hy,hz,ymin,zmin)
                call setlctab_CompDom(Ageom,lctabnmx,lctabnmy,lctabrgx,&
                                     lctabrgy,npx,npy,myidx,myidy)
              endif
              call getlcmnum_CompDom(Ageom,lcgrid)
              Nxlocal = lcgrid(1) 
              if(npy.gt.1) then
                Nylocal = lcgrid(2) + 2
              else
                Nylocal = lcgrid(2) 
              endif
              if(npx.gt.1) then
                Nzlocal = lcgrid(3) + 2
              else
                Nzlocal = lcgrid(3) 
              endif
              call getlcrange_CompDom(Ageom,lcrange)
              deallocate(chgdens)
              allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
              call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,npy)

              ! pass particles to local space domain on new processor.
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                      Nplcmax,lcrange)
              else
                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                 Nplcmax,lcrange)
              endif
              ! assign new 'Nplocal' local particles on each processor.
              call setnpt_BeamBunch(Bpts,Nplocal)

              ! deposit particles onto grid to obtain charge density.
              call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,&
                                    Ageom,grid2d,chgdens,Flagbc,Perdlen)
            endif
            ibal = ibal + 1
!end load balance.
!-------------------------------------------------------------------
            if(npx.gt.1) then
              nzlcr = Nzlocal-2
              kadd = 1
            else
              nzlcr = Nzlocal
              kadd = 0
            endif
            if(npy.gt.1) then
              nylcr = Nylocal-2
              jadd = 1
            else
              nylcr = Nylocal
              jadd = 0
            endif
           
!-------------------------------------------------------------------------
! solve 3D Poisson's equation
            if(flagsc.eq.1) then

            if(Flagbc.eq.1) then
              ! solve Poisson's equation using 3D isolated boundary condition.
!comment out just for test purpose 8/15/06
               call update3O_FieldQuant(Potential,chgdens,Ageom,&
               grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.2) then
              ! solve Poisson's equation using 2D isolated 1D periodic 
              ! boundary condition.
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              call update2O1P_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.3) then
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update3_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr,&
              besscoef,bessnorm,gml,modth,nmod)
            else if(Flagbc.eq.4) then
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update4_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr,Perdlen)
            else if(Flagbc.eq.5) then
              ! solve Poisson's equation using 2D rectangular pipe, 1D open
              ! boundary condition
              if(myidy.eq.(npy-1)) nylcr = nylcr-1
              call update5_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.6) then
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update6_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.8) then
              ! solve 2D Poisson's equation using 2D isolated boundary condition.
               call update7_FieldQuant(Potential,chgdens,Ageom,&
               grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else
              print*,"no such boundary condition type!!!"
              stop
            endif

            else
              Potential%FieldQ = 0.0
            endif
            !for sc calculation
            !---------------------------------
             if(flagdecomp.ne.1) then
               print*,"under development...."
               stop
             else
             endif


            !end sc for substep
            endif
!----------------------------

              !includes wakefield
              if(flagwake.eq.1) then
                !hzwake = (range(6)-range(5))*1.0000001/(Nz-1) !avoid over index
                !no need to avoid over index since range already has 2 extra grids.
                hzwake = (range(6)-range(5))/(Nz-1) 
                xwakelc = 0.0
                ywakelc = 0.0
                denszlc = 0.0
                xwakez = 0.0
                ywakez = 0.0
                densz = 0.0
                !linear interpolation
                do ipt = 1, Nplocal
                  iizz = (Bpts%Pts1(5,ipt)-range(5))/hzwake + 1
                  iizz1 = iizz + 1
                  ab = ((range(5)-Bpts%Pts1(5,ipt))+iizz*hzwake)/hzwake
                  xwakelc(iizz) = xwakelc(iizz) + Bpts%Pts1(1,ipt)*ab
                  xwakelc(iizz1) = xwakelc(iizz1) + Bpts%Pts1(1,ipt)*(1.0-ab)
                  ywakelc(iizz) = ywakelc(iizz) + Bpts%Pts1(3,ipt)*ab
                  ywakelc(iizz1) = ywakelc(iizz1) + Bpts%Pts1(3,ipt)*(1.0-ab)
                  denszlc(iizz) = denszlc(iizz) + ab
                  denszlc(iizz1) = denszlc(iizz1) + 1.0 -ab
                enddo
                call MPI_ALLREDUCE(denszlc,densz,Nz,MPI_DOUBLE_PRECISION,&
                                   MPI_SUM,comm2d,ierr)
                call MPI_ALLREDUCE(xwakelc,xwakez,Nz,MPI_DOUBLE_PRECISION,&
                                   MPI_SUM,comm2d,ierr)
                call MPI_ALLREDUCE(ywakelc,ywakez,Nz,MPI_DOUBLE_PRECISION,&
                                   MPI_SUM,comm2d,ierr)

                !hzwake divided by gamma is due to fact that the particle z coordinate is in the beam
                !frame instead of the lab frame (during the T -> Z conversion).
                hzwake = hzwake/(-Bpts%refptcl(6))
                !get the line charge density along z
                do kz = 1, Nz
                  recvdensz(kz,1) = densz(kz)*Bpts%Current/Scfreq/Np/(hzwake)*&
                                    Bpts%Charge/abs(Bpts%Charge)
                enddo

                !due to the edge deposition
                do kz = 1, Nz
                  if(densz(kz).ne.0.0) then
                    xwakez(kz) = xwakez(kz)/densz(kz)
                    ywakez(kz) = ywakez(kz)/densz(kz)
                  else
                    xwakez(kz) = 0.0
                    ywakez(kz) = 0.0
                  endif
                  recvdensz(kz,2) = ywakez(kz)
                enddo

                !calculate the longitudinal and transverse wakefield from
                !the line charge density and analytical wake function
                call wakefield_FieldQuant(Nz,xwakez,ywakez,recvdensz,exwake,eywake,ezwake,&
                     hzwake,aawk,ggwk,lengwk,flagbtw)
                exwake = scwk*exwake
                eywake = scwk*eywake
                ezwake = scwk*ezwake

              endif

              flagcsr = 0
!---------------------------------------------------
!---------------------------------------------------
              if(flagcsr.eq.1) then
                hzwake = (range(6)-range(5))/(Nz-1)
                denszlc = 0.0
                densz = 0.0
                !linear interpolation
                do ipt = 1, Nplocal
                  iizz = (Bpts%Pts1(5,ipt)-range(5))/hzwake + 1
                  iizz1 = iizz + 1
                  ab = ((range(5)-Bpts%Pts1(5,ipt))+iizz*hzwake)/hzwake
                  denszlc(iizz) = denszlc(iizz) + ab
                  denszlc(iizz1) = denszlc(iizz1) + 1.0 -ab
                enddo
                call MPI_ALLREDUCE(denszlc,densz,Nz,MPI_DOUBLE_PRECISION,&
                                   MPI_SUM,comm2d,ierr)
                !hzwake divided by gamma is due to fact that the particle z coordinate is in the beam
                !frame instead of the lab frame (during the T -> Z conversion).
                hzwake = hzwake/(-Bpts%refptcl(6))
                !get the line charge density along z
                do kz = 1, Nz
                  densz(kz) = densz(kz)*Bpts%Current/Scfreq/Np/(hzwake)*&
                                    Bpts%Charge/abs(Bpts%Charge)
                enddo
 
                gam = -Bpts%refptcl(6)

! The following smooth the density distribution using different type
! of filters.
                    ldsg = 0
                    nlsg = 64
                    nrsg = 64
                    npsg = 129
                    msg = 1
                    csg = 0.0

                    ans=densz
                    call filterIP_FieldQuant(Nz,ans,densz,denszp,hzwake)

                    ans=denszp

                    call filterIP_FieldQuant(Nz,ans,denszp,denszpp,hzwake)

                    if(bitype.eq.4) then
                      zwkmin = range(5)/(-Bpts%refptcl(6)) + (z-zbleng)
                      bendlen = blength !inside the bend
                    else
                      zwkmin = range(5)/(-Bpts%refptcl(6)) + (z-zbleng+blengthold)
                      bendlen = blengthold !out of the bend
                    endif

                    !This includes both transient and steady state csr wake
                    ezwake = 0.0
!new csrwake including transition effects
                    call csrwakeTr3_FieldQuant(Nz,r0,zwkmin,hzwake,&
                                   bendlen,densz,denszp,denszpp,ezwake)
                    if(myid.eq.0) then
                      do iz = 1, Nz
                        write(1,777)iz*1.0,hzwake,densz(iz),denszp(iz),&
                                    denszpp(iz),ezwake(iz)
                      enddo
                      call flush(1)
                    endif
777                 format(6(1x,e16.7))

! steady state csr wake
                ssll = r0*(z/r0)**3/24 !???

                exwake = 0.0
                eywake = 0.0
              endif

            else
            endif
200         continue
! end space charge field calcualtion for curr > 0 or flagcol = 1
!-------------------------------------------------------------------
            t_shrink = t_shrink + elapsedtime_Timer(t3)
            call starttime_Timer(t4)

!-------------------------------------------------------------------
! use linear map or nonlinear Lorentz integrator to advance particles.
!------------------------------------------------------
            !use transfer map
            if(Flagmap.eq.1) then
              if((flagwake.eq.1) .or. (flagcsr.eq.1)) then
                if(flagdecomp.eq.1) then
                  !Potential%FieldQ = 0.0d0 !test wakefield
                  call kick1wake_BeamBunch(Bpts,tau2,Nxlocal,Nylocal,Nzlocal,&
                   Potential%FieldQ,Ageom,grid2d,Flagbc,Perdlen,exwake,eywake,&
                   ezwake,Nz,npx,npy,flagcoll)
                else
                endif
              else
                if(Flagbc.eq.7) then
                  call map2_BeamBunch(Bpts,z,tau2,Nplocal,Nx,Ny,xrad,yrad,&
                   Flagbc,flagcoll)
                elseif(flagdecomp.eq.1) then
                  call map2_BeamBunch(Bpts,tau2,Nxlocal,Nylocal,Nzlocal,&
                   Potential%FieldQ,Ageom,grid2d,Flagbc,Perdlen,flagcoll)
                else
                endif
              endif
              if(bitype.ne.4) then
                if(bitype.ne.-14) then
                call map1_BeamBunch(Bpts,Blnelem(i),z,tau1,bitype)
                endif
              else
                call starttime_Timer(t6)
                if(dparam(4).le.300) then !TRANSPORT
                if(BcurrImp.gt.0.0) then 
                  !this is due to the fact that after space-charge, particle coordinates
                  !go to Impact unit
                  do ipt = 1, Nplocal
                    rcpgammai = 1.0/(-Bpts%Pts1(6,ipt)+gamma0)
                    betai = sqrt(1.0-rcpgammai*rcpgammai*(1+Bpts%Pts1(2,ipt)**2+ &
                                 Bpts%Pts1(4,ipt)**2) )
                    Bpts%Pts1(1,ipt) = Bpts%Pts1(1,ipt)*Scxl
                    Bpts%Pts1(3,ipt) = Bpts%Pts1(3,ipt)*Scxl
                    Bpts%Pts1(5,ipt) = -gamma0*betai*Bpts%Pts1(5,ipt)*Scxl
                  enddo
                endif
                do ipt = 1, Nplocal
                  gamn = gamma0 - Bpts%Pts1(6,ipt)
                  gambetz = sqrt(gamn**2-1.0-Bpts%Pts1(2,ipt)**2-&
                               Bpts%Pts1(4,ipt)**2)
                  ptarry(2) = Bpts%Pts1(2,ipt)/gambetz
                  ptarry(4) = Bpts%Pts1(4,ipt)/gambetz
                  ptarry(5) = Bpts%Pts1(5,ipt)/gamma0
                  ptarry(6) = -Bpts%Pts1(6,ipt)/beta0/gambet - &
                            (Bpts%Pts1(7,ipt)-qm0)/qm0
                  Bpts%Pts1(2,ipt) = ptarry(2)
                  Bpts%Pts1(4,ipt) = ptarry(4)
                  Bpts%Pts1(5,ipt) = ptarry(5)
                  Bpts%Pts1(6,ipt) = ptarry(6)
                enddo
                call Sector_Dipole(tau1,beta0,hd0,hd1,Bpts%Pts1,&
                                   Nplocal,qm0)
                z = z + tau1
                angz = angz + tau1*hd0
                else !MaryLie
                  jslice = j
                  ihalf=2
                  call get_sbendmap(pp,xmh,h,jslice,nslices,slfrac,ihalf)
                  if(ntaysym.eq.1)then
                    call brkts(h)
                    call eval(Bpts%Pts1,Nplocal,xmh,norder,ntrace)
                  else
                    call canx(xmh,h,norder)
                    call rearr
!                    call evalsr(Bpts%Pts1,Nplocal,xmh,norder,ntrace)
                    tmpPts(1:6,1:Nplocal)=Bpts%Pts1(1:6,1:Nplocal)
                    call evalsr(tmpPts,Nplocal,xmh,norder,ntrace)
                    Bpts%Pts1(1:6,1:Nplocal)=tmpPts(1:6,1:Nplocal)
                  endif
                  z = z + tau1
                endif
                t_guardexch = t_guardexch + elapsedtime_Timer(t6)
              endif
            else
!------------------------------------------------------
            !use Lorentz integrator
              if((flagwake.eq.1) .or. (flagcsr.eq.1)) then
                if(flagdecomp.eq.1) then
                  call kick2wake_BeamBunch(Bpts,Blnelem(i),z,tau2,Nxlocal,Nylocal,&
                           Nzlocal,Potential%FieldQ,Ageom,grid2d,Flagbc,Flagerr,&
                           exwake,eywake,ezwake,Nz,npx,npy,flagcoll)
                else
                  !call kick2wakeglb_BeamBunch(Bpts,Blnelem(i),z,tau2,Nx,Ny,&
                  !         Nz,glpot,Ageom,grid2d,Flagbc,Flagerr,&
                  !         exwake,eywake,ezwake,Nz,npx,npy,flagcoll)
                endif
              else
                if(flagdecomp.eq.1) then
                  call map2_BeamBunch(Bpts,Blnelem(i),z,tau2,Nxlocal,Nylocal,&
                           Nzlocal,Potential%FieldQ,Ageom,grid2d,Flagbc,Flagerr,flagcoll)
                else
                  !call kick2glb_BeamBunch(Bpts,Blnelem(i),z,tau2,Nx,Ny,&
                  !         Nz,glpot,Ageom,grid2d,Flagbc,Flagerr,flagcoll)
                endif
              endif
              if(bitype.ne.4) then
                call map1_BeamBunch(Bpts,z,tau1)
              else
                if(dparam(4).le.300) then !Transport

                if(BcurrImp.gt.0.0) then
                  !this is due to the fact that after space-charge, particle coordinates
                  !go to Impact unit
                  do ipt = 1, Nplocal
                    rcpgammai = 1.0/(-Bpts%Pts1(6,ipt)+gamma0)
                    betai = sqrt(1.0-rcpgammai*rcpgammai*(1+Bpts%Pts1(2,ipt)**2+ &
                                 Bpts%Pts1(4,ipt)**2) )
                    Bpts%Pts1(1,ipt) = Bpts%Pts1(1,ipt)*Scxl
                    Bpts%Pts1(3,ipt) = Bpts%Pts1(3,ipt)*Scxl
                    Bpts%Pts1(5,ipt) = -gamma0*betai*Bpts%Pts1(5,ipt)*Scxl
                  enddo
                endif
                do ipt = 1, Nplocal
                  gamn = gamma0 - Bpts%Pts1(6,ipt)
                  gambetz = sqrt(gamn**2-1.0-Bpts%Pts1(2,ipt)**2-&
                               Bpts%Pts1(4,ipt)**2)
                  ptarry(2) = Bpts%Pts1(2,ipt)/gambetz
                  ptarry(4) = Bpts%Pts1(4,ipt)/gambetz
                  ptarry(5) = Bpts%Pts1(5,ipt)/gamma0
                  ptarry(6) = -Bpts%Pts1(6,ipt)/beta0/gambet - &
                            (Bpts%Pts1(7,ipt)-qm0)/qm0
                  Bpts%Pts1(2,ipt) = ptarry(2)
                  Bpts%Pts1(4,ipt) = ptarry(4)
                  Bpts%Pts1(5,ipt) = ptarry(5)
                  Bpts%Pts1(6,ipt) = ptarry(6)
                enddo
                call Sector_Dipole(tau1,beta0,hd0,hd1,Bpts%Pts1,&
                                   Nplocal,qm0)
                z = z + tau1
                angz = angz + tau1*hd0

                else !MaryLie

                  jslice = j
                  ihalf=2
                  call get_sbendmap(pp,xmh,h,jslice,nslices,slfrac,ihalf)
                  if(ntaysym.eq.1)then
                    call brkts(h)
                    call eval(Bpts%Pts1,Nplocal,xmh,norder,ntrace)
                  else
                    call canx(xmh,h,norder)
                    call rearr
!                    call evalsr(Bpts%Pts1,Nplocal,xmh,norder,ntrace)
                    tmpPts(1:6,1:Nplocal)=Bpts%Pts1(1:6,1:Nplocal)
                    call evalsr(tmpPts,Nplocal,xmh,norder,ntrace)
                    Bpts%Pts1(1:6,1:Nplocal)=tmpPts(1:6,1:Nplocal)
                  endif
                  z = z + tau1
                endif
              endif
            endif

            if((Flagerr.eq.1).and.(Flagmap.eq.1)) then
              call geomerrT_BeamBunch(Bpts,Blnelem(i)) 
            end if

            t_guardsum = t_guardsum + elapsedtime_Timer(t4)
            call starttime_Timer(t5)

            nstep = nstep + 1

        !Output after nstep
          !Special diagnostic for symplectic SC solver added:
               if(Flagbc.eq.7) then
                  call map2_BeamBunch(Bpts,z,tau2,Nplocal,Nx,Ny,xrad,yrad,&
                   1,flagcoll)
               endif
          !Standard diagnostics:
          if(Flagdiag.eq.1) then
          ! <<<<<<<<<<<<<<<<<<<< Kilean <<<<<<<<<<<<<<<<<<<<<<<
          !  call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
          !else
          !  call diagnostic2_Output(Bpts,z,nchrg,nptlist0)
            call diagnostic1_Output(z,Bpts,nchrg,[Bpts%Npt])
          else
            call diagnostic2_Output(Bpts,z,nchrg,[Bpts%Npt])
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          endif

!--------------------------------------------------------------------------
          !end of nstep in an element
            t_boundint = t_boundint + elapsedtime_Timer(t5)
          end do
!          !shift back the energy centroid
!          if(bitype.eq.4 .and. dparam(4).gt.300) then
!            Bpts%refptcl(6) = Bpts%refptcl(6) - pz0
!            do ipt = 1, Nplocal
!              Bpts%Pts1(6,ipt) = Bpts%Pts1(6,ipt) + pz0
!            enddo
!          endif
            if(bitype.eq.4 .and. dparam(4).le.300) then
              call Bpol_Dipole(hd0,hB,tanphiB,tanphiBb,hd1,&
                                 psi2,Bpts%Pts1,angB,Nplocal)
              !convert back to the Impact coordinate
              do ipt = 1, Nplocal
                ptarry(1) = Bpts%Pts1(1,ipt)
                ptarry(2) = Bpts%Pts1(2,ipt)
                ptarry(3) = Bpts%Pts1(3,ipt)
                ptarry(4) = Bpts%Pts1(4,ipt)
                ptarry(5) = Bpts%Pts1(5,ipt)
                ptarry(6) = Bpts%Pts1(6,ipt)
                Bpts%Pts1(6,ipt) = -(ptarry(6)+(Bpts%Pts1(7,ipt)-qm0)/qm0)*beta0*gambet
                gamn = gamma0 - Bpts%Pts1(6,ipt)
                Bpts%Pts1(1,ipt) = ptarry(1)/Scxl
                Bpts%Pts1(3,ipt) = ptarry(3)/Scxl
                Bpts%Pts1(5,ipt) = -ptarry(5)/(Scxl*beta0)
                gambetz = sqrt((gamn**2-1)/(1+ptarry(2)**2+ptarry(4)**2))
                Bpts%Pts1(2,ipt) = ptarry(2)*gambetz
                Bpts%Pts1(4,ipt) = ptarry(4)*gambetz
              enddo
              if(Flagdiag.eq.1) then
              ! <<<<<<<<<<<<<<<<<<<< Kilean <<<<<<<<<<<<<<<<<<<<<<<
              !  call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
              !else
              !  call diagnostic2_Output(Bpts,z,nchrg,nptlist0)
                call diagnostic1_Output(z,Bpts,nchrg,[Bpts%Npt])
              else
                call diagnostic2_Output(Bpts,z,nchrg,[Bpts%Npt])
              !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              endif
            endif
!------------------------------------------------------------------------------
          else !//bend magnet using Time domain
            call getparam_BeamLineElem(Blnelem(i),dparam)
            if(dparam(4).gt.100.0) then !file id
            else

            !//go to T frame
            tg = Bpts%refptcl(5)
            call getparam_BeamLineElem(Blnelem(i),4,rfile)
            nfile = int(rfile + 0.1)
            call read3_Data(nfile) !input the geometry information of bend
            !Bheit = Fcoef(2)/Scxl !vertical bend height (wrong one)
            Bpts%refptcl(6) = -Fcoef(2)
            !//go to T frame
            call convZT_BeamBunch(Bpts)
            !//loop through bnseg steps
            gam = sqrt(1.0+Bpts%refptcl(6)**2)
            !//normalized t (omega t)
            tau2 = 2*pi*Scfreq*blength/(Clight*Bpts%refptcl(6)/gam)/bnseg
            tv = 0.0
            call diagnosticT_Output(tv,Bpts)
            Flagbctmp = 1
            zz = 0.0
            do j = 1, 2*bnseg
              call drifthalfT_BeamBunch(Bpts,tv,tau2)
              ptref = Bpts%refptcl
              !//go to the local "ptref" coordinates
              call rotto_BeamBunch(Bpts,ptref,ptrange)
              if(BcurrImp.lt.1.0e-10) goto 400
              call update_CompDom(Ageom,ptrange,grid2d,Flagbctmp)
              call getlcmnum_CompDom(Ageom,lcgrid)
              Nxlocal = lcgrid(1)
              if(npy.gt.1) then
                Nylocal = lcgrid(2) + 2
              else
                Nylocal = lcgrid(2)
              endif
              if(npx.gt.1) then
                Nzlocal = lcgrid(3) + 2
              else
                Nzlocal = lcgrid(3)
              endif
              call getlcrange_CompDom(Ageom,lcrange)
              call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                Nplcmax,lcrange)
              ! assign new 'Nplocal' local particles on each processor.
              call setnpt_BeamBunch(Bpts,Nplocal)
              call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,&
                                npy)
              deallocate(chgdens)
              allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
              ! deposit particles onto grid to obtain charge density.
              call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,Ageom,&
                                  grid2d,chgdens,Flagbc,Perdlen)
              ibal = ibal + 1

              if(npx.gt.1) then
                 nzlcr = Nzlocal-2
              else
                 nzlcr = Nzlocal
              endif
              if(npy.gt.1) then
                nylcr = Nylocal-2
              else
                nylcr = Nylocal
              endif
 
              call update3O_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
400           continue
              !//kick particles in the rotated local coordinates.
              call kickT_BeamBunch(Bpts,Blnelem(i),tv,tau2,Nxlocal,Nylocal,&
                      Nzlocal,Potential%FieldQ,Ageom,grid2d,Flagbc,Flagerr)
              call rotback_BeamBunch(Bpts,ptref)
              call drifthalfT_BeamBunch(Bpts,tv,tau2)
              tv = tv + tau2
              call diagnosticT_Output(tv,Bpts)
              itmp = 50 + j
555           format(9(1x,e14.6))
              vref = sqrt(Bpts%refptcl(2)**2+Bpts%refptcl(6)**2) &
                /sqrt(1.0+Bpts%refptcl(2)**2+Bpts%refptcl(6)**2)
              zz = zz + vref*tau2*Scxl
              if(zz.gt.blength) then
                exit
              endif
            enddo
            tg = tg+tv
            call convTZ_BeamBunch(Bpts,tg) !//go back to Z frame
            z = z + blength

            endif
          endif  !//end bend magnets
          bitypeold = bitype
          blengthold = blength
          zbleng = zbleng + blength
!output after each element
          call diagnostic1_Output(z,Bpts,nchrg,[Bpts%Npt])
        !-------------------------------------------------
        ! <<<<<<<<<<<<<<<<<<<<<<<<< Kilean <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        call write_lost_pData(lost_pdata,nlost,bitype,i)
        nlost = 0
        pipe_override = .false.
        pipeID = 1
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        enddo !end loop of N beam line elements

        if(myid.eq.0) print*,"iturn: ",iturn
        if(Flagdiag.eq.1) then
          ! <<<<<<<<<<<<<<<<<<<< Kilean <<<<<<<<<<<<<<<<<<<<<<<
          !  call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
          !else
          !  call diagnostic2_Output(Bpts,z,nchrg,nptlist0)
            call diagnostic1_Output(z,Bpts,nchrg,[Bpts%Npt])
          else
            call diagnostic2_Output(Bpts,z,nchrg,[Bpts%Npt])
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        endif
      enddo !end loop of Nturn 
!-------------------


! final output.
        call MPI_BARRIER(comm2d,ierr)
        t_integ = t_integ + elapsedtime_Timer(t0)
        call showtime_Timer()

        deallocate(lctabnmx,lctabnmy)
        deallocate(lctabrgx,lctabrgy)
        deallocate(temptab)
        deallocate(chgdens)
        deallocate(tmpptcs)
        if(Flagbc.eq.3) then
          deallocate(besscoef)
          deallocate(bessnorm)
          deallocate(gml)
          deallocate(pydisp)
          deallocate(modth)
        endif
        !deallocate(ztable)
        !deallocate(zdisp)
        deallocate(denszlc)
        deallocate(densz)
        deallocate(denszp)
        deallocate(denszpp)
        deallocate(csg)
        deallocate(ans)
        deallocate(xwakelc)
        deallocate(xwakez)
        deallocate(ywakelc)
        deallocate(ywakez)
        deallocate(recvdensz)
        deallocate(sendensz)
        deallocate(exwake)
        deallocate(eywake)
        deallocate(ezwake)
        deallocate(tmppot)
        deallocate(recvtmppot)
        deallocate(gltmppot)
        deallocate(glpot)
        deallocate(tmpPts)

        end subroutine run_AccSimulator

        subroutine restart_AccSimulator()
        implicit none
        include 'mpif.h'
        integer :: i,j,bnseg,bmpstp,bitype
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,npx,npy,&
                   totnp,ierr,nylcr,nzlcr
        integer :: nbal,ibal,nstep,ifile
        integer :: tmpfile,nfile
        integer, dimension(3) :: lcgrid
        integer, allocatable, dimension(:) :: lctabnmx,lctabnmy
        integer, allocatable, dimension(:,:,:) :: temptab
        double precision :: z0,z,tau1,tau2,blength,t0
        double precision, allocatable, dimension(:,:) :: lctabrgx, lctabrgy
        double precision, dimension(6) :: lcrange, range, ptrange
        double precision, dimension(3) :: msize
        double precision :: hy,hz,ymin,zmin,piperad,zedge
        double precision :: tmp1,tmp2,tmp3,tmp4,rfile
        double precision, allocatable, dimension(:,:,:) :: chgdens
        double precision, allocatable, dimension(:,:,:) :: besscoef
        double precision, allocatable, dimension(:,:) :: bessnorm,gml
        integer, allocatable, dimension(:) :: modth,pydisp
        integer :: nmod,k
        !double precision :: sumtest, sumtest2, sumtest3
        double precision :: piperad2
        integer :: flagcoll,flagtmp,flag1056
  
        flag1056 = 0
        flagcoll = 1
!-------------------------------------------------------------------
! prepare initial parameters, allocate temporary array.
        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid2d,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid2d,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)

        call starttime_Timer(t0)

        allocate(lctabnmx(0:npx-1))
        allocate(lctabnmy(0:npy-1))
        allocate(lctabrgx(2,0:npx-1))
        allocate(lctabrgy(2,0:npy-1))
        allocate(temptab(2,0:npx-1,0:npy-1))

        nbal = 5
        ibal = ibalend
        nstep = nstepend
        z = zend

        if(Flagdiag.eq.1) then
          ! <<<<<<<<<<<<<<<<<<<< Kilean <<<<<<<<<<<<<<<<<<<<<<<
          !  call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
          !else
          !  call diagnostic2_Output(Bpts,z,nchrg,nptlist0)
            call diagnostic1_Output(z,Bpts,nchrg,[Bpts%Npt])
          else
            call diagnostic2_Output(Bpts,z,nchrg,[Bpts%Npt])
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        endif

        allocate(chgdens(1,1,1))
!-------------------------------------------------------------------
! prepare for round pipe, open longitudinal
        if(Flagbc.eq.3) then
          allocate(pydisp(0:npy-1))
          call getlctabnm_CompDom(Ageom,temptab)
          pydisp(0) = 0
          do i = 1, npy-1
            pydisp(i) = pydisp(i-1) + temptab(2,0,i-1)
          enddo
          call getlcmnum_CompDom(Ageom,lcgrid)
          Nxlocal = lcgrid(1)
          if(npcol.gt.1) then
            Nylocal = lcgrid(2) + 2
          else
            Nylocal = lcgrid(2)
          endif
          if(nprow.gt.1) then
            Nzlocal = lcgrid(3) + 2
          else
            Nzlocal = lcgrid(3)
          endif

          if(nprow.gt.1) then
            nzlcr = Nzlocal-2
          else
            nzlcr = Nzlocal
          endif
          if(npcol.gt.1) then
            nylcr = Nylocal-2
          else
            nylcr = Nylocal
          endif
          if(myidy.eq.(npcol-1)) nylcr = nylcr - 1
          allocate(modth(nylcr))
          do i = 1, nylcr
            modth(i) = (pydisp(myidy)+i-2+1)/2
          enddo
          if(myidy.eq.0) then
            modth(1) = 0
            modth(2) = (ny-1)/2
            nmod = modth(nylcr) - modth(1) + 2
          else
          nmod = modth(nylcr) - modth(1) + 1
          endif
          allocate(besscoef(lcgrid(1),lcgrid(1),nmod))
          allocate(bessnorm(lcgrid(1),nmod))
          allocate(gml(lcgrid(1),nmod))

          call getmsize_CompDom(Ageom,msize)
          call Bessprep_Bessel(msize(1),lcgrid(1),nylcr,nmod,modth,&
                        besscoef,bessnorm,gml)
        endif

!-------------------------------------------------------------------
! start looping through 'Nblem' beam line elements.
          call getparam_BeamLineElem(Blnelem(iend),blength,bnseg,bmpstp,&
                                     bitype)
          call getradius_BeamLineElem(Blnelem(iend),piperad,piperad2)
          nfile = 0
          tau1 = 0.0
          if(bitype.ge.0) tau1 = 0.5*blength/bnseg
          tau2 = 2.0*tau1

!-------------------------------------------------------------------
! read in the on axis E field for rf cavities.
          if(bitype.gt.100) then
            call getparam_BeamLineElem(Blnelem(iend),5,rfile)
            nfile = int(rfile + 0.1)
            tmpfile = 0
            ifile = nfile
            if(ifile.ne.tmpfile)then
              !for linear map integrator
              if(Flagmap.eq.1) then
                !read in Ez, Ez', Ez'' on the axis
                call read1_Data(ifile)
              else
                !read in Er, Ez, H_theta on r-z grid 
                !call read2_Data(ifile)
                !read in Fourier coefficients of Ez on the axis
                call read3_Data(ifile)
              endif
              tmpfile=ifile
            endif
          endif

!-------------------------------------------------------------------
! loop through 'bnseg' numerical segments in each beam element
! using 2 step symplectic integeration (ie. leap frog).
          zedge = z
          call setparam_BeamLineElem(Blnelem(iend),1,zedge)
          do j = jend+1, bnseg
            if((Flagerr.eq.1).and.(Flagmap.eq.1)) then
              call geomerrL_BeamBunch(Bpts,Blnelem(iend)) 
            end if

!-------------------------------------------------------------------
! use linear map or nonlinear Lorentz integrator to advance particles.
            ! spatial drift.
            !linear map integrator
            if(Flagmap.eq.1) then
              call map1_BeamBunch(Bpts,Blnelem(i),z,tau1,bitype)
            else
              call map1_BeamBunch(Bpts,z,tau1)
            endif

!-------------------------------------------------------------------
! escape the space charge calculation for 0 current case
            if(BcurrImp.lt.1.0e-30) goto 200
! escape the PIC-based space charge calculation when the symplectic
! spectral space charge solver is used
            if(Flagbc.eq.7) goto 200

            call conv1st_BeamBunch(Bpts,tau2,Nplocal,Np,ptrange,&
                                   Flagbc,Perdlen,piperad,piperad2)

            !fix the global range for sub-cycle of space charge potential.
            if(Flagsubstep.eq.1) then
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                ptrange(1) = 0.0
                ptrange(2) = piperad
                ptrange(3) = 0.0
                ptrange(4) = 4*asin(1.0)
              else
                ptrange(1) = -piperad
                ptrange(2) = piperad
                ptrange(3) = -piperad2
                ptrange(4) = piperad2
              endif
              ptrange(5) = -Perdlen/2
              ptrange(6) = Perdlen/2
            endif

            ! get new boundary from the range of beam particles.
            if(Flagbc.eq.4) then
            else
              call update_CompDom(Ageom,ptrange,grid2d,Flagbc)
            endif
            call getlcmnum_CompDom(Ageom,lcgrid)
            Nxlocal = lcgrid(1) 
            if(npy.gt.1) then
              Nylocal = lcgrid(2) + 2
            else
              Nylocal = lcgrid(2) 
            endif
            if(npx.gt.1) then
              Nzlocal = lcgrid(3) + 2
            else
              Nzlocal = lcgrid(3)
            endif
            deallocate(chgdens)
            allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
            call getlcrange_CompDom(Ageom,lcrange)
            call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,&
                                npy)

            if(totnp.gt.1) then
              ! pass particles to local space domain on new processor.
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                      Nplcmax,lcrange)
              else
                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                Nplcmax,lcrange)
              endif
            endif
            ! assign new 'Nplocal' local particles on each processor.
            call setnpt_BeamBunch(Bpts,Nplocal)

            ! deposit particles onto grid to obtain charge density.
            call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,Ageom,&
                                  grid2d,chgdens,Flagbc,Perdlen)

!-------------------------------------------------------------------
! start load balance.
            if((mod(ibal,nbal).eq.0).and.(totnp.gt.1)) then
              call MPI_BARRIER(comm2d,ierr)

              call getlctabnm_CompDom(Ageom,temptab)
              lctabnmx(0:npx-1) = temptab(1,0:npx-1,0)
              lctabnmy(0:npy-1) = temptab(2,0,0:npy-1)
              call getmsize_CompDom(Ageom,msize)
              hy = msize(2)
              hz = msize(3) 
              call getrange_CompDom(Ageom,range)
              ymin = range(3)
              zmin = range(5)
              if(Flagbc.eq.3) then
                call balance_CompDom(chgdens,lctabnmx,&
                lctabrgx,npx,npy,commrow,commcol, &
                lcgrid(1),lcgrid(2),lcgrid(3),Ny,Nz,hz,zmin)
                call setlctab_CompDom(Ageom,lctabnmx,lctabrgx,&
                                     npx,npy,myidx,myidy)
              else
                call balance_CompDom(chgdens,lctabnmx,&
                lctabnmy,lctabrgx,lctabrgy,npx,npy,commrow,commcol, &
                lcgrid(1),lcgrid(2),lcgrid(3),Ny,Nz,hy,hz,ymin,zmin)
                call setlctab_CompDom(Ageom,lctabnmx,lctabnmy,lctabrgx,&
                                     lctabrgy,npx,npy,myidx,myidy)
              endif
              call getlcmnum_CompDom(Ageom,lcgrid)
              Nxlocal = lcgrid(1) 
              if(npy.gt.1) then
                Nylocal = lcgrid(2) + 2
              else
                Nylocal = lcgrid(2) 
              endif
              if(npx.gt.1) then
                Nzlocal = lcgrid(3) + 2
              else
                Nzlocal = lcgrid(3) 
              endif
              call getlcrange_CompDom(Ageom,lcrange)
              deallocate(chgdens)
              allocate(chgdens(Nxlocal,Nylocal,Nzlocal))
              call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,npy)

              ! pass particles to local space domain on new processor.
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                      Nplcmax,lcrange)
              else
                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                 Nplcmax,lcrange)
              endif
              ! assign new 'Nplocal' local particles on each processor.
              call setnpt_BeamBunch(Bpts,Nplocal)

              ! deposit particles onto grid to obtain charge density.
              call charge_BeamBunch(Bpts,Nplocal,Nxlocal,Nylocal,Nzlocal,&
                                    Ageom,grid2d,chgdens,Flagbc,Perdlen)
            endif
            ibal = ibal + 1
!end load balance.
!-------------------------------------------------------------------

            if(npx.gt.1) then
              nzlcr = Nzlocal-2
            else
              nzlcr = Nzlocal
            endif
            if(npy.gt.1) then
              nylcr = Nylocal-2
            else
              nylcr = Nylocal
            endif
           
!-------------------------------------------------------------------------
! solve 3D Poisson's equation
            if(Flagbc.eq.1) then
              ! solve Poisson's equation using 3D isolated boundary condition.
              call update3O_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.2) then
              ! solve Poisson's equation using 2D isolated 1D periodic 
              ! boundary condition.
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              call update2O1P_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.3) then
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update3_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr,&
              besscoef,bessnorm,gml,modth,nmod)
            else if(Flagbc.eq.4) then
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update4_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr,Perdlen)
            else if(Flagbc.eq.5) then
              ! solve Poisson's equation using 2D rectangular pipe, 1D open
              ! boundary condition
              if(myidy.eq.(npy-1)) nylcr = nylcr-1
              call update5_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.6) then
              if(myidx.eq.(npx-1)) nzlcr = nzlcr - 1
              if(myidy.eq.(npy-1)) nylcr = nylcr - 1
              call update6_FieldQuant(Potential,chgdens,Ageom,&
              grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else if(Flagbc.eq.8) then
              ! solve 2D Poisson's equation using 2D isolated boundary condition.
               call update7_FieldQuant(Potential,chgdens,Ageom,& 
               grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr)
            else
              print*,"no such boundary condition type!!!"
              stop
            endif

            call cvbkforth1st_BeamBunch(Bpts)
            if(totnp.gt.1) then
              if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
                call ptsmv5r_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                 Nplcmax,lcrange)
              else
                call ptsmv2_Ptclmger(Bpts%Pts1,Nplocal,grid2d,Pdim,&
                                 Nplcmax,lcrange)
              endif
            endif
            call setnpt_BeamBunch(Bpts,Nplocal)

200         continue
! end space charge field calcualtion.
!-------------------------------------------------------------------

!-------------------------------------------------------------------
! use linear map or nonlinear Lorentz integrator to advance particles.
            if(Flagmap.eq.1) then
            ! kick particles in velocity space.
              !call map2_BeamBunch(Bpts,tau2,Nxlocal,Nylocal,Nzlocal,&
              !     Potential%FieldQ,Ageom,grid2d,Flagbc,Perdlen,flagcoll)
              !  if(flagdecomp.eq.1) then
               if(Flagbc.eq.7) then
                  call map2_BeamBunch(Bpts,z,tau2,Nplocal,Nx,Ny,xrad,yrad,&
                   Flagbc,flagcoll)
               else
                  call map2_BeamBunch(Bpts,tau2,Nxlocal,Nylocal,Nzlocal,&
                   Potential%FieldQ,Ageom,grid2d,Flagbc,Perdlen,flagcoll)
               endif
              !  else
              !    call kick1glb_BeamBunch(Bpts,tau2,Nx,Ny,Nz,&
              !     glpot,Ageom,grid2d,Flagbc,Perdlen,flagcoll)
              !  endif
                call map1_BeamBunch(Bpts,Blnelem(i),z,tau1,bitype)

            else
              call map2_BeamBunch(Bpts,Blnelem(iend),z,tau2,Nxlocal,Nylocal,&
                   Nzlocal,Potential%FieldQ,Ageom,grid2d,Flagbc,Flagerr,&
                           flagcoll)
              !    call kick2glb_BeamBunch(Bpts,Blnelem(i),z,tau2,Nx,Ny,&
              !             Nz,glpot,Ageom,grid2d,Flagbc,Flagerr,flagcoll)
              call map1_BeamBunch(Bpts,z,tau1)
            endif

            if(Flagerr.eq.1) then
              call geomerrT_BeamBunch(Bpts,Blnelem(iend)) 
            end if

            if(Flagdiag.eq.1) then
          ! <<<<<<<<<<<<<<<<<<<< Kilean <<<<<<<<<<<<<<<<<<<<<<<
          !  call diagnostic1_Output(z,Bpts,nchrg,nptlist0)
          !else
          !  call diagnostic2_Output(Bpts,z,nchrg,nptlist0)
            call diagnostic1_Output(z,Bpts,nchrg,[Bpts%Npt])
          else
            call diagnostic2_Output(Bpts,z,nchrg,[Bpts%Npt])
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            endif

            nstep = nstep + 1

          end do

          zend = z

        call MPI_BARRIER(comm2d,ierr)

        deallocate(lctabnmx,lctabnmy)
        deallocate(lctabrgx,lctabrgy)
        deallocate(temptab)
        deallocate(chgdens)
        if(Flagbc.eq.3) then
          deallocate(besscoef)
          deallocate(bessnorm)
          deallocate(gml)
          deallocate(pydisp)
          deallocate(modth)
        endif

        end subroutine restart_AccSimulator

        subroutine destruct_AccSimulator(time)
        implicit none
        include 'mpif.h'
        double precision :: time
 
        call destruct_Data()
        call destruct_BeamBunch(Bpts)
        call destruct_FieldQuant(Potential)
        call destruct_CompDom(Ageom)

        deallocate(nptlist0)
        deallocate(currlist0)
        deallocate(qmcclist0)
 
        call end_Output(time)

        end subroutine destruct_AccSimulator

      end module AccSimulatorclass
