!----------------------------------------------------------------
! (c) Copyright, 2003 by the Regents of the University of California.
! Distributionclass: Initial distribution of charged beam bunch class in 
!                    Beam module of APPLICATION layer.
! Version: 2.0
! Author: Ji Qiang, Robert Ryne, LBNL, 7/25/03
! Description: This class defines initial distributions for the charged 
!              particle beam bunch information in the accelerator.
! Comments: we have added three attributes to each particle:
!           x,px,y,py,t,pt,charge/mass,charge weight,id
!----------------------------------------------------------------
      module Distributionclass
        use Pgrid2dclass
        use CompDomclass
        use BeamBunchclass      
        use Timerclass
        use NumConstclass
        use PhysConstclass
        !<<<<<<<<<<<<<<< begin distIOTA type(Kilean)  <<<<<<<<<<<<<<<<<<
        use Multipoleclass
        type, private :: distIOTA_class
          double precision :: t,c,beta,betap,emittance
          contains
            procedure :: init => distIOTA_init
            procedure :: getH => distIOTA_getH
            procedure :: genP_waterbag => distIOTA_genP_waterbag
            procedure :: genP_gaussian => distIOTA_genP_gaussian
            procedure :: secant_method => distIOTA_secant_method
        end type
        private :: distIOTA_init, distIOTA_getH, &
                   distIOTA_genP_waterbag, distIOTA_genP_gaussian, &
                   distIOTA_secant_method
        !>>>>>>>>>>>>>>>>>>>>> end distIOTA tpye >>>>>>>>>>>>>>>>>>>>>>>
      contains
        ! sample the particles with intial distribution.
        subroutine sample_Dist(this,distparam,nparam,flagdist,geom,grid,Flagbc,&
                               nchrg,nptlist,qmcclist,currlist)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nparam,Flagbc,nchrg
        double precision, dimension(nparam) :: distparam
        double precision, dimension(nchrg) :: qmcclist,currlist
        integer*8, dimension(nchrg) :: nptlist
        type (BeamBunch), intent(inout) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: flagdist
        integer :: myid, myidx, myidy,seedsize,i,isize
        !integer seedarray(1)
        integer, allocatable, dimension(:) :: seedarray
        real rancheck
        integer :: nslice

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(100001+myid)*(myid+7)
!        call random_seed(put=seedarray(1:1))
!        write(6,*)'seedarray=',seedarray

        call random_seed(SIZE=seedsize)
        allocate(seedarray(seedsize))
        do i = 1, seedsize
!          seedarray(i) = (1000+5*myid)*(myid+7)+i-1 !//original one.
          !//2nd group
!          seedarray(i) = (2000+5*myid)*(myid+7)+i-1
          !//3rd group
!used in previous simulator before 08/15/07
!          seedarray(i) = (3000+5*myid)*(myid+7)+i-1
!          seedarray(i) = (3000+5*myid)*(myid+2000000)+i-1
!          seedarray(i) = (300000+5*myid)*(myid+2000000)+i-1
          seedarray(i) = (3000+5*myid)*(myid+20000)+i-1
        enddo
        call random_seed(PUT=seedarray)
        call random_number(rancheck)
        !the following is added new for 2nd random group ....
        do i = 1, 3000
          call random_number(rancheck)
        enddo
!        write(6,*)'myid,rancheck=',seedarray,myid,rancheck

        nslice = distparam(16)+0.1

        if(flagdist.eq.1) then
          call Uniform_Dist(this,nparam,distparam,geom,grid)
        else if(flagdist.eq.2) then
!          call Gauss1_Dist(this,nparam,distparam,geom,grid)
          call Gauss3_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.3) then
          call Waterbag_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.4) then
          call Semigauss_Dist(this,nparam,distparam,geom,grid)
        else if(flagdist.eq.5) then
          call KV3d_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.6) then
          call regen_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.7) then
          call GaussGamma_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.8) then
          call WaterGamma_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.9) then
          call KVGamma_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.10) then
          call regen2_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.11) then
          call Gauss4new_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.12) then
          call Gauss7_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.14) then
          call Regen7_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.15) then
          call GaussDouble_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.16) then
          call WaterbagMC_Dist(this,nparam,distparam,grid,0,nchrg,&
                               nptlist,qmcclist,currlist)
        else if(flagdist.eq.17) then
          call GaussMC_Dist(this,nparam,distparam,grid,0,nchrg,&
                               nptlist,qmcclist,currlist)
        else if(flagdist.eq.22) then
          !here the repopulation parameter can be controlled
          !by external parameters from distparam.
          call readElegant_Dist(this,nparam,distparam,geom,grid,Flagbc,nslice)
          !call readElegantDB_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.23) then
          call read_Dist(this)
        else if(flagdist.eq.-23) then
          call read_Dist_binary(this,nparam,distparam)
        else if(flagdist.eq.24) then
          call readElegant2_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.25) then
          call readElegantRot_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.26) then
          call readElegantNcor_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.27) then
          !call readElegantlaser_Dist(this,nparam,distparam,geom,grid,Flagbc)
          call readElegantlaser2_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.28) then
          call readElegantNcorMod_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.29) then
          call readElegantNcorMod3_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.30) then
          call readElegantNcorMod4_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.33) then
          call Gauss3ldrd_Dist(this,nparam,distparam,grid,0)
          !for test the Franklin version purpose 6/25/08
          !call Gauss3ldrd4_Dist(this,nparam,distparam,grid,0)
        else if(flagdist.eq.34) then
          call readElegantCor_Dist(this,nparam,distparam,geom,grid,Flagbc,nslice)
        else if(flagdist.eq.35) then
          call readimpt_Dist(this,nparam,distparam,geom,grid,Flagbc)
        else if(flagdist.eq.38) then
          call Gauss3dSoblcls_Dist(this,nparam,distparam,grid)
        else if(flagdist.eq.39) then
          call readMLI_Dist(this,nparam,distparam,grid)
        !<<<<<<<<<<<<<<<<<< Kilean <<<<<<<<<<<<<<<<<<<<<<<
        else if(flagdist.eq.81) then
          call distIOTA_waterbag(this,nparam,distparam)
        else if(flagdist.eq.82) then
          call distIOTA_gaussian(this,nparam,distparam)
        !>>>>>>>>>>>>>>>>>> Kilean >>>>>>>>>>>>>>>>>>>>>>>
        else
          print*,"Initial distribution not available!!"
          stop
        endif

        deallocate(seedarray)

        end subroutine sample_Dist
       
        subroutine Gauss1_Dist(this,nparam,distparam,geom,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(6,2) :: a
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        double precision :: r1, r2
        double precision,allocatable,dimension(:,:) :: ptstmp
        integer :: totnp,npy,npx,myid,myidy,myidx,comm2d, &
                   commcol,commrow,ierr
        integer :: avgpts,numpts0,i,ii,i0,j,jj,intvsamp,pid
        double precision :: t0
        double precision, allocatable, dimension(:,:) :: x1,x2,x3,x4,x5,x6
        integer :: npttmp,intv

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        avgpts = this%Npt/(npx*npy)

        call getlcrange_CompDom(geom,lcrange)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        numpts0 = 0

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
! The performance inside this loop might be improved due to
! a lot of subroutine call in this loop.    
        npttmp = this%Npt/intv
        allocate(x1(2,npttmp))
        allocate(x2(2,npttmp))
        allocate(x3(2,npttmp))
        allocate(x4(2,npttmp))
        allocate(x5(2,npttmp))
        allocate(x6(2,npttmp))
        call normVec(x1,npttmp)
        call normVec(x2,npttmp)
        call normVec(x3,npttmp)
        call normVec(x4,npttmp)
        call normVec(x5,npttmp)
        call normVec(x6,npttmp)
        !print*,"sig1: ",sig1,sig2,sig3,sig4,sig5,sig6
        do ii = 1, npttmp
          !x-px:
!          call normdv(x1)
!          call normdv(x2) 
!         Correct Gaussian distribution.
          a(1,1) = xmu1 + sig1*x1(1,ii)/sq12
          a(2,1) = xmu2 + sig2*(-muxpx*x1(1,ii)/sq12+x1(2,ii))
          a(1,2) = xmu1 + sig1*x2(1,ii)/sq12
          a(2,2) = xmu2 + sig2*(-muxpx*x2(1,ii)/sq12+x2(2,ii))
!         Rob's Gaussian distribution.
          !a(1,1) = xmu1 + sig1*x1(1)
          !a(2,1) = xmu2 + sig2*(muxpx*x1(1)+sq12*x1(2))
          !a(1,2) = xmu1 + sig1*x2(1)
          !a(2,2) = xmu2 + sig2*(muxpx*x2(1)+sq12*x2(2))
          !y-py
!          call normdv(x1)
!          call normdv(x2) 
!         Correct Gaussian distribution.
          a(3,1) = xmu3 + sig3*x3(1,ii)/sq34
          a(4,1) = xmu4 + sig4*(-muypy*x3(1,ii)/sq34+x3(2,ii))
          a(3,2) = xmu3 + sig3*x4(1,ii)/sq34
          a(4,2) = xmu4 + sig4*(-muypy*x4(1,ii)/sq34+x4(2,ii))
!         Rob's Gaussian distribution.
          !a(3,1) = xmu3 + sig3*x1(1)
          !a(4,1) = xmu4 + sig4*(muypy*x1(1)+sq34*x1(2))
          !a(3,2) = xmu3 + sig3*x2(1)
          !a(4,2) = xmu4 + sig4*(muypy*x2(1)+sq34*x2(2))
          !z-pz
!          call normdv(x1)
!          call normdv(x2) 
!         Correct Gaussian distribution.
          a(5,1) = xmu5 + sig5*x5(1,ii)/sq56
          a(6,1) = xmu6 + sig6*(-muzpz*x5(1,ii)/sq56+x5(2,ii))
          a(5,2) = xmu5 + sig5*x6(1,ii)/sq56
          a(6,2) = xmu6 + sig6*(-muzpz*x6(1,ii)/sq56+x6(2,ii))
!         Rob's Gaussian distribution.
          !a(5,1) = xmu5 + sig5*x1(1)
          !a(6,1) = xmu6 + sig6*(muzpz*x1(1)+sq56*x1(2))
          !a(5,2) = xmu5 + sig5*x2(1)
          !a(6,2) = xmu6 + sig6*(muzpz*x2(1)+sq56*x2(2))

          do jj = 1, intv

            pid = (ii-1)*2+jj
          if((npx.gt.1).and.(npy.gt.1)) then

          if((myidx.eq.0).and.(myidy.eq.0)) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidx.eq.0) then
            if((a(5,jj).le.lcrange(6)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidx.eq.(npx-1)) then
            if((a(5,jj).gt.lcrange(5)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidy.eq.0) then
            if((a(3,jj).le.lcrange(4)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidy.eq.(npy-1)) then
            if((a(3,jj).gt.lcrange(3)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          else
            if( ((a(3,jj).gt.lcrange(3)).and.(a(3,jj).le.lcrange(4))) &
               .and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(6,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 6
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(6,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 6
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          endif

          else if(npx.gt.1) then
            if(myidx.eq.0) then
              if(a(5,jj).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            else if(myidx.eq.(npx-1)) then
              if(a(5,jj).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            else 
              if((a(5,jj).gt.lcrange(5)).and.(a(5,jj).le.lcrange(6))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            endif
          else if(npy.gt.1) then
            if(myidy.eq.0) then
              if(a(3,jj).le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            else if(myidy.eq.(npy-1)) then
              if(a(3,jj).gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            else
              if((a(3,jj).gt.lcrange(3)).and.(a(3,jj).le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  this%Pts1(j,numpts0) = a(j,jj)
                enddo

                if(mod(numpts0,avgpts).eq.0) then
                  allocate(ptstmp(6,numpts0))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      ptstmp(j,i0) = this%Pts1(j,i0)
                    enddo
                  enddo
                  deallocate(this%Pts1)
                  allocate(this%Pts1(6,numpts0+avgpts))
                  do i0 = 1, numpts0
                    do j = 1, 6
                      this%Pts1(j,i0) = ptstmp(j,i0)
                    enddo
                  enddo
                  deallocate(ptstmp)
                endif
              endif
            endif
          else
            numpts0 = numpts0 + 1
            do j = 1, 6
              this%Pts1(j,numpts0) = a(j,jj)
            enddo
          endif

          enddo
        enddo

        deallocate(x1)
        deallocate(x2)
        deallocate(x3)
        deallocate(x4)
        deallocate(x5)
        deallocate(x6)
          
!        call MPI_BARRIER(comm2d,ierr)
        allocate(ptstmp(9,numpts0))
        do i0 = 1, numpts0
          do j = 1, 6
            ptstmp(j,i0) = this%Pts1(j,i0)
          enddo
        enddo
        deallocate(this%Pts1)
        allocate(this%Pts1(9,numpts0))
        do i0 = 1, numpts0
          do j = 1, 6
            this%Pts1(j,i0) = ptstmp(j,i0)
          enddo
!            this%Pts1(5,i0) = -ptstmp(5,i0)
!            this%Pts1(6,i0) = -ptstmp(6,i0)
        enddo
        deallocate(ptstmp)

        this%Nptlocal = numpts0
!        print*,"numpts0: ",numpts0

!        call MPI_BARRIER(comm2d,ierr)
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss1_Dist

        ! sample the particles with intial TRUNCATED GAUSS distribution 
        ! using rejection method. 
        subroutine Gauss2_Dist(this,nparam,distparam,geom,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r,r1,r2,x0,y0,x1,x2
        double precision :: zlcmin,zlcmax,ylcmin,ylcmax,fvalue,z0
        double precision :: tmp1,tmp2,xrange,pxrange,yrange,pyrange,zrange,pzrange
        integer :: totnp,npy,npx
        integer :: avgpts,numpts,isamz,isamy
        integer :: myid,myidx,myidy
!        integer seedarray(1)
        double precision :: t0,x11

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        seedarray(1)=(1021+myid)*(myid+7)
!        seedarray(1)=(1121+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        xrange = 4.0*sig1/sqrt(1.0-muxpx*muxpx)
        pxrange = 4.0*sig2/sqrt(1.0-muxpx*muxpx)
        yrange = 4.0*sig3/sqrt(1.0-muypy*muypy)
        pyrange = 4.0*sig4/sqrt(1.0-muypy*muypy)
        zrange = 4.0*sig5/sqrt(1.0-muzpz*muzpz)
        pzrange = 4.0*sig6/sqrt(1.0-muzpz*muzpz)

        call getlcrange_CompDom(geom,lcrange)
        zlcmin = lcrange(5)
        zlcmax = lcrange(6)
        ylcmin = lcrange(3)
        ylcmax = lcrange(4)

        ! scale z and y for sampling purpose.
        zlcmin = zlcmin*sqrt(1.0-muzpz*muzpz)/sig5
        zlcmax = zlcmax*sqrt(1.0-muzpz*muzpz)/sig5
        ylcmin = ylcmin*sqrt(1.0-muypy*muypy)/sig3
        ylcmax = ylcmax*sqrt(1.0-muypy*muypy)/sig3
        if(zlcmax.le.0) then
          z0 = zlcmax
        else if(zlcmin.ge.0) then
          z0 = zlcmin
        else if((zlcmax.gt.0).and.(zlcmin.lt.0)) then
          z0 = 0.0
        else
        endif
        if(ylcmax.le.0) then
          y0 = ylcmax
        else if(ylcmin.ge.0) then
          y0 = ylcmin
        else if((ylcmax.gt.0).and.(ylcmin.lt.0)) then
          y0 = 0.0
        else
        endif
!        write(6,1234)myid,myidx,myidy,zlcmin,zlcmax,ylcmin,ylcmax,z0,y0
 1234   format(3i4,1x,6(1pe10.3,1x))
 
        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
        numpts = 0
        isamz = 0
        isamy = 0
        do 
          ! rejection sample.
10        call random_number(r)
          r1 = zlcmin + r*(zlcmax-zlcmin)
          fvalue = exp(-0.5*(r1*r1-z0*z0))
          call random_number(r2)
          isamz = isamz + 1
          if(r2.gt.fvalue) goto 10
          x1 = r1 
20        call random_number(r)
          r1 = ylcmin + r*(ylcmax-ylcmin)
          fvalue = exp(-0.5*(r1*r1-y0*y0))
          call random_number(r2)
          isamy = isamy + 1
          if(r2.gt.fvalue) goto 20
          x2 = r1 

          numpts = numpts + 1
          if(numpts.gt.avgpts) exit
!          tmp1 = sig5*x1
!          tmp2 = sig3*x2
!          if((abs(tmp1).gt.zrange).or.(abs(tmp2).gt.yrange)) goto 10
          !z-y
          this%Pts1(5,numpts) = xmu5 + sig5*x1/rootz
          this%Pts1(3,numpts) = xmu3 + sig3*x2/rooty
          !rob's distribution:
          !this%Pts1(5,numpts) = xmu5 + sig5*x1
          !this%Pts1(3,numpts) = xmu3 + sig3*x2
          !pz-py
30        call normdv(gs)
          tmp1 = sig6*(-muzpz*x1/rootz+gs(1))
          tmp2 = sig4*(-muypy*x2/rooty+gs(2))
          if((abs(tmp1).gt.pzrange).or.(abs(tmp2).gt.pyrange)) goto 30
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x1/rootz+gs(1))
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x2/rooty+gs(2))
          !rob's distribution:
          !this%Pts1(6,numpts) = xmu6 + sig6*(muzpz*x1+rootz*gs(1))
          !this%Pts1(4,numpts) = xmu4 + sig4*(muypy*x2+rooty*gs(2))
          !x-px
40        call normdv(gs)
          tmp1 = sig1*gs(1)/rootx
          tmp2 = sig2*(-muxpx*gs(1)/rootx+gs(2))
          if((abs(tmp1).gt.xrange).or.(abs(tmp2).gt.pxrange)) goto 40
          this%Pts1(1,numpts) = xmu1 + sig1*gs(1)/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*gs(1)/rootx+gs(2))
          !rob's distribution.
          !this%Pts1(1,numpts) = xmu1 + sig1*gs(1)
          !this%Pts1(2,numpts) = xmu2 + sig2*(muxpx*gs(1)+rootx*gs(2))
        enddo
          
        this%Nptlocal = avgpts
       
!        print*,avgpts,isamz,isamy

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss2_Dist

        subroutine Gauss3_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp
!        integer seedarray(1)
        double precision :: t0,x11,pid

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)

        avgpts = this%Npt/(npx*npy)

        if(mod(avgpts,10).ne.0) then
          print*,"The number of particles has to be an integer multiple of 10Nprocs"
          stop
        endif
        
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

!        allocate(x1(2,avgpts))
!        allocate(x2(2,avgpts))
!        allocate(x3(2,avgpts))
!        call normVec(x1,avgpts)
!        call normVec(x2,avgpts)
!        call normVec(x3,avgpts)
        
        intvsamp = 10
        allocate(x1(2,intvsamp))
        allocate(x2(2,intvsamp))
        allocate(x3(2,intvsamp))

        do j = 1, avgpts/intvsamp
          call normVec(x1,intvsamp)
          call normVec(x2,intvsamp)
          call normVec(x3,intvsamp)
          do k = 1, intvsamp
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq56
            this%Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,k)/sq56+x3(2,k))
          enddo
        enddo
          
        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts

        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss3_Dist

        subroutine Gauss4_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, dimension(2) :: x1 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i
!        integer seedarray(1)
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1,&
                            ddx,ddy,ddz,ddx1,ddy1,ddz1

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0

        !dd = 0.000005
        !dd = 0.000010
        ddx = 0.00002
        !ddy = 2.0 ! for test2
        ddy = 0.0003 ! for viz
        !ddy = 0.0005
        !ddx = 0.00002
        !ddy = 0.00003
        !ddx = 0.00001
        !ddy = 0.00001
        ddz = 0.0
        !ddx1 = 0.000001
        !ddy1 = 0.000001
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0
        numpts = 0
        do
          !x-px:
10        call random_number(r)
!          r1 = -10*sig1 + r*20*sig1
          r1 = -20*sig1 + r*40*sig1
          call random_number(r)
!          r2 = -10*sig2 + r*20*sig2
          r2 = -20*sig2 + r*40*sig2
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*r2*muxpx/sig1/sig2+r2*r2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+0.1*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+0.1*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+ddx*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*(r1/sig1)**3+&
                   ddx1*(r1/sig1)**5)*muxpx/sig1/sig2+&
                   (r2+ddx*(r1/sig1)**3+ddx1*(r1/sig1)**5)**2/sig2/sig2)/2)
          call random_number(r3)
          if(r3.gt.fvalue) goto 10
          xr1 = r1
          px1 = r2
          !y-py
20        call random_number(r)
          r1 = -20*sig3 + r*40*sig3
          call random_number(r)
          r2 = -20*sig4 + r*40*sig4
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*r2*muypy/sig3/sig4+r2*r2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-0.0001*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-0.0001*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-ddy*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
! used for previous 3 case studies.
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*(r1/sig3)**3-&
!                   ddy1*(r1/sig3)**5)*muypy/sig3/sig4+&
!                   (r2-ddy*(r1/sig3)**3-ddy1*(r1/sig3)**5)**2/sig4/sig4)/2)
!test 1
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3-r1 &
!                   )*muypy/sig3/sig4+&
!                   (ddy*(r2/sig4)**3-r1)**2/sig4/sig4)/2)
!test 2
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3+r2 &
!                   )*muypy/sig3/sig4+&
!                   (ddy*(r2/sig4)**3+r2)**2/sig4/sig4)/2)
!test 3
          fvalue = exp(-((r1-ddy*(r2/sig4)**3)**2/sig3/sig3+2*&
                   (r1-ddy*(r2/sig4)**3)*r2*muypy/sig3/sig4+&
                   r2**2/sig4/sig4)/2)
          call random_number(r3)
          if(r3.gt.fvalue) goto 20
          y1 = r1
          py1 = r2
          numpts = numpts + 1
!          this%Pts1(1,numpts) = xmu1 + xr1*16.8*1.75
!          this%Pts1(2,numpts) = xmu2 + px1/2.5/1.75
!          this%Pts1(3,numpts) = xmu3 + y1*18.7*1.75
!          this%Pts1(4,numpts) = xmu4 + py1/3.9/1.75
          this%Pts1(1,numpts) = xmu1 + xr1
          this%Pts1(2,numpts) = xmu2 + px1
          this%Pts1(3,numpts) = xmu3 + y1
          this%Pts1(4,numpts) = xmu4 + py1
          !z-pz
          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(5,numpts) = xmu5 + sig5*x1(1)/sq56
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x1(1)/sq56+x1(2))
          if(numpts.ge.avgpts) exit
        enddo
        print*,"numpts: ",numpts
          
        this%Nptlocal = avgpts
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss4_Dist

        subroutine Gauss4new_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        double precision, dimension(3,3) :: rr
        double precision, dimension(3) :: frac
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,pid
        integer seedarray(1)
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy
        double precision :: r11x,r22x,r12x
        double precision :: r11y,r22y,r12y
        double precision :: r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1
        double precision, allocatable, dimension(:,:) :: x1 
        double precision, allocatable, dimension(:) :: ranumx,ranumy,ranum2x,ranum2y 
        integer :: isamx,isamy,iranx,irany,intvsamp

        call starttime_Timer(t0)

        frac = 1.0
        rr = 1.0
        rr(3,:) = 0.0

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)
        !print*,myid,x11

        avgpts = this%Npt/(npx*npy)
        if(mod(avgpts,10).ne.0) then
          print*,"The number of particles has to be an integer multiple of 10Nprocs" 
          stop
        endif
 
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0

        factx = frac(1)
        facty = frac(2)
        factz = frac(3)

        !ddx = 0.0002
        !ddy = 0.0002
        !ddx = 0.00002 !for test2 and 3
        !ddx = 0.0000 !from T3E
        ddx = 0.00004
        !ddy = 2.0 !for test2
        !ddy = 0.0003 !for test3
        !ddy = 0.0001 !from T3E
        ddy = 0.0003 
        ddz = 0.0
        !ddx1 = 0.0002 ! for emitdistort2
        !ddy1 = 0.0002 ! for emitdistort2
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0
        numpts = 0

       
        !intvsamp = avgpts
        intvsamp = 10
        allocate(ranum2x(2*intvsamp))
        allocate(ranum2y(2*intvsamp))
        allocate(ranumx(intvsamp))
        allocate(ranumy(intvsamp))
        allocate(x1(2,intvsamp))
        call normVec(x1,intvsamp)

        isamx = 0
        isamy = 0
        do
          !x-px:
10        continue 
          isamx = isamx + 1
          if(mod(isamx-1,intvsamp).eq.0) then
            call random_number(ranum2x)
            call random_number(ranumx)
          endif
          iranx = 2*mod(isamx-1,intvsamp)
!          call random_number(r)
          r1 = -20*sig1 + ranum2x(iranx+1)*40*sig1
!          call random_number(r)
          r2 = -20*sig2 + ranum2x(iranx+2)*40*sig2
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*r2*muxpx/sig1/sig2+r2*r2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+0.1*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+0.1*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+ddx*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*(r1/sig1)**3+&
                   ddx1*(r1/sig1)**5)*muxpx/sig1/sig2+&
                   (r2+ddx*(r1/sig1)**3+ddx1*(r1/sig1)**5)**2/sig2/sig2)/2)
!          call random_number(r3)
          r3 = ranumx(iranx/2+1)
          if(r3.gt.fvalue) goto 10
          xr1 = r1
          px1 = r2
          !y-py
20        continue
          isamy = isamy + 1
          if(mod(isamy-1,intvsamp).eq.0) then
            call random_number(ranum2y)
            call random_number(ranumy)
          endif
          irany = 2*mod(isamy-1,intvsamp)
!          call random_number(r)
          r1 = -20*sig3 + ranum2y(irany+1)*40*sig3
!          call random_number(r)
          r2 = -20*sig4 + ranum2y(irany+2)*40*sig4
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*r2*muypy/sig3/sig4+r2*r2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-0.0001*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-0.0001*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-ddy*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
!used for the previous 3 case studies.
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*(r1/sig3)**3-&
!                   ddy1*(r1/sig3)**5)*muypy/sig3/sig4+&
!                   (r2-ddy*(r1/sig3)**3-ddy1*(r1/sig3)**5)**2/sig4/sig4)/2)
!test 1
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3-r1 &
!                   )*muypy/sig3/sig4+&
!                   (ddy*(r2/sig4)**3-r1)**2/sig4/sig4)/2)
!test 2
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3+r2 &
!                   )*muypy/sig3/sig4+&
!                   (ddy*(r2/sig4)**3+r2)**2/sig4/sig4)/2)
!test 3
          fvalue = exp(-((r1-ddy*(r2/sig4)**3-ddy1*(r2/sig4)**5)**2/sig3/sig3+2*&
                   (r1-ddy*(r2/sig4)**3-ddy1*(r2/sig4)**5)*r2*muypy/sig3/sig4+&
                   r2**2/sig4/sig4)/2)

          r3 = ranumy(irany/2+1)
!          call random_number(r3)
          if(r3.gt.fvalue) goto 20
          y1 = r1
          py1 = r2
          numpts = numpts + 1
          this%Pts1(1,numpts) = xmu1 + xr1*factx
          this%Pts1(2,numpts) = xmu2 + px1*factx
          this%Pts1(3,numpts) = xmu3 + y1*facty
          this%Pts1(4,numpts) = xmu4 + py1*facty
          !z-pz
!          call normdv(x1)
          this%Pts1(5,numpts) = xmu5 + sig5*x1(1,numpts)/sq56*factz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x1(1,numpts)/sq56+x1(2,numpts))*factz
          if(numpts.ge.avgpts) exit
        enddo
        
        deallocate(ranum2x)
        deallocate(ranum2y)
        deallocate(ranumx)
        deallocate(ranumy)
        deallocate(x1)
          
        this%Nptlocal = avgpts

        r11x = rr(1,1)
        r22x = rr(2,1)
        r12x = rr(3,1) 
        r11y = rr(1,2)
        r22y = rr(2,2)
        r12y = rr(3,2) 

        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xmu1 + xxx
          this%Pts1(2,i) = xmu2 + yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xmu3 + xxx
          this%Pts1(4,i) = xmu4 + yyy
        enddo

        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        end subroutine Gauss4new_Dist

        subroutine Gauss5_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, dimension(2) :: x1 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i
!        integer seedarray(1)
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy
        double precision :: al0x,al1x,ga0x,ga1x,b0x,b1x,r11x,r22x,r12x
        double precision :: al0y,al1y,ga0y,ga1y,b0y,b1y,r11y,r22y,r12y
        double precision :: al0z,al1z,ga0z,ga1z,b0z,b1z,r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0

        !used for the original comparison at beginning
        factx = 0.93894
        facty = 0.72253
        !factx = 0.971047
        !facty = 0.793620
        !used for 3rd order comparison at quad 20
        !factx = 0.95776371
        !facty = 0.85342924
        ! used for 5th order comparison at quad 20
        !factx = 0.8466795
        !facty = 0.8491518
        factz = 1.0
        !dd = 0.000005
        !dd = 0.000010
        !used for the original comparison at beginning
        ddx = 0.0002
        ddy = 0.0002
        !used for 3rd order comparison at quad 20
        !ddx = 0.00002
        !ddy = 0.00003
        ! used for 5th order comparison at quad 20
        !ddx = 0.00001
        !ddy = 0.00001
        !ddx = 0.0
        !ddy = 0.0
        ddz = 0.0
        !ddx1 = 0.000001
        !ddy1 = 0.000001
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0
        numpts = 0

!test 2
        ddx = 0.00002
        ddy = 2.0
        factx = 0.95618 
        facty = 6.55394
        do
          !x-px:
10        call random_number(r)
          !used for the original comparison at beginning
!          r1 = -10*sig1 + r*20*sig1
          ! used for 5th order comparison at quad 20
          r1 = -20*sig1 + r*40*sig1
          call random_number(r)
          !used for the original comparison at beginning
!          r2 = -10*sig2 + r*20*sig2
          ! used for 5th order comparison at quad 20
          r2 = -20*sig2 + r*40*sig2
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*r2*muxpx/sig1/sig2+r2*r2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+0.1*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+0.1*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
!          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*r1*r1*r1/sig1/sig1/sig1)*muxpx/sig1/sig2+&
!                        (r2+ddx*r1*r1*r1/sig1/sig1/sig1)**2/sig2/sig2)/2)
          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*(r1/sig1)**3+&
                   ddx1*(r1/sig1)**5)*muxpx/sig1/sig2+&
                   (r2+ddx*(r1/sig1)**3+ddx1*(r1/sig1)**5)**2/sig2/sig2)/2)
          call random_number(r3)
          if(r3.gt.fvalue) goto 10
          xr1 = r1
          px1 = r2
          !y-py
20        call random_number(r)
          r1 = -20*sig3 + r*40*sig3
          call random_number(r)
          r2 = -20*sig4 + r*40*sig4
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*r2*muypy/sig3/sig4+r2*r2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-0.0001*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-0.0001*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*r1*r1*r1/sig3/sig3/sig3)*muypy/sig3/sig4+&
!                         (r2-ddy*r1*r1*r1/sig3/sig3/sig3)**2/sig4/sig4)/2)
!used for the previous 3 case studies.
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(r2-ddy*(r1/sig3)**3-&
!                   ddy1*(r1/sig3)**5)*muypy/sig3/sig4+&
!                   (r2-ddy*(r1/sig3)**3-ddy1*(r1/sig3)**5)**2/sig4/sig4)/2)
!test 1
!          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3-r1 &
!                   )*muypy/sig3/sig4+&
!                   (ddy*(r2/sig4)**3-r1)**2/sig4/sig4)/2)
!test 2
          fvalue = exp(-(r1*r1/sig3/sig3+2*r1*(ddy*(r2/sig4)**3+r2 &
                   )*muypy/sig3/sig4+&
                   (ddy*(r2/sig4)**3+r2)**2/sig4/sig4)/2)

          call random_number(r3)
          if(r3.gt.fvalue) goto 20
          y1 = r1
          py1 = r2
          numpts = numpts + 1
!          this%Pts1(1,numpts) = xmu1 + xr1*16.8*1.75
!          this%Pts1(2,numpts) = xmu2 + px1/2.5/1.75
!          this%Pts1(3,numpts) = xmu3 + y1*18.7*1.75
!          this%Pts1(4,numpts) = xmu4 + py1/3.9/1.75
          this%Pts1(1,numpts) = xmu1 + xr1*factx
          this%Pts1(2,numpts) = xmu2 + px1*factx
          this%Pts1(3,numpts) = xmu3 + y1*facty
          this%Pts1(4,numpts) = xmu4 + py1*facty
          !z-pz
          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(5,numpts) = xmu5 + sig5*x1(1)/sq56*factz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x1(1)/sq56+x1(2))*factz
          if(numpts.ge.avgpts) exit
        enddo
          
        this%Nptlocal = avgpts
        !al0x = 1.16489
        !al1x = 1.56003
        !b0x = 2.123973
        !b1x = 1.977976
        al0x = 1.41627
        !al1x = 1.91224
        al1x = 1.66392
        b0x = 3.215056
        !b1x = 2.779062
        b1x = 2.237308
! test 2
        al0x = 1.16954
        al1x = 1.59331
        b0x = 2.12649
        b1x = 1.96820
        ga0x = (1.0+al0x*al0x)/b0x*0.1363
        ga1x = (1.0+al1x*al1x)/b1x*0.1363
        r11x = sqrt(ga1x/ga0x)
        r22x = 1.0/r11x
        r12x = (al1x-al0x)/sqrt(ga1x*ga0x) 
        al0y = -1.41952
        !al1y = -1.94419
        al1y = -1.65727
        b0y = 3.236013
        !b1y = 2.318362
        b1y = 2.2658332 
! test 2
        al0y = -1.63889
        b0y = 4.58248
        al1y = -0.77426
        b1y = 101.5378
        ga0y = (1.0+al0y*al0y)/b0y*0.1363
        ga1y = (1.0+al1y*al1y)/b1y*0.1363
        r11y = sqrt(ga1y/ga0y)
        r22y = 1.0/r11y
        r12y = (al1y-al0y)/sqrt(ga1y*ga0y) 
!        r11x = 1.0
!        r12x = 0.0
!        r22x = 1.0
!        r11y = 1.0
!        r12y = 0.0
!        r22y = 1.0
        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xmu1 + xxx
          this%Pts1(2,i) = xmu2 + yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xmu3 + xxx
          this%Pts1(4,i) = xmu4 + yyy
        enddo
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss5_Dist

        subroutine Gauss6_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, dimension(2) :: x1 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i
!        integer seedarray(1)
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy
        double precision :: al0x,al1x,ga0x,ga1x,b0x,b1x,r11x,r22x,r12x
        double precision :: al0y,al1y,ga0y,ga1y,b0y,b1y,r11y,r22y,r12y
        double precision :: al0z,al1z,ga0z,ga1z,b0z,b1z,r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1
        double precision, dimension(3) :: frac,al0,ga0,epson0,al1,ga1,epson1
        double precision, dimension(3,3) :: rr

        call starttime_Timer(t0)

        frac = 1.0
        rr = 1.0
        rr(3,:) = 0.0

        call Waterbag_Dist(this,nparam,distparam,grid,0)
        call gammaepson_Dist(this,al0,ga0,epson0)
        !print*,"a10: ",al0
        !print*,"epson0: ",epson0
        call Gaussdistort_Dist(this,nparam,distparam,grid,rr,frac,0)
        call gammaepson_Dist(this,al1,ga1,epson1)
        !print*,"a11: ",al1
        !print*,"epson1: ",epson1

        do i = 1, 3
          frac(i) = sqrt(epson0(i)/epson1(i))
          rr(1,i) = sqrt(ga1(i)/ga0(i))
          rr(2,i) = 1.0/rr(1,i)
          rr(3,i) = (al1(i)-al0(i))/sqrt(ga1(i)*ga0(i))
        enddo
        !no distort for z-pz yet.
        !frac(3) = 1.0

        call Gaussdistort_Dist(this,nparam,distparam,grid,rr,frac,1)

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss6_Dist

        subroutine Gaussdistort_Dist(this,nparam,distparam,grid,rr,frac,resamp)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,resamp
        double precision, dimension(nparam) :: distparam
        double precision, dimension(3,3) :: rr
        double precision, dimension(3) :: frac
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1,z1,pz1
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy
        double precision :: r11x,r22x,r12x
        double precision :: r11y,r22y,r12y
        double precision :: r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1,xx

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call random_number(x11)
        !print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        if(resamp.eq.1) then
          call Waterbag_Dist(this,nparam,distparam,grid,1)
        endif

        ddx = -4.0e2 !for alphax > 0.0
        ddx = -8.0e2 !for alphax > 0.0
        !ddy = -3.0e2 
        ddy = 1.0e7  !for alphay < 0.0
        ddy = 3.0e7  !for alphay < 0.0
        !ddy = 6.0e7  !for alphay < 0.0 for the purpose of display
        !ddz = -5.0e-3
        !ddz = 5.0e11 !for alphaz < 0.0
        ddz = 0.0
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0
        do i = 1, avgpts
          xx = this%Pts1(1,i) - xmu1
          this%Pts1(2,i) = this%Pts1(2,i) + ddx*xx**3
          !for alpha > 0
          !xx = this%Pts1(3,i) - xmu3
          !this%Pts1(4,i) = this%Pts1(4,i) + ddy*xx**3
          !for alpha < 0
          xx = this%Pts1(4,i) - xmu4
          this%Pts1(3,i) = this%Pts1(3,i) + ddy*xx**3
          !xx = this%Pts1(5,i) - xmu5
          !this%Pts1(6,i) = this%Pts1(6,i) + ddz*xx**3
          xx = this%Pts1(6,i) - xmu6
          this%Pts1(5,i) = this%Pts1(5,i) + ddz*xx**3
        enddo

        factx = frac(1)
        facty = frac(2)
        factz = frac(3)

        do i = 1, avgpts
          this%Pts1(1,i) = this%Pts1(1,i)*factx
          this%Pts1(2,i) = this%Pts1(2,i)*factx
          this%Pts1(3,i) = this%Pts1(3,i)*facty
          this%Pts1(4,i) = this%Pts1(4,i)*facty
          this%Pts1(5,i) = this%Pts1(5,i)*factz
          this%Pts1(6,i) = this%Pts1(6,i)*factz
        enddo
        
        this%Nptlocal = avgpts

        r11x = rr(1,1)
        r22x = rr(2,1)
        r12x = rr(3,1) 
        r11y = rr(1,2)
        r22y = rr(2,2)
        r12y = rr(3,2) 
        r11z = rr(1,3)
        r22z = rr(2,3)
        r12z = rr(3,3)

        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xxx
          this%Pts1(2,i) = yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xxx
          this%Pts1(4,i) = yyy
          xxx = this%Pts1(5,i)*r11z + this%Pts1(6,i)*r12z
          yyy = this%Pts1(6,i)*r22z
          this%Pts1(5,i) = xxx
          this%Pts1(6,i) = yyy
        enddo
       
        end subroutine Gaussdistort_Dist

        subroutine Gaussdistortold_Dist(this,nparam,distparam,grid,rr,frac)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        double precision, dimension(3,3) :: rr
        double precision, dimension(3) :: frac
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i
!        integer seedarray(1)
        double precision :: t0,x11,r1,r2,r3,fvalue,r,xr1,px1,y1,py1,z1,pz1
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy
        double precision :: r11x,r22x,r12x
        double precision :: r11y,r22y,r12y
        double precision :: r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1
        double precision, allocatable, dimension(:) :: ranumx,ranumy,&
                          ranum2x,ranum2y,ranumz,ranum2z 
        integer :: isamx,isamy,iranx,irany,intvsamp,isamz,iranz

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)
        !print*,myid,x11

        avgpts = this%Npt/(npx*npy)
        if(mod(avgpts,10).ne.0) then
          print*,"The number of particles has to be an integer multiple of 10Nprocs" 
          stop
        endif
 
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        !allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0

        factx = frac(1)
        facty = frac(2)
        factz = frac(3)

        !ddx = 0.0002
        !ddy = 0.0002
        !ddx = 0.00002 !for test2 and 3
        !ddx = 0.0000 !from T3E
        ddx = 0.00004
        !ddy = 2.0 !for test2
        !ddy = 0.0003 !for test3
        !ddy = 0.0001 !from T3E
        ddy = 0.0003 
        !ddz = 0.3
        !ddz = 0.0005
        ddz = 3.0005
        !ddx1 = 0.0002 ! for emitdistort2
        !ddy1 = 0.0002 ! for emitdistort2
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0
        numpts = 0
       
        !intvsamp = avgpts
        intvsamp = 10
        allocate(ranum2x(2*intvsamp))
        allocate(ranum2y(2*intvsamp))
        allocate(ranum2z(2*intvsamp))
        allocate(ranumx(intvsamp))
        allocate(ranumy(intvsamp))
        allocate(ranumz(intvsamp))

        isamx = 0
        isamy = 0
        isamz = 0
        do
          !x-px:
10        continue 
          isamx = isamx + 1
          if(mod(isamx-1,intvsamp).eq.0) then
            call random_number(ranum2x)
            call random_number(ranumx)
          endif
          iranx = 2*mod(isamx-1,intvsamp)
          r1 = -20*sig1 + ranum2x(iranx+1)*40*sig1
          r2 = -20*sig2 + ranum2x(iranx+2)*40*sig2
          fvalue = exp(-(r1*r1/sig1/sig1+2*r1*(r2+ddx*(r1/sig1)**3+&
                   ddx1*(r1/sig1)**5)*muxpx/sig1/sig2+&
                   (r2+ddx*(r1/sig1)**3+ddx1*(r1/sig1)**5)**2/sig2/sig2)/2)
          r3 = ranumx(iranx/2+1)
          if(r3.gt.fvalue) goto 10
          xr1 = r1
          px1 = r2
          !y-py
20        continue
          isamy = isamy + 1
          if(mod(isamy-1,intvsamp).eq.0) then
            call random_number(ranum2y)
            call random_number(ranumy)
          endif
          irany = 2*mod(isamy-1,intvsamp)
          r1 = -20*sig3 + ranum2y(irany+1)*40*sig3
          r2 = -20*sig4 + ranum2y(irany+2)*40*sig4
          fvalue = exp(-((r1-ddy*(r2/sig4)**3-ddy1*(r2/sig4)**5)**2/sig3/sig3+2*&
                   (r1-ddy*(r2/sig4)**3-ddy1*(r2/sig4)**5)*r2*muypy/sig3/sig4+&
                   r2**2/sig4/sig4)/2)

          r3 = ranumy(irany/2+1)
          if(r3.gt.fvalue) goto 20
          y1 = r1
          py1 = r2
          !z-pz
30        continue
          isamz = isamz + 1
          if(mod(isamz-1,intvsamp).eq.0) then
            call random_number(ranum2z)
            call random_number(ranumz)
          endif
          iranz = 2*mod(isamz-1,intvsamp)
          r1 = -20*sig5 + ranum2z(iranz+1)*40*sig5
          r2 = -20*sig6 + ranum2z(iranz+2)*40*sig6
          fvalue = exp(-((r1-ddz*(r2/sig5)**3-ddz1*(r2/sig6)**5)**2/sig5/sig5+2*&
                   (r1-ddz*(r2/sig6)**3-ddz1*(r2/sig6)**5)*r2*muzpz/sig5/sig6+&
                   r2**2/sig6/sig6)/2)
          r3 = ranumz(iranz/2+1)
          if(r3.gt.fvalue) goto 30
          z1 = r1
          pz1 = r2

          numpts = numpts + 1
          this%Pts1(1,numpts) = xmu1 + xr1*factx
          this%Pts1(2,numpts) = xmu2 + px1*factx
          this%Pts1(3,numpts) = xmu3 + y1*facty
          this%Pts1(4,numpts) = xmu4 + py1*facty
          this%Pts1(5,numpts) = xmu5 + z1*factz
          this%Pts1(6,numpts) = xmu6 + pz1*factz

          if(numpts.ge.avgpts) exit
        enddo
        
        deallocate(ranum2x)
        deallocate(ranum2y)
        deallocate(ranum2z)
        deallocate(ranumx)
        deallocate(ranumy)
        deallocate(ranumz)
          
        this%Nptlocal = avgpts

        r11x = rr(1,1)
        r22x = rr(2,1)
        r12x = rr(3,1) 
        r11y = rr(1,2)
        r22y = rr(2,2)
        r12y = rr(3,2) 
        r11z = rr(1,3)
        r22z = rr(2,3)
        r12z = rr(3,3)

        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xmu1 + xxx
          this%Pts1(2,i) = xmu2 + yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xmu3 + xxx
          this%Pts1(4,i) = xmu4 + yyy
          xxx = this%Pts1(5,i)*r11z + this%Pts1(6,i)*r12z
          yyy = this%Pts1(6,i)*r22z
          this%Pts1(5,i) = xmu5 + xxx
          this%Pts1(6,i) = xmu6 + yyy
        enddo
       
        end subroutine Gaussdistortold_Dist

        subroutine gammaepson_Dist(this,CSalpha,CSgamma,CSepson)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(in) :: this
        double precision, dimension(3), intent(out) :: CSalpha,CSgamma,CSepson
        integer :: innp
        integer*8 :: nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,yrms
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0
        double precision ::sqx,sqpx,sqy,sqpy,sqz,sqpz
        integer :: i,ierr
        double precision:: qmc,xl,xt,pi
        double precision, dimension(15) :: tmplc,tmpgl
        double precision, dimension(3)  :: CSbeta

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg

        innp = this%Nptlocal
        nptot = this%Npt

        den1 = 1.0/float(nptot)
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
        pi = 2*asin(1.0)

        do i = 1, innp
          x0lc = x0lc + this%Pts1(1,i)
          sqsum1local = sqsum1local + this%Pts1(1,i)*this%Pts1(1,i)
          xpxlocal = xpxlocal + this%Pts1(1,i)*this%Pts1(2,i)
          px0lc = px0lc + this%Pts1(2,i)
          sqsum2local = sqsum2local + this%Pts1(2,i)*this%Pts1(2,i)
          y0lc = y0lc + this%Pts1(3,i)
          sqsum3local = sqsum3local + this%Pts1(3,i)*this%Pts1(3,i)
          ypylocal = ypylocal + this%Pts1(3,i)*this%Pts1(4,i)
          py0lc = py0lc + this%Pts1(4,i)
          sqsum4local = sqsum4local + this%Pts1(4,i)*this%Pts1(4,i)
          z0lc = z0lc + this%Pts1(5,i)
          sqsum5local = sqsum5local + this%Pts1(5,i)*this%Pts1(5,i)

          zpzlocal = zpzlocal + this%Pts1(5,i)*this%Pts1(6,i)
          pz0lc = pz0lc + this%Pts1(6,i)
          sqsum6local = sqsum6local + this%Pts1(6,i)*this%Pts1(6,i)
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

        epsx2 = (sqsum1*sqsum2-xpx*xpx)
        epsy2 = (sqsum3*sqsum4-ypy*ypy)
        epsz2 = (sqsum5*sqsum6-zpz*zpz)
        epx = sqrt(max(epsx2,0.0d0))
        epy = sqrt(max(epsy2,0.0d0))
        epz = sqrt(max(epsz2,0.0d0))
        xrms = sqrt(sqsum1)
        yrms = sqrt(sqsum3)
        zrms = sqrt(sqsum5)

        CSalpha(1) = -xpx/epx
        CSalpha(2) = -ypy/epy
        CSalpha(3) = -zpz/epz
        CSbeta(1) = (xrms*xl)**2/(epx*xl)
        CSbeta(2) = (yrms*xl)**2/(epy*xl)
        CSbeta(3) = (zrms*xt)**2/(epz*qmc*xt)
        
        CSgamma(:) = (1.0 + CSalpha(:)*CSalpha(:))/CSbeta(:)*xl
        CSgamma(3) = (1.0 + CSalpha(3)*CSalpha(3))/CSbeta(3)/qmc*xt

        CSepson(1) = epx*xl
        CSepson(2) = epy*xl
        CSepson(3) = epz*qmc*xt

        end subroutine gammaepson_Dist

        subroutine normdv(y)
        implicit none
        include 'mpif.h'
        double precision, dimension(2), intent(out) :: y
        double precision :: twopi,x1,x2,epsilon

        epsilon = 1.0e-18

        twopi = 4.0*asin(1.0)
        call random_number(x2)
10      call random_number(x1)
!        x1 = 0.5
!10      x2 = 0.6
        if(x1.eq.0.0) goto 10
!        if(x1.eq.0.0) x1 = epsilon
        y(1) = sqrt(-2.0*log(x1))*cos(twopi*x2)
        y(2) = sqrt(-2.0*log(x1))*sin(twopi*x2)

        end subroutine normdv

        subroutine normVec(y,num)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: num
        double precision, dimension(2,num), intent(out) :: y
        double precision :: twopi,epsilon
        double precision, dimension(num) :: x1,x2
        integer :: i

        epsilon = 1.0e-18

        twopi = 4.0*asin(1.0)
        call random_number(x2)
        call random_number(x1)
        do i = 1, num
          if(x1(i).eq.0.0) x1(i) = epsilon
          y(1,i) = sqrt(-2.0*log(x1(i)))*cos(twopi*x2(i))
          y(2,i) = sqrt(-2.0*log(x1(i)))*sin(twopi*x2(i))
        enddo

        end subroutine normVec

        ! sample the particles with intial distribution 
        ! using rejection method. 
        subroutine Waterbag_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx
        integer :: avgpts,numpts,isamz,isamy
        integer :: myid,myidx,myidy,iran,intvsamp,pid,j
!        integer seedarray(2)
        double precision :: t0,x11
        double precision, allocatable, dimension(:) :: ranum6

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        seedarray(2)=(101+2*myid)*(myid+4)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray)
        call random_number(x11)

        avgpts = this%Npt/(npx*npy)
        !if(mod(avgpts,10).ne.0) then
        !  print*,"The number of particles has to be an integer multiple of 10Nprocs" 
        !  stop
        !endif
 
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif
        numpts = 0
        isamz = 0
        isamy = 0
        intvsamp = avgpts
        !intvsamp = 10
        allocate(ranum6(6*intvsamp))

        do 
          ! rejection sample.
10        continue 
          isamz = isamz + 1
          if(mod(isamz-1,intvsamp).eq.0) then
            call random_number(ranum6)
          endif
          iran = 6*mod(isamz-1,intvsamp)
          r1 = 2.0*ranum6(iran+1)-1.0
          r2 = 2.0*ranum6(iran+2)-1.0
          r3 = 2.0*ranum6(iran+3)-1.0
          r4 = 2.0*ranum6(iran+4)-1.0
          r5 = 2.0*ranum6(iran+5)-1.0
          r6 = 2.0*ranum6(iran+6)-1.0
          if(r1**2+r2**2+r3**2+r4**2+r5**2+r6**2.gt.1.0) goto 10
          isamy = isamy + 1
          numpts = numpts + 1
          if(numpts.gt.avgpts) exit
!x-px:
          x1 = r1*sqrt(8.0)
          x2 = r2*sqrt(8.0)
          !Correct transformation.
          this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          x3 = r3*sqrt(8.0)
          x4 = r4*sqrt(8.0)
          !correct transformation
          this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          x5 = r5*sqrt(8.0)
          x6 = r6*sqrt(8.0)
          !correct transformation
          this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          !Rob's transformation
          !this%Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
          !this%Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
        enddo

        deallocate(ranum6)
          
        this%Nptlocal = avgpts
       
        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Waterbag_Dist

        subroutine KV3d_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,j,pid
!        integer seedarray(1)
        double precision :: t0,x11,twopi

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
        twopi = 4*asin(1.0)

        do numpts = 1, avgpts
          call random_number(r1)
          call random_number(r2)
          call random_number(r3)
          r4 = sqrt(r1)
          r5 = sqrt(1.0-r1)
          r2 = r2*twopi
          r3 = r3*twopi
          x1 = 2*r4*cos(r2)
          x2 = 2*r4*sin(r2)
          x3 = 2*r5*cos(r3)
          x4 = 2*r5*sin(r3)
!x-px:
          !Correct transformation.
          this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation.
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          !correct transformation
          this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          call random_number(r5)
          r5 = 2*r5 - 1.0
          call random_number(r6)
          r6 = 2*r6 - 1.0
          x5 = r5*sqrt(3.0)
          x6 = r6*sqrt(3.0)
          !correct transformation
!          this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
!          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          !Rob's transformation
          !this%Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
          !this%Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
          !this distribution uses chirp rotation and uncorrelated energy spread 
          this%Pts1(5,numpts) = xmu5 + sig5*x5*cos(muzpz)-sig6*x6*sin(muzpz)
          this%Pts1(6,numpts) = xmu6 + sig6*x6*cos(muzpz)+sig5*x5*sin(muzpz)
        enddo
          
        this%Nptlocal = avgpts
        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine KV3d_Dist

! This sampling does not work properly if one-dimensional PE is 1
        subroutine Uniform_Dist(this,nparam,distparam,geom,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(6,1) :: a
        double precision, dimension(2) :: x1, x2
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        double precision :: r1, r2
        integer :: totnp,npy,npx,myid,myidy,myidx,comm2d, &
                   commcol,commrow,ierr
        integer :: avgpts,numpts0,i,ii,i0,j,jj,pid
        double precision :: t0

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        avgpts = this%Npt/(npx*npy)

        call getlcrange_CompDom(geom,lcrange)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        numpts0 = 0

        ! initial allocate 'avgpts' particles on each processor.
        this%Pts1 = 0.0
    
        do ii = 1, this%Npt
          call random_number(r1)
          r1 = (2*r1 - 1.0)*sqrt(3.0)
          call random_number(r2)
          r2 = (2*r2 - 1.0)*sqrt(3.0)
          a(1,1) = xmu1 + sig1*r1/sq12
          !a(2,1) = xmu2 + sig2*(-muxpx*r2/sq12+r2)
          a(2,1) = xmu2 + sig2*(-muxpx*r1/sq12+r2)
          call random_number(r1)
          r1 = (2*r1 - 1.0)*sqrt(3.0)
          call random_number(r2)
          r2 = (2*r2 - 1.0)*sqrt(3.0)
          a(3,1) = xmu3 + sig3*r1/sq34
          !a(4,1) = xmu4 + sig4*(-muypy*r2/sq34+r2)
          a(4,1) = xmu4 + sig4*(-muypy*r1/sq34+r2)
          call random_number(r1)
          r1 = (2*r1 - 1.0)*sqrt(3.0)
          call random_number(r2)
          r2 = (2*r2 - 1.0)*sqrt(3.0)
          a(5,1) = xmu5 + sig5*r1/sq56
          !a(6,1) = xmu6 + sig6*(-muzpz*r2/sq56+r2)
          a(6,1) = xmu6 + sig6*(-muzpz*r1/sq56+r2)

          do jj = 1, 1

          if(totnp.ne.1) then

          if((myidx.eq.0).and.(myidy.eq.0)) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidx.eq.0) then
            if((a(5,jj).le.lcrange(6)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidx.eq.(npx-1)) then
            if((a(5,jj).gt.lcrange(5)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidy.eq.0) then
            if((a(3,jj).le.lcrange(4)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          elseif(myidy.eq.(npy-1)) then
            if((a(3,jj).gt.lcrange(3)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          else
            if( ((a(3,jj).gt.lcrange(3)).and.(a(3,jj).le.lcrange(4))) &
               .and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo

            endif
          endif

          else
            numpts0 = numpts0 + 1
            do j = 1, 6
              this%Pts1(j,numpts0) = a(j,jj)
            enddo
          endif

          enddo
        enddo
          
!        call MPI_BARRIER(comm2d,ierr)

        this%Nptlocal = numpts0
       
!        call MPI_BARRIER(comm2d,ierr)
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)
        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo


        end subroutine Uniform_Dist

! This sampling does not work properly if one-dimensional PE is 1
        subroutine Semigauss_Dist(this,nparam,distparam,geom,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy,yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(6,2) :: a
        double precision, dimension(3) :: x1, x2
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        double precision :: r1, r2, r, r3
        double precision,allocatable,dimension(:,:) :: ptstmp
        integer :: totnp,npy,npx,myid,myidy,myidx,comm2d, &
                   commcol,commrow,ierr
        integer :: avgpts,numpts0,i,ii,i0,j,jj,pid
        double precision :: t0

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        avgpts = this%Npt/(npx*npy)

        call getlcrange_CompDom(geom,lcrange)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        numpts0 = 0

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
! The performance inside this loop might be improved due to
! a lot of subroutine call in this loop.    
        do ii = 1, this%Npt
          ! rejection sample.
10        call random_number(r1)
          call random_number(r2)
          call random_number(r3)
          r1 = 2.0*r1-1.0
          r2 = 2.0*r2-1.0
          r3 = 2.0*r3-1.0
          if(r1**2+r2**2+r3**2.gt.1.0) goto 10
          x2(1) = r1
          x2(2) = r2
          x2(3) = r3
          call normdv2(x1)

          !x-px:
!         Correct Gaussian distribution.
          a(1,1) = xmu1 + sig1*x2(1)/sq12*sqrt(5.0)
          a(2,1) = xmu2 + sig2*(-muxpx*x2(1)/sq12+x1(1))
!         Rob's Gaussian distribution.
          !a(1,1) = xmu1 + sig1*x2(1)*sqrt(5.0)
          !a(2,1) = xmu2 + sig2*(muxpx*x2(1)+sq12*x1(1))
          !y-py
!         Correct Gaussian distribution.
          a(3,1) = xmu3 + sig3*x2(2)/sq34*sqrt(5.0)
          a(4,1) = xmu4 + sig4*(-muypy*x2(2)/sq34+x1(2))
!         Rob's Gaussian distribution.
          !a(3,1) = xmu3 + sig3*x2(2)*sqrt(5.0)
          !a(4,1) = xmu4 + sig4*(muypy*x2(2)+sq34*x1(2))
          !z-pz
!         Correct Gaussian distribution.
          a(5,1) = xmu5 + sig5*x2(3)/sq56*sqrt(5.0)
          a(6,1) = xmu6 + sig6*(-muzpz*x2(3)/sq56+x1(3))
!         Rob's Gaussian distribution.
          !a(5,1) = xmu5 + sig5*x2(3)*sqrt(5.0)
          !a(6,1) = xmu6 + sig6*(muzpz*x2(3)+sq56*x1(3))

          do jj = 1, 1

          if(totnp.ne.1) then

          if((myidx.eq.0).and.(myidy.eq.0)) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).le.lcrange(4))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).gt.lcrange(5)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
            if((a(5,jj).le.lcrange(6)).and.(a(3,jj).gt.lcrange(3))) &
            then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidx.eq.0) then
            if((a(5,jj).le.lcrange(6)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidx.eq.(npx-1)) then
            if((a(5,jj).gt.lcrange(5)).and.((a(3,jj).gt.lcrange(3))&
               .and.(a(3,jj).le.lcrange(4))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidy.eq.0) then
            if((a(3,jj).le.lcrange(4)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          elseif(myidy.eq.(npy-1)) then
            if((a(3,jj).gt.lcrange(3)).and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          else
            if( ((a(3,jj).gt.lcrange(3)).and.(a(3,jj).le.lcrange(4))) &
               .and.((a(5,jj).gt.lcrange(5))&
               .and.(a(5,jj).le.lcrange(6))) )  then
              numpts0 = numpts0 + 1
              do j = 1, 6
                this%Pts1(j,numpts0) = a(j,jj)
              enddo
              this%Pts1(7,numpts0) = this%Charge/this%mass
              this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
              this%Pts1(9,numpts0) = ii

              if(mod(numpts0,avgpts).eq.0) then
                allocate(ptstmp(9,numpts0))
                do i0 = 1, numpts0
                  do j = 1, 9
                    ptstmp(j,i0) = this%Pts1(j,i0)
                  enddo
                enddo
                deallocate(this%Pts1)
                allocate(this%Pts1(9,numpts0+avgpts))
                do i0 = 1, numpts0
                  do j = 1, 9
                    this%Pts1(j,i0) = ptstmp(j,i0)
                  enddo
                enddo
                deallocate(ptstmp)
              endif
            endif
          endif

          else
            numpts0 = numpts0 + 1
            do j = 1, 6
              this%Pts1(j,numpts0) = a(j,jj)
            enddo
            this%Pts1(7,numpts0) = this%Charge/this%mass
            this%Pts1(8,numpts0) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge) 
            this%Pts1(9,numpts0) = ii

          endif

          enddo
        enddo
          
!        call MPI_BARRIER(comm2d,ierr)
        allocate(ptstmp(9,numpts0))
        do i0 = 1, numpts0
          do j = 1, 9
            ptstmp(j,i0) = this%Pts1(j,i0)
          enddo
        enddo
        deallocate(this%Pts1)
        allocate(this%Pts1(9,numpts0))
        do i0 = 1, numpts0
          do j = 1, 9
            this%Pts1(j,i0) = ptstmp(j,i0)
          enddo
        enddo
        deallocate(ptstmp)

        this%Nptlocal = numpts0

!        call MPI_BARRIER(comm2d,ierr)
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Semigauss_Dist

        subroutine normdv2(y)
        implicit none
        include 'mpif.h'
        double precision, dimension(3), intent(out) :: y
        double precision :: sumtmp,x
        integer :: i

        sumtmp = 0.0
        do i = 1, 12
          call random_number(x)
          sumtmp = sumtmp + x
        enddo
        y(1) = sumtmp - 6.0

        sumtmp = 0.0
        do i = 1, 12
          call random_number(x)
          sumtmp = sumtmp + x
        enddo
        y(2) = sumtmp - 6.0

        sumtmp = 0.0
        do i = 1, 12
          call random_number(x)
          sumtmp = sumtmp + x
        enddo
        y(3) = sumtmp - 6.0

        end subroutine normdv2

        subroutine regenold_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')
        read(12,*)inipts,kenergy,synangle
!        allocate(Ptcl(Pdim,inipts))
        allocate(Ptcl(9,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:6,i)
          sumx = sumx + Ptcl(1,i)
          sumx2 = sumx2 + Ptcl(1,i)*Ptcl(1,i)
          sumy = sumy + Ptcl(3,i)
          sumy2 = sumy2 + Ptcl(3,i)*Ptcl(3,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        do j = 1, inipts
          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/100.0/xl)*yscale + xmu3
          gamma = 1.0+Ptcl(6,j)*1.0e6/mccq
          gammabet = sqrt(gamma*gamma-1.0)
          Ptcl(2,j) = (Ptcl(2,j)*gammabet)*pxscale + xmu2
          Ptcl(4,j) = (Ptcl(4,j)*gammabet)*pyscale + xmu4
!          betaz = sqrt((1.0-1.0/gamma/gamma)/(1+Ptcl(2,j)**2+Ptcl(4,j)**2))
!          Ptcl(2,j) = (Ptcl(2,j)*gamma*betaz)*pxscale + xmu2
!          Ptcl(4,j) = (Ptcl(4,j)*gamma*betaz)*pyscale + xmu4
!          Ptcl(5,j) = (Ptcl(5,j)-synangle)*2.0*asin(1.0)/180.0
          !unit in rad
          Ptcl(5,j) = (Ptcl(5,j)-synangle*2.0*asin(1.0)/180.0)*zscale+xmu5  
          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6 
          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
        enddo

        numpts0 = 0
        if(totnp.eq.1) then
          numpts0 = inipts
        else if(npx.eq.1) then
 
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if(myidy.eq.0) then
              if(thi.le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidy.eq.npy-1) then
              if(thi.gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (thi.gt.lcrange(3)).and.(thi.le.lcrange(4)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else if(npy.eq.1) then
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if(myidx.eq.0) then
              if(a(5).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidx.eq.npx-1) then
              if(a(5).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (a(5).gt.lcrange(5)).and.(a(5).le.lcrange(6)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if((myidx.eq.0).and.(myidy.eq.0)) then
              if((a(5).le.lcrange(6)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
              if((a(5).gt.lcrange(5)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
              if((a(5).gt.lcrange(5)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
              if((a(5).le.lcrange(6)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.0) then
              if((a(5).le.lcrange(6)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.(npx-1)) then
              if((a(5).gt.lcrange(5)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.0) then
              if((thi.le.lcrange(4)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.(npy-1)) then
              if((thi.gt.lcrange(3)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( ((thi.gt.lcrange(3)).and.(thi.le.lcrange(4))) &
                 .and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        endif

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
!            do k = 1, 6
!              call random_number(r)
!              r = (2.0*r-1.0)*0.001
!              ii = j+nset*(i-1)
!              this%Pts1(k,ii) = Ptcl(k,i)*(1.0+r)
!            enddo
              call random_number(r)
              !r = (2.0*r-1.0)*0.003*xmax
              !r = (2.0*r-1.0)*0.01*xmax
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.05*xmax
                r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.04*xmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(1,ii) = Ptcl(3,i)+r
              this%Pts1(1,ii) = Ptcl(1,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.04*pxmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(2,ii) = Ptcl(4,i)+r
              this%Pts1(2,ii) = Ptcl(2,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*ymax
                !r = (2.0*r-1.0)*0.04*ymax
                r = (2.0*r-1.0)*0.02*ymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(3,ii) = Ptcl(1,i)+r
              this%Pts1(3,ii) = Ptcl(3,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*pymax
                !r = (2.0*r-1.0)*0.04*pymax
                r = (2.0*r-1.0)*0.02*pymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(4,ii) = Ptcl(2,i)+r
              this%Pts1(4,ii) = Ptcl(4,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*zmax
                !r = (2.0*r-1.0)*0.04*zmax
                r = (2.0*r-1.0)*0.005*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pzmax
                !r = (2.0*r-1.0)*0.04*pzmax
                r = (2.0*r-1.0)*0.002*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(2,i)
        enddo


        deallocate(Ptcl)

        do j = 1, this%Nptlocal
          pid = j + myid*totnp
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        end subroutine regenold_Dist

        subroutine regen_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(6) :: tmptcl

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')
        read(12,*)inipts,kenergy,synangle

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          allocate(Ptcl(6,avgpts))

        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do j = 1, inipts
          read(12,*)tmptcl(1:6)
          sumx = sumx + tmptcl(1)
          sumx2 = sumx2 + tmptcl(1)*tmptcl(1)
          sumy = sumy + tmptcl(3)
          sumy2 = sumy2 + tmptcl(3)*tmptcl(3)
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:6,i) = tmptcl(1:6)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        do j = 1, avgpts
          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/100.0/xl)*yscale + xmu3
          gamma = 1.0+Ptcl(6,j)*1.0e6/mccq
          gammabet = sqrt(gamma*gamma-1.0)
          Ptcl(2,j) = (Ptcl(2,j)*gammabet)*pxscale + xmu2
          Ptcl(4,j) = (Ptcl(4,j)*gammabet)*pyscale + xmu4
!          betaz = sqrt((1.0-1.0/gamma/gamma)/(1+Ptcl(2,j)**2+Ptcl(4,j)**2))
!          Ptcl(2,j) = (Ptcl(2,j)*gamma*betaz)*pxscale + xmu2
!          Ptcl(4,j) = (Ptcl(4,j)*gamma*betaz)*pyscale + xmu4
!          Ptcl(5,j) = (Ptcl(5,j)-synangle)*2.0*asin(1.0)/180.0
          !unit in rad
          Ptcl(5,j) = (Ptcl(5,j)-synangle*2.0*asin(1.0)/180.0)*zscale+xmu5  
          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6 
          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
!            do k = 1, 6
!              call random_number(r)
!              r = (2.0*r-1.0)*0.001
!              ii = j+nset*(i-1)
!              this%Pts1(k,ii) = Ptcl(k,i)*(1.0+r)
!            enddo
              call random_number(r)
              !r = (2.0*r-1.0)*0.003*xmax
              !r = (2.0*r-1.0)*0.01*xmax
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.05*xmax
                r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.04*xmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(1,ii) = Ptcl(3,i)+r
              this%Pts1(1,ii) = Ptcl(1,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.04*pxmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(2,ii) = Ptcl(4,i)+r
              this%Pts1(2,ii) = Ptcl(2,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*ymax
                !r = (2.0*r-1.0)*0.04*ymax
                r = (2.0*r-1.0)*0.02*ymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(3,ii) = Ptcl(1,i)+r
              this%Pts1(3,ii) = Ptcl(3,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*pymax
                !r = (2.0*r-1.0)*0.04*pymax
                r = (2.0*r-1.0)*0.02*pymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
!              this%Pts1(4,ii) = Ptcl(2,i)+r
              this%Pts1(4,ii) = Ptcl(4,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*zmax
                !r = (2.0*r-1.0)*0.04*zmax
                r = (2.0*r-1.0)*0.005*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pzmax
                !r = (2.0*r-1.0)*0.04*pzmax
                r = (2.0*r-1.0)*0.002*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(2,i)
        enddo


        deallocate(Ptcl)

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        end subroutine regen_Dist

        subroutine GaussGamma_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4
        double precision :: sq12,sq34
        double precision, dimension(2) :: x1 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,inipts,j,pid
!        integer seedarray(1)
        double precision :: t0,x11
        double precision :: alphaz,alphapz,lambdaz,lambdapz,pi,synangle,&
        kenergy,mccq,tmp1,tmp2,Emin,Emax,phimin,phimax,gamma0,fe,femax,fphi,&
        fphimax,e1,phi1,r1,r2,r3,r4

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        alphaz = distparam(15)
        lambdaz = distparam(16)
        alphapz = distparam(17)
        lambdapz = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)

        open(unit=12,file='partcl.data',status='old')
        read(12,*)inipts,kenergy,synangle
        close(12)
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      
        pi = 2*asin(1.0)
        synangle = synangle*pi/180.0

        Emin = 0.0
        !Emax = 6.9
        Emax = 7.5
        phimin = -1.8
        phimax = 0.6

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
        numpts = 0
        print*,"kenergy: ",kenergy,synangle,avgpts,alphaz,lambdaz,&
               alphapz,lambdapz
        tmp1 = Emax - alphaz/lambdaz
        tmp2 = phimin + alphapz/lambdapz
        do
          ! rejection sample.
10        call random_number(r1)
          r1 = Emin + (Emax-Emin)*r1
          femax = (Emax - tmp1)**alphaz*exp(-lambdaz*(Emax-tmp1))
          fe = (Emax - r1)**alphaz*exp(-lambdaz*(Emax-r1))/femax
          call random_number(r3)
          if(r3.gt.fe) goto 10
          e1 = r1

20        call random_number(r2)
          r2 = phimin + (phimax-phimin)*r2
          fphimax = (tmp2-phimin)**alphapz*exp(-lambdapz*(tmp2-phimin))
          fphi = (r2-phimin)**alphapz*exp(-lambdapz*(r2-phimin))/fphimax
          call random_number(r4)
          if(r4.gt.fphi) goto 20
          phi1 = r2

          numpts = numpts + 1
          if(numpts.gt.avgpts) exit

          !x-px:
          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(1,numpts) = xmu1 + sig1*x1(1)/sq12
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1(1)/sq12+x1(2))
          !y-py
          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(3,numpts) = xmu3 + sig3*x1(1)/sq34
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x1(1)/sq34+x1(2))
          !z-pz
          this%Pts1(5,numpts) = xmu5 + (phi1-synangle)
          this%Pts1(6,numpts) = xmu6 + gamma0 - (1+e1*1.0e6/mccq)
        enddo  

        this%Nptlocal = avgpts
        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine GaussGamma_Dist

        subroutine WaterGamma_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4
        double precision :: sq12,sq34
        double precision, dimension(2) :: x1 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,inipts,j,pid
!        integer seedarray(1)
        double precision :: t0,x11
        double precision :: alphaz,alphapz,lambdaz,lambdapz,pi,synangle,&
        kenergy,mccq,tmp1,tmp2,Emin,Emax,phimin,phimax,gamma0,fe,femax,fphi,&
        fphimax,e1,phi1,r1,r2,r3,r4,r5,r6,xx1,xx2,xx3,xx4

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        alphaz = distparam(15)
        lambdaz = distparam(16)
        alphapz = distparam(17)
        lambdapz = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)
        tmp1 = Emax - alphaz/lambdaz
        tmp2 = phimin + alphapz/lambdapz
        do
          ! rejection sample.
10        call random_number(r1)
          r1 = Emin + (Emax-Emin)*r1
          femax = (Emax - tmp1)**alphaz*exp(-lambdaz*(Emax-tmp1))
          fe = (Emax - r1)**alphaz*exp(-lambdaz*(Emax-r1))/femax
          call random_number(r3)
          if(r3.gt.fe) goto 10
          e1 = r1

20        call random_number(r2)
          r2 = phimin + (phimax-phimin)*r2
          fphimax = (tmp2-phimin)**alphapz*exp(-lambdapz*(tmp2-phimin))
          fphi = (r2-phimin)**alphapz*exp(-lambdapz*(r2-phimin))/fphimax
          call random_number(r4)
          if(r4.gt.fphi) goto 20
          phi1 = r2

30        call random_number(r1)
          call random_number(r2)
          call random_number(r3)
          call random_number(r4)
          call random_number(r5)
          call random_number(r6)
          r1 = 2.0*r1-1.0
          r2 = 2.0*r2-1.0
          r3 = 2.0*r3-1.0
          r4 = 2.0*r4-1.0
          r5 = 2.0*r5-1.0
          r6 = 2.0*r6-1.0
          if(r1**2+r2**2+r3**2+r4**2+r5**2+r6**2.gt.1.0) goto 30
!x-px:
          xx1 = r1*sqrt(8.0)
          xx2 = r2*sqrt(8.0)
          xx3 = r3*sqrt(8.0)
          xx4 = r4*sqrt(8.0)

          numpts = numpts + 1
          if(numpts.gt.avgpts) exit

          !x-px:
!          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(1,numpts) = xmu1 + sig1*xx1/sq12
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*xx1/sq12+xx2)
          !y-py
!          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(3,numpts) = xmu3 + sig3*xx3/sq34
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*xx3/sq34+xx4)
          !z-pz
          this%Pts1(5,numpts) = xmu5 + (phi1-synangle)
          this%Pts1(6,numpts) = xmu6 + gamma0 - (1+e1*1.0e6/mccq)
        enddo  

        this%Nptlocal = avgpts
        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine WaterGamma_Dist

        subroutine KVGamma_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4
        double precision :: sq12,sq34
        double precision, dimension(2) :: x1 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,inipts,j,pid
!        integer seedarray(1)
        double precision :: t0,x11
        double precision :: alphaz,alphapz,lambdaz,lambdapz,pi,synangle,&
        kenergy,mccq,tmp1,tmp2,Emin,Emax,phimin,phimax,gamma0,fe,femax,fphi,&
        fphimax,e1,phi1,r1,r2,r3,r4,r5,xx1,xx2,xx3,xx4,twopi

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        alphaz = distparam(15)
        lambdaz = distparam(16)
        alphapz = distparam(17)
        lambdapz = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)

        open(unit=12,file='partcl.data',status='old')
        read(12,*)inipts,kenergy,synangle
        close(12)
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      
        pi = 2*asin(1.0)
        synangle = synangle*pi/180.0

        Emin = 0.0
        !Emax = 6.9
        Emax = 7.5
        phimin = -1.8
        phimax = 0.6

        avgpts = this%Npt/(npx*npy)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)

        ! initial allocate 'avgpts' particles on each processor.
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
        numpts = 0
        tmp1 = Emax - alphaz/lambdaz
        tmp2 = phimin + alphapz/lambdapz
        twopi = 4*asin(1.0)
        do
          ! rejection sample.
10        call random_number(r1)
          r1 = Emin + (Emax-Emin)*r1
          femax = (Emax - tmp1)**alphaz*exp(-lambdaz*(Emax-tmp1))
          fe = (Emax - r1)**alphaz*exp(-lambdaz*(Emax-r1))/femax
          call random_number(r3)
          if(r3.gt.fe) goto 10
          e1 = r1

20        call random_number(r2)
          r2 = phimin + (phimax-phimin)*r2
          fphimax = (tmp2-phimin)**alphapz*exp(-lambdapz*(tmp2-phimin))
          fphi = (r2-phimin)**alphapz*exp(-lambdapz*(r2-phimin))/fphimax
          call random_number(r4)
          if(r4.gt.fphi) goto 20
          phi1 = r2

          numpts = numpts + 1
          if(numpts.gt.avgpts) exit

          call random_number(r1)
          call random_number(r2)
          call random_number(r3)
          r4 = sqrt(r1)
          r5 = sqrt(1.0-r1)
          r2 = r2*twopi
          r3 = r3*twopi
          xx1 = 2*r4*cos(r2)
          xx2 = 2*r4*sin(r2)
          xx3 = 2*r5*cos(r3)
          xx4 = 2*r5*sin(r3)

          !x-px:
!          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(1,numpts) = xmu1 + sig1*xx1/sq12
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*xx1/sq12+xx2)
          !y-py
!          call normdv(x1)
!         Correct Gaussian distribution.
          this%Pts1(3,numpts) = xmu3 + sig3*xx3/sq34
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*xx3/sq34+xx4)
          !z-pz
          this%Pts1(5,numpts) = xmu5 + (phi1-synangle)
          this%Pts1(6,numpts) = xmu6 + gamma0 - (1+e1*1.0e6/mccq)
        enddo  

        this%Nptlocal = avgpts
        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine KVGamma_Dist

        subroutine regen2_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')
        read(12,*)inipts,kenergy,synangle
!        allocate(Ptcl(Pdim,inipts))
        allocate(Ptcl(9,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:6,i)
          sumx = sumx + Ptcl(1,i)
          sumx2 = sumx2 + Ptcl(1,i)*Ptcl(1,i)
          sumy = sumy + Ptcl(3,i)
          sumy2 = sumy2 + Ptcl(3,i)*Ptcl(3,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        do j = 1, inipts
          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/100.0/xl)*yscale + xmu3
          gamma = 1.0+Ptcl(6,j)*1.0e6/mccq
          gammabet = sqrt(gamma*gamma-1.0)
          Ptcl(2,j) = (Ptcl(2,j)*gammabet)*pxscale + xmu2
          Ptcl(4,j) = (Ptcl(4,j)*gammabet)*pyscale + xmu4
!          betaz = sqrt((1.0-1.0/gamma/gamma)/(1+Ptcl(2,j)**2+Ptcl(4,j)**2))
!          Ptcl(2,j) = (Ptcl(2,j)*gamma*betaz)*pxscale + xmu2
!          Ptcl(4,j) = (Ptcl(4,j)*gamma*betaz)*pyscale + xmu4
!          Ptcl(5,j) = (Ptcl(5,j)-synangle)*2.0*asin(1.0)/180.0
          !unit in rad
          Ptcl(5,j) = (Ptcl(5,j)-synangle*2.0*asin(1.0)/180.0)*zscale+xmu5  
          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6 
          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
        enddo

        numpts0 = 0
        if(totnp.eq.1) then
          numpts0 = inipts
        else if(npx.eq.1) then
 
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if(myidy.eq.0) then
              if(thi.le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidy.eq.npy-1) then
              if(thi.gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (thi.gt.lcrange(3)).and.(thi.le.lcrange(4)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else if(npy.eq.1) then
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if(myidx.eq.0) then
              if(a(5).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidx.eq.npx-1) then
              if(a(5).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (a(5).gt.lcrange(5)).and.(a(5).le.lcrange(6)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if((myidx.eq.0).and.(myidy.eq.0)) then
              if((a(5).le.lcrange(6)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
              if((a(5).gt.lcrange(5)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
              if((a(5).gt.lcrange(5)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
              if((a(5).le.lcrange(6)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.0) then
              if((a(5).le.lcrange(6)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.(npx-1)) then
              if((a(5).gt.lcrange(5)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.0) then
              if((thi.le.lcrange(4)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.(npy-1)) then
              if((thi.gt.lcrange(3)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( ((thi.gt.lcrange(3)).and.(thi.le.lcrange(4))) &
                 .and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        endif

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
!            do k = 1, 6
!              call random_number(r)
!              r = (2.0*r-1.0)*0.001
!              ii = j+nset*(i-1)
!              this%Pts1(k,ii) = Ptcl(k,i)*(1.0+r)
!            enddo
              call random_number(r)
              !r = (2.0*r-1.0)*0.003*xmax
              !r = (2.0*r-1.0)*0.01*xmax
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.05*xmax
                !r = (2.0*r-1.0)*0.015*xmax !old one
                r = (2.0*r-1.0)*0.04*xmax !tt59
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(1,ii) = Ptcl(3,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax !old one
                r = (2.0*r-1.0)*0.04*pxmax !tt59
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(2,ii) = Ptcl(4,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*ymax
                r = (2.0*r-1.0)*0.04*ymax !tt59
                !r = (2.0*r-1.0)*0.02*ymax !old one
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(3,ii) = Ptcl(1,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*pymax
                r = (2.0*r-1.0)*0.04*pymax !tt59
                !r = (2.0*r-1.0)*0.02*pymax !old one
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(4,ii) = Ptcl(2,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*zmax
                r = (2.0*r-1.0)*0.04*zmax !tt59
                !r = (2.0*r-1.0)*0.005*zmax !old one
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pzmax
                r = (2.0*r-1.0)*0.04*pzmax !tt59
                !r = (2.0*r-1.0)*0.002*pzmax !old one
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo

          do j = 1,1
              ii = j+nset*(i-1)
              !tt35
              !this%Pts1(1,ii) = this%Pts1(2,ii)*1.2
              !this%Pts1(2,ii) = this%Pts1(2,ii)*2.0
              !this%Pts1(3,ii) = this%Pts1(2,ii)*1.2
              !this%Pts1(4,ii) = this%Pts1(4,ii)*2.0
              !tt37
              !this%Pts1(1,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(2,ii)*0.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt38
              !this%Pts1(1,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(3,ii) = this%Pts1(2,ii)*0.0
              !this%Pts1(3,ii) = 0.0
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.5
              !tt39
              !this%Pts1(1,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt32,tt40, nominal
              !tt41
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.0
              !this%Pts1(2,ii) = this%Pts1(2,ii)*2.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.0
              !this%Pts1(4,ii) = this%Pts1(4,ii)*2.0
              !tt42
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt43
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*2.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5
              !this%Pts1(4,ii) = this%Pts1(4,ii)*2.0
              !tt44 nominal case with centroid offset for off-energy particles
              !this%Pts1(3,ii) = this%Pts1(3,ii) + 0.00025/xl
              !tt45
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5 + 0.001/xl
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt46
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5 + 0.001/xl
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5 + 0.001/xl
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt47
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5 + 0.001/xl 
              !this%Pts1(2,ii) = this%Pts1(2,ii)*2.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5 + 0.001/xl
              !this%Pts1(4,ii) = this%Pts1(4,ii)*2.0
              !tt48
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.1 + 0.001/xl 
              !this%Pts1(2,ii) = this%Pts1(2,ii)*2.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.1 + 0.001/xl
              !this%Pts1(4,ii) = this%Pts1(4,ii)*2.0
              !tt49
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5 + 0.00025/xl 
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5 
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt50 new partcl.data, different percentage off-energy particle
              !tt51,tt52 new partcl.data, Quad scan
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5  
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5 
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt53 new partcl.data
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5  
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5 
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.5
              call random_number(r)
              this%Pts1(5,ii) = 2*pi*r - synangle*pi/180.0
              !this%Pts1(6,ii) = gamma0 - (1+0.075e6/mccq) 
              call random_number(r)
              this%Pts1(6,ii) = gamma0 - (1+(6.3+0.3*r)*1.0e6/mccq) 
              !this%Pts1(6,ii) = gamma0 - (1+(6.3+0.4*r)*1.0e6/mccq) 
              !this%Pts1(6,ii) = gamma0 - (1+(6.0+0.7*r)*1.0e6/mccq) 
              !ttmatch using nominal output only, no off-energy particles
              !tt53 new partcl.data, no off-energy particle
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.5  
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.0
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.5 
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.0
              !tt54 new partcl.data, no off-energy particle
              !this%Pts1(1,ii) = this%Pts1(1,ii)*0.5  
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.5
              !this%Pts1(3,ii) = this%Pts1(3,ii)*0.5 
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.5
              !tt55 new partcl.data, no off-energy particle,conserve RMS
              !this%Pts1(1,ii) = this%Pts1(1,ii)*1.17  
              !this%Pts1(2,ii) = this%Pts1(2,ii)*1.17
              !this%Pts1(3,ii) = this%Pts1(3,ii)*1.17 
              !this%Pts1(4,ii) = this%Pts1(4,ii)*1.17
              !tt56 new partcl.data, no off-energy particle,conserve RMS
              !this%Pts1(1,ii) = this%Pts1(1,ii)*0.56  
              !this%Pts1(2,ii) = this%Pts1(2,ii)*0.56
              !this%Pts1(3,ii) = this%Pts1(3,ii)*0.56 
              !this%Pts1(4,ii) = this%Pts1(4,ii)*0.56
              !tt57 test steering magnets, the other same as tt56
              !tt58 test steering magnets, no off-energy pt, matched beam.
              !tt59, 1% and 5% off-energy, big mismatch,new repopulation radii
!              this%Pts1(1,ii) = this%Pts1(1,ii)*2.00  
!              this%Pts1(2,ii) = this%Pts1(2,ii)*2.00
!              this%Pts1(3,ii) = this%Pts1(3,ii)*2.00 
!              this%Pts1(4,ii) = this%Pts1(4,ii)*2.00
              !tt60, 5% off-energy, big mismatch,new repopulation radii:400kev
          enddo

          sumx = sumx + this%Pts1(2,i)
        enddo


        deallocate(Ptcl)

        do j = 1, this%Nptlocal
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = j
        enddo

        end subroutine regen2_Dist

        subroutine Gauss7_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(3) :: al0,ga0,epson0
        double precision :: t0

        call starttime_Timer(t0)

        call Waterbag_Dist(this,nparam,distparam,grid,0)
        call gammaepson_Dist(this,al0,ga0,epson0)
        call Distort_Dist(this,al0,ga0,epson0,grid)

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss7_Dist

        subroutine Distort_Dist(this,al0,ga0,epson0,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(3) :: al0,ga0,epson0
        integer :: avgpts
        integer :: i,ierr,j,pid
        integer*8 :: nptot
        double precision :: t0
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy,vtmp
        double precision :: r11x,r22x,r12x
        double precision :: r11y,r22y,r12y
        double precision :: r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1,xx
        double precision, dimension(3,3) :: rr
        double precision, dimension(3) :: frac,CSbeta,CSalpha,CSgamma,CSepson
        double precision, dimension(6) :: ptctmp
        double precision, dimension(15) :: tmplc,tmpgl
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,yrms
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0
        double precision ::sqx,sqpx,sqy,sqpy,sqz,sqpz
        double precision:: qmc,xl,xt,pi
        integer :: totnp,npy,npx,myidy,myid,myidx

        call starttime_Timer(t0)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg
        nptot = this%Npt
        avgpts = this%Nptlocal

        ddx = -4.0e2 !for alphax > 0.0
        !ddx = -8.0e2 !for alphax > 0.0
        ddx = -4.0e2 !for alphax > 0.0
        !ddy = -3.0e2 
        ddy = 1.0e7  !for alphay < 0.0
        ddy = 3.0e7  !for alphay < 0.0
        ddy = 2.0e7  !for alphay < 0.0 for vis
        !ddy = -4.0e2 !for alphay > 0.0 for visulization 
        !ddy = 6.0e7  !for alphay < 0.0 for the purpose of display
        !ddz = -5.0e-3
        !ddz = 5.0e11 !for alphaz < 0.0
        ddz = 0.0
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0

        den1 = 1.0/float(nptot)
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
        pi = 2*asin(1.0)

        do i = 1, avgpts
          ptctmp(1) = this%Pts1(1,i)
          vtmp = this%Pts1(1,i) 
          ptctmp(2) = this%Pts1(2,i) + ddx*vtmp**3

          !for alpha > 0
          !ptctmp(3) = this%Pts1(3,i)
          !vtmp = this%Pts1(3,i) 
          !ptctmp(4) = this%Pts1(4,i) + ddy*vtmp**3
          !for alpha < 0
          vtmp = this%Pts1(4,i) 
          ptctmp(3) = this%Pts1(3,i) + ddy*vtmp**3
          ptctmp(4) = this%Pts1(4,i)
         
          !ptctmp(5) = this%Pts1(5,i)
          !vtmp = this%Pts1(5,i) 
          !ptctmp(6) = this%Pts1(6,i) + ddz*vtmp**3
          vtmp = this%Pts1(6,i)
          ptctmp(5) = this%Pts1(5,i) + ddz*vtmp**3
          ptctmp(6) = this%Pts1(6,i)

          x0lc = x0lc + ptctmp(1)
          sqsum1local = sqsum1local + ptctmp(1)*ptctmp(1)
          xpxlocal = xpxlocal + ptctmp(1)*ptctmp(2)
          px0lc = px0lc + ptctmp(2)
          sqsum2local = sqsum2local + ptctmp(2)*ptctmp(2)
          y0lc = y0lc + ptctmp(3)
          sqsum3local = sqsum3local + ptctmp(3)*ptctmp(3)
          ypylocal = ypylocal + ptctmp(3)*ptctmp(4)
          py0lc = py0lc + ptctmp(4)
          sqsum4local = sqsum4local + ptctmp(4)*ptctmp(4)
          z0lc = z0lc + ptctmp(5)
          sqsum5local = sqsum5local + ptctmp(5)*ptctmp(5)
          zpzlocal = zpzlocal + ptctmp(5)*ptctmp(6)
          pz0lc = pz0lc + ptctmp(6)
          sqsum6local = sqsum6local + ptctmp(6)*ptctmp(6)
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

        epsx2 = (sqsum1*sqsum2-xpx*xpx)
        epsy2 = (sqsum3*sqsum4-ypy*ypy)
        epsz2 = (sqsum5*sqsum6-zpz*zpz)
        epx = sqrt(max(epsx2,0.0d0))
        epy = sqrt(max(epsy2,0.0d0))
        epz = sqrt(max(epsz2,0.0d0))
        xrms = sqrt(sqsum1)
        yrms = sqrt(sqsum3)
        zrms = sqrt(sqsum5)

        CSalpha(1) = -xpx/epx
        CSalpha(2) = -ypy/epy
        CSalpha(3) = -zpz/epz
        CSbeta(1) = (xrms*xl)**2/(epx*xl)
        CSbeta(2) = (yrms*xl)**2/(epy*xl)
        CSbeta(3) = (zrms*xt)**2/(epz*qmc*xt)
        CSgamma(:) = (1.0 + CSalpha(:)*CSalpha(:))/CSbeta(:)*xl
        CSgamma(3) = (1.0 + CSalpha(3)*CSalpha(3))/CSbeta(3)/qmc*xt

        CSepson(1) = epx*xl
        CSepson(2) = epy*xl
        CSepson(3) = epz*qmc*xt

        do i = 1, 3
          frac(i) = sqrt(epson0(i)/CSepson(i))
          rr(1,i) = sqrt(CSgamma(i)/ga0(i))
          rr(2,i) = 1.0/rr(1,i)
          rr(3,i) = (CSalpha(i)-al0(i))/sqrt(CSgamma(i)*ga0(i))
        enddo

        do i = 1, avgpts
          vtmp = this%Pts1(1,i) 
          this%Pts1(2,i) = this%Pts1(2,i) + ddx*vtmp**3
          !for alpha > 0
          !vtmp = this%Pts1(3,i)
          !this%Pts1(4,i) = this%Pts1(4,i) + ddy*vtmp**3
          !for alpha < 0
          vtmp = this%Pts1(4,i)
          this%Pts1(3,i) = this%Pts1(3,i) + ddy*vtmp**3
          !vtmp = this%Pts1(5,i)
          !this%Pts1(6,i) = this%Pts1(6,i) + ddz*vtmp**3
          vtmp = this%Pts1(6,i)
          this%Pts1(5,i) = this%Pts1(5,i) + ddz*vtmp**3
        enddo

        factx = frac(1)
        facty = frac(2)
        factz = frac(3)

        do i = 1, avgpts
          this%Pts1(1,i) = this%Pts1(1,i)*factx
          this%Pts1(2,i) = this%Pts1(2,i)*factx
          this%Pts1(3,i) = this%Pts1(3,i)*facty
          this%Pts1(4,i) = this%Pts1(4,i)*facty
          this%Pts1(5,i) = this%Pts1(5,i)*factz
          this%Pts1(6,i) = this%Pts1(6,i)*factz
        enddo
        
        r11x = rr(1,1)
        r22x = rr(2,1)
        r12x = rr(3,1) 
        r11y = rr(1,2)
        r22y = rr(2,2)
        r12y = rr(3,2) 
        r11z = rr(1,3)
        r22z = rr(2,3)
        r12z = rr(3,3)

        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xxx
          this%Pts1(2,i) = yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xxx
          this%Pts1(4,i) = yyy
          xxx = this%Pts1(5,i)*r11z + this%Pts1(6,i)*r12z
          yyy = this%Pts1(6,i)*r22z
          this%Pts1(5,i) = xxx
          this%Pts1(6,i) = yyy
        enddo

        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo
       
        end subroutine Distort_Dist

        subroutine Regen7_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        double precision, dimension(3) :: al0,ga0,epson0
        double precision :: t0

        call starttime_Timer(t0)

        call Waterbag_Dist(this,nparam,distparam,grid,0)
        call gammaepson_Dist(this,al0,ga0,epson0)
        call regendstort_Dist(this,nparam,distparam,geom,grid,Flagbc)
        call DistReg_Dist(this,al0,ga0,epson0,grid)

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Regen7_Dist

        subroutine DistReg_Dist(this,al0,ga0,epson0,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        double precision, dimension(3) :: al0,ga0,epson0
        type (Pgrid2d), intent(in) :: grid
        integer :: avgpts
        integer :: i,ierr,j,pid
        integer*8 :: nptot
        double precision :: t0
        double precision :: ddx,ddy,ddz,factx,facty,factz,xxx,yyy,vtmp
        double precision :: r11x,r22x,r12x
        double precision :: r11y,r22y,r12y
        double precision :: r11z,r22z,r12z
        double precision :: ddx1,ddy1,ddz1,xx
        double precision, dimension(3,3) :: rr
        double precision, dimension(3) :: frac,CSbeta,CSalpha,CSgamma,CSepson
        double precision, dimension(6) :: ptctmp
        double precision, dimension(15) :: tmplc,tmpgl
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,yrms
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0
        double precision ::sqx,sqpx,sqy,sqpy,sqz,sqpz
        double precision:: qmc,xl,xt,pi
        integer :: totnp,npx,npy,myid,myidx,myidy

        call starttime_Timer(t0)

        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        qmc = this%Mass/1.0e6
        xl = Scxl
        xt = Rad2deg
        nptot = this%Npt
        avgpts = this%Nptlocal

        ddx = 0.0 
        ddy = 0.0 
        ddz = 0.0
        ddx1 = 0.0
        ddy1 = 0.0
        ddz1 = 0.0

        den1 = 1.0/float(nptot)
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
        pi = 2*asin(1.0)

        do i = 1, avgpts
          ptctmp(1) = this%Pts1(1,i)
          vtmp = this%Pts1(1,i) 
          ptctmp(2) = this%Pts1(2,i) + ddx*vtmp**3

          !for alpha > 0
          !ptctmp(3) = this%Pts1(3,i)
          !vtmp = this%Pts1(3,i) 
          !ptctmp(4) = this%Pts1(4,i) + ddy*vtmp**3
          !for alpha < 0
          vtmp = this%Pts1(4,i) 
          ptctmp(3) = this%Pts1(3,i) + ddy*vtmp**3
          ptctmp(4) = this%Pts1(4,i)
         
          !ptctmp(5) = this%Pts1(5,i)
          !vtmp = this%Pts1(5,i) 
          !ptctmp(6) = this%Pts1(6,i) + ddz*vtmp**3
          vtmp = this%Pts1(6,i)
          ptctmp(5) = this%Pts1(5,i) + ddz*vtmp**3
          ptctmp(6) = this%Pts1(6,i)

          x0lc = x0lc + ptctmp(1)
          sqsum1local = sqsum1local + ptctmp(1)*ptctmp(1)
          xpxlocal = xpxlocal + ptctmp(1)*ptctmp(2)
          px0lc = px0lc + ptctmp(2)
          sqsum2local = sqsum2local + ptctmp(2)*ptctmp(2)
          y0lc = y0lc + ptctmp(3)
          sqsum3local = sqsum3local + ptctmp(3)*ptctmp(3)
          ypylocal = ypylocal + ptctmp(3)*ptctmp(4)
          py0lc = py0lc + ptctmp(4)
          sqsum4local = sqsum4local + ptctmp(4)*ptctmp(4)
          z0lc = z0lc + ptctmp(5)
          sqsum5local = sqsum5local + ptctmp(5)*ptctmp(5)
          zpzlocal = zpzlocal + ptctmp(5)*ptctmp(6)
          pz0lc = pz0lc + ptctmp(6)
          sqsum6local = sqsum6local + ptctmp(6)*ptctmp(6)
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

        epsx2 = (sqsum1*sqsum2-xpx*xpx)
        epsy2 = (sqsum3*sqsum4-ypy*ypy)
        epsz2 = (sqsum5*sqsum6-zpz*zpz)
        epx = sqrt(max(epsx2,0.0d0))
        epy = sqrt(max(epsy2,0.0d0))
        epz = sqrt(max(epsz2,0.0d0))
        xrms = sqrt(sqsum1)
        yrms = sqrt(sqsum3)
        zrms = sqrt(sqsum5)

        CSalpha(1) = -xpx/epx
        CSalpha(2) = -ypy/epy
        CSalpha(3) = -zpz/epz
        CSbeta(1) = (xrms*xl)**2/(epx*xl)
        CSbeta(2) = (yrms*xl)**2/(epy*xl)
        CSbeta(3) = (zrms*xt)**2/(epz*qmc*xt)
        CSgamma(:) = (1.0 + CSalpha(:)*CSalpha(:))/CSbeta(:)*xl
        CSgamma(3) = (1.0 + CSalpha(3)*CSalpha(3))/CSbeta(3)/qmc*xt

        CSepson(1) = epx*xl
        CSepson(2) = epy*xl
        CSepson(3) = epz*qmc*xt

        do i = 1, 3
          frac(i) = sqrt(epson0(i)/CSepson(i))
          rr(1,i) = sqrt(CSgamma(i)/ga0(i))
          rr(2,i) = 1.0/rr(1,i)
          rr(3,i) = (CSalpha(i)-al0(i))/sqrt(CSgamma(i)*ga0(i))
        enddo

        do i = 1, avgpts
          vtmp = this%Pts1(1,i) 
          this%Pts1(2,i) = this%Pts1(2,i) + ddx*vtmp**3
          !for alpha > 0
          !vtmp = this%Pts1(3,i)
          !this%Pts1(4,i) = this%Pts1(4,i) + ddy*vtmp**3
          !for alpha < 0
          vtmp = this%Pts1(4,i)
          this%Pts1(3,i) = this%Pts1(3,i) + ddy*vtmp**3
          !vtmp = this%Pts1(5,i)
          !this%Pts1(6,i) = this%Pts1(6,i) + ddz*vtmp**3
          vtmp = this%Pts1(6,i)
          this%Pts1(5,i) = this%Pts1(5,i) + ddz*vtmp**3
        enddo

        factx = frac(1)
        facty = frac(2)
        factz = frac(3)

        do i = 1, avgpts
          this%Pts1(1,i) = this%Pts1(1,i)*factx
          this%Pts1(2,i) = this%Pts1(2,i)*factx
          this%Pts1(3,i) = this%Pts1(3,i)*facty
          this%Pts1(4,i) = this%Pts1(4,i)*facty
          this%Pts1(5,i) = this%Pts1(5,i)*factz
          this%Pts1(6,i) = this%Pts1(6,i)*factz
        enddo
        
        r11x = rr(1,1)
        r22x = rr(2,1)
        r12x = rr(3,1) 
        r11y = rr(1,2)
        r22y = rr(2,2)
        r12y = rr(3,2) 
        r11z = rr(1,3)
        r22z = rr(2,3)
        r12z = rr(3,3)

        do i = 1, avgpts
          xxx = this%Pts1(1,i)*r11x + this%Pts1(2,i)*r12x
          yyy = this%Pts1(2,i)*r22x
          this%Pts1(1,i) = xxx
          this%Pts1(2,i) = yyy
          xxx = this%Pts1(3,i)*r11y + this%Pts1(4,i)*r12y
          yyy = this%Pts1(4,i)*r22y
          this%Pts1(3,i) = xxx
          this%Pts1(4,i) = yyy
          xxx = this%Pts1(5,i)*r11z + this%Pts1(6,i)*r12z
          yyy = this%Pts1(6,i)*r22z
          this%Pts1(5,i) = xxx
          this%Pts1(6,i) = yyy
        enddo

        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo
       
        end subroutine DistReg_Dist

        subroutine regendstort_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,ikeep
        double precision, dimension(6) :: lcrange,a,tmptcl
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')
        read(12,*)inipts,kenergy,synangle
!        allocate(Ptcl(Pdim,inipts))
        allocate(Ptcl(9,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        ikeep = 0
        do i = 1, inipts
          read(12,*)tmptcl(1:6)
          if(abs(tmptcl(6)-kenergy/1.0e6).lt.1.0) then
            ikeep = ikeep + 1
            Ptcl(1:6,ikeep) = tmptcl(1:6)
            sumx = sumx + Ptcl(1,ikeep)
            sumx2 = sumx2 + Ptcl(1,ikeep)*Ptcl(1,ikeep)
            sumy = sumy + Ptcl(3,ikeep)
            sumy2 = sumy2 + Ptcl(3,ikeep)*Ptcl(3,ikeep)
          else
          endif
        enddo
        inipts = ikeep
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        gamma0 = 1.0+kenergy/mccq      !2.5 MeV at beginning of DTL.
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        do j = 1, inipts
          Ptcl(1,j) = (Ptcl(1,j)/100.0/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/100.0/xl)*yscale + xmu3
          gamma = 1.0+Ptcl(6,j)*1.0e6/mccq
          gammabet = sqrt(gamma*gamma-1.0)
          Ptcl(2,j) = (Ptcl(2,j)*gammabet)*pxscale + xmu2
          Ptcl(4,j) = (Ptcl(4,j)*gammabet)*pyscale + xmu4
!          betaz = sqrt((1.0-1.0/gamma/gamma)/(1+Ptcl(2,j)**2+Ptcl(4,j)**2))
!          Ptcl(2,j) = (Ptcl(2,j)*gamma*betaz)*pxscale + xmu2
!          Ptcl(4,j) = (Ptcl(4,j)*gamma*betaz)*pyscale + xmu4
!          Ptcl(5,j) = (Ptcl(5,j)-synangle)*2.0*asin(1.0)/180.0
          !unit in rad
          Ptcl(5,j) = (Ptcl(5,j)-synangle*2.0*asin(1.0)/180.0)*zscale+xmu5  
          Ptcl(6,j) = (-gamma + gamma0)*pzscale + xmu6 
          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
        enddo

        numpts0 = 0
        if(totnp.eq.1) then
          numpts0 = inipts
        else if(npx.eq.1) then
 
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if(myidy.eq.0) then
              if(thi.le.lcrange(4)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidy.eq.npy-1) then
              if(thi.gt.lcrange(3)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (thi.gt.lcrange(3)).and.(thi.le.lcrange(4)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else if(npy.eq.1) then
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if(myidx.eq.0) then
              if(a(5).le.lcrange(6)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else if(myidx.eq.npx-1) then
              if(a(5).gt.lcrange(5)) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( (a(5).gt.lcrange(5)).and.(a(5).le.lcrange(6)) ) then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        else
          do i = 1, inipts
            a(1:6) = Ptcl(1:6,i)
            if((Flagbc.eq.3).or.(Flagbc.eq.4)) then
              ri = sqrt(a(1)*a(1)+a(3)*a(3))
              if(a(1).gt.0.0) then
                if(a(3).gt.0.0) then
                  thi = asin(a(3)/ri)
                else
                  thi = 2*pi+asin(a(3)/ri)
                endif
              else
                thi = pi - asin(a(3)/ri)
              endif
            else
              thi = a(3)
            endif
            if((myidx.eq.0).and.(myidy.eq.0)) then
              if((a(5).le.lcrange(6)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.0)) then
              if((a(5).gt.lcrange(5)).and.(thi.le.lcrange(4))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.(npx-1)).and.(myidy.eq.(npy-1))) then
              if((a(5).gt.lcrange(5)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif((myidx.eq.0).and.(myidy.eq.(npy-1))) then
              if((a(5).le.lcrange(6)).and.(thi.gt.lcrange(3))) &
              then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.0) then
              if((a(5).le.lcrange(6)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidx.eq.(npx-1)) then
              if((a(5).gt.lcrange(5)).and.((thi.gt.lcrange(3))&
                 .and.(thi.le.lcrange(4))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.0) then
              if((thi.le.lcrange(4)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            elseif(myidy.eq.(npy-1)) then
              if((thi.gt.lcrange(3)).and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            else
              if( ((thi.gt.lcrange(3)).and.(thi.le.lcrange(4))) &
                 .and.((a(5).gt.lcrange(5))&
                 .and.(a(5).le.lcrange(6))) )  then
                numpts0 = numpts0 + 1
                do j = 1, 6
                  Ptcl(j,numpts0) = a(j)
                enddo
              endif
            endif
          enddo
        endif

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        deallocate(this%Pts1)
        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
!            do k = 1, 6
!              call random_number(r)
!              r = (2.0*r-1.0)*0.001
!              ii = j+nset*(i-1)
!              this%Pts1(k,ii) = Ptcl(k,i)*(1.0+r)
!            enddo
              call random_number(r)
              !r = (2.0*r-1.0)*0.003*xmax
              !r = (2.0*r-1.0)*0.01*xmax
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.05*xmax
                r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.04*xmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(1,ii) = Ptcl(3,i)+r
!              this%Pts1(1,ii) = Ptcl(1,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.04*pxmax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(2,ii) = Ptcl(4,i)+r
!              this%Pts1(2,ii) = Ptcl(2,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*ymax
                !r = (2.0*r-1.0)*0.04*ymax
                r = (2.0*r-1.0)*0.02*ymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(3,ii) = Ptcl(1,i)+r
!              this%Pts1(3,ii) = Ptcl(3,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.01*pymax
                !r = (2.0*r-1.0)*0.04*pymax
                r = (2.0*r-1.0)*0.02*pymax
              endif
              ii = j+nset*(i-1)
! for LEDA only
              this%Pts1(4,ii) = Ptcl(2,i)+r
!              this%Pts1(4,ii) = Ptcl(4,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*zmax
                !r = (2.0*r-1.0)*0.04*zmax
                r = (2.0*r-1.0)*0.005*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pzmax
                !r = (2.0*r-1.0)*0.04*pzmax
                r = (2.0*r-1.0)*0.002*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(2,i)
        enddo


        deallocate(Ptcl)

        end subroutine regendstort_Dist

        subroutine GaussDouble_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp,pid
!        integer seedarray(1)
        double precision :: t0,x11,Xfrac,Yfrac
        double precision :: sig1a,sig2a,sig3a,sig4a,sig5a,sig6a
        double precision :: sig1b,sig2b,sig3b,sig4b,sig5b,sig6b
        integer :: num1,num2

        call starttime_Timer(t0)

        !Xfrac = 5.0
        !Yfrac = 0.01
        !Xfrac = 4.0
        !Yfrac = 0.0
        !Xfrac = 4.0
        !Yfrac = 0.001
        !Xfrac = 4.0
        !Yfrac = 0.0002
        !Xfrac = 4.0
        !Yfrac = 0.01
        !Xfrac = 4.0
        !Yfrac = 0.05
        !Xfrac = 4.0
        !Yfrac = 0.005
        !Xfrac = 2.0
        !Yfrac = 0.01
        Xfrac = 2.0
        Yfrac = 0.05

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)

        avgpts = this%Npt/(npx*npy)

        if(mod(avgpts,10).ne.0) then
          print*,"The number of particles has to be an integer multiple of 10Nprocs"
          stop
        endif
        
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sig1a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig1
        sig2a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig2
        sig3a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig3
        sig4a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig4
        sig5a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig5
        sig6a = sqrt((1+Yfrac)/(1+Yfrac*Xfrac*Xfrac))*sig6

        sig1b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig1
        sig2b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig2
        sig3b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig3
        sig4b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig4
        sig5b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig5
        sig6b = sqrt((Xfrac*Xfrac+Xfrac*Xfrac*Yfrac)/(1+Yfrac*Xfrac*Xfrac))&
                *sig6

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

!        allocate(x1(2,avgpts))
!        allocate(x2(2,avgpts))
!        allocate(x3(2,avgpts))
!        call normVec(x1,avgpts)
!        call normVec(x2,avgpts)
!        call normVec(x3,avgpts)
        
        intvsamp = 10000 !used in the previous simulation
!        intvsamp = avgpts
        allocate(x1(2,intvsamp))
        allocate(x2(2,intvsamp))
        allocate(x3(2,intvsamp))

!        num1 = intvsamp*(1-Yfrac) + 1
        num2 = intvsamp*(Yfrac/(1.0+Yfrac))
        num1 = intvsamp - num2


        do j = 1, avgpts/intvsamp
          call normVec(x1,intvsamp)
          call normVec(x2,intvsamp)
          call normVec(x3,intvsamp)
          do k = 1, num1
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1a*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2a*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3a*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4a*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(5,i) = xmu5 + sig5a*x3(1,k)/sq56
            this%Pts1(6,i) = xmu6 + sig6a*(-muzpz*x3(1,k)/sq56+x3(2,k))
          enddo
          do k = num1+1, intvsamp
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1b*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2b*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3b*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4b*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(5,i) = xmu5 + sig5b*x3(1,k)/sq56
            this%Pts1(6,i) = xmu6 + sig6b*(-muzpz*x3(1,k)/sq56+x3(2,k))
          enddo
        enddo
          
        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts
        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine GaussDouble_Dist

        ! sample the particles with intial distribution 
        ! using rejection method for multi-charge state. 
        subroutine WaterbagMC_Dist(this,nparam,distparam,grid,flagalloc,nchrg,&
                                   nptlist,qmcclist,currlist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc,nchrg
        double precision, dimension(nparam) :: distparam
        double precision, dimension(nchrg) :: qmcclist,currlist
        integer*8, dimension(nchrg) :: nptlist
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision, dimension(6) :: lcrange 
        double precision, dimension(2) :: gs
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: rootx,rooty,rootz,r1,r2,x1,x2
        double precision :: r3,r4,r5,r6,x3,x4,x5,x6
        integer :: totnp,npy,npx
        integer :: avgpts,numpts,isamz,isamy
        integer :: myid,myidx,myidy,iran,intvsamp,pid,j,i,ii,kk
!        integer seedarray(2)
        double precision :: t0,x11
        double precision, allocatable, dimension(:) :: ranum6

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        seedarray(2)=(101+2*myid)*(myid+4)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray)
        call random_number(x11)

        avgpts = this%Npt/(npx*npy)
        !if(mod(avgpts,10).ne.0) then
        !  print*,"The number of particles has to be an integer multiple of 10Nprocs" 
        !  stop
        !endif
 
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        rootx=sqrt(1.-muxpx*muxpx)
        rooty=sqrt(1.-muypy*muypy)
        rootz=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif
        numpts = 0
        isamz = 0
        isamy = 0
        intvsamp = avgpts
        !intvsamp = 10
        allocate(ranum6(6*intvsamp))

        do 
          ! rejection sample.
10        continue 
          isamz = isamz + 1
          if(mod(isamz-1,intvsamp).eq.0) then
            call random_number(ranum6)
          endif
          iran = 6*mod(isamz-1,intvsamp)
          r1 = 2.0*ranum6(iran+1)-1.0
          r2 = 2.0*ranum6(iran+2)-1.0
          r3 = 2.0*ranum6(iran+3)-1.0
          r4 = 2.0*ranum6(iran+4)-1.0
          r5 = 2.0*ranum6(iran+5)-1.0
          r6 = 2.0*ranum6(iran+6)-1.0
          if(r1**2+r2**2+r3**2+r4**2+r5**2+r6**2.gt.1.0) goto 10
          isamy = isamy + 1
          numpts = numpts + 1
          if(numpts.gt.avgpts) exit
!x-px:
          x1 = r1*sqrt(8.0)
          x2 = r2*sqrt(8.0)
          !Correct transformation.
          this%Pts1(1,numpts) = xmu1 + sig1*x1/rootx
          this%Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
          !Rob's transformation
          !this%Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
          !this%Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
!y-py:
          x3 = r3*sqrt(8.0)
          x4 = r4*sqrt(8.0)
          !correct transformation
          this%Pts1(3,numpts) = xmu3 + sig3*x3/rooty
          this%Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
          !Rob's transformation
          !this%Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
          !this%Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
!t-pt:
          x5 = r5*sqrt(8.0)
          x6 = r6*sqrt(8.0)
          !correct transformation
          this%Pts1(5,numpts) = xmu5 + sig5*x5/rootz
          this%Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
          !Rob's transformation
          !this%Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
          !this%Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
        enddo

        deallocate(ranum6)
          
        this%Nptlocal = avgpts
       
        
        ii = 0
        avgpts = sum(nptlist)/totnp
        do i = 1, nchrg
          kk = nptlist(i)/totnp
          do j = 1, kk
            ii = ii + 1
            this%Pts1(7,ii) = qmcclist(i)
            this%Pts1(8,ii) = currlist(i)/Scfreq/nptlist(i)*qmcclist(i)/abs(qmcclist(i))
            pid = ii + myid*totnp
            this%Pts1(9,ii) = pid
          enddo
        enddo

!        do j = 1, avgpts
!          pid = j + myid*totnp
!          this%Pts1(7,j) = this%Charge/this%mass
!          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
!          this%Pts1(9,j) = pid
!        enddo

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine WaterbagMC_Dist

        subroutine GaussMC_Dist(this,nparam,distparam,grid,flagalloc,nchrg,&
                                   nptlist,qmcclist,currlist)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc,nchrg
        double precision, dimension(nparam) :: distparam
        double precision, dimension(nchrg) :: qmcclist,currlist
        integer*8, dimension(nchrg) :: nptlist
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp,pid,ii,kk
!        integer seedarray(1)
        double precision :: t0,x11

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)

        avgpts = this%Npt/(npx*npy)

!        if(mod(avgpts,10).ne.0) then
!          print*,"The number of particles has to be an integer multiple of 10Nprocs"
!          stop
!        endif
        
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

!        allocate(x1(2,avgpts))
!        allocate(x2(2,avgpts))
!        allocate(x3(2,avgpts))
!        call normVec(x1,avgpts)
!        call normVec(x2,avgpts)
!        call normVec(x3,avgpts)
        
!        intvsamp = 10
        intvsamp = avgpts
        allocate(x1(2,intvsamp))
        allocate(x2(2,intvsamp))
        allocate(x3(2,intvsamp))

        do j = 1, avgpts/intvsamp
          call normVec(x1,intvsamp)
          call normVec(x2,intvsamp)
          call normVec(x3,intvsamp)
          do k = 1, intvsamp
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq56
            this%Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,k)/sq56+x3(2,k))
          enddo
        enddo
          
        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts

        ii = 0
        avgpts = sum(nptlist)/totnp
        do i = 1, nchrg
          kk = nptlist(i)/totnp
          do j = 1, kk
            ii = ii + 1
            this%Pts1(7,ii) = qmcclist(i)
            this%Pts1(8,ii) = currlist(i)/Scfreq/nptlist(i)*qmcclist(i)/abs(qmcclist(i))
            pid = ii + myid*totnp
            this%Pts1(9,ii) = pid
          enddo
        enddo
 
!        do j = 1, avgpts
!          pid = j + myid*totnp
!          this%Pts1(7,j) = this%Charge/this%mass
!          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
!          this%Pts1(9,j) = pid
!        enddo
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine GaussMC_Dist

        subroutine readElegantold_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        real*8 :: beta0,sumeng
        integer :: jlow,jhigh,nleft,avgpts

        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')
        read(12,*)inipts
!        allocate(Ptcl(Pdim,inipts))
        allocate(Ptcl(7,inipts))
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        sumeng = 0.0
        do i = 1, inipts
          read(12,*)Ptcl(1:7,i)
          sumx = sumx + Ptcl(2,i)
          sumx2 = sumx2 + Ptcl(2,i)*Ptcl(2,i)
          sumy = sumy + Ptcl(4,i)
          sumy2 = sumy2 + Ptcl(4,i)*Ptcl(4,i)
          sumeng = sumeng + Ptcl(7,i)
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        sumeng = sumeng/inipts
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+sumeng*1.0e6/mccq  
        gamma0 = -this%refptcl(6) 
        beta0 = sqrt(1.0-1./gamma0**2)
        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0

        avgpts = inipts/totnp
        nleft = inipts - avgpts*totnp
        if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
        else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
        endif

        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0

        do j = 1, inipts
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
            i = j - jlow + 1
            this%Pts1(1,i) = Ptcl(2,j)*xscale 
            this%Pts1(3,i) = Ptcl(4,j)*yscale 
            this%Pts1(5,i) = Ptcl(6,j)*zscale 
            gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
            this%Pts1(2,i) = Ptcl(3,j)*gammabet*pxscale + xmu2
            this%Pts1(4,i) = Ptcl(5,j)*gammabet*pyscale + xmu4
            !this%Pts1(6,i) = gammabet*pzscale + xmu6
            this%Pts1(6,i) = gamma0 - sqrt(1.0 + Ptcl(7,j)**2) 
            this%Pts1(9,i) = j
          endif
        enddo

        this%Nptlocal = avgpts

        do i = 1, this%Nptlocal
          this%Pts1(1,i) = this%Pts1(1,i)/Scxl + xmu1
          !this%Pts1(2,i) = 0.0
          this%Pts1(3,i) = this%Pts1(3,i)/Scxl + xmu3
          !this%Pts1(4,i) = 0.0
          this%Pts1(5,i) = this%Pts1(5,i)*Scfreq*2*Pi + xmu5
          !!this%Pts1(6,i) = xmu6
        enddo

        deallocate(Ptcl)

        do j = 1, this%Nptlocal
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
        enddo

        end subroutine readElegantold_Dist

        subroutine readElegantold2_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.abs(tmptcl(6))) zmax = abs(tmptcl(6)) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                r = (2.0*r-1.0)*0.125*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                r = (2.0*r-1.0)*0.125*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                r = (2.0*r-1.0)*0.02*zmax !used in the 1st 1B simulation
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(6,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(7,i)+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo


        deallocate(Ptcl)


        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        this%Npt = nset*inipts

        end subroutine readElegantold2_Dist

        subroutine readElegant_Dist(this,nparam,distparam,geom,grid,Flagbc,nslice)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc,nslice
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        !integer, parameter :: nslice = 100 !used in the 2nd 1B
        !integer, parameter :: nslice = 50 !used in the 3rd 1B
        !integer, parameter :: nslice = 25
        integer :: iz
        real*8, dimension(nslice) :: a11lc,a11,a12lc,a12,a22lc,a22,y1lc,&
                                     y1,y2lc,y2,aa,bb,rk,gam1,zz1
        real*8 :: hz,zmin,detaabb
        integer :: iizz 
        real*8 :: ziizz1,gamz1,ziizz1b,gamz1b
        real*8 :: ziizz2,gamz2,ziizz2b,gamz2b
        !half box size of repopulation
        real*8 :: xrp,pxrp,yrp,pyrp,zrp,pzrp

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        xrp = distparam(2)
        pxrp = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yrp = distparam(9)
        pyrp = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        pzrp = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        zmin = 1.0e10
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.tmptcl(6)) zmax = tmptcl(6) 
          if(zmin.ge.tmptcl(6)) zmin = tmptcl(6) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)

          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        zmin = zmin*Scfreq*2*Pi*zscale
        !avoid index overflow
        zmin = zmin - (zmax-zmin)*0.5e-5
        zmax = zmax + (zmax-zmin)*0.5e-5
        hz = (zmax-zmin)/nslice

        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        a11lc = 0.0
        a12lc = 0.0
        a22lc = 0.0
        y1lc = 0.0
        y2lc = 0.0
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
          iz = (Ptcl(6,j)-zmin)/hz + 1
          a11lc(iz) = a11lc(iz) + Ptcl(6,j)**2
          a12lc(iz) = a12lc(iz) + Ptcl(6,j)
          a22lc(iz) = a22lc(iz) + 1.0
          y1lc(iz) = y1lc(iz) + Ptcl(7,j)
          y2lc(iz) = y2lc(iz) + Ptcl(6,j)*Ptcl(7,j)
        enddo

        call MPI_ALLREDUCE(a11lc,a11,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a12lc,a12,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a22lc,a22,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y1lc,y1,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y2lc,y2,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        do iz = 1, nslice
          if(a22(iz).gt.0.0) then
            detaabb = a22(iz)*a11(iz)-a12(iz)**2
            aa(iz) = (a11(iz)*y1(iz)-a12(iz)*y2(iz))/detaabb
            bb(iz) = (a22(iz)*y2(iz)-a12(iz)*y1(iz))/detaabb
          else
            aa(iz) = 0.0
            bb(iz) = 0.0
          endif
        enddo

!-----------------------------------------------
!the following lines are used for the piecewise continuous of the fitting of E-phi
      iizz = nslice/2
      ziizz1 = zmin+(iizz-1)*hz
      gamz1 = aa(iizz)+bb(iizz)*ziizz1
      zz1(iizz) = ziizz1
      gam1(iizz) = gamz1
      rk(iizz) = bb(iizz)
      do i = iizz-1,1,-1
        rk(i) = (y2(i)-ziizz1*y1(i)+ziizz1*gamz1*a22(i)-gamz1*a12(i))/ &
                (a11(i)+ziizz1*ziizz1*a22(i)-2*ziizz1*a12(i))
        ziizz1b = zmin+(i-1)*hz
        gamz1b = gamz1 + rk(i)*(ziizz1b-ziizz1)
        gamz1 = gamz1b
        ziizz1 = ziizz1b
        zz1(i) = ziizz1
        gam1(i) = gamz1
      enddo
      ziizz2 = zmin+iizz*hz
      gamz2 = aa(iizz)+bb(iizz)*ziizz2
      do i = iizz+1,nslice
        rk(i) = (y2(i)-ziizz2*y1(i)+ziizz2*gamz2*a22(i)-gamz2*a12(i))/ &
                (a11(i)+ziizz2*ziizz2*a22(i)-2*ziizz2*a12(i))
        ziizz2b = zmin+i*hz
        gamz2b = gamz2 + rk(i)*(ziizz2b-ziizz2)
        zz1(i) = ziizz2
        gam1(i) = gamz2
        gamz2 = gamz2b
        ziizz2 = ziizz2b
      enddo
!-----------------------------------------------

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                !r = (2.0*r-1.0)*0.125*xmax
                r = (2.0*r-1.0)*xrp*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*0.03*pxmax
                 r = (2.0*r-1.0)*pxrp*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                !r = (2.0*r-1.0)*0.125*ymax
                r = (2.0*r-1.0)*yrp*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                r = (2.0*r-1.0)*pyrp*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                !r = (2.0*r-1.0)*0.02*zmax 
                r = (2.0*r-1.0)*hz 
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(6,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               !r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0) !15keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 2*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/15 !2keV is rms
               !r = 5*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/15 !5keV is rms
                r = pzrp*(2.0*r-1.0)*sqrt(3.0d0) !5keV is rms
              endif
              ii = j+nset*(i-1)
              iz = (this%Pts1(5,ii)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              if(nset.eq.1) then
                this%Pts1(6,ii) = Ptcl(7,i)+r
              else
                !this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+r
                this%Pts1(6,ii) = gam1(iz)+rk(iz)*&
                                (this%Pts1(5,ii)-zz1(iz))+r
              endif
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        !print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)

        !print*,"Nplocal: ",this%Nptlocal

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        this%Npt = nset*inipts

        end subroutine readElegant_Dist

        subroutine readElegantDB_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        !integer, parameter :: nslice = 100 !used in the 2nd 1B
        integer, parameter :: nslice = 50 !used in the 3rd 1B
        !integer, parameter :: nslice = 25
        integer :: iz
        real*8, dimension(nslice) :: a11lc,a11,a12lc,a12,a22lc,a22,y1lc,&
                                     y1,y2lc,y2,aa,bb
        real*8 :: hz,zmin,detaabb
        real*8, allocatable, dimension(:) :: rrdd
        real*8 :: aa0,aaseed
        parameter (aa0 = 3.d0**33 + 100.d0)
        integer :: nrand,ij

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          !print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        zmin = 1.0e10
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.tmptcl(6)) zmax = tmptcl(6) 
          if(zmin.ge.tmptcl(6)) zmin = tmptcl(6) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)

          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        !print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        zmin = zmin*Scfreq*2*Pi*zscale
        !avoid index overflow
        zmin = zmin - (zmax-zmin)*0.5e-5
        zmax = zmax + (zmax-zmin)*0.5e-5
        hz = (zmax-zmin)/nslice

        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        !print*,"xmax: ",xmax,ymax,zmin,zmax,pxmax,pymax,pzmax
        a11lc = 0.0
        a12lc = 0.0
        a22lc = 0.0
        y1lc = 0.0
        y2lc = 0.0
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
          iz = (Ptcl(6,j)-zmin)/hz + 1
          a11lc(iz) = a11lc(iz) + Ptcl(6,j)**2
          a12lc(iz) = a12lc(iz) + Ptcl(6,j)
          a22lc(iz) = a22lc(iz) + 1.0
          y1lc(iz) = y1lc(iz) + Ptcl(7,j)
          y2lc(iz) = y2lc(iz) + Ptcl(6,j)*Ptcl(7,j)
        enddo

        call MPI_ALLREDUCE(a11lc,a11,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a12lc,a12,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a22lc,a22,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y1lc,y1,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y2lc,y2,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        !print*,"a22: ",sum(a22)
        do iz = 1, nslice
          if(a22(iz).lt.1.0e-8) print*,"a22: ",iz,a22(iz)
        enddo
        do iz = 1, nslice
          if(a22(iz).gt.0.0) then
            detaabb = a22(iz)*a11(iz)-a12(iz)**2
            aa(iz) = (a11(iz)*y1(iz)-a12(iz)*y2(iz))/detaabb
            bb(iz) = (a22(iz)*y2(iz)-a12(iz)*y1(iz))/detaabb
          else
            aa(iz) = 0.0
            bb(iz) = 0.0
          endif
          if(myid.eq.0) print*,"isilce: ",iz,detaabb,aa(iz),bb(iz),a22(iz)
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        nrand = 6*nset*numpts0
        allocate(rrdd(nrand))
        aaseed = aa0 + 53*nrand*myid
        call bcnrand(nrand,aaseed,rrdd)
 
        sumx = 0.0
        ij = 0
        do i = 1, numpts0
          do j = 1, nset
              !call random_number(r)
              ij = ij + 1
              r = rrdd(ij)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                r = (2.0*r-1.0)*0.125*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              !call random_number(r)
              ij = ij + 1
              r = rrdd(ij)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              !call random_number(r)
              ij = ij + 1
              r = rrdd(ij)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                r = (2.0*r-1.0)*0.125*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              !call random_number(r)
              ij = ij + 1
              r = rrdd(ij)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              !call random_number(r)
              ij = ij + 1
              r = rrdd(ij)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                !r = (2.0*r-1.0)*0.02*zmax 
                r = (2.0*r-1.0)*hz 
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(6,i)+r
              !call random_number(r)
              ij = ij + 1
              r = rrdd(ij)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               !r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0) !15keV is rms
               r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 2*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/15 !15keV is rms
              endif
              if(i.eq.1 .and. j.eq.1) print*,"r, ",r,ij,aaseed
              if(i.eq.100 .and. j.eq.1) print*,"r, ",r,ij
              ii = j+nset*(i-1)
              iz = (this%Pts1(5,ii)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              !this%Pts1(6,ii) = Ptcl(7,i)+r
              this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        !print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)
        deallocate(rrdd)

        !print*,"Nplocal: ",this%Nptlocal

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        this%Npt = nset*inipts

        end subroutine readElegantDB_Dist

        !read in an initial distribution with format from ImpactT
        subroutine read_Dist(this)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer :: i,j,jlow,jhigh,avgpts,myid,nproc,ierr,nleft
        integer*8 :: nptot
        double precision, dimension(9) :: tmptcl
        double precision :: sum1,sum2
 
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
 
        open(unit=12,file='partcl.data',status='old')
 
        sum1 = 0.0
        sum2 = 0.0
          read(12,*)nptot
          avgpts = nptot/nproc
          nleft = nptot - avgpts*nproc
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
          !jlow = myid*avgpts + 1
          !jhigh = (myid+1)*avgpts
          print*,"avgpts, jlow, and jhigh: ",avgpts,jlow,jhigh
          do j = 1, nptot
            read(12,*)tmptcl(1:9)
            sum1 = sum1 + tmptcl(1)
            sum2 = sum2 + tmptcl(3)
            if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              this%Pts1(1:9,i) = tmptcl(1:9)
            endif
          enddo
          !<<<<<<<< commenting out unnecessary printing(kilean) <<<<<
          !print*,"sumx1,sumy1: ",sum1/nptot,sum2/nptot
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
          close(12)
 
          this%Nptlocal = avgpts
 
        end subroutine read_Dist
        
        
        !<<<<<<<<<< Binary pData read in ( Kilean ) <<<<<<<<<<
        !read in binary distribution of ImpactT format
        subroutine read_Dist_binary(this,nparam,distparam)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam
        double precision, dimension(nparam) :: distparam
        integer :: i,j,jlow,jhigh,avgpts,myid,nproc,ierr,nleft,fID
        integer*8 :: nptot
        double precision, dimension(9) :: tmptcl
        double precision :: sum1,sum2
 
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
        fID = int(distparam(1))
 
        open(unit=fID,form='unformatted',action='read')
 
        sum1 = 0.0
        sum2 = 0.0
        nptot = this%Npt
        avgpts = nptot/nproc
        nleft = nptot - avgpts*nproc
        if(myid.lt.nleft) then
          avgpts = avgpts+1
          jlow = myid*avgpts + 1
          jhigh = (myid+1)*avgpts
        else
          jlow = myid*avgpts + 1 + nleft
          jhigh = (myid+1)*avgpts + nleft
        endif
        allocate(this%Pts1(9,avgpts))
        this%Pts1 = 0.0
        do j = 1, nptot
          read(fID)tmptcl(1:9)
          sum1 = sum1 + tmptcl(1)
          sum2 = sum2 + tmptcl(3)
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
            i = j - jlow + 1
            this%Pts1(1:9,i) = tmptcl(1:9)
          endif
        enddo
        close(fID)
 
        this%Nptlocal = avgpts
        
        end subroutine read_Dist_binary
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        
        !read particle distribution from Elegant output
        subroutine readElegant2_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(8) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz,phi0lc

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        phi0lc = this%refptcl(5)
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(8,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do j = 1, inipts
          read(12,*)tmptcl(1:8)
          !read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(1))) xmax = abs(tmptcl(1)) 
          if(ymax.le.abs(tmptcl(3))) ymax = abs(tmptcl(3)) 
          if(zmax.le.abs(tmptcl(5))) zmax = abs(tmptcl(5)) 
          gammabet = tmptcl(6)/sqrt(1.d0+tmptcl(2)**2+tmptcl(4)**2)
          tmppx = tmptcl(2)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(4)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(6)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(1)
          sumx2 = sumx2 + tmptcl(1)*tmptcl(1)
          !sumy = sumy + tmptcl(3)
          !sumy2 = sumy2 + tmptcl(3)*tmptcl(3)
          sumy = sumy + tmptcl(5)
          sumy2 = sumy2 + tmptcl(5)*tmptcl(5)
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:8,i) = tmptcl(1:8)
              !Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        print*,"sumx2: ",sumx2,sumx/inipts,sumy2,sumy/inipts
        close(12)
        call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        print*,"xmax: ",xmax,ymax,zmax,pxmax,pymax,pzmax
        do j = 1, avgpts
          Ptcl(1,j) = (Ptcl(1,j)/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/xl)*yscale + xmu3
          !Ptcl(5,j) = (Ptcl(5,j)*Scfreq*2*Pi-phi0lc)*zscale + xmu5
          Ptcl(5,j) = (Ptcl(5,j)-sumy/inipts)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(6,j)/sqrt(1.d0+Ptcl(2,j)**2+Ptcl(4,j)**2)
          Ptcl(2,j) = Ptcl(2,j)*gammabet*pxscale + xmu2
          Ptcl(4,j) = Ptcl(4,j)*gammabet*pyscale + xmu4
          Ptcl(6,j) = (gamma0 - sqrt(1.0 + Ptcl(6,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                r = (2.0*r-1.0)*0.03*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(1,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(2,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                r = (2.0*r-1.0)*0.03*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(3,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(4,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.005*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                r = (2.0*r-1.0)*0.01*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                !r = (2.0*r-1.0)*0.004*pzmax
                r = (2.0*r-1.0)*0.01*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal

        deallocate(Ptcl)

        print*,"Nplocal: ",this%Nptlocal

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        end subroutine readElegant2_Dist

        subroutine readElegantRot_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pilc
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz,xtmp,pxtmp,ytmp,pytmp,theta1,theta2,&
                  theta1i,theta2i,dtheta1,dtheta2,rr1,rr2
        integer :: Nrot,ipt,ntmp

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pilc = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts,Nrot

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.abs(tmptcl(6))) zmax = abs(tmptcl(6)) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        dtheta1 = 2*pilc/Nrot
        dtheta2 = 2*pilc/Nrot
        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        print*,"xmax: ",xmax,ymax,zmax,pxmax,pymax,pzmax
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
        enddo

        numpts0 = avgpts 

        nset = nptot/(inipts*Nrot)
        nremain = nptot - nset*inipts*Nrot
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0*Nrot

        allocate(this%Pts1(9,nset*numpts0*Nrot))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                r = (2.0*r-1.0)*0.03*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                r = (2.0*r-1.0)*0.03*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                r = (2.0*r-1.0)*0.001*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(6,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               r = (2.0*r-1.0)*0.001*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(7,i)+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)

        print*,"Nplocal: ",this%Nptlocal,nset,Nrot

        if(Nrot.gt.1) then
          ntmp = numpts0*nset
          ipt = ntmp
          do j = 1, ntmp
            if(this%Pts1(1,j).gt.0.0d0) then
              theta1 = atan(this%Pts1(3,j)/this%Pts1(1,j)) 
            else if(this%Pts1(1,j).lt.0.0d0) then
              theta1 = pilc+atan(this%Pts1(3,j)/this%Pts1(1,j)) 
            else
              if(this%Pts1(3,j).gt.0.0d0) then
                theta1 = pilc/2
              else
                theta1 = 3.0d0*pilc/2
              endif
            endif
            if(this%Pts1(2,j).gt.0.0d0) then
              theta2 = atan(this%Pts1(4,j)/this%Pts1(2,j)) 
            else if(this%Pts1(2,j).lt.0.0d0) then
              theta2 = pilc+atan(this%Pts1(4,j)/this%Pts1(2,j)) 
            else
              if(this%Pts1(4,j).gt.0.0d0) then
                theta2 = pilc/2
              else
                theta2 = 3.0d0*pilc/2
              endif
            endif
            rr1 = sqrt(this%Pts1(1,j)**2+this%Pts1(3,j)**2)
            rr2 = sqrt(this%Pts1(2,j)**2+this%Pts1(4,j)**2)
            do i = 2, Nrot
              theta1i = theta1 + (i-1)*dtheta1
              theta2i = theta2 + (i-1)*dtheta2
              xtmp = rr1*cos(theta1i) 
              pxtmp = rr2*cos(theta2i) 
              ytmp = rr1*sin(theta1i) 
              pytmp = rr2*sin(theta2i) 
              ipt = ipt + 1 
              this%Pts1(1,ipt) = xtmp
              this%Pts1(2,ipt) = pxtmp
              this%Pts1(3,ipt) = ytmp
              this%Pts1(4,ipt) = pytmp
              this%Pts1(5,ipt) = this%Pts1(5,j)
              this%Pts1(6,ipt) = this%Pts1(6,j)
            enddo
          enddo
        endif
        print*,"ipt: ",ipt,dtheta1,dtheta2

        jlow = (jlow-1)*nset*Nrot
        do j = 1, this%Nptlocal
          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        end subroutine readElegantRot_Dist

        !remove the transverse-longitudinal correlation by assuming a given
        !transverse Gaussian distribution.
        subroutine readElegantNcorold_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2
        double precision  :: sigx,sigpx,muxpx,sigy,&
        sigpy,muypy,sigz,sigpz,muzpz

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale
 
        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.abs(tmptcl(6))) zmax = abs(tmptcl(6)) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        print*,"xmax: ",xmax,ymax,zmax,pxmax,pymax,pzmax
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(x1(2,nset*numpts0))
        allocate(x2(2,nset*numpts0))
 
        call normVec(x1,nset*numpts0)
        call normVec(x2,nset*numpts0)

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                r = (2.0*r-1.0)*0.03*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              !Gaussian distribution
              this%Pts1(1,ii) = xmu1 + sig1*x1(1,ii)/sq12
              this%Pts1(2,ii) = xmu2 + sig2*(-muxpx*x1(1,ii)/sq12+x1(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                r = (2.0*r-1.0)*0.03*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              !transverse Gaussian distribution
              this%Pts1(3,ii) = xmu3 + sig3*x2(1,ii)/sq34
              this%Pts1(4,ii) = xmu4 + sig4*(-muypy*x2(1,ii)/sq34+x2(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                r = (2.0*r-1.0)*0.001*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(6,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               r = (2.0*r-1.0)*0.016*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(7,i)+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)

        print*,"Nplocal: ",this%Nptlocal

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        this%Npt = nset*inipts

        deallocate(x1)
        deallocate(x2)

        end subroutine readElegantNcorold_Dist

        subroutine readElegantNcorold2_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        !integer, parameter :: nslice = 100 !used in the 2nd 1B
        !integer, parameter :: nslice = 50 !used in the 3rd 1B
        !integer, parameter :: nslice = 25
        integer, parameter :: nslice = 1
        integer :: iz
        real*8, dimension(nslice) :: a11lc,a11,a12lc,a12,a22lc,a22,y1lc,&
                                     y1,y2lc,y2,aa,bb
        real*8 :: hz,zmin,detaabb
        double precision, allocatable, dimension(:,:) :: x1,x2
        double precision  :: sigx,sigpx,muxpx,sigy,&
        sigpy,muypy,sigz,sigpz,muzpz
        real*8 :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        real*8, dimension(2) :: yy2

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale
 
        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          !print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        zmin = 1.0e10
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.tmptcl(6)) zmax = tmptcl(6) 
          if(zmin.ge.tmptcl(6)) zmin = tmptcl(6) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)

          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        !print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        zmin = zmin*Scfreq*2*Pi*zscale
        !avoid index overflow
        zmin = zmin - (zmax-zmin)*0.5e-5
        zmax = zmax + (zmax-zmin)*0.5e-5
        hz = (zmax-zmin)/nslice

        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        !print*,"xmax: ",xmax,ymax,zmin,zmax,pxmax,pymax,pzmax
        a11lc = 0.0
        a12lc = 0.0
        a22lc = 0.0
        y1lc = 0.0
        y2lc = 0.0
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
          iz = (Ptcl(6,j)-zmin)/hz + 1
          a11lc(iz) = a11lc(iz) + Ptcl(6,j)**2
          a12lc(iz) = a12lc(iz) + Ptcl(6,j)
          a22lc(iz) = a22lc(iz) + 1.0
          y1lc(iz) = y1lc(iz) + Ptcl(7,j)
          y2lc(iz) = y2lc(iz) + Ptcl(6,j)*Ptcl(7,j)
        enddo

        call MPI_ALLREDUCE(a11lc,a11,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a12lc,a12,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a22lc,a22,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y1lc,y1,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y2lc,y2,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        !print*,"a22: ",sum(a22)
        do iz = 1, nslice
          if(a22(iz).lt.1.0e-8) print*,"a22: ",iz,a22(iz)
        enddo
        do iz = 1, nslice
          if(a22(iz).gt.0.0) then
            detaabb = a22(iz)*a11(iz)-a12(iz)**2
            aa(iz) = (a11(iz)*y1(iz)-a12(iz)*y2(iz))/detaabb
            bb(iz) = (a22(iz)*y2(iz)-a12(iz)*y1(iz))/detaabb
          else
            aa(iz) = 0.0
            bb(iz) = 0.0
          endif
          if(myid.eq.0) print*,"isilce: ",i,detaabb,aa(iz),bb(iz),a22(iz)
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(x1(2,nset*numpts0))
        allocate(x2(2,nset*numpts0))
 
        call normVec(x1,nset*numpts0)
        call normVec(x2,nset*numpts0)

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                r = (2.0*r-1.0)*0.125*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              !Gaussian distribution
              this%Pts1(1,ii) = xmu1 + sig1*x1(1,ii)/sq12
              this%Pts1(2,ii) = xmu2 + sig2*(-muxpx*x1(1,ii)/sq12+x1(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                r = (2.0*r-1.0)*0.125*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              !transverse Gaussian distribution
              this%Pts1(3,ii) = xmu3 + sig3*x2(1,ii)/sq34
              this%Pts1(4,ii) = xmu4 + sig4*(-muypy*x2(1,ii)/sq34+x2(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                !r = (2.0*r-1.0)*0.02*zmax 
                !r = (2.0*r-1.0)*hz
                r = (2.0*r-1.0)*hz/2
              endif
              ii = j+nset*(i-1)
              !this%Pts1(5,ii) = Ptcl(6,i)+r
              this%Pts1(5,ii) = r
              !call random_number(r)
              call normdv(yy2)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               !r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0) !15keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               r = yy2(1)*0.029354/15*3 !3keV is rms Gaussian
              endif
              ii = j+nset*(i-1)
              iz = (this%Pts1(5,ii)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              !this%Pts1(6,ii) = Ptcl(7,i)+r
              !this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+r
              this%Pts1(6,ii) = r !longitudinal uniform
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        !print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)

        !print*,"Nplocal: ",this%Nptlocal

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        this%Npt = nset*inipts

        end subroutine readElegantNcorold2_Dist

        subroutine readElegantlaser_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        !integer, parameter :: nslice = 100 !used in the 2nd 1B
        integer, parameter :: nslice = 50 !used in the 3rd 1B
        !integer, parameter :: nslice = 25
        integer :: iz
        real*8, dimension(nslice) :: a11lc,a11,a12lc,a12,a22lc,a22,y1lc,&
                                     y1,y2lc,y2,aa,bb
        real*8 :: hz,zmin,detaabb,dee

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          !print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        zmin = 1.0e10
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.tmptcl(6)) zmax = tmptcl(6) 
          if(zmin.ge.tmptcl(6)) zmin = tmptcl(6) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)

          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        !print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        zmin = zmin*Scfreq*2*Pi*zscale
        !avoid index overflow
        zmin = zmin - (zmax-zmin)*0.5e-5
        zmax = zmax + (zmax-zmin)*0.5e-5
        hz = (zmax-zmin)/nslice

        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        !print*,"xmax: ",xmax,ymax,zmin,zmax,pxmax,pymax,pzmax
        a11lc = 0.0
        a12lc = 0.0
        a22lc = 0.0
        y1lc = 0.0
        y2lc = 0.0
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
          iz = (Ptcl(6,j)-zmin)/hz + 1
          a11lc(iz) = a11lc(iz) + Ptcl(6,j)**2
          a12lc(iz) = a12lc(iz) + Ptcl(6,j)
          a22lc(iz) = a22lc(iz) + 1.0
          y1lc(iz) = y1lc(iz) + Ptcl(7,j)
          y2lc(iz) = y2lc(iz) + Ptcl(6,j)*Ptcl(7,j)
        enddo

        call MPI_ALLREDUCE(a11lc,a11,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a12lc,a12,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a22lc,a22,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y1lc,y1,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y2lc,y2,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        !print*,"a22: ",sum(a22)
        do iz = 1, nslice
          if(a22(iz).lt.1.0e-8) print*,"a22: ",iz,a22(iz)
        enddo
        do iz = 1, nslice
          if(a22(iz).gt.0.0) then
            detaabb = a22(iz)*a11(iz)-a12(iz)**2
            aa(iz) = (a11(iz)*y1(iz)-a12(iz)*y2(iz))/detaabb
            bb(iz) = (a22(iz)*y2(iz)-a12(iz)*y1(iz))/detaabb
          else
            aa(iz) = 0.0
            bb(iz) = 0.0
          endif
          if(myid.eq.0) print*,"isilce: ",i,detaabb,aa(iz),bb(iz),a22(iz)
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                r = (2.0*r-1.0)*0.125*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                r = (2.0*r-1.0)*0.125*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                !r = (2.0*r-1.0)*0.02*zmax 
                r = (2.0*r-1.0)*hz 
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(6,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               !r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0) !15keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/15/10 !1keV /10 is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*5.0 !1keV /10 is rms
               r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*7.0 !1keV /10 is rms
              endif
              ii = j+nset*(i-1)

              iz = (Ptcl(6,i)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              dee = Ptcl(7,i)-(aa(iz)+bb(iz)*Ptcl(6,i))

              iz = (this%Pts1(5,ii)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              !this%Pts1(6,ii) = Ptcl(7,i)+r

              !if r=0, there is no repopulation in dE
              !if dee = 0, uniform resampling in dE
              this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+dee+r

          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        !print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)

        !print*,"Nplocal: ",this%Nptlocal

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        this%Npt = nset*inipts

        end subroutine readElegantlaser_Dist

        subroutine readElegantlaser2_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        !integer, parameter :: nslice = 100 !used in the 2nd 1B
        integer, parameter :: nslice = 50 !used in the 3rd 1B
        !integer, parameter :: nslice = 25
        integer :: iz
        real*8, dimension(nslice) :: a11lc,a11,a12lc,a12,a22lc,a22,y1lc,&
                                     y1,y2lc,y2,aa,bb
        real*8 :: hz,zmin,detaabb,dee
        real*8, allocatable, dimension(:) :: rrdd
        real*8 :: aa0,aaseed
        parameter (aa0 = 3.d0**33 + 100.d0)
        integer :: nrand,ij

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
!        call random_number(xx)
!        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          !print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        zmin = 1.0e10
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.tmptcl(6)) zmax = tmptcl(6) 
          if(zmin.ge.tmptcl(6)) zmin = tmptcl(6) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)

          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        !print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        zmin = zmin*Scfreq*2*Pi*zscale
        !avoid index overflow
        zmin = zmin - (zmax-zmin)*0.5e-5
        zmax = zmax + (zmax-zmin)*0.5e-5
        hz = (zmax-zmin)/nslice

        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        !print*,"xmax: ",xmax,ymax,zmin,zmax,pxmax,pymax,pzmax
        a11lc = 0.0
        a12lc = 0.0
        a22lc = 0.0
        y1lc = 0.0
        y2lc = 0.0
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
          iz = (Ptcl(6,j)-zmin)/hz + 1
          a11lc(iz) = a11lc(iz) + Ptcl(6,j)**2
          a12lc(iz) = a12lc(iz) + Ptcl(6,j)
          a22lc(iz) = a22lc(iz) + 1.0
          y1lc(iz) = y1lc(iz) + Ptcl(7,j)
          y2lc(iz) = y2lc(iz) + Ptcl(6,j)*Ptcl(7,j)
        enddo

        call MPI_ALLREDUCE(a11lc,a11,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a12lc,a12,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a22lc,a22,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y1lc,y1,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y2lc,y2,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        !print*,"a22: ",sum(a22)
        do iz = 1, nslice
          if(a22(iz).lt.1.0e-8) print*,"a22: ",iz,a22(iz)
        enddo
        do iz = 1, nslice
          if(a22(iz).gt.0.0) then
            detaabb = a22(iz)*a11(iz)-a12(iz)**2
            aa(iz) = (a11(iz)*y1(iz)-a12(iz)*y2(iz))/detaabb
            bb(iz) = (a22(iz)*y2(iz)-a12(iz)*y1(iz))/detaabb
          else
            aa(iz) = 0.0
            bb(iz) = 0.0
          endif
          if(myid.eq.0) print*,"isilce: ",i,detaabb,aa(iz),bb(iz),a22(iz)
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        nrand = 6*nset*numpts0
        allocate(rrdd(nrand))
        aaseed = aa0 + 53*nrand*myid
        call bcnrand(nrand,aaseed,rrdd)

        sumx = 0.0
        ij = 0
        do i = 1, numpts0
          do j = 1, nset
              ij = ij + 1
              !call random_number(r)
              r = rrdd(ij)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                r = (2.0*r-1.0)*0.125*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              !call random_number(r)
              ij = ij + 1
              r = rrdd(ij)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              !call random_number(r)
              ij = ij + 1
              r = rrdd(ij) 
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                r = (2.0*r-1.0)*0.125*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              !call random_number(r)
              ij = ij + 1
              r = rrdd(ij) 
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              !call random_number(r)
              ij = ij + 1
              r = rrdd(ij) 
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                !r = (2.0*r-1.0)*0.02*zmax 
                r = (2.0*r-1.0)*hz 
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(6,i)+r
              !call random_number(r)
              ij = ij + 1
              r = rrdd(ij) 
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               !r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0) !15keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/15/10 !1keV /10 is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*5.0 !1keV /10 is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*7.0 !1keV /10 is rms
              endif
              ii = j+nset*(i-1)

              iz = (Ptcl(6,i)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              dee = Ptcl(7,i)-(aa(iz)+bb(iz)*Ptcl(6,i))

              iz = (this%Pts1(5,ii)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              !this%Pts1(6,ii) = Ptcl(7,i)+r

              !if r=0, there is no repopulation in dE
              !if dee = 0, uniform resampling in dE
              this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+dee+r

          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        !print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)
        deallocate(rrdd)

        !print*,"Nplocal: ",this%Nptlocal

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        this%Npt = nset*inipts

        end subroutine readElegantlaser2_Dist

        subroutine readElegantNcor_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        !integer, parameter :: nslice = 100 !used in the 2nd 1B
        integer, parameter :: nslice = 50 !used in the 3rd 1B
        !integer, parameter :: nslice = 25
        integer :: iz
        real*8, dimension(nslice) :: a11lc,a11,a12lc,a12,a22lc,a22,y1lc,&
                                     y1,y2lc,y2,aa,bb
        real*8 :: hz,zmin,detaabb
        double precision, allocatable, dimension(:,:) :: x1,x2
        double precision  :: sigx,sigpx,muxpx,sigy,&
        sigpy,muypy,sigz,sigpz,muzpz
        real*8 :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale
 
        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          !print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        zmin = 1.0e10
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.tmptcl(6)) zmax = tmptcl(6) 
          if(zmin.ge.tmptcl(6)) zmin = tmptcl(6) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)

          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        !print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        zmin = zmin*Scfreq*2*Pi*zscale
        !avoid index overflow
        zmin = zmin - (zmax-zmin)*0.5e-5
        zmax = zmax + (zmax-zmin)*0.5e-5
        hz = (zmax-zmin)/nslice

        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        !print*,"xmax: ",xmax,ymax,zmin,zmax,pxmax,pymax,pzmax
        a11lc = 0.0
        a12lc = 0.0
        a22lc = 0.0
        y1lc = 0.0
        y2lc = 0.0
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
          iz = (Ptcl(6,j)-zmin)/hz + 1
          a11lc(iz) = a11lc(iz) + Ptcl(6,j)**2
          a12lc(iz) = a12lc(iz) + Ptcl(6,j)
          a22lc(iz) = a22lc(iz) + 1.0
          y1lc(iz) = y1lc(iz) + Ptcl(7,j)
          y2lc(iz) = y2lc(iz) + Ptcl(6,j)*Ptcl(7,j)
        enddo

        call MPI_ALLREDUCE(a11lc,a11,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a12lc,a12,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a22lc,a22,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y1lc,y1,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y2lc,y2,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        !print*,"a22: ",sum(a22)
        do iz = 1, nslice
          if(a22(iz).lt.1.0e-8) print*,"a22: ",iz,a22(iz)
        enddo
        do iz = 1, nslice
          if(a22(iz).gt.0.0) then
            detaabb = a22(iz)*a11(iz)-a12(iz)**2
            aa(iz) = (a11(iz)*y1(iz)-a12(iz)*y2(iz))/detaabb
            bb(iz) = (a22(iz)*y2(iz)-a12(iz)*y1(iz))/detaabb
          else
            aa(iz) = 0.0
            bb(iz) = 0.0
          endif
          if(myid.eq.0) print*,"isilce: ",i,detaabb,aa(iz),bb(iz),a22(iz)
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(x1(2,nset*numpts0))
        allocate(x2(2,nset*numpts0))
 
        call normVec(x1,nset*numpts0)
        call normVec(x2,nset*numpts0)

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                r = (2.0*r-1.0)*0.125*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              !Gaussian distribution
              this%Pts1(1,ii) = xmu1 + sig1*x1(1,ii)/sq12
              this%Pts1(2,ii) = xmu2 + sig2*(-muxpx*x1(1,ii)/sq12+x1(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                r = (2.0*r-1.0)*0.125*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              !transverse Gaussian distribution
              this%Pts1(3,ii) = xmu3 + sig3*x2(1,ii)/sq34
              this%Pts1(4,ii) = xmu4 + sig4*(-muypy*x2(1,ii)/sq34+x2(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                !r = (2.0*r-1.0)*0.02*zmax 
                r = (2.0*r-1.0)*hz 
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(6,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               !r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
               r = (2.0*r-1.0)*0.029354*sqrt(3.0d0) !15keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
              endif
              ii = j+nset*(i-1)
              iz = (this%Pts1(5,ii)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              !this%Pts1(6,ii) = Ptcl(7,i)+r
              this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        !print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)

        !print*,"Nplocal: ",this%Nptlocal

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        this%Npt = nset*inipts

        deallocate(x1)
        deallocate(x2)

        end subroutine readElegantNcor_Dist

        subroutine readElegantNcorMod_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        !integer, parameter :: nslice = 100 !used in the 2nd 1B
        integer, parameter :: nslice = 50 !used in the 3rd 1B
        !integer, parameter :: nslice = 25
        integer :: iz
        real*8, dimension(nslice) :: a11lc,a11,a12lc,a12,a22lc,a22,y1lc,&
                                     y1,y2lc,y2,aa,bb
        real*8 :: hz,zmin,detaabb
        double precision, allocatable, dimension(:,:) :: x1,x2
        double precision  :: sigx,sigpx,muxpx,sigy,&
        sigpy,muypy,sigz,sigpz,muzpz
        real*8 :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        real*8 :: zz,zbunch,rkk,sig1z,sig3z,ampz

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale
 
        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          !print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        zmin = 1.0e10
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.tmptcl(6)) zmax = tmptcl(6) 
          if(zmin.ge.tmptcl(6)) zmin = tmptcl(6) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)

          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        !print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        zmin = zmin*Scfreq*2*Pi*zscale
        !avoid index overflow
        zmin = zmin - (zmax-zmin)*0.5e-5
        zmax = zmax + (zmax-zmin)*0.5e-5
        hz = (zmax-zmin)/nslice

        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        !print*,"xmax: ",xmax,ymax,zmin,zmax,pxmax,pymax,pzmax
        a11lc = 0.0
        a12lc = 0.0
        a22lc = 0.0
        y1lc = 0.0
        y2lc = 0.0
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
          iz = (Ptcl(6,j)-zmin)/hz + 1
          a11lc(iz) = a11lc(iz) + Ptcl(6,j)**2
          a12lc(iz) = a12lc(iz) + Ptcl(6,j)
          a22lc(iz) = a22lc(iz) + 1.0
          y1lc(iz) = y1lc(iz) + Ptcl(7,j)
          y2lc(iz) = y2lc(iz) + Ptcl(6,j)*Ptcl(7,j)
        enddo

        call MPI_ALLREDUCE(a11lc,a11,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a12lc,a12,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a22lc,a22,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y1lc,y1,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y2lc,y2,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        !print*,"a22: ",sum(a22)
        do iz = 1, nslice
          if(a22(iz).lt.1.0e-8) print*,"a22: ",iz,a22(iz)
        enddo
        do iz = 1, nslice
          if(a22(iz).gt.0.0) then
            detaabb = a22(iz)*a11(iz)-a12(iz)**2
            aa(iz) = (a11(iz)*y1(iz)-a12(iz)*y2(iz))/detaabb
            bb(iz) = (a22(iz)*y2(iz)-a12(iz)*y1(iz))/detaabb
          else
            aa(iz) = 0.0
            bb(iz) = 0.0
          endif
          if(myid.eq.0) print*,"isilce: ",i,detaabb,aa(iz),bb(iz),a22(iz)
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(x1(2,nset*numpts0))
        allocate(x2(2,nset*numpts0))
 
        call normVec(x1,nset*numpts0)
        call normVec(x2,nset*numpts0)

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                r = (2.0*r-1.0)*0.125*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              !Gaussian distribution
              this%Pts1(1,ii) = xmu1 + sig1*x1(1,ii)/sq12
              this%Pts1(2,ii) = xmu2 + sig2*(-muxpx*x1(1,ii)/sq12+x1(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                r = (2.0*r-1.0)*0.125*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              !transverse Gaussian distribution
              this%Pts1(3,ii) = xmu3 + sig3*x2(1,ii)/sq34
              this%Pts1(4,ii) = xmu4 + sig4*(-muypy*x2(1,ii)/sq34+x2(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                !r = (2.0*r-1.0)*0.02*zmax 
                r = (2.0*r-1.0)*hz 
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(6,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               !r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0) !15keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*2.0d0 !2keV is rms
              endif
              ii = j+nset*(i-1)
              iz = (this%Pts1(5,ii)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              !this%Pts1(6,ii) = Ptcl(7,i)+r
              this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        !print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)

        !print*,"Nplocal: ",this%Nptlocal

        zmin = zmin - hz
        zmax = zmax + hz
        zbunch = zmax-zmin
        rkk = 2*asin(1.0d0)*2/zbunch
        rkk = 2*2*asin(1.0d0)*2/zbunch
        rkk = 4*2*asin(1.0d0)*2/zbunch
        rkk = 512*2*asin(1.0d0)*2/zbunch
        ampz = 0.2d0
        ampz = 0.0d0
        ampz = 0.2d0

        if(myid.eq.0) then
          print*,"zmin: ",zmin,zmax,zbunch,sig1,sig3,sq12,sq34
        endif
        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
         
          zz = this%Pts1(5,j)-zmin
          sig1z = sig1*(1.0+ampz*sin(rkk*zz))
          sig3z = sig3*(1.0+ampz*sin(rkk*zz))
          this%Pts1(1,j) = xmu1 + sig1z*x1(1,j)/sq12
          this%Pts1(3,j) = xmu3 + sig3z*x2(1,j)/sq34

          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        deallocate(x1)
        deallocate(x2)

        this%Npt = nset*inipts

        end subroutine readElegantNcorMod_Dist

        !longitudinal uniform density distribution with linear chirp
        subroutine readElegantNcorMod2_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        !integer, parameter :: nslice = 100 !used in the 2nd 1B
        !integer, parameter :: nslice = 50 !used in the 3rd 1B
        !integer, parameter :: nslice = 25
        integer, parameter :: nslice = 1 !used in the 3rd 1B
        integer :: iz
        real*8, dimension(nslice) :: a11lc,a11,a12lc,a12,a22lc,a22,y1lc,&
                                     y1,y2lc,y2,aa,bb
        real*8 :: hz,zmin,detaabb
        double precision, allocatable, dimension(:,:) :: x1,x2
        double precision  :: sigx,sigpx,muxpx,sigy,&
        sigpy,muypy,sigz,sigpz,muzpz
        real*8 :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        real*8 :: zz,zbunch,rkk,sig1z,sig3z,ampz

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale
 
        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          !print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        zmin = 1.0e10
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.tmptcl(6)) zmax = tmptcl(6) 
          if(zmin.ge.tmptcl(6)) zmin = tmptcl(6) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)

          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        !print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        zmin = zmin*Scfreq*2*Pi*zscale
        !avoid index overflow
        zmin = zmin - (zmax-zmin)*0.5e-5
        zmax = zmax + (zmax-zmin)*0.5e-5
        hz = (zmax-zmin)/nslice

        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        !print*,"xmax: ",xmax,ymax,zmin,zmax,pxmax,pymax,pzmax
        a11lc = 0.0
        a12lc = 0.0
        a22lc = 0.0
        y1lc = 0.0
        y2lc = 0.0
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
          iz = (Ptcl(6,j)-zmin)/hz + 1
          a11lc(iz) = a11lc(iz) + Ptcl(6,j)**2
          a12lc(iz) = a12lc(iz) + Ptcl(6,j)
          a22lc(iz) = a22lc(iz) + 1.0
          y1lc(iz) = y1lc(iz) + Ptcl(7,j)
          y2lc(iz) = y2lc(iz) + Ptcl(6,j)*Ptcl(7,j)
        enddo

        call MPI_ALLREDUCE(a11lc,a11,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a12lc,a12,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a22lc,a22,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y1lc,y1,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y2lc,y2,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        !print*,"a22: ",sum(a22)
        do iz = 1, nslice
          if(a22(iz).lt.1.0e-8) print*,"a22: ",iz,a22(iz)
        enddo
        do iz = 1, nslice
          if(a22(iz).gt.0.0) then
            detaabb = a22(iz)*a11(iz)-a12(iz)**2
            aa(iz) = (a11(iz)*y1(iz)-a12(iz)*y2(iz))/detaabb
            bb(iz) = (a22(iz)*y2(iz)-a12(iz)*y1(iz))/detaabb
          else
            aa(iz) = 0.0
            bb(iz) = 0.0
          endif
          if(myid.eq.0) print*,"isilce: ",i,detaabb,aa(iz),bb(iz),a22(iz)
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(x1(2,nset*numpts0))
        allocate(x2(2,nset*numpts0))
 
        call normVec(x1,nset*numpts0)
        call normVec(x2,nset*numpts0)

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                r = (2.0*r-1.0)*0.125*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              !Gaussian distribution
              this%Pts1(1,ii) = xmu1 + sig1*x1(1,ii)/sq12
              this%Pts1(2,ii) = xmu2 + sig2*(-muxpx*x1(1,ii)/sq12+x1(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                r = (2.0*r-1.0)*0.125*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              !transverse Gaussian distribution
              this%Pts1(3,ii) = xmu3 + sig3*x2(1,ii)/sq34
              this%Pts1(4,ii) = xmu4 + sig4*(-muypy*x2(1,ii)/sq34+x2(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                !r = (2.0*r-1.0)*0.02*zmax 
                r = (2.0*r-1.0)*hz 
              endif
              ii = j+nset*(i-1)
              !this%Pts1(5,ii) = Ptcl(6,i)+r
              this%Pts1(5,ii) = r/2 !uniform longiudinal distribution from -hz/2 to hz/2
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               !r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0) !15keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*2.0d0 !2keV is rms
              endif
              ii = j+nset*(i-1)
              iz = (this%Pts1(5,ii)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              !this%Pts1(6,ii) = Ptcl(7,i)+r
              this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        !print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)

        !print*,"Nplocal: ",this%Nptlocal

        zmin = zmin - hz
        zmax = zmax + hz
        zbunch = zmax-zmin
        rkk = 2*asin(1.0d0)*2/zbunch
        rkk = 2*2*asin(1.0d0)*2/zbunch
        rkk = 4*2*asin(1.0d0)*2/zbunch
        rkk = 512*2*asin(1.0d0)*2/zbunch
        ampz = 0.2d0
        ampz = 0.0d0
        ampz = 0.2d0
        ampz = 0.0d0

        if(myid.eq.0) then
          print*,"zmin: ",zmin,zmax,zbunch,sig1,sig3,sq12,sq34
        endif
        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
         
          zz = this%Pts1(5,j)-zmin
          sig1z = sig1*(1.0+ampz*sin(rkk*zz))
          sig3z = sig3*(1.0+ampz*sin(rkk*zz))
          this%Pts1(1,j) = xmu1 + sig1z*x1(1,j)/sq12
          this%Pts1(3,j) = xmu3 + sig3z*x2(1,j)/sq34

          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        deallocate(x1)
        deallocate(x2)

        this%Npt = nset*inipts

        end subroutine readElegantNcorMod2_Dist

        !longitudinal uniform density distribution but with nonlinear chirp
        subroutine readElegantNcorMod3_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        !integer, parameter :: nslice = 100 !used in the 2nd 1B
        integer, parameter :: nslice = 50 !used in the 3rd 1B
        !integer, parameter :: nslice = 25
        !integer, parameter :: nslice = 1 !used in the 3rd 1B
        integer :: iz
        real*8, dimension(nslice) :: a11lc,a11,a12lc,a12,a22lc,a22,y1lc,&
                                     y1,y2lc,y2,aa,bb
        real*8 :: hz,zmin,detaabb
        double precision, allocatable, dimension(:,:) :: x1,x2
        double precision  :: sigx,sigpx,muxpx,sigy,&
        sigpy,muypy,sigz,sigpz,muzpz
        real*8 :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        real*8 :: zz,zbunch,rkk,sig1z,sig3z,ampz,zsize

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale
 
        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          !print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        zmin = 1.0e10
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.tmptcl(6)) zmax = tmptcl(6) 
          if(zmin.ge.tmptcl(6)) zmin = tmptcl(6) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)

          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        !print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        zmin = zmin*Scfreq*2*Pi*zscale
        !avoid index overflow
        zmin = zmin - (zmax-zmin)*0.5e-5
        zmax = zmax + (zmax-zmin)*0.5e-5
        hz = (zmax-zmin)/nslice
        zsize = (zmax-zmin)

        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        !print*,"xmax: ",xmax,ymax,zmin,zmax,pxmax,pymax,pzmax
        a11lc = 0.0
        a12lc = 0.0
        a22lc = 0.0
        y1lc = 0.0
        y2lc = 0.0
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
          iz = (Ptcl(6,j)-zmin)/hz + 1
          a11lc(iz) = a11lc(iz) + Ptcl(6,j)**2
          a12lc(iz) = a12lc(iz) + Ptcl(6,j)
          a22lc(iz) = a22lc(iz) + 1.0
          y1lc(iz) = y1lc(iz) + Ptcl(7,j)
          y2lc(iz) = y2lc(iz) + Ptcl(6,j)*Ptcl(7,j)
        enddo

        call MPI_ALLREDUCE(a11lc,a11,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a12lc,a12,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a22lc,a22,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y1lc,y1,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y2lc,y2,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        !print*,"a22: ",sum(a22)
        do iz = 1, nslice
          if(a22(iz).lt.1.0e-8) print*,"a22: ",iz,a22(iz)
        enddo
        do iz = 1, nslice
          if(a22(iz).gt.0.0) then
            detaabb = a22(iz)*a11(iz)-a12(iz)**2
            aa(iz) = (a11(iz)*y1(iz)-a12(iz)*y2(iz))/detaabb
            bb(iz) = (a22(iz)*y2(iz)-a12(iz)*y1(iz))/detaabb
          else
            aa(iz) = 0.0
            bb(iz) = 0.0
          endif
          if(myid.eq.0) print*,"isilce: ",i,detaabb,aa(iz),bb(iz),a22(iz)
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(x1(2,nset*numpts0))
        allocate(x2(2,nset*numpts0))
 
        call normVec(x1,nset*numpts0)
        call normVec(x2,nset*numpts0)

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                r = (2.0*r-1.0)*0.125*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              !Gaussian distribution
              this%Pts1(1,ii) = xmu1 + sig1*x1(1,ii)/sq12
              this%Pts1(2,ii) = xmu2 + sig2*(-muxpx*x1(1,ii)/sq12+x1(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                r = (2.0*r-1.0)*0.125*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              !transverse Gaussian distribution
              this%Pts1(3,ii) = xmu3 + sig3*x2(1,ii)/sq34
              this%Pts1(4,ii) = xmu4 + sig4*(-muypy*x2(1,ii)/sq34+x2(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                !r = (2.0*r-1.0)*0.02*zmax 
                !r = (2.0*r-1.0)*hz 
                r = (2.0*r-1.0)*zsize 
              endif
              ii = j+nset*(i-1)
              !this%Pts1(5,ii) = Ptcl(6,i)+r
              this%Pts1(5,ii) = r/2 !uniform longiudinal distribution from -hz/2 to hz/2
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               !r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0) !15keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*2.0d0 !2keV is rms
              endif
              ii = j+nset*(i-1)
              iz = (this%Pts1(5,ii)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              !this%Pts1(6,ii) = Ptcl(7,i)+r
              this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        !print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)

        !print*,"Nplocal: ",this%Nptlocal

        zmin = zmin - hz
        zmax = zmax + hz
        zbunch = zmax-zmin
        rkk = 2*asin(1.0d0)*2/zbunch
        rkk = 2*2*asin(1.0d0)*2/zbunch
        rkk = 4*2*asin(1.0d0)*2/zbunch
        rkk = 512*2*asin(1.0d0)*2/zbunch
        rkk = 256*2*asin(1.0d0)*2/zbunch
        ampz = 0.2d0
        ampz = 0.0d0
        ampz = 0.2d0
        ampz = 0.0d0
        ampz = 0.2d0

        if(myid.eq.0) then
          print*,"zmin: ",zmin,zmax,zbunch,sig1,sig3,sq12,sq34
        endif
        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
         
          zz = this%Pts1(5,j)-zmin
          sig1z = sig1*(1.0+ampz*sin(rkk*zz))
          sig3z = sig3*(1.0+ampz*sin(rkk*zz))
          this%Pts1(1,j) = xmu1 + sig1z*x1(1,j)/sq12
          this%Pts1(3,j) = xmu3 + sig3z*x2(1,j)/sq34

          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        deallocate(x1)
        deallocate(x2)

        this%Npt = nset*inipts

        end subroutine readElegantNcorMod3_Dist

        !longitudinal orginal particle distribution
        subroutine readElegantNcorMod4_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        !integer, parameter :: nslice = 100 !used in the 2nd 1B
        !integer, parameter :: nslice = 50 !used in the 3rd 1B
        integer, parameter :: nslice = 25
        integer :: iz
        real*8, dimension(nslice) :: a11lc,a11,a12lc,a12,a22lc,a22,y1lc,&
                                     y1,y2lc,y2,aa,bb,rk,gam1,zz1
        real*8 :: hz,zmin,detaabb
        double precision, allocatable, dimension(:,:) :: x1,x2
        double precision  :: sigx,sigpx,muxpx,sigy,&
        sigpy,muypy,sigz,sigpz,muzpz
        real*8 :: sig1,sig2,sig3,sig4,sig5,sig6,sq12,sq34,sq56
        real*8 :: zz,zbunch,rkk,sig1z,sig3z,ampz
        integer :: iizz
        real*8 :: ziizz1,gamz1,ziizz1b,gamz1b
        real*8 :: ziizz2,gamz2,ziizz2b,gamz2b

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale
 
        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          !print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        zmin = 1.0e10
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.tmptcl(6)) zmax = tmptcl(6) 
          if(zmin.ge.tmptcl(6)) zmin = tmptcl(6) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)

          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        !print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        zmin = zmin*Scfreq*2*Pi*zscale
        !avoid index overflow
        zmin = zmin - (zmax-zmin)*0.5e-5
        zmax = zmax + (zmax-zmin)*0.5e-5
        hz = (zmax-zmin)/nslice

        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        !print*,"xmax: ",xmax,ymax,zmin,zmax,pxmax,pymax,pzmax
        a11lc = 0.0
        a12lc = 0.0
        a22lc = 0.0
        y1lc = 0.0
        y2lc = 0.0
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
          iz = (Ptcl(6,j)-zmin)/hz + 1
          a11lc(iz) = a11lc(iz) + Ptcl(6,j)**2
          a12lc(iz) = a12lc(iz) + Ptcl(6,j)
          a22lc(iz) = a22lc(iz) + 1.0
          y1lc(iz) = y1lc(iz) + Ptcl(7,j)
          y2lc(iz) = y2lc(iz) + Ptcl(6,j)*Ptcl(7,j)
        enddo

        call MPI_ALLREDUCE(a11lc,a11,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a12lc,a12,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a22lc,a22,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y1lc,y1,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y2lc,y2,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        !print*,"a22: ",sum(a22)
        do iz = 1, nslice
          if(a22(iz).lt.1.0e-8) print*,"a22: ",iz,a22(iz)
        enddo
        do iz = 1, nslice
          if(a22(iz).gt.0.0) then
            detaabb = a22(iz)*a11(iz)-a12(iz)**2
            aa(iz) = (a11(iz)*y1(iz)-a12(iz)*y2(iz))/detaabb
            bb(iz) = (a22(iz)*y2(iz)-a12(iz)*y1(iz))/detaabb
          else
            aa(iz) = 0.0
            bb(iz) = 0.0
          endif
          if(myid.eq.0) print*,"isilce: ",i,detaabb,aa(iz),bb(iz),a22(iz)
        enddo

      iizz = nslice/2
      ziizz1 = zmin+(iizz-1)*hz
      gamz1 = aa(iizz)+bb(iizz)*ziizz1
      zz1(iizz) = ziizz1
      gam1(iizz) = gamz1
      rk(iizz) = bb(iizz)
      do i = iizz-1,1,-1
        rk(i) = (y2(i)-ziizz1*y1(i)+ziizz1*gamz1*a22(i)-gamz1*a12(i))/ &
                (a11(i)+ziizz1*ziizz1*a22(i)-2*ziizz1*a12(i))
        ziizz1b = zmin+(i-1)*hz
        gamz1b = gamz1 + rk(i)*(ziizz1b-ziizz1)
        gamz1 = gamz1b
        ziizz1 = ziizz1b
        zz1(i) = ziizz1
        gam1(i) = gamz1
      enddo
      ziizz2 = zmin+iizz*hz
      gamz2 = aa(iizz)+bb(iizz)*ziizz2
      do i = iizz+1,nslice
        rk(i) = (y2(i)-ziizz2*y1(i)+ziizz2*gamz2*a22(i)-gamz2*a12(i))/ &
                (a11(i)+ziizz2*ziizz2*a22(i)-2*ziizz2*a12(i))
        ziizz2b = zmin+i*hz
        gamz2b = gamz2 + rk(i)*(ziizz2b-ziizz2)
        zz1(i) = ziizz2
        gam1(i) = gamz2
        gamz2 = gamz2b
        ziizz2 = ziizz2b
      enddo


        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(x1(2,nset*numpts0))
        allocate(x2(2,nset*numpts0))
 
        call normVec(x1,nset*numpts0)
        call normVec(x2,nset*numpts0)

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                r = (2.0*r-1.0)*0.125*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              !Gaussian distribution
              this%Pts1(1,ii) = xmu1 + sig1*x1(1,ii)/sq12
              this%Pts1(2,ii) = xmu2 + sig2*(-muxpx*x1(1,ii)/sq12+x1(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                r = (2.0*r-1.0)*0.125*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              !transverse Gaussian distribution
              this%Pts1(3,ii) = xmu3 + sig3*x2(1,ii)/sq34
              this%Pts1(4,ii) = xmu4 + sig4*(-muypy*x2(1,ii)/sq34+x2(2,ii))

              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                !r = (2.0*r-1.0)*0.02*zmax 
                r = (2.0*r-1.0)*hz 
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(6,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               !r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0) !15keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*2.0d0 !2keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*4.0d0 !4keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*7.5d0 !7.5keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*15.0d0 !7.5keV is rms
               r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/15*7.0d0 !7.5keV is rms
              endif
              ii = j+nset*(i-1)
              iz = (this%Pts1(5,ii)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              !this%Pts1(6,ii) = Ptcl(7,i)+r
              !this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+r
              this%Pts1(6,ii) = gam1(iz)+rk(iz)*&
                                (this%Pts1(5,ii)-zz1(iz))+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        !print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)

        !print*,"Nplocal: ",this%Nptlocal

        zmin = zmin - hz
        zmax = zmax + hz
        zbunch = zmax-zmin
        rkk = 2*asin(1.0d0)*2/zbunch
        rkk = 2*2*asin(1.0d0)*2/zbunch
        rkk = 4*2*asin(1.0d0)*2/zbunch
        rkk = 512*2*asin(1.0d0)*2/zbunch
        !rkk = 1024*2*asin(1.0d0)*2/zbunch
        !rkk = 256*2*asin(1.0d0)*2/zbunch
        ampz = 0.2d0
        ampz = 0.0d0
        ampz = 0.2d0
        ampz = 0.0d0

        if(myid.eq.0) then
          print*,"zmin: ",zmin,zmax,zbunch,sig1,sig3,sq12,sq34
        endif
        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
         
          zz = this%Pts1(5,j)-zmin
          sig1z = sig1*(1.0+ampz*sin(rkk*zz))
          sig3z = sig3*(1.0+ampz*sin(rkk*zz))
          this%Pts1(1,j) = xmu1 + sig1z*x1(1,j)/sq12
          this%Pts1(3,j) = xmu3 + sig3z*x2(1,j)/sq34

          pid = j + jlow
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo

        deallocate(x1)
        deallocate(x2)

        this%Npt = nset*inipts

        end subroutine readElegantNcorMod4_Dist

subroutine bcnrand (n, a, x)

!   This routine generates a sequence of IEEE 64-bit floating-point pseudorandom
!   numbers in the range (0, 1), based on the recently discovered class of 
!   normal numbers described in the paper "Random Generators and Normal Numbers"
!   by DHB and Richard Crandall, available at http://www.nersc.gov/~dhbailey.

!   The sequence generated is (x_k), where x_k = first 53 binary digits of the
!   binary expansion of alpha_{2,3} beginning at binary position 53*(k-1) + a,
!   normalized by 2^53.  The parameter a is the starting point of the sequence,
!   or in other words the seed of the pseudorandom sequence.  To obtain the 
!   maximum period, a should be set to at least 3^33 = 5.5590606e15.  The value
!   of a must not exceed 2^53 = 9.0071992e15.  When a is set in this range, the
!   period of the generator is 2x3^32 = 3.7060404e15.  In general, the period
!   of the generator is 2x3^(p-1), where 3^p is the largest power of 3 less 
!   than a.  However the generated sequence will not match the binary digits of
!   alpha_{2,3} if the range of sequence indices, namely [a, a+53*n], includes
!   a power of 3.

!   The bcnrand routine is designed for simple parallelization.  For example,
!   in an MPI program, suppose that kpr is the processor number and npr is the
!   number of processors.  Then the line

!   call bcnrand (n/npr, a + 53*n/npr*kpr, x)

!   generates on each processor a section of length n/npr.  In this way, the
!   npr processors collectively have the same n-long sequence (provided that
!   n is divisible by npr) as is generated on a single processor system by
!   means of the line

!   call bcnrand (n, a, x)

!   On IBM and Apple system use the compiler flags -O3 -qstrict, or results 
!   will not be right.  Also, on these systems and others with a fused
!   multiply-add instruction, see !> comments in subroutines ddmuldd and dddivd
!   below for a simple code change that significantly improves performance.

!   David H. Bailey    2002-08-10

implicit none
integer i, ib, j, k, n
real*8 a, aa, d1, d2, d3, dd1(2), dd2(2), dd3(2), p2(53), p3(34), &
  p3i, t53, x(n)
!external expm2

d1 = 1.d0

do i = 1, 53
  d1 = 2.d0 * d1
  p2(i) = d1
enddo

t53 = p2(53)

!   Check input parameters.

if (a > t53) then
  write (6, 1) a
1 format ('bcnrand: error - second argument exceeds 2^53'/ &
  'value =',1p,d25.15)
  goto 200
endif

d1 = 1.d0

do i = 1, 34
  d1 = 3.d0 * d1
  p3(i) = d1
  if (p3(i) > a) goto 100
enddo

100 continue

p3i = p3(i-1)

if (n > 2.d0 * p3i / 3.d0) then
  write (6, 2) n, 2.d0 * p3i / 3.d0
2 format ('bcnrand:  warning - number of elements exceeds period'/ &
  'n, period =',i12,f18.0)
endif

!   Calculate starting element.

d2 = expm2 (a - p3i, p2, p3i)
d3 = aint (0.5d0 * p3i)
call ddmuldd (d2, d3, dd1)
call dddivd (dd1, p3i, dd2)
d1 = aint (dd2(1))
call ddmuldd (d1, p3i, dd2)
call ddsub (dd1, dd2, dd3)
d1 = dd3(1)
x(1) = d1 / p3i

!   Calculate successive members of sequence.

do i = 2, n
  dd1(1) = t53 * d1
  dd1(2) = 0.d0
  call dddivd (dd1, p3i, dd2)
  d2 = aint (dd2(1))
  call ddmuldd (p3i, d2, dd2)
  call ddsub (dd1, dd2, dd3)
  d1 = dd3(1)
  if (d1 < 0.d0) d1 = d1 + p3i
  x(i) = d1 / p3i
enddo

200 continue

return
end subroutine bcnrand

function expm2 (p, p2, am)

!   expm2 = 2^p mod am.  p2 is a table with powers of 2, i.e., p2(i) = 2^i.
!   This routine uses a left-to-right binary exponentiation scheme.

implicit none
integer i
real*8 am, d1, d2, dd1(2), dd2(2), dd3(2), ddm(2), expm2, p, p1, p2(53), pt1, r

do i = 1, 53
  if (p2(i) > p) goto 100
enddo

100 continue

p1 = p
pt1 = p2(i-1)
r = 1.d0
ddm(1) = am
ddm(2) = 0.d0

110 continue

if (p1 .ge. pt1) then
! r = mod (2.d0 * r, am)
  call ddmuldd (2.d0, r, dd1)
  if (dd1(1) > am) then
    call ddsub (dd1, ddm, dd2)
    dd1(1) = dd2(1)
    dd1(2) = dd2(2)
  endif
  r = dd1(1)
  p1 = p1 - pt1
endif
pt1 = 0.5d0 * pt1
if (pt1 .ge. 1.d0) then
!  r = mod (r * r, am)
  call ddmuldd (r, r, dd1)
  call dddivd (dd1, am, dd2)
  d2 = aint (dd2(1))
  call ddmuldd (am, d2, dd2)
  call ddsub (dd1, dd2, dd3)
  r = dd3(1)
  if (r < 0.d0) r = r + am
  goto 110
endif

expm2 = r

return
end function expm2

subroutine ddmuldd (da, db, ddc)

!   This subroutine computes ddc = da x db.

implicit none
real*8 a1, a2, b1, b2, cona, conb, da, db, ddc(2), split, s1, s2
parameter (split = 134217729.d0)

!>
!   On systems with a fused multiply add, such as IBM systems, it is faster to
!   uncomment the next two lines and comment out the following lines until !>.
!   On other systems, do the opposite.

s1 = da * db
s2 = da * db - s1

!   This splits da and db into high-order and low-order words.

! cona = da * split
! conb = db * split
! a1 = cona - (cona - da)
! b1 = conb - (conb - db)
! a2 = da - a1
! b2 = db - b1
! s1 = da * db
! s2 = (((a1 * b1 - s1) + a1 * b2) + a2 * b1) + a2 * b2
!>
ddc(1) = s1
ddc(2) = s2

return
end subroutine

subroutine dddivd (dda, db, ddc)

!   This routine divides the DD number A by the DP number B to yield
!   the DD quotient C.  

implicit none
real*8 dda(2), db, ddc(2)
real*8 a1, a2, b1, b2, cona, conb, e, split, t1, t2, t11, t12, t21, t22
parameter (split = 134217729.d0)

!   Compute a DP approximation to the quotient.

t1 = dda(1) / db
!>
!   On systems with a fused multiply add, such as IBM systems, it is faster to
!   uncomment the next two lines and comment out the following lines until !>.
!   On other systems, do the opposite.

t12 = t1 * db
t22 = t1 * db - t12

!   This splits t1 and db into high-order and low-order words.

! cona = t1 * split
! conb = db * split
! a1 = cona - (cona - t1)
! b1 = conb - (conb - db)
! a2 = t1 - a1
! b2 = db - b1
! t12 = t1 * db
! t22 = (((a1 * b1 - t12) + a1 * b2) + a2 * b1) + a2 * b2
!>
!   Compute dda - (t12, t22) using Knuth's trick.

t11 = dda(1) - t12
e = t11 - dda(1)
t21 = ((-t12 - e) + (dda(1) - (t11 - e))) + dda(2) - t22

!   Compute high-order word of (t11, t21) and divide by db.

t2 = (t11 + t21) / db

!   The result is t1 + t2, after normalization.

ddc(1) = t1 + t2
ddc(2) = t2 - (ddc(1) - t1)
return
end subroutine

subroutine ddsub (dda, ddb, ddc)

!   This subroutine computes ddc = dda - ddb.

implicit none
real*8 dda(2), ddb(2), ddc(2)
real*8 e, t1, t2

!   Compute dda + ddb using Knuth's trick.

t1 = dda(1) - ddb(1)
e = t1 - dda(1)
t2 = ((-ddb(1) - e) + (dda(1) - (t1 - e))) + dda(2) - ddb(2)

!   The result is t1 + t2, after normalization.

ddc(1) = t1 + t2
ddc(2) = t2 - (ddc(1) - t1)
return
end subroutine

        !use the formulae from Sasha for longitudinal
        !phase space sampling (including 2nd and 3rd )
        subroutine Gauss3ldrd_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp
!        integer seedarray(1)
        double precision :: t0,x11,zz
! for testing April 4 distribution. 
        real*8 :: r1,r2,r4,x1tmp,twopi,pid

        call starttime_Timer(t0)

        twopi = 4*asin(1.0d0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)
        print*,myid,x11,this%Npt

        avgpts = this%Npt/(npx*npy)

        !if(mod(avgpts,10).ne.0) then
        !  print*,"The number of particles has to be an integer multiple of 10Nprocs"
        !  stop
        !endif
        
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        sig5 = sigz*zscale
        sig6 = sigpz*pzscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

!        allocate(x1(2,avgpts))
!        allocate(x2(2,avgpts))
!        allocate(x3(2,avgpts))
!        call normVec(x1,avgpts)
!        call normVec(x2,avgpts)
!        call normVec(x3,avgpts)
        
        !intvsamp = 10
        intvsamp = 1
        allocate(x1(2,intvsamp))
        allocate(x2(2,intvsamp))
        allocate(x3(2,intvsamp))

        do j = 1, avgpts/intvsamp
          call normVec(x1,intvsamp)
          call normVec(x2,intvsamp)
          call normVec(x3,intvsamp)
          do k = 1, intvsamp
            !x-px:
!            call normdv(x1)
!           Correct Gaussian distribution.
            i = (j-1)*intvsamp + k
            this%Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
            this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
            !y-py
!            call normdv(x1)
!           Correct Gaussian distribution.
            this%Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
            this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
            !z-pz
!            call normdv(x1)
!           Correct Gaussian distribution.
!comment out for testing April 4 distribution
!            call random_number(zz)
!            zz = 2*zz-1.0
!            this%Pts1(5,i) = xmu5 + sig5*zz*sqrt(3.0)
!            !convert into ps unit
!            zz = this%Pts1(5,i)*1000/(4*asin(1.0d0)*1.3)
!            this%Pts1(6,i) = xmu6 + sig6*x3(2,k)-&
!             (13.975*zz-1.79146*zz**2+1.43235*zz**3) 

!for testing April 4 distribution------------
            call random_number(r1)
            call random_number(r2)
            r4 = sqrt(r1)
            r2 = r2*twopi
            x1tmp = 2*r4*cos(r2)
            this%Pts1(5,i) = xmu5 + sig5*x1tmp/sq12
! 05/06/08 Gaussian current profile
!            this%Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq12
            !convert into ps unit
            zz = this%Pts1(5,i)*1000/(4*asin(1.0d0)*1.3)
!            this%Pts1(6,i) = xmu6+sig6*x3(2,k)+0.0154033*zz-&
!                             0.000330052*zz*zz
!for testing April 28 distribution..
            this%Pts1(6,i) = xmu6+sig6*x3(2,k)+0.00210534*zz+&
                             0.000190623*zz*zz
!for test 6/25/08
!            this%Pts1(6,i) = xmu6
!---------------------------------------------

          enddo
        enddo
          
        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts

        print*,"avgpts: ",this%Nptlocal
        do j = 1, avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = j + myid*avgpts
        enddo
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss3ldrd_Dist

        !read in an initial distribution and do repopulation.
        !Here, the T-Pt is fitted by a piecewise linear function
        !to capture the global correlated energy spread. The uncorrelated
        !energy spread can be added to the top of original energy spread.
        !By this way, the original T-dPt correlation is included. 
        subroutine readElegantCor_Dist(this,nparam,distparam,geom,grid,Flagbc,nslice)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc,nslice
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sumx2,sumy2,xmax,pxmax,ymax,pymax,zmax,pzmax
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(7) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz
        !integer, parameter :: nslice = 100 !used in the 2nd 1B
        !integer, parameter :: nslice = 50 !used in the 3rd 1B
        !integer, parameter :: nslice = 25
        integer :: iz
        real*8, dimension(nslice) :: a11lc,a11,a12lc,a12,a22lc,a22,y1lc,&
                                     y1,y2lc,y2,aa,bb,rk,gam1,zz1
        real*8 :: hz,zmin,detaabb
        integer :: iizz 
        real*8 :: ziizz1,gamz1,ziizz1b,gamz1b
        real*8 :: ziizz2,gamz2,ziizz2b,gamz2b
        !half box size of repopulation
        real*8 :: xrp,pxrp,yrp,pyrp,zrp,pzrp
        real*8 :: dee

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        xrp = distparam(2)
        pxrp = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yrp = distparam(9)
        pyrp = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        pzrp = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          !print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(7,avgpts))

        xmax = 0.0
        pxmax = 0.0
        ymax = 0.0
        pymax = 0.0
        zmax = 0.0
        pzmax = 0.0
        sumx2 = 0.0
        sumy2 = 0.0
        sumx = 0.0
        sumy = 0.0
        zmin = 1.0e10
        do j = 1, inipts
          read(12,*)tmptcl(1:7)
          if(xmax.le.abs(tmptcl(2))) xmax = abs(tmptcl(2)) 
          if(ymax.le.abs(tmptcl(4))) ymax = abs(tmptcl(4)) 
          if(zmax.le.tmptcl(6)) zmax = tmptcl(6) 
          if(zmin.ge.tmptcl(6)) zmin = tmptcl(6) 
          gammabet = tmptcl(7)/sqrt(1.d0+tmptcl(3)**2+tmptcl(5)**2)
          tmppx = tmptcl(3)*gammabet
          if(pxmax.le.abs(tmppx)) pxmax = abs(tmppx)
          tmppy = tmptcl(5)*gammabet
          if(pymax.le.abs(tmppy)) pymax = abs(tmppy)
          tmppz = gamma0 - sqrt(1.0+tmptcl(7)**2)
          if(pzmax.le.abs(tmppz)) pzmax = abs(tmppz)

          sumx = sumx + tmptcl(2)
          sumx2 = sumx2 + tmptcl(2)*tmptcl(2)
          sumy = sumy + tmptcl(4)
          sumy2 = sumy2 + tmptcl(4)*tmptcl(4)

          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:7,i) = tmptcl(1:7)
          endif
        enddo
        sumx2 = sqrt(sumx2/inipts-(sumx/inipts)**2)
        sumy2 = sqrt(sumy2/inipts-(sumy/inipts)**2)
        !print*,"sumx2: ",sumx2,sumy2
        close(12)
        !call MPI_BARRIER(comm2d,ierr)

        xmax = xmax/xl*xscale
        ymax = ymax/xl*yscale
        zmax = zmax*Scfreq*2*Pi*zscale
        zmin = zmin*Scfreq*2*Pi*zscale
        !avoid index overflow
        zmin = zmin - (zmax-zmin)*0.5e-5
        zmax = zmax + (zmax-zmin)*0.5e-5
        hz = (zmax-zmin)/nslice

        pxmax = pxmax*pxscale
        pymax = pymax*pyscale
        pzmax = pzmax*pzscale
        !print*,"xmax: ",xmax,ymax,zmin,zmax,pxmax,pymax,pzmax
        a11lc = 0.0
        a12lc = 0.0
        a22lc = 0.0
        y1lc = 0.0
        y2lc = 0.0
        do j = 1, avgpts
          Ptcl(2,j) = (Ptcl(2,j)/xl)*xscale + xmu1
          Ptcl(4,j) = (Ptcl(4,j)/xl)*yscale + xmu3
          Ptcl(6,j) = Ptcl(6,j)*Scfreq*2*Pi*zscale + xmu5
          gammabet = Ptcl(7,j)/sqrt(1.d0+Ptcl(3,j)**2+Ptcl(5,j)**2)
          Ptcl(3,j) = Ptcl(3,j)*gammabet*pxscale + xmu2
          Ptcl(5,j) = Ptcl(5,j)*gammabet*pyscale + xmu4
          Ptcl(7,j) = (gamma0 - sqrt(1.0 + Ptcl(7,j)**2))*pzscale + xmu6 
!          if(xmax.le.abs(Ptcl(1,j))) xmax = abs(Ptcl(1,j)) 
!          if(pxmax.le.abs(Ptcl(2,j))) pxmax = abs(Ptcl(2,j)) 
!          if(ymax.le.abs(Ptcl(3,j))) ymax = abs(Ptcl(3,j)) 
!          if(pymax.le.abs(Ptcl(4,j))) pymax = abs(Ptcl(4,j)) 
!          if(zmax.le.abs(Ptcl(5,j))) zmax = abs(Ptcl(5,j)) 
!          if(pzmax.le.abs(Ptcl(6,j))) pzmax = abs(Ptcl(6,j)) 
          iz = (Ptcl(6,j)-zmin)/hz + 1
          a11lc(iz) = a11lc(iz) + Ptcl(6,j)**2
          a12lc(iz) = a12lc(iz) + Ptcl(6,j)
          a22lc(iz) = a22lc(iz) + 1.0
          y1lc(iz) = y1lc(iz) + Ptcl(7,j)
          y2lc(iz) = y2lc(iz) + Ptcl(6,j)*Ptcl(7,j)
        enddo

        call MPI_ALLREDUCE(a11lc,a11,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a12lc,a12,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(a22lc,a22,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y1lc,y1,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(y2lc,y2,nslice,MPI_DOUBLE_PRECISION,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        !print*,"a22: ",sum(a22)
        do iz = 1, nslice
          if(a22(iz).lt.1.0e-8) print*,"a22: ",iz,a22(iz)
        enddo
        do iz = 1, nslice
          if(a22(iz).gt.0.0) then
            detaabb = a22(iz)*a11(iz)-a12(iz)**2
            aa(iz) = (a11(iz)*y1(iz)-a12(iz)*y2(iz))/detaabb
            bb(iz) = (a22(iz)*y2(iz)-a12(iz)*y1(iz))/detaabb
          else
            aa(iz) = 0.0
            bb(iz) = 0.0
          endif
          if(myid.eq.0) print*,"isilce: ",i,detaabb,aa(iz),bb(iz),a22(iz)
        enddo

!-----------------------------------------------
!the following lines are used for the piecewise continuous of the fitting of E-phi
      iizz = nslice/2
      ziizz1 = zmin+(iizz-1)*hz
      gamz1 = aa(iizz)+bb(iizz)*ziizz1
      zz1(iizz) = ziizz1
      gam1(iizz) = gamz1
      rk(iizz) = bb(iizz)
      do i = iizz-1,1,-1
        rk(i) = (y2(i)-ziizz1*y1(i)+ziizz1*gamz1*a22(i)-gamz1*a12(i))/ &
                (a11(i)+ziizz1*ziizz1*a22(i)-2*ziizz1*a12(i))
        ziizz1b = zmin+(i-1)*hz
        gamz1b = gamz1 + rk(i)*(ziizz1b-ziizz1)
        gamz1 = gamz1b
        ziizz1 = ziizz1b
        zz1(i) = ziizz1
        gam1(i) = gamz1
      enddo
      ziizz2 = zmin+iizz*hz
      gamz2 = aa(iizz)+bb(iizz)*ziizz2
      do i = iizz+1,nslice
        rk(i) = (y2(i)-ziizz2*y1(i)+ziizz2*gamz2*a22(i)-gamz2*a12(i))/ &
                (a11(i)+ziizz2*ziizz2*a22(i)-2*ziizz2*a12(i))
        ziizz2b = zmin+i*hz
        gamz2b = gamz2 + rk(i)*(ziizz2b-ziizz2)
        zz1(i) = ziizz2
        gam1(i) = gamz2
        gamz2 = gamz2b
        ziizz2 = ziizz2b
      enddo
!-----------------------------------------------

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                !r = (2.0*r-1.0)*0.03*xmax
                !r = (2.0*r-1.0)*0.14*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*1.0*xmax
                !r = (2.0*r-1.0)*0.5*xmax
                !r = (2.0*r-1.0)*0.1*xmax
                !r = (2.0*r-1.0)*0.03*xmax !used in the 1st 1B 
                !r = (2.0*r-1.0)*0.125*xmax
                r = (2.0*r-1.0)*xrp*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(2,i)+r
              !this%Pts1(1,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                !r = (2.0*r-1.0)*0.03*pxmax
                ! r = (2.0*r-1.0)*0.14*pxmax
                ! r = (2.0*r-1.0)*0.0*pxmax
                ! r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*1.0*pxmax
                 !r = (2.0*r-1.0)*0.5*pxmax
                 !r = (2.0*r-1.0)*0.1*pxmax
                 !r = (2.0*r-1.0)*0.03*pxmax
                 r = (2.0*r-1.0)*pxrp*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(3,i)+r
              !this%Pts1(2,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                !r = (2.0*r-1.0)*0.03*ymax
                !r = (2.0*r-1.0)*0.14*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*1.0*ymax
                !r = (2.0*r-1.0)*0.5*ymax
                !r = (2.0*r-1.0)*0.1*ymax
                !r = (2.0*r-1.0)*0.03*ymax !used in the 1st 1B
                !r = (2.0*r-1.0)*0.125*ymax
                r = (2.0*r-1.0)*yrp*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(4,i)+r
              !this%Pts1(3,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                !r = (2.0*r-1.0)*0.14*pymax
                !r = (2.0*r-1.0)*0.0*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*1.0*pymax
                !r = (2.0*r-1.0)*0.5*pymax
                !r = (2.0*r-1.0)*0.1*pymax
                !r = (2.0*r-1.0)*0.03*pymax
                r = (2.0*r-1.0)*pyrp*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(5,i)+r
              !this%Pts1(4,ii) = r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else  
              !  r = (2.0*r-1.0)*0.0080*zmax  
                 !r = (2.0*r-1.0)*0.0159*zmax
                 !r = (2.0*r-1.0)*0.0*zmax
                 !r = (2.0*r-1.0)*0.05*zmax
                !r = (2.0*r-1.0)*0.004*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                !r = (2.0*r-1.0)*0.24*zmax
                !r = (2.0*r-1.0)*0.001*zmax !used in the 1st 1B simulation
                !r = (2.0*r-1.0)*0.02*zmax 
                r = (2.0*r-1.0)*hz 
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(6,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                ! r = (2.0*r-1.0)*0.0159*pzmax
                ! r = (2.0*r-1.0)*0.0*pzmax
                !r = (2.0*r-1.0)*0.01*pzmax
               !r = (2.0*r-1.0)*0.24*pzmax
               !r = (2.0*r-1.0)*0.001*pzmax
               !r = (2.0*r-1.0)*0.005*pzmax
               !r = (2.0*r-1.0)*0.017*pzmax
               !r = (2.0*r-1.0)*0.015*pzmax
               !r = (2.0*r-1.0)*0.016*pzmax !used in 1st 1B simulation
               !r = (2.0*r-1.0)*0.014*pzmax
               !r = (2.0*r-1.0)*0.010*pzmax
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0) !15keV is rms
               !r = (2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 0.0*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/2 !7.5keV is rms
               !r = 2*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/15 !2keV is rms
               !r = 5*(2.0*r-1.0)*0.029354*sqrt(3.0d0)/15 !5keV is rms
                r = pzrp*(2.0*r-1.0)*sqrt(3.0d0) !5keV is rms
              endif
              !ii = j+nset*(i-1)
              !iz = (this%Pts1(5,ii)-zmin)/hz + 1
              !if(iz.le.0) iz = 1
              !if(iz.gt.nslice) iz = nslice
              !if(nset.eq.1) then
              !  this%Pts1(6,ii) = Ptcl(7,i)+r
              !else
              !  !this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+r
              !  this%Pts1(6,ii) = gam1(iz)+rk(iz)*&
              !                  (this%Pts1(5,ii)-zz1(iz))+r
              !endif

              ii = j+nset*(i-1)

              iz = (Ptcl(6,i)-zmin)/hz + 1
              if(iz.le.0) iz = 1
              if(iz.gt.nslice) iz = nslice
              !dee = Ptcl(7,i)-(aa(iz)+bb(iz)*Ptcl(6,i))
              dee = Ptcl(7,i)-(gam1(iz)+rk(iz)*(Ptcl(6,i)-zz1(iz)))
 
              if(nset.eq.1) then
                this%Pts1(6,ii) = Ptcl(7,i)+r
              else
                !use dee will help to keep the Pz-z correlation in the
                !original distribution. Hence, the uncorrelated energy
                !spread is no longer fully determined by r. The function
                !of r is actually like a box. 
                iz = (this%Pts1(5,ii)-zmin)/hz + 1
                if(iz.le.0) iz = 1
                if(iz.gt.nslice) iz = nslice
 
                !if r=0, there is no repopulation in dE
                !if dee = 0, uniform resampling in dE
                !this%Pts1(6,ii) = aa(iz)+bb(iz)*this%Pts1(5,ii)+dee+r
                this%Pts1(6,ii) = gam1(iz)+rk(iz)*&
                                (this%Pts1(5,ii)-zz1(iz))+dee+r
              endif 
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        !print*,"sumxnew: ",sumx/this%Nptlocal, this%Pts1(6,1)  

        deallocate(Ptcl)

        !print*,"Nplocal: ",this%Nptlocal

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = j + jlow
        enddo

        this%Npt = nset*inipts

        end subroutine readElegantCor_Dist

        !read particle distribution from ImpactT output
        subroutine readimpt_Dist(this,nparam,distparam,geom,grid,Flagbc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,Flagbc
        double precision, dimension(nparam) :: distparam
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer*8 :: nptot
        integer :: ierr,i,j,k,ii,nset,nremain,numpts0
        integer :: myid,myidy,myidx,totnp,npy,npx,comm2d,commcol,&
                   commrow,inipts,pid
        double precision, dimension(6) :: lcrange,a
        double precision, allocatable, dimension(:,:) :: Ptcl
        double precision :: r,xl,gamma0,gamma,synangle,betaz
        double precision :: sum1,sum2
!        integer seedarray(1)
        double precision :: xx,mccq,kenergy,gammabet
        double precision :: xscale,xmu1,xmu2,yscale,xmu3,xmu4,zscale,&
        xmu5,xmu6,sumx,sumy,pxscale,pyscale,pzscale,ri,thi,pi
        integer :: jlow,jhigh,avgpts,nleft
        real*8, dimension(6) :: tmptcl
        real*8 :: tmppx,tmppy,tmppz,phi0lc
        real*8 :: xmax,ymax,zmax,pxmax,pymax,pzmax

        xl = Scxl
        mccq = this%Mass
        !gamma0 = 1.0+kenergy/mccq      
        phi0lc = this%refptcl(5)
        gamma0 = -this%refptcl(6)
        pi = 2*asin(1.0)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        nptot = this%Npt

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        call getsize_Pgrid2d(grid,totnp,npy,npx)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

!        seedarray(1)=(1021+myid)*(myid+7)
!        call random_seed(PUT=seedarray(1:1))
        call random_number(xx)
        write(6,*)myid,xx

        call getlcrange_CompDom(geom,lcrange)

        open(unit=12,file='partcl.data',status='old')

        read(12,*)inipts

          avgpts = inipts/totnp
          nleft = inipts - avgpts*totnp
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          print*,"avgpts: ",avgpts,myid,inipts,totnp
          allocate(Ptcl(6,avgpts))

        xmax = 0.0
        ymax = 0.0
        zmax = 0.0
        pxmax = 0.0
        pymax = 0.0
        pzmax = 0.0

        sum1 = 0.0d0
        sum2 = 0.0d0
        do j = 1, inipts
          read(12,*)tmptcl(1:6)
          sum1 = sum1 + tmptcl(5)
          gamma = sqrt(1.0+tmptcl(2)**2+tmptcl(4)**2+tmptcl(6)**2)
          sum2 = sum2 + gamma
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              Ptcl(1:6,i) = tmptcl(1:6)
          endif
        enddo
        sum1 = sum1/inipts
        sum2 = sum2/inipts
        print*,"mean z location and energy: ",sum1, sum2
        close(12)
        call MPI_BARRIER(comm2d,ierr)
        this%refptcl(6) = -sum2 

        do j = 1, avgpts
          Ptcl(1,j) = (Ptcl(1,j)/xl)*xscale + xmu1
          Ptcl(3,j) = (Ptcl(3,j)/xl)*yscale + xmu3
          Ptcl(5,j) = ((Ptcl(5,j)-sum1)/(xl*sqrt(1.d0-1.d0/sum2**2)))*zscale + xmu5
          gamma = sqrt(1.0+Ptcl(2,j)**2+Ptcl(4,j)**2+Ptcl(6,j)**2)
          Ptcl(6,j) = (sum2-gamma)*pzscale + xmu6
          Ptcl(2,j) = Ptcl(2,j)*pxscale + xmu2
          Ptcl(4,j) = Ptcl(4,j)*pyscale + xmu4
        enddo

        numpts0 = avgpts 

        nset = nptot/inipts
        nremain = nptot - nset*inipts
        if(nremain.ne.0) then
          if(myid.eq.0) then
            print*,"The total number of particle is not ",nptot
          endif
        endif

        this%Nptlocal = nset*numpts0

        allocate(this%Pts1(9,nset*numpts0))
        sumx = 0.0
        do i = 1, numpts0
          do j = 1, nset
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*xmax
                r = (2.0*r-1.0)*0.03*xmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(1,ii) = Ptcl(1,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.015*pxmax
                r = (2.0*r-1.0)*0.03*pxmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(2,ii) = Ptcl(2,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*ymax
                r = (2.0*r-1.0)*0.03*ymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(3,ii) = Ptcl(3,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.02*pymax
                r = (2.0*r-1.0)*0.03*pymax
              endif
              ii = j+nset*(i-1)
              this%Pts1(4,ii) = Ptcl(4,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.005*zmax
                !r = (2.0*r-1.0)*0.01*zmax
                r = (2.0*r-1.0)*0.01*zmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(5,ii) = Ptcl(5,i)+r
              call random_number(r)
              if(nset.eq.1) then
                r = 0.0
              else
                !r = (2.0*r-1.0)*0.002*pzmax
                !r = (2.0*r-1.0)*0.004*pzmax
                r = (2.0*r-1.0)*0.01*pzmax
              endif
              ii = j+nset*(i-1)
              this%Pts1(6,ii) = Ptcl(6,i)+r
          enddo
          sumx = sumx + this%Pts1(1,i)
        enddo

        print*,"sumxnew: ",sumx/this%Nptlocal

        deallocate(Ptcl)

        print*,"Nplocal: ",this%Nptlocal

        jlow = (jlow-1)*nset
        do j = 1, this%Nptlocal
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = j + jlow
        enddo

        end subroutine readimpt_Dist

        subroutine Gauss3ldrd4_Dist(this,nparam,distparam,grid,flagalloc)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam,flagalloc
        double precision, dimension(nparam) :: distparam
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3 
        integer :: totnp,npy,npx
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp,pid
!        integer seedarray(1)
        double precision :: t0,x11,zz
! for testing April 4 distribution.
        real*8 :: r1,r2,r4,x1tmp,twopi
        real*8 :: r,fvalue,zlcmin,zlcmax,ampmod,kk
        real*8 :: x33,cs,ss

        call starttime_Timer(t0)

        twopi = 4*asin(1.0d0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)
        print*,myid,x11

        avgpts = this%Npt/(npx*npy)

        !if(mod(avgpts,10).ne.0) then
        !  print*,"The number of particles has to be an integer multiple of 10Nprocs"
        !  stop
        !endif
        
        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale
        !sig5 = sigz*zscale
        !sig6 = sigpz*pzscale
        sig5 = sigz
        sig6 = sigpz

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        if(flagalloc.eq.1) then
          this%Pts1 = 0.0
        else
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
        endif

!        allocate(x1(2,avgpts))
!        allocate(x2(2,avgpts))
!        allocate(x3(2,avgpts))
!        call normVec(x1,avgpts)
!        call normVec(x2,avgpts)
!        call normVec(x3,avgpts)
        
        intvsamp = 1
        allocate(x1(2,intvsamp))
        allocate(x2(2,intvsamp))
        allocate(x3(2,intvsamp))

        zlcmin = -sig5*sqrt(3.0d0)
        zlcmax = sig5*sqrt(3.0d0)
        ampmod = zscale
        kk = pzscale

        i = 0
        k = 1
        do
          ! rejection sample.
10        call random_number(r)
          r1 = zlcmin + r*(zlcmax-zlcmin)
          fvalue = (1.0d0+ampmod*sin(kk*r1))/(1.0d0+ampmod)
          call random_number(r2)
          if(r2.gt.fvalue) goto 10
          zz = r1
          i = i + 1
          if(i.gt.avgpts) exit

!uniform transverse cylinder 6/21/08
          call random_number(r)
          x11 = sig1*sqrt(r)
          x33 = sig3*sqrt(r)
          call random_number(r)
          cs = cos(r*twopi)
          ss = sin(r*twopi)
          this%Pts1(1,i) = xmu1 + x11*cs
          this%Pts1(2,i) = xmu2 
          this%Pts1(3,i) = xmu3 + x33*ss
          this%Pts1(4,i) = xmu4 

!          call normVec(x1,intvsamp)
!          call normVec(x2,intvsamp)
!          call normVec(x3,intvsamp)
!          !x-px:
!!         call normdv(x1)
!!         Correct Gaussian distribution.
!          this%Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
!          this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
!          !y-py
!!         call normdv(x1)
!!         Correct Gaussian distribution.
!          this%Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
!          this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))

          this%Pts1(5,i) = xmu5 + zz
          this%Pts1(6,i) = xmu6+sig6*x3(2,k)
!          !this linear chirp before chicane is for benchmark purpose 
!          this%Pts1(6,i) = xmu6+sig6*x3(2,k)+0.55188-118.614*zz
        enddo


!comment out for testing current modulization 6/9/08
!        do j = 1, avgpts/intvsamp
!          call normVec(x1,intvsamp)
!          call normVec(x2,intvsamp)
!          call normVec(x3,intvsamp)
!          do k = 1, intvsamp
!            !x-px:
!!            call normdv(x1)
!!           Correct Gaussian distribution.
!            i = (j-1)*intvsamp + k
!            this%Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
!            this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
!            !y-py
!!            call normdv(x1)
!!           Correct Gaussian distribution.
!            this%Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
!            this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
!            !z-pz
!!            call normdv(x1)
!!           Correct Gaussian distribution.
!!comment out for testing April 4 distribution
!            call random_number(zz)
!            zz = 2*zz-1.0
!            this%Pts1(5,i) = xmu5 + sig5*zz*sqrt(3.0)
!!            !convert into ps unit
!!            zz = this%Pts1(5,i)*1000/(4*asin(1.0d0)*1.3)
!!            this%Pts1(6,i) = xmu6 + sig6*x3(2,k)-&
!!             (13.975*zz-1.79146*zz**2+1.43235*zz**3) 
!
!!for testing April 4 distribution------------
!!            call random_number(r1)
!!            call random_number(r2)
!!            r4 = sqrt(r1)
!!            r2 = r2*twopi
!!            x1tmp = 2*r4*cos(r2)
!!            this%Pts1(5,i) = xmu5 + sig5*x1tmp/sq12
!! 05/06/08 Gaussian current profile
!!            this%Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq12
!            !convert into ps unit
!            zz = this%Pts1(5,i)*1000/(4*asin(1.0d0)*1.3)
!!            this%Pts1(6,i) = xmu6+sig6*x3(2,k)+0.0154033*zz-&
!!                             0.000330052*zz*zz
!!for testing April 28 distribution..
!            this%Pts1(6,i) = xmu6+sig6*x3(2,k)+0.00210534*zz+&
!                             0.000190623*zz*zz
!!for 6/6 test
!            zz = this%Pts1(5,i)
!            this%Pts1(6,i) = xmu6+sig6*x3(2,k)+0.55188-118.614*zz
!!---------------------------------------------
!
!          enddo
!        enddo
          
        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts

        do j = 1, avgpts
          pid = j + myid*avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = pid
        enddo
       
        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss3ldrd4_Dist

        subroutine Gauss3dSoblcls_Dist(this,nparam,distparam,grid)
        implicit none 
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this 
        integer, intent(in) :: nparam 
        double precision, dimension(nparam) :: distparam 
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        double precision :: sig1,sig2,sig3,sig4,sig5,sig6
        double precision :: sq12,sq34,sq56
        double precision, allocatable, dimension(:,:) :: x1,x2,x3 
        integer :: totnp,npy,npx,npt
        integer :: avgpts,numpts
        integer :: myid,myidx,myidy,i,j,k,intvsamp
!        integer seedarray(1)
        double precision :: t0,x11
        integer :: ierr,ist,ilow,ihigh,nleft
        real*8 :: rkz,engamp

        call starttime_Timer(t0)

        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6)
        xmu2 = distparam(7)
        sigy = distparam(8)
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11)
        pyscale = distparam(12)
        xmu3 = distparam(13)
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        !here zscale is used as modulation wavelength for lcls study
        zscale = distparam(18)
        !here pzscale is used as modulation amplitude for lcls study
        pzscale = distparam(19)
        xmu5 = distparam(20)
        xmu6 = distparam(21)

        !normalized modulation frequency
        rkz = Clight/zscale/Scfreq
        !energy modulation amplitude
        engamp = pzscale

        call getsize_Pgrid2d(grid,totnp,npy,npx)

        call getpost_Pgrid2d(grid,myid,myidy,myidx)
!        seedarray(1)=(1001+myid)*(myid+7)
!        write(6,*)'seedarray=',seedarray
!        call random_seed(PUT=seedarray(1:1))
        call random_number(x11)

        avgpts = this%Npt/totnp
        nleft = this%Npt - avgpts*totnp
        if(myid.lt.nleft) then
            avgpts = avgpts+1
            ilow = myid*avgpts + 1
            ihigh = (myid+1)*avgpts
        else
            ilow = myid*avgpts + 1 + nleft
            ihigh = (myid+1)*avgpts + nleft
        endif
        ist = 0

        sig1 = sigx*xscale
        sig2 = sigpx*pxscale
        sig3 = sigy*yscale
        sig4 = sigpy*pyscale

        sq12=sqrt(1.-muxpx*muxpx)
        sq34=sqrt(1.-muypy*muypy)
        sq56=sqrt(1.-muzpz*muzpz)

        ! initial allocate 'avgpts' particles on each processor.
        !if(flagalloc.eq.1) then
        !  Pts1 = 0.0
        !else
        !  allocate(Pts1(6,avgpts))
        !  Pts1 = 0.0
        !endif
        allocate(this%Pts1(6,avgpts))

        !intvsamp = 10
        intvsamp = avgpts
        allocate(x1(2,Npt))
        allocate(x2(2,Npt))
        allocate(x3(2,Npt))

!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        do j = 1, 1
          call normVecSob(x1,x2,x3,Npt,ist)
!          print*,"pass Sob: ",sum(x1),sum(x2),sum(x3),ilow,ihigh,myidy,myid
          do k = 1, Npt
            if(k.gt.ilow .and. k.le.ihigh) then
              i = k-ilow + 1
              !x-px:
              this%Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
              this%Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
              !y-py
              this%Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
              this%Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
              !z-pz
              this%Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq56
              this%Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,k)/sq56+x3(2,k)) + &
                               engamp*sin(rkz*this%Pts1(5,i))
            endif
          enddo
        enddo

        deallocate(x1)
        deallocate(x2)
        deallocate(x3)

        this%Nptlocal = avgpts
 
        do j = 1, avgpts
          this%Pts1(7,j) = this%Charge/this%mass
          this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
          this%Pts1(9,j) = ilow + j - 1
        enddo 

        t_kvdist = t_kvdist + elapsedtime_Timer(t0)

        end subroutine Gauss3dSoblcls_Dist

        subroutine readMLI_Dist(this,nparam,distparam,grid)
        implicit none
        include 'mpif.h'
        type (BeamBunch), intent(inout) :: this
        integer, intent(in) :: nparam 
        double precision, dimension(nparam) :: distparam 
        type (Pgrid2d), intent(in) :: grid
        double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
        sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
        double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
        integer :: i,j,jlow,jhigh,avgpts,myid,nproc,ierr,nleft
        integer*8 :: nptot
        double precision, dimension(6) :: tmptcl
        double precision :: sum1,sum2
        integer*8 :: pid
        integer :: nset
 
        sigx = distparam(1)
        sigpx = distparam(2)
        muxpx = distparam(3)
        xscale = distparam(4)
        pxscale = distparam(5)
        xmu1 = distparam(6) 
        xmu2 = distparam(7) 
        sigy = distparam(8) 
        sigpy = distparam(9)
        muypy = distparam(10)
        yscale = distparam(11) 
        pyscale = distparam(12)
        xmu3 = distparam(13) 
        xmu4 = distparam(14)
        sigz = distparam(15)
        sigpz = distparam(16)
        muzpz = distparam(17)
        zscale = distparam(18)
        pzscale = distparam(19)
        xmu5 = distparam(20) 
        xmu6 = distparam(21) 

        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
 
        open(unit=12,file='partcl.data',status='old')
 
        sum1 = 0.0
        sum2 = 0.0
          read(12,*)nptot
          avgpts = nptot/nproc
          nleft = nptot - avgpts*nproc
          if(myid.lt.nleft) then
            avgpts = avgpts+1
            jlow = myid*avgpts + 1
            jhigh = (myid+1)*avgpts
          else
            jlow = myid*avgpts + 1 + nleft
            jhigh = (myid+1)*avgpts + nleft
          endif
          allocate(this%Pts1(9,avgpts))
          this%Pts1 = 0.0
          !jlow = myid*avgpts + 1
          !jhigh = (myid+1)*avgpts
          print*,"avgpts, jlow, and jhigh: ",avgpts,jlow,jhigh
          do j = 1, nptot
            read(12,*)tmptcl(1:6)
            sum1 = sum1 + tmptcl(1)
            sum2 = sum2 + tmptcl(3)
            if( (j.ge.jlow).and.(j.le.jhigh) ) then
              i = j - jlow + 1
              this%Pts1(1:6,i) = tmptcl(1:6)
            endif
          enddo
          print*,"sumx1,sumy1: ",sum1/nptot,sum2/nptot
 
          close(12)
 
          this%Nptlocal = avgpts

!          nset = 1
!          jlow = (jlow-1)*nset
          do j = 1, this%Nptlocal
            pid = j + jlow-1
            this%Pts1(1,j) = (this%Pts1(1,j)*1.0d0/Scxl)*xscale+xmu1
            this%Pts1(2,j) = this%Pts1(2,j)*pxscale+xmu2
            this%Pts1(3,j) = (this%Pts1(3,j)*1.0d0/Scxl)*yscale+xmu3
            this%Pts1(4,j) = this%Pts1(4,j)*pyscale+xmu4
            this%Pts1(5,j) = this%Pts1(5,j)*zscale+xmu5
            this%Pts1(6,j) = (this%Pts1(6,j)*0.838338d0)*pzscale+xmu6
            this%Pts1(7,j) = this%Charge/this%mass
            this%Pts1(8,j) = this%Current/Scfreq/this%Npt*this%Charge/abs(this%Charge)
            this%Pts1(9,j) = pid
          enddo

        end subroutine readMLI_Dist

        subroutine normVecSob(y1,y2,y3,num,iseed)
        implicit none
        include 'mpif.h'
        integer :: num
        double precision, dimension(2,num) :: y1,y2,y3
        double precision :: twopi,epsilon
        integer :: i,ndim,iseed
        real*8, dimension(6) :: xtmp 

        epsilon = 1.0d-18
        ndim = 6

        twopi = 4.0d0*dasin(1.0d0)
        !call sobseqV(-2,xtmp)
        !iseed = 303
        call sobseqV(-2,xtmp,iseed)
        !print*,"sobseq: ",0,xtmp
        do i = 1, num
          call sobseqV(ndim,xtmp,iseed)
          !print*,"sobseq: ",i,xtmp
          if(xtmp(1).eq.0.0) xtmp(1) = epsilon
          if(xtmp(3).eq.0.0) xtmp(3) = epsilon
          if(xtmp(5).eq.0.0) xtmp(5) = epsilon
          y1(1,i) = sqrt(-2.0*log(xtmp(1)))*cos(twopi*xtmp(2))
          y1(2,i) = sqrt(-2.0*log(xtmp(1)))*sin(twopi*xtmp(2))
          y2(1,i) = sqrt(-2.0*log(xtmp(3)))*cos(twopi*xtmp(4))
          y2(2,i) = sqrt(-2.0*log(xtmp(3)))*sin(twopi*xtmp(4))
          y3(1,i) = sqrt(-2.0*log(xtmp(5)))*cos(twopi*xtmp(6))
          y3(2,i) = sqrt(-2.0*log(xtmp(5)))*sin(twopi*xtmp(6))
          !write(2,*)y1(:,i),y2(:,i),y3(:,i)
        enddo

        end subroutine normVecSob

      SUBROUTINE sobseqV(n,x,iseed)
      INTEGER n,MAXBIT,MAXDIM,iseed
      REAL*8 x(*)
      PARAMETER (MAXBIT=30,MAXDIM=6)
      INTEGER i,im,in,ipp,j,k,l,ip(MAXDIM),iu(MAXDIM,MAXBIT),&
      iv(MAXBIT*MAXDIM),ix(MAXDIM),mdeg(MAXDIM)
      REAL*8 fac
      SAVE ip,mdeg,ix,iv,in,fac
      EQUIVALENCE (iv,iu)
      DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
      DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
      if (n.lt.0) then
        do 14 k=1,MAXDIM
          do 11 j=1,mdeg(k)
            iu(k,j)=iu(k,j)*2**(MAXBIT-j)
11        continue
          do 13 j=mdeg(k)+1,MAXBIT
            ipp=ip(k)
            i=iu(k,j-mdeg(k))
            i=ieor(i,i/2**mdeg(k))
            do 12 l=mdeg(k)-1,1,-1
              if(iand(ipp,1).ne.0)i=ieor(i,iu(k,j-l))
              ipp=ipp/2
12          continue
            iu(k,j)=i
13        continue
14      continue
        fac=1./2.**MAXBIT
        in=iseed
      else
        im=in
        do 15 j=1,MAXBIT
          if(iand(im,1).eq.0)goto 1
          im=im/2
15      continue
        pause 'MAXBIT too small in sobseq'
1       im=(j-1)*MAXDIM
        do 16 k=1,min(n,MAXDIM)
          ix(k)=ieor(ix(k),iv(im+k))
          x(k)=ix(k)*fac
16      continue
        in=in+1
      endif
      return
      END subroutine sobseqV

!<<<<<<<<<<<<<<   begin : distIOTA member routines(Kilean) <<<<<<<<<<<<<
  subroutine distIOTA_init(self,distParam)
    implicit none
    class(distIOTA_class) :: self
    double precision, intent(in) :: distParam(2)
    integer,parameter :: t_=1,c_=2
    self%t = distParam(t_)
    self%c = distParam(c_)
  end subroutine distIOTA_init
    
  double precision function distIOTA_getH(self, xHat)
    implicit none
    class(distIOTA_class) :: self
    double precision, intent(in):: xHat(4)
    double complex  :: zeta
    integer, parameter :: x_=1,y_=3
    
    zeta = dcmplx(xHat(x_),xHat(y_))/self%c
    distIOTA_getH = 0.5d0*sum(xHat*xHat) + self%c*self%c*self%t*&
                    real(zeta/croot(zeta)*carcsin(zeta))
    return
  end function distIOTA_getH
  
  subroutine distIOTA_genP_waterbag(self,BB,beta,betap,emittanceCut)
    !generate particle in synergia unit
    implicit none
    include 'mpif.h'
    class(distIOTA_class) :: self
    type(BeamBunch),intent(inout) :: BB
    double precision, intent(in) :: beta,betap,emittanceCut
    integer :: nproc,avgpts,nleft,myid,i
    double precision  :: xHat(4),xMax,U(2),newH,bg,trialH,newHxy,p,p_theta
    integer, parameter :: x_=1,px_=2,y_=3,py_=4
    
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,i)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,i)
    avgpts = BB%Npt/nproc
    nleft = BB%Npt - avgpts*nproc
    if(myid.lt.nleft) then
      avgpts = avgpts+1
    endif
    !if(allocated(BB%Pts1)) deallocate(BB%Pts1)
    allocate(BB%Pts1(1:9,avgpts))
    BB%Pts1 = 0d0
    BB%Nptlocal = avgpts
    bg = sqrt(BB%refptcl(6)**2 - 1d0)
    do i=1,avgpts
      xHat = 0d0
      do
        call random_number(p_theta)
        if(p_theta>0) exit
      enddo
      newH = emittanceCut*sqrt(p_theta)
      xMax = self%secant_method(sqrt(newH), newH) 
      do
        call random_number(U)
        xHat(x_) = 2d0*(0.5d0-U(1))*xMax
        xHat(y_) = 2d0*(0.5d0-U(2))*xMax
        newHxy = self%getH(xHat)
        if(newHxy < newH) exit
      enddo
      p = sqrt(2d0*(newH-newHxy))
      call random_number(p_theta)
      p_theta = 2.0*pi*p_theta
      xHat(px_) = p*cos(p_theta)
      xHat(py_) = p*sin(p_theta)
      BB%Pts1(x_,i) = xHat(x_)*sqrt(beta)/Scxl
      BB%Pts1(y_,i) = xHat(y_)*sqrt(beta)/Scxl
      BB%Pts1(px_,i) = (xHat(px_) + 0.5d0*betap*xHat(x_))/sqrt(beta)*bg
      BB%Pts1(py_,i) = (xHat(py_) + 0.5d0*betap*xHat(y_))/sqrt(beta)*bg
      BB%Pts1(9,i) = i
    enddo
    BB%Pts1(7,:) = BB%Charge/BB%mass
    BB%Pts1(8,:) = BB%Current/Scfreq/BB%Npt*BB%Charge/abs(BB%Charge)
    if(myid.lt.nleft) then
      BB%Pts1(9,i) = BB%Pts1(9,i) + myid*avgpts
    else
      BB%Pts1(9,i) = BB%Pts1(9,i) + myid*avgpts + nleft
    endif
  end subroutine distIOTA_genP_waterbag
  
  subroutine distIOTA_waterbag(this,nparam,distparam)
    implicit none
    type (BeamBunch), intent(inout) :: this
    integer, intent(in) :: nparam
    double precision, dimension(nparam) :: distparam
    type (distIOTA_class) :: IOTA
    integer,parameter :: beta_=3,betap_=4,emittance_=5
    
    call IOTA%init(distParam(1:2))
    call IOTA%genP_waterbag(this,distparam(beta_),distparam(betap_),distparam(emittance_))
  end subroutine distIOTA_waterbag
  
  double precision function distIOTA_secant_method(self,x0,emittance)
    !secant method to find a maximum x0 of given emittance
    implicit none
    class(distIOTA_class) :: self
    double precision, intent(in) :: x0,emittance
    double precision, parameter  :: tol=1.48d-8
    integer,parameter :: maxiter = 50
    integer :: i
    double precision :: p0,p1,p,q0,q1
        
    p0 = x0
    if (x0>=0) then
      p1 = x0*1.0001d0 + 1d-4
    else
      p1 = x0*1.0001d0 - 1d-4
    endif
    q0 = emittance - self%getH([0d0,0d0,p0,0d0]) !func eval
    q1 = emittance - self%getH([0d0,0d0,p1,0d0])
    do i=1,maxiter
      if (q1==q0) then
        if (p1/=p0) then
          print*, 'Msg from [Distributionclass::distIOTA_secant_method] -> toloerance reached'
          distIOTA_secant_method = 0.5d0*(p1+p0)
          return
        endif
      else
        p = p1 - q1*(p1 - p0)/(q1 - q0)
      endif
      if(abs(p - p1) < tol) then
        distIOTA_secant_method = p
        return
      endif
      p0 = p1
      q0 = q1
      p1 = p
      q1 = emittance - self%getH([0d0,0d0,p1,0d0])
    enddo
    print*, 'Msg from [Distributionclass::distIOTA_secant_method] -> Failed to converge after maxiterations'
    distIOTA_secant_method = p
  end function distIOTA_secant_method
  
  subroutine distIOTA_genP_gaussian(self,BB,beta,betap,emittance,cutoff)
    !generate particle in synergia unit
    implicit none
    include 'mpif.h'
    class(distIOTA_class) :: self
    type(BeamBunch),intent(inout) :: BB
    double precision, intent(in) :: beta,betap,emittance,cutoff
    integer :: nproc,avgpts,nleft,myid,i
    double precision  :: xHat(4),xMax,U(2),newH,bg,trialH,newHxy,p,p_theta
    integer, parameter :: x_=1,px_=2,y_=3,py_=4
    
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,i)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,i)
    avgpts = BB%Npt/nproc
    nleft = BB%Npt - avgpts*nproc
    if(myid.lt.nleft) then
      avgpts = avgpts+1
    endif
    !if(allocated(BB%Pts1)) deallocate(BB%Pts1)
    allocate(BB%Pts1(1:9,avgpts))
    BB%Pts1 = 0d0
    BB%Nptlocal = avgpts
    bg = sqrt(BB%refptcl(6)**2 - 1d0)
    do i=1,avgpts
      xHat = 0d0
      do
        do
          call random_number(U)
          if(U(1)*U(2)>0) exit
        enddo
        trialH = -log(U(1)*U(2))
        if (trialH < cutoff) exit
      enddo
      newH = emittance*trialH
      xMax = self%secant_method(sqrt(newH), newH) 
      do
        call random_number(U)
        xHat(x_) = 2d0*(0.5d0-U(1))*xMax
        xHat(y_) = 2d0*(0.5d0-U(2))*xMax
        newHxy = self%getH(xHat)
        if(newHxy < newH) exit
      enddo
      p = sqrt(2d0*(newH-newHxy))
      call random_number(p_theta)
      p_theta = 2.0*pi*p_theta
      xHat(px_) = p*cos(p_theta)
      xHat(py_) = p*sin(p_theta)
      BB%Pts1(x_,i) = xHat(x_)*sqrt(beta)/Scxl
      BB%Pts1(y_,i) = xHat(y_)*sqrt(beta)/Scxl
      BB%Pts1(px_,i) = (xHat(px_) + 0.5d0*betap*xHat(x_))/sqrt(beta)*bg
      BB%Pts1(py_,i) = (xHat(py_) + 0.5d0*betap*xHat(y_))/sqrt(beta)*bg
      BB%Pts1(9,i) = i
    enddo
    BB%Pts1(7,:) = BB%Charge/BB%mass
    BB%Pts1(8,:) = BB%Current/Scfreq/BB%Npt*BB%Charge/abs(BB%Charge)
    if(myid.lt.nleft) then
      BB%Pts1(9,i) = BB%Pts1(9,i) + myid*avgpts
    else
      BB%Pts1(9,i) = BB%Pts1(9,i) + myid*avgpts + nleft
    endif
  end subroutine distIOTA_genP_gaussian

  subroutine distIOTA_gaussian(this,nparam,distparam)
    implicit none
    type (BeamBunch), intent(inout) :: this
    integer, intent(in) :: nparam
    double precision, dimension(nparam) :: distparam
    type (distIOTA_class) :: IOTA
    integer,parameter :: beta_=3,betap_=4,emittance_=5,cutoff_=6
    
    call IOTA%init(distParam(1:2))
    call IOTA%genP_gaussian(this,distparam(beta_),distparam(betap_),distparam(emittance_),distparam(cutoff_))
  end subroutine distIOTA_gaussian
  
  !>>>>>>>>>>>>>>> end : distIOTA member routines(Kilean) >>>>>>>>>>>>>>>
      end module Distributionclass
