!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! must be included in order of dependence
#include "HMC_MODULES/modulesRHMC.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program hmc
      use options
      use gaugemodule
      use condmod
      use gammas
      use paulimodule
      use indices
      implicit none
      integer j
      double precision num

      print *,"Start program Dirac HMC"

      if (.false.) then
        print *,"Start tests"
        theta=0
        GAUGETYPE=2 ! 1=compact 2=non-compact link field
        MDW=one ! domain wall height
        call setGammas
        call setPauliMatrices
        call setIndices
        call coef(u,theta)
        call test()
        print *,"End tests"
        stop
      end if

      call init()
      do j=1,Naux
        print *,"make field ",j," of ",Naux
        call makeGaugeField(.false.)
        call measureCondensateInstance()
      end do
!      call  averageCondensate();

      call finalise()
      end program hmc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine init()
      use timer
      use countmod
      use rvmodule
      use options
      use gammas
      use paulimodule
      use indices
      use gaugefield
      use zolomodule
      use dwcoeffs
      use condmod
      implicit none
      real(prc) tav
      type(zolotarev) :: zolo
      real(prc) lmin,lmax

      call initRVs(.false.,.false.,0) ! uses a time based seed initialiser

      oc_idx=0
      outer_count=0
      ic_idx=0
      inner_count=0

      call init_timer()

      HMCtype=2 ! 1:HMC-DW, 2: HMC-OL, 3:RHMC-DW

      Nferms=1 ! hmmm ... take this out, beta controls this too
      GAUGETYPE=2 ! 1=compact 2=non-compact link field
      MDW=one ! domain wall height

      MTYPE=3 ! mass term type
      DWkernel=1 ! 1 for Shamir, 2 for Wilson, 3 for OWilson
      OLTYPE=1 ! 1 is direct calculation, 2 is for DW (K)-type
      baremass=0.02d0
      gbeta=1.0

      HMC_etime=0.2 ! average total time for HMC trajectory
      HMC_dt=0.02  ! time step size for HMC trajectory
      HMC_tsmax=1000 ! maximum no of timesteps before forcing end of HMC tranjectory
      QUENCHED=.false.

      Naux=10 ! no of auxiliary fields to be generated
      Nswp=5 ! how many trajectories for each each field

      open(unit=11,file="caseInput.dat",status='unknown')
      read(11,*) HMC_dt,gbeta,MDW,baremass,MTYPE,tsav,DWkernel,QUENCHED
      close(11)
      HMC_etime=tsav*HMC_dt
      open(unit=11,file="runtimeInput.dat",status='unknown')
      read(11,*) Naux,Nswp
      close(11)

      write(102,*) "Ns:",Ns,"Nt:",Nt,"Ls:",Ls
      write(102,*) "dwkernel:",DWkernel
      write(102,*) "mtype:",MTYPE
      write(102,*) "gbeta:",gbeta
      write(102,*) "mass:",baremass
      write(102,*) "MD length:",HMC_etime,"Time Step:",HMC_dt
      write(102,*) "Max steps:",HMC_tsmax
      write(102,*) "QUENCHED:",QUENCHED

      call setGammas
      call setPauliMatrices
#ifdef TWODIMENSIONS
      call setIndices2d
#else
      call setIndices
#endif

      theta=0
!      call setGRVs(3*Nv,theta)
      call setRealRVs(3*Nv,theta)
      theta=theta-0.5d0
      write(150,*) theta
      call coef(u,theta)

      omega=1.0
      if (dwkernel.eq.2) then
        call setHTcoeffs(Ls,SRF)
      else if (dwkernel.eq.3) then
        lmin=5d-2 ; lmax=5d0
        call setZoloCoeffs(Ls,SRF,lmin,lmax)
        call setZolo(lmin,lmax,Ls,zolo)
        call getRoots(zolo)
        omega=one/zolo%roots
        print *,"omega:",omega
      end if

      return
      end subroutine init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine finalise()
      implicit none

       return
       end subroutine finalise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine test()
      use testDequivmod
      implicit none

       call convKDDW4()
!       call equivOLmtype()
!       call equivDOL_DDW()

       return
       end subroutine test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

