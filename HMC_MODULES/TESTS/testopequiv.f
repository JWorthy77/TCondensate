      module testDequivmod
      use pacc
      use arraysizes
      use numbers
      use options
      use gammas
c      use basicdiracopsmod
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testOperatorEquivalence()
      implicit none

      print *,"testOperatorEquivalence"
      call equivDOL_DDW()

      return
      end subroutine testOperatorEquivalence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine convKDDW4()
      use overlap
      use rvmodule
      use gaugefield
      use ratfuncs
      use dwcoeffs
      use zolomodule
      use domainwallmod
      implicit none
      type(sgnratfunc) :: SRF
      complex(prc) RR(Nv,4)
      complex(prc),dimension(Nv,4) ::  HM1test,HM3test
      complex(prc),dimension(Nv,4) ::  ZM1test,ZM3test
      real ev(1)
      real(prc) lmin,lmax,lmass
      type(zolotarev) :: zolo

      call setRVs(Nv*4,RR)
      lmass=one/20
      lmin=1d-2 ; lmax=5d0
      dwkernel=3 

      omega=1.0

      MTYPE=1
      call KDDW4(RR,HM1test,u,.false.,lmass)
      MTYPE=3
      call KDDW4(RR,HM3test,u,.false.,lmass)
      if (dwkernel.eq.3) then
        call setZolo(5d-2,5d0,Ls,zolo)
        call getRoots(zolo)
        omega=one/zolo%roots
        print *,"omega:",omega
      end if
      MTYPE=1
      call KDDW4(RR,ZM1test,u,.false.,lmass)
      MTYPE=3
      call KDDW4(RR,ZM3test,u,.false.,lmass)

       print *,maxval(abs(HM1test)),
     &           maxval(abs(HM3test))
       print *,maxval(abs(ZM1test)),
     &           maxval(abs(ZM3test))
!      print *,"No reason why they should be identical"
!      print *,"Condensate should be identical"

      return
      end subroutine convKDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine equivOLmtype()
      use overlap
      use rvmodule
      use gaugefield
      use ratfuncs
!      use spectra
!      use kernelspectrarange
      implicit none
!     use 12x12 mesh, and create quenched gauge field
      logical ESTEIGS
      character(len=80) :: fname
      type(sgnratfunc) :: SRF
      complex(prc) R(Nv,4)
      complex(prc),dimension(Nv,4,100) ::  HM1test,HM3test,HM4test
      complex(prc),dimension(Nv,4,100) ::  HM5test,HM6test
      complex(prc),dimension(Nv,4,100) ::  ZM1test,ZM3test,ZM4test
      complex(prc),dimension(Nv,4,100) ::  ZM5test,ZM6test
      real ev(1)
      real(prc) lmin,lmax,lmass
      integer j,jmin,jmax,shft

      call setRVs(Nv*4,R)
      lmass=one/20
      jmin=6 ; jmax=20 ; shft=2
      lmin=1d-2 ; lmax=5d0
      dwkernel=1 ! no implementation of Shamir OL here!
      do j=jmin,jmax,shft
        call setHTcoeffs(j,SRF)
        MTYPE=1
        call DOLMW(R,HM1test(:,:,j),u,.false.,lmass,SRF)
        MTYPE=3
        call DOLMW(R,HM3test(:,:,j),u,.false.,lmass,SRF)
        MTYPE=4
        call DOLMW(R,HM4test(:,:,j),u,.false.,lmass,SRF)
        MTYPE=5
        call DOLMW(R,HM5test(:,:,j),u,.false.,lmass,SRF)
        call setZoloCoeffs(j,SRF,lmin,lmax)
        MTYPE=1
        call DOLMW(R,ZM1test(:,:,j),u,.false.,lmass,SRF)
        MTYPE=3
        call DOLMW(R,ZM3test(:,:,j),u,.false.,lmass,SRF)
        MTYPE=4
        call DOLMW(R,ZM4test(:,:,j),u,.false.,lmass,SRF)
        MTYPE=5
        call DOLMW(R,ZM5test(:,:,j),u,.false.,lmass,SRF)

       print *,j,maxval(abs(HM1test(:,:,j))),
     &           maxval(abs(HM3test(:,:,j))),
     &           maxval(abs(HM4test(:,:,j))),
     &           maxval(abs(HM5test(:,:,j)))
       print *,j,maxval(abs(ZM1test(:,:,j))),
     &           maxval(abs(ZM3test(:,:,j))),
     &           maxval(abs(ZM4test(:,:,j))),
     &           maxval(abs(ZM5test(:,:,j)))
      end do
      print *,"No reason why they should be identical"
      print *,"Condensate should be identical"

      return
      end subroutine equivOLmtype
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine equivDOL_DDW()
      use overlap
      use dwcoeffs
      use zolomodule
      use domainwallmod
      use rvmodule
      use gaugefield
      use ratfuncs
      use kernelspectrarange
      implicit none
      complex(prc),dimension(Nv,4) :: RR,DR1,DR2,TMP
      real(prc) :: err
      type(sgnratfunc) :: SRF
      type(zolotarev) :: zolo
      integer i1,i2
      real(prc) lmin,lmax
      real(prc),dimension(Nv,3) :: thetastore
      complex(prc) :: R5(Nv,4,Ls)
      complex(prc) :: DR5(Nv,4,Ls)
      complex(prc) :: TMP5(Nv,4,Ls)

!      call setHTcoeffs(Ls,SRF)
      call getMinMaxHEigs(lmin,lmax)
      lmin=5e-2 ; lmax=10.0
      call setZoloCoeffs(Ls,SRF,lmin,lmax)
      call setZolo(lmin,lmax,Ls,zolo)
      call getRoots(zolo)
      omega=one/zolo%roots
      print *,"roots:",zolo%roots
      print *,"omega:",omega

      call setRVs(Nv*4,RR)
      call setGRVs(Nv*3,thetastore)
!      theta=0
      theta=thetastore
      call coef(u,theta)

      print *,""
      print *,"TEST INVERSION OF OVERLAP WILSON FORMULATIONS"
      print *,""

      baremass=0.05
      MTYPE=1
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"DAGGER:    false"
      call DOLMW(RR,DR1,u,.false.,baremass,SRF)
      call IDOLMW(DR1,DR2,u,.false.,baremass,SRF)
      print *,"DOL.IDOL:",maxval(abs(DR2-RR))
      call IDOLMW(RR,DR1,u,.false.,baremass,SRF)
      call DOLMW(DR1,DR2,u,.false.,baremass,SRF)
      print *,"IDOL.DOL:",maxval(abs(DR2-RR))

      baremass=0.05
      MTYPE=1
      print *,""
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"DAGGER:    true"
!      print *,"dwkernel:",dwkernel
      call DOLMW(RR,DR1,u,.true.,baremass,SRF)
      call IDOLMW(DR1,DR2,u,.true.,baremass,SRF)
      print *,"DOL.IDOL:",maxval(abs(DR2-RR))
      call IDOLMW(RR,DR1,u,.true.,baremass,SRF)
      call DOLMW(DR1,DR2,u,.true.,baremass,SRF)
      print *,"IDOL.DOL:",maxval(abs(DR2-RR))

      baremass=0.05
      MTYPE=3
      print *,""
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"DAGGER:    false"
      call DOLMW(RR,DR1,u,.false.,baremass,SRF)
      call IDOLMW(DR1,DR2,u,.false.,baremass,SRF)
      print *,"DOL.IDOL:",maxval(abs(DR2-RR))
      call IDOLMW(RR,DR1,u,.false.,baremass,SRF)
      call DOLMW(DR1,DR2,u,.false.,baremass,SRF)
      print *,"IDOL.DOL:",maxval(abs(DR2-RR))

      baremass=0.05
      MTYPE=3
      print *,""
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"DAGGER:    true"
!      print *,"dwkernel:",dwkernel
      call DOLMW(RR,DR1,u,.true.,baremass,SRF)
      call IDOLMW(DR1,DR2,u,.true.,baremass,SRF)
      print *,"DOL.IDOL:",maxval(abs(DR2-RR))
      call IDOLMW(RR,DR1,u,.true.,baremass,SRF)
      call DOLMW(DR1,DR2,u,.true.,baremass,SRF)
      print *,"IDOL.DOL:",maxval(abs(DR2-RR))


      print *,""
      print *,""
      print *,"TEST INVERSION OF KDDW4 WILSON FORMULATIONS"
      print *,""

      baremass=0.05d0
      MTYPE=1
      dwkernel=2
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"dwkernel:",dwkernel
      print *,"DAGGER:    false"
      call KDDW4(RR,DR1,u,.false.,baremass)
      call IKDDW4(DR1,DR2,u,.false.,baremass)
      print *,"KDDW4.IKDDW4:",maxval(abs(DR2-RR))
      call IKDDW4(RR,DR1,u,.false.,baremass)
      call KDDW4(DR1,DR2,u,.false.,baremass)
      print *,"IDOL.DOL:",maxval(abs(DR2-RR))

      baremass=0.05d0
      MTYPE=1
      dwkernel=2
      print *,""
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"dwkernel:",dwkernel
      print *,"DAGGER:    true"
      call KDDW4(RR,DR1,u,.true.,baremass)
      call IKDDW4(DR1,DR2,u,.true.,baremass)
      print *,"KDDW4.IKDDW4:",maxval(abs(DR2-RR))
      call IKDDW4(RR,DR1,u,.true.,baremass)
      call KDDW4(DR1,DR2,u,.true.,baremass)
      print *,"IKDDW4.KDDW4:",maxval(abs(DR2-RR))

      baremass=0.05d0
      MTYPE=1
      dwkernel=3
      print *,""
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"dwkernel:",dwkernel
      print *,"DAGGER:    false"
      call KDDW4(RR,DR1,u,.false.,baremass)
      call IKDDW4(DR1,DR2,u,.false.,baremass)
      print *,"KDDW4.IKDDW4:",maxval(abs(DR2-RR))
      call IKDDW4(RR,DR1,u,.false.,baremass)
      call KDDW4(DR1,DR2,u,.false.,baremass)
      print *,"IKDDW4.KDDW4:",maxval(abs(DR2-RR))

      baremass=0.05d0
      MTYPE=1
      dwkernel=3
      print *,""
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"dwkernel:",dwkernel
      print *,"DAGGER:    true"
      call KDDW4(RR,DR1,u,.true.,baremass)
      call IKDDW4(DR1,DR2,u,.true.,baremass)
      print *,"KDDW4.IKDDW4:",maxval(abs(DR2-RR))
      call IKDDW4(RR,DR1,u,.true.,baremass)
      call KDDW4(DR1,DR2,u,.true.,baremass)
      print *,"IKDDW4.KDDW4:",maxval(abs(DR2-RR))

      baremass=0.05d0
      MTYPE=3
      dwkernel=3
      print *,""
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"dwkernel:",dwkernel
      print *,"DAGGER:    false"
      call KDDW4(RR,DR1,u,.false.,baremass)
      call IKDDW4(DR1,DR2,u,.false.,baremass)
      print *,"KDDW4.IKDDW4:",maxval(abs(DR2-RR))
      call IKDDW4(RR,DR1,u,.false.,baremass)
      call KDDW4(DR1,DR2,u,.false.,baremass)
      print *,"IKDDW4.KDDW4:",maxval(abs(DR2-RR))

      baremass=0.05d0
      MTYPE=3
      dwkernel=3
      print *,""
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"dwkernel:",dwkernel
      print *,"DAGGER:    true"
      call KDDW4(RR,DR1,u,.true.,baremass)
      call IKDDW4(DR1,DR2,u,.true.,baremass)
      print *,"KDDW4.IKDDW4:",maxval(abs(DR2-RR))
      call IKDDW4(RR,DR1,u,.true.,baremass)
      call KDDW4(DR1,DR2,u,.true.,baremass)
      print *,"IKDDW4.KDDW4:",maxval(abs(DR2-RR))

      print *,""
      print *,"TEST EQUIVALENCE KDDW4 AND OVERLAP WILSON FORMULATIONS"
      print *,""


      baremass=0.05d0
      MTYPE=1
      dwkernel=3
      print *,""
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"dwkernel:",dwkernel
      print *,"DAGGER:    false"
      call KDDW4(RR,DR1,u,.false.,baremass)
      call DOLMW(RR,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOL-KDDW4 ',err
      call IKDDW4(RR,DR1,u,.false.,baremass)
      call IDOLMW(RR,DR2,u,.false.,baremass,SRF)
!      call IDOLop(RR,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'IDOL-IKDDW4 ',err

      baremass=0.05d0
      MTYPE=1
      dwkernel=3
      print *,""
      print *,"baremass:",baremass
      print *,"MTYPE:",MTYPE
      print *,"dwkernel:",dwkernel
      print *,"DAGGER:    true"
      call KDDW4(RR,DR1,u,.true.,baremass)
      call DOLMW(RR,DR2,u,.true.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOL-KDDW4 ',err
      call IKDDW4(RR,DR1,u,.true.,baremass)
      call IDOLMW(RR,DR2,u,.true.,baremass,SRF)
!      call IDOLop(RR,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'IDOL-IKDDW4 ',err


      stop

      print *,""
      print *,"TEST EQUIVALENCE OF OVERLAP AND DOMAIN WALL FORMULATIONS"
      print *,""

      MTYPE=1
      do i2=1,2
        dwkernel=i2 ! 1:Shamir 2:Wilson
        print *,""
        print *,"MTYPE:",MTYPE,"DWKERNEL:",dwkernel
        print *,""
        call KDDW4(RR,DR1,u,.false.,baremass)
!        call DOLop(RR,DR2,u,.false.,baremass,SRF)
        err=maxval(abs(DR1-DR2))
        print *,'DOL-KDDW4:',err
        call KDDW4(RR,DR1,u,.true.,baremass)
!        call DOLop(RR,DR2,u,.true.,baremass,SRF)
        err=maxval(abs(DR1-DR2))
        print *,'DOLdag-KDDW4dag ',err
      end do

      return

      dwkernel=1
      MTYPE=3
      call KDDW4(RR,DR1,u,.false.,baremass)
      MTYPE=4
!      call DOLop(RR,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOL-KDDW4 M4/M3 Shamir:',err
      MTYPE=3
      call KDDW4(RR,DR1,u,.true.,baremass)
      MTYPE=4
!      call DOLop(RR,DR2,u,.true.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOLdag-KDDW4dag M4/M3 Shamir:',err

      dwkernel=2
      MTYPE=3
      call KDDW4(RR,DR1,u,.false.,baremass)
      MTYPE=4
!      call DOLop(RR,DR2,u,.false.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOL-KDDW4 M4/M3 Wilson:',err
      MTYPE=3
      call KDDW4(RR,DR1,u,.true.,baremass)
      MTYPE=4
!      call DOLop(RR,DR2,u,.true.,baremass,SRF)
      err=maxval(abs(DR1-DR2))
      print *,'DOLdag-KDDW4dag M4/M3 Wilson:',err

      return
      end subroutine equivDOL_DDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module testDequivmod
