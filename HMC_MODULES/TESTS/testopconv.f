      module testDconvmod
!     contains routines to test the Ls convergence of polynomial and
!     rational approximations to Dirac and related operators
      use arraysizes
      use numbers
      use options
      use rvmodule
      use ratfuncs
      use gaugefield
      use basicdiracopsmod
      use overlapmoduledev
      implicit none
      logical,parameter :: VBS_TDC=.true. 
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testOperatorConvergence()
      use gaugemodule
      use rvmodule
      implicit none
      character(len=80) fname

      MDW=one ! domain wall height
      MTYPE=3 ! mass term type
      dwkernel=2 ! 2-Wilson
      GAUGETYPE=2 ! 2-non compact
      baremass=0.05
      gbeta=2.0 

      call makeQuenchedAuxField()

!      call testOLConvergence()
      call testKDDW4Convergence()
c      call testDDWConvergence()

      return
      end subroutine testOperatorConvergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testOLConvergence()
      implicit none
      character(len=80) fname

      MTYPE=1 ; dwkernel=2 ; fname="convOLWM1.dat"
      call convOL(fname) 

      MTYPE=1 ; dwkernel=1 ; fname="convOLSM1.dat"
      call convOL(fname) 

      MTYPE=3 ; dwkernel=2 ; fname="convOLWM3.dat"
      call convOL(fname) 

      MTYPE=3 ; dwkernel=1 ; fname="convOLSM3.dat"
      call convOL(fname) 

      return
      end subroutine testOLConvergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine convOL(fname)
      use spectra
      use kernelspectrarange
      implicit none
!     use 12x12 mesh, and create quenched gauge field
      logical ESTEIGS
      character(len=80) :: fname
      type(sgnratfunc) :: SRF
      complex(prc) R(Nv,4)
      complex(prc),dimension(Nv,4,100) ::  HTtest,Zolotest
      real ev(1)
      real(prc) lmin,lmax
      integer j,jmin,jmax,shft

      print *,'test Ls convergence of overlap operators'
!      call setRVs(Nv*4,R)
      open(unit=12,file="R12x12.dat",status='unknown',form='formatted')
      read(12,*) R
      close(12)
      call getMinMaxHEigs(lmin,lmax)
      print *,"lmin:",lmin
      print *,"lmax:",lmax

      jmin=6 ; jmax=36 ; shft=2
      do j=jmin,jmax,shft
        call setHTcoeffs(j,SRF)
        call DOLop(R,HTtest(:,:,j),u,.false.,one/50,SRF)
        call setZoloCoeffs(j,SRF,lmin,lmax)
        call DOLop(R,Zolotest(:,:,j),u,.false.,one/50,SRF)
       print *,j,maxval(abs(HTtest(:,:,j))),maxval(abs(Zolotest(:,:,j)))
      end do

      open(unit=11,file=fname,status='unknown',form='formatted')

        do j=jmin+shft,jmax,shft
          print *,j,maxval(abs(HTtest(:,:,j)-HTtest(:,:,j-shft))),
     &              maxval(abs(Zolotest(:,:,j)-Zolotest(:,:,j-shft))),
     &              maxval(abs(Zolotest(:,:,j)-HTtest(:,:,j)))
          write(11,*) j,maxval(abs(HTtest(:,:,j)-HTtest(:,:,j-shft))),
     &              maxval(abs(Zolotest(:,:,j)-Zolotest(:,:,j-shft))),
     &              maxval(abs(Zolotest(:,:,j)-HTtest(:,:,j)))
        end do

      close(11)

      return
      end subroutine convOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testKDDW4Convergence()
      implicit none
      character(len=80) fstub

      MTYPE=1 ; dwkernel=2 ; fstub="WM1_"
      call convKDDW4(fstub) 

      MTYPE=1 ; dwkernel=1 ; fstub="SM1_"
      call convKDDW4(fstub) 

      MTYPE=3 ; dwkernel=2 ; fstub="WM3_"
      call convKDDW4(fstub) 

      MTYPE=3 ; dwkernel=1 ; fstub="SM3_"
      call convKDDW4(fstub) 

c      MTYPE=1 ; dwkernel=3 ; fstub="WZM1_"
c      call convKDDW4(fstub) 

c      MTYPE=3 ; dwkernel=3 ; fstub="WZM3_"
c      call convKDDW4(fstub) 

      return
      end subroutine testKDDW4Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine convKDDW4(fstub)
      use domainwallmod
      use dwcoeffs
      use zolomodule
      use spectra
      use kernelspectrarange
      use IOmodule
      implicit none
!     use 12x12 mesh, and create quenched gauge field
      logical ESTEIGS
      character(len=80) :: fstub,fname
      type(zolotarev) :: zolo
      type(sgnratfunc) :: SRF
      complex(prc) RR(Nv,4)
      complex(prc),dimension(Nv,4) ::  DWtest
      real ev(1)
      real(prc) lmin,lmax
      integer nc,i
      character(len=3) cnum

      print *,'test Ls convergence of overlap operators'
!      call setRVs(Nv*4,R)
      open(unit=12,file="R12x12.dat",status='unknown',form='formatted')
      read(12,*) RR
      close(12)

      call getMinMaxHEigs(lmin,lmax)
      call setZoloCoeffs(Ls,SRF,lmin,lmax)
      call setZolo(lmin,lmax,Ls,zolo)
      call getRoots(zolo)
      omega=one/zolo%roots
c      print *,"lmin:",lmin
c      print *,"lmax:",lmax

      cnum=itoa(Ls)
c      nc=len(trim(cnum))
c      do i=nc+1,3
c        cnum="0"//cnum
c      enddo
      print *,"cnum: ",cnum
      print *,"fstub: ",fstub
      
      fname=trim(fstub)//trim(cnum)//".dat"
      print *,"fname: ",fname

      call KDDW4(RR,DWtest,u,.false.,one/50)

      open(unit=11,file=fname,status='unknown',form='formatted')
      write(11,*) DWtest
      close(11)

      return
      end subroutine convKDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDDWConvergence()
      implicit none
      character(len=80) fstub

      print *,"need to think this through!"
      stop

      MTYPE=1 ; dwkernel=1 ; fstub="DDW_WM1_"
      call convDDW(fstub) 

!      MTYPE=1 ; dwkernel=1 ; fstub="DDW_SM1_"
!      call convDDW(fstub) 

      MTYPE=1 ; dwkernel=3 ; fstub="DDW_WZM1_"
      call convDDW(fstub) 

!      MTYPE=3 ; dwkernel=2 ; fstub="DDW_WM3_"
!      call convDDW(fstub) 

!      MTYPE=3 ; dwkernel=1 ; fstub="DDW_SM3_"
!      call convDDW(fstub) 

c      MTYPE=1 ; dwkernel=3 ; fstub="WZM1_"
c      call convKDDW4(fstub) 

c      MTYPE=3 ; dwkernel=3 ; fstub="WZM3_"
c      call convKDDW4(fstub) 

      return
      end subroutine testDDWConvergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine convDDW(fstub)
      use domainwallmod
      use dwcoeffs
      use zolomodule
      use spectra
      use kernelspectrarange
      use IOmodule
      implicit none
!     use 12x12 mesh, and create quenched gauge field
      logical ESTEIGS
      character(len=80) :: fstub,fname
      type(zolotarev) :: zolo
      type(sgnratfunc) :: SRF
      complex(prc) RR(Nv,4,Ls)
      complex(prc),dimension(Nv,4,Ls) ::  DWtest
      real ev(1)
      real(prc) lmin,lmax
      integer nc,i
      character(len=3) cnum

      print *,'test Ls convergence of overlap operators'
c      call setRVs(Nv*4*Ls,RR)
      open(unit=12,file="R5_12.dat",status='unknown',form='formatted')
c      write(12,*) RR
      read(12,*) RR
      close(12)

      call getMinMaxHEigs(lmin,lmax)
      call setZoloCoeffs(Ls,SRF,lmin,lmax)
      call setZolo(lmin,lmax,Ls,zolo)
      call getRoots(zolo)
      omega=one/zolo%roots
c      print *,"lmin:",lmin
c      print *,"lmax:",lmax

      cnum=itoa(Ls)
c      nc=len(trim(cnum))
c      do i=nc+1,3
c        cnum="0"//cnum
c      enddo
      print *,"cnum: ",cnum
      print *,"fstub: ",fstub
      
      fname=trim(fstub)//trim(cnum)//".dat"
      print *,"fname: ",fname

      call DDW(RR,DWtest,u,.false.,one/50)

      open(unit=11,file=fname,status='unknown',form='formatted')
      write(11,*) DWtest
      close(11)

      return
      end subroutine convDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makeQuenchedAuxField()
      use gaugemodule
      use rvmodule
      implicit none

      seed=0
      seed(1) = 167868904
      seed(2) = 32712
      seed(3) = seed(1)
      seed(4) = seed(3)
      print *,"seed:",seed
      call random_seed(put=seed)

      if (GAUGETYPE.eq.1) then
        call makeQuenchedCosineThirringField()
      elseif (GAUGETYPE.eq.2) then
        call makeQuenchedGaussianThirringField()
      endif

      return
      end subroutine makeQuenchedAuxField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module testDconvmod
