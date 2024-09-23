      module kernelspectrarange
      use pacc
      use arraysizes
      use numbers
      use utilsmod
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcKernelRange()
      use gaugefield
      use options
      use IOmodule
      use statsmod
      implicit none
      integer icf,idx
      character(len=80) fname 
      real(prc) lmax,lmin,lmaxav,lmaxsd,lminav,lminsd,cnav,cnsd
      real(prc),allocatable,dimension(:) :: vlmax,vlmin,vcn
      integer fstart,fstop,fskip,Nf
      logical evalMax,evalMin
      integer Nmax,Nmin,MAXMIN
      logical FEXISTS

      open(unit=31,file='kernelopts.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) fstart,fstop,fskip,evalMax,Nmax,evalMin,Nmin,dwkernel
     &,GAUGETYPE
        close(31)
        print *,"calcKernelRange"
        print *,"fstart,fstop,fskip:",fstart,fstop,fskip
        print *,"evalMax,evalMin:",evalMax,evalMin
        print *,"Nmax,Nmin:",Nmax,Nmin
        print *,"kernel:",dwkernel
        print *,"gaugetype:",GAUGETYPE
      else
        print *,"kernel options file not found"
        stop
      endif
      close(31)

      lmax=0
      lmin=0
      Nf=(fstop-fstart)/fskip+1
      allocate(vlmax(Nf),vlmin(Nf),vcn(Nf))
      idx=1
      do icf=fstart,fstop,fskip
        open(unit=11,file='DetailsExtrema.dat',access='append',
     &           status='unknown',form='formatted')
        call readThetaFromFile(icf,theta)
        call coef(u,theta)
  
        if (evalMax) then
          MAXMIN=1
          call estimateKernelExtrema(Nmax,MAXMIN,lmax)
        endif

        if (evalMin) then
          MAXMIN=-1
          call estimateKernelExtrema(Nmin,MAXMIN,lmin)
        endif

        vlmax(idx)=lmax
        vlmin(idx)=lmin
        vcn(idx)=lmax/lmin
        idx=idx+1

        write(11,'(I6,F5.1,F9.3,E12.4E2)') icf,-MDW,lmax,lmin
        close(11)
      end do

      call calcVarReal(Nf,vlmax,lmaxav,lmaxsd) 
      call calcVarReal(Nf,vlmin,lminav,lminsd) 
      call calcVarReal(Nf,vcn,cnav,cnsd) 

      open(unit=11,file='KernelExtrema.dat',access='append',
     &           status='unknown',form='formatted')
        write(11,*) "Nterms,av,max,min,sd,err"
        write(11,*) "Max:"
        write(11,*) Nf,lmaxav,maxval(vlmax),minval(vlmax),lmaxsd,lmaxsd/
     &                                                 sqrt(real(Nf))
        write(11,*) "Min:"
        write(11,*) Nf,lminav,maxval(vlmin),minval(vlmin),lminsd,lminsd/
     &                                                 sqrt(real(Nf))
        write(11,*) "Condition number:"
        write(11,*) Nf,cnav,maxval(vcn),minval(vcn),cnsd,cnsd/
     &                                                 sqrt(real(Nf))
      close(11)

      return
      end subroutine calcKernelRange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcQuenchedKernelRange()
      use gaugefield
      use gaugemodule
      use options
      use IOmodule
      use statsmod
      implicit none
      integer icf
      character(len=80) fname 
      real(prc) lmax,lmin,lmaxav,lmaxsd,lminav,lminsd,cnav,cnsd
      real(prc),allocatable,dimension(:) :: vlmax,vlmin,vcn
      integer Nf
      integer Nmax,Nmin,MAXMIN
      logical FEXISTS

      open(unit=31,file='qkernelopts.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) Nf,Nmax,Nmin,dwkernel,GAUGETYPE,gbeta
        close(31)
        print *,"calcQuencedKernelRange"
        print *,"Nf:",Nf
        print *,"Nmax,Nmin:",Nmax,Nmin
        print *,"kernel:",dwkernel
        print *,"gaugetype:",GAUGETYPE
        print *,"gbeta:",gbeta
      else
        print *,"kernel options file not found"
        stop
      endif
      close(31)

      lmax=0
      lmin=0
      allocate(vlmax(Nf),vlmin(Nf),vcn(Nf))
      do icf=1,Nf
        open(unit=11,file='DetailsExtrema.dat',access='append',
     &           status='unknown',form='formatted')

        if (GAUGETYPE.eq.1) then
          call makeQuenchedCosineThirringField()
        elseif (GAUGETYPE.eq.2) then
          call makeQuenchedGaussianThirringField()
        endif
        print *,"made quenched field"
        if (Nmax > 0) then
          MAXMIN=1
          call estimateKernelExtrema(Nmax,MAXMIN,lmax)
        endif

        if (Nmin > 0) then
          MAXMIN=-1
          call estimateKernelExtrema(Nmin,MAXMIN,lmin)
        endif

        vlmax(icf)=lmax
        vlmin(icf)=lmin
        vcn(icf)=lmax/lmin

        write(11,*) icf,-MDW,lmax,lmin
        close(11)
      end do

      print *,"calculate statistics"

      call calcVarReal(Nf,vlmax,lmaxav,lmaxsd) 
      call calcVarReal(Nf,vlmin,lminav,lminsd) 
      call calcVarReal(Nf,vcn,cnav,cnsd) 

      open(unit=11,file='QuenchedKernelExtrema.dat',access='append',
     &           status='unknown',form='formatted')
        write(11,*) "Nterms,av,max,min,sd,err"
        write(11,*) "Max:"
        write(11,*) Nf,lmaxav,maxval(vlmax),minval(vlmax),lmaxsd,lmaxsd/
     &                                                 sqrt(real(Nf))
        write(11,*) "Min:"
        write(11,*) Nf,lminav,maxval(vlmin),minval(vlmin),lminsd,lminsd/
     &                                                 sqrt(real(Nf))
        write(11,*) "Condition number:"
        write(11,*) Nf,cnav,maxval(vcn),minval(vcn),cnsd,cnsd/
     &                                                 sqrt(real(Nf))
      close(11)

      return
      end subroutine calcQuenchedKernelRange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine estimateKernelExtrema(Nmax,MAXMIN,eig)
!     the basic power method
      use options
      use rvmodule
      use gaugefield
      use WilsonDirac
      use WilsonExtraMod
      implicit none
      integer Nmax,MAXMIN
      real(prc) eig
      complex(prc),dimension(Nv,Ndc) :: R,TMP,DR
      integer i

      call setRVs(Nv*4,R)
      call normalise(R)
      open(unit=81,file='DetailsEigs.dat',access='append',
     &           status='unknown',form='formatted')
      do i=1,Nmax
        if (MAXMIN.eq.1) then
          call Hkernel(R,DR,u,.false.,-MDW)
        elseif (MAXMIN.eq.-1) then
          call IHkernel(R,DR,u,.false.,-MDW)
        endif
        R=DR
        eig=mag(R)
        R=R/eig
        print *,i,eig
        write(81,*) i,eig
      end do
      close(81)
      if (MAXMIN.eq.-1) then
        eig=one/eig
      endif

      return
      end subroutine estimateKernelExtrema
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getMinMaxHeigs(lmin,lmax)
      implicit none
      real(prc) lmin,lmax
 
      call estimateKernelExtrema(20,1,lmax)
      call estimateKernelExtrema(20,-1,lmin)
      return
      end subroutine getMinMaxHeigs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      real(prc) function mag(V)
c      implicit none
c      complex(prc) :: V(Nv*Ndc)
c      real(prc) :: rV(Nv*Ndc),iV(Nv*Ndc)
c
c      rV=real(V,prc)
c      iV=dimag(V)
c      mag=dot_product(rV,rV)+dot_product(iV,iV)
c      mag=sqrt(mag)
c      return
c      end function mag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine normalise(V)
c      implicit none
c      complex(prc) :: V(Nv*Ndc)
c      real(prc) :: nrm
c
c      nrm=mag(V)
c      V=V/nrm
c      return
c      end subroutine normalise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      real(prc) function magv(N,V)
c      implicit none
c      integer N
c      complex(prc) :: V(N)
c      real(prc) :: rV(N),iV(N)
c
c      rV=real(V,prc)
c      iV=dimag(V)
c      magv=dot_product(rV,rV)+dot_product(iV,iV)
c      magv=sqrt(magv)
c      return
c      end function magv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine normalisev(N,V)
c      implicit none
c      integer N
c      complex(prc) :: V(N)
c      real(prc) :: nrm
c
c      nrm=magv(N,V)
c      V=V/nrm
c      return
c      end subroutine normalisev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module kernelspectrarange
