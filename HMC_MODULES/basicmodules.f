!!!#define TWODIMENSIONS

! declare usage of Intel Fortran
!#define IFORTRAN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module pacc
      implicit none
      integer,parameter :: prc = selected_real_kind(15,307) 
      integer,parameter :: zprc = selected_real_kind(15,307) 
c      integer,parameter :: prc = selected_real_kind(33,4931) 
c      real(prc),parameter :: resid=1q-8
c      real(prc),parameter :: resid=1q-10
c      real(prc),parameter :: resid=1q-12
c      real(prc),parameter :: resid=1q-14
      real(prc),parameter :: resid=1q-20
c      real(prc),parameter :: resid=1q-20
      end module pacc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module arraysizes
      implicit none
      integer,parameter :: Ns = 4
      integer,parameter :: Nt = 4
#ifdef TWODIMENSIONS
      integer,parameter :: NDS = 1 ! no spatial dimensions
      integer,parameter :: NDT = 2 ! no total dimensions
#else
      integer,parameter :: NDS = 2 ! no spatial dimensions
      integer,parameter :: NDT = 3 ! no total dimensions
#endif
      integer,parameter :: Nv = Ns**NDS*Nt
      integer,parameter :: Ndc = 4 ! no dirac components (1 for staggered, 2 for 2-comp, 4 for 4-comp)
      end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module numbers
      use pacc
      implicit none
      real(prc),parameter :: zero = 0.0q0
      real(prc),parameter :: one = 1.0q0
      real(prc),parameter :: two = 2.0q0
      real(prc),parameter :: half = one/two
      complex(prc),parameter :: cone = cmplx(one,zero,prc)
      complex(prc),parameter :: czero = cmplx(zero,zero,prc)
      complex(prc),parameter :: ctwo = cmplx(two,zero,prc)
      complex(prc),parameter :: zi = cmplx(zero,one,prc)
      real(prc),parameter :: pi = two*two*atan(one)
      real(prc),parameter :: rt3 = sqrt(3.0q0)
      real(prc),parameter :: rt2 = sqrt(2.0q0)
      end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module options
      use pacc
      use numbers
      implicit none

!     debug options
      integer,parameter :: VERBOSE=1 ! 0,1,2
      logical :: DUPLICATE

!     gauge field options
      integer GAUGETYPE ! 1=QED3, 2=Thirring3, 3=Gross-Neveau
      real(prc) :: gbeta

!     HMC parameters
      real(prc) :: HMC_dt,HMC_etime
      integer :: HMC_tsmax,tsav
      logical :: QUENCHED

!     dirac options
      integer :: HMCtype ! 1 for HMC, 2 for RHMC
      integer Nferms ! set by HMCtype 
      integer DWkernel ! 1 for Shamir, 2 for Wilson, Mobius implementation incomplete
      integer MTYPE ! 1,2,3,4 - 2 is alt form of 3 (5 not implemented)
      integer OLTYPE ! 1,2 - 1 is direct calculation, 2 is domain wall calculation (K-type)
      integer,parameter :: RFTYPE=1 ! 1,2 - 1 is partial fraction, 2 is multiplicative 
      real(prc) :: baremass ! bare mass
      real(prc) :: MDW ! domain wall height/overlap mass
      integer,parameter :: Npf=8 ! no of terms for partial fraction DOL formulation, redundant
!      integer,parameter :: Ls=2*Npf+1 ! for partial fraction reconstruction
#ifndef LSVALUE
#define LSVALUE 8
#endif
      integer,parameter :: Ls=LSVALUE ! for domain wall construction

!     runtime options
      integer Naux ! generate Nauxilary fields
      integer Nswp ! use Nswp HMC sweeps for each aux field

!      real(prc),parameter :: amob=one ! am for Mobius Dm=(a+b).Dw.(2+(a-b)Dw)^-1 
!      real(prc),parameter :: bmob=zero ! b=1, Dw, b=0, Shamir (with a = 1) 

!     measurement options
c      integer,parameter :: Nnoise=20 ! number of loops for noisy estimator
      contains
      subroutine printOptions()
      implicit none

      print *,"GAUGETYPE:",GAUGETYPE
      print *,"gbeta:",gbeta

      print *,"Nferms:",Nferms
      print *,"DWkernel:",DWkernel
      print *,"MTYPE:",MTYPE 
      print *,"OLTYPE:",OLTYPE 
      print *,"baremass:",baremass
      print *,"MDW(dom wall height):",MDW
      print *,"Npf:",Npf
      print *,"Ls:",Ls
      return
      end subroutine printOptions
      end module options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module countmod
      implicit none

      logical COUNT_IS_STARTED
      integer oc_idx,ic_idx 
      integer outer_count
      integer inner_count

      end module countmod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module timer
      implicit none
      character(len=10) :: sdate,stime,szone,edate
      integer,dimension(8) :: svals,evals
      real tstart,tend ! cpu seconds

      contains

      subroutine init_timer()
      implicit none

      call date_and_time(sdate,stime,szone,svals)
      print *,"Date:",sdate,"time:",stime
      call cpu_time(tstart)
      return 
      end subroutine init_timer
   

      end module timer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module paulimodule
      use pacc
      implicit none
      complex(prc) :: sigPauli(2,2,3)
      integer :: pauliin(3,2)
      complex(prc) :: paulival(3,2)
   
      contains

      subroutine setPauliMatrices
      use numbers
      implicit none

      sigPauli=czero
      sigPauli(1,2,1)=cone
      sigPauli(2,1,1)=cone
      pauliin(1,1)=2
      pauliin(1,2)=1
      paulival(1,1)=cone
      paulival(1,2)=cone
      sigPauli(1,2,2)=-zi
      sigPauli(2,1,2)=zi
      pauliin(2,1)=2
      pauliin(2,2)=1
      paulival(2,1)=-zi
      paulival(2,2)=zi
      sigPauli(1,1,3)=cone
      sigPauli(2,2,3)=-cone
      pauliin(3,1)=1
      pauliin(3,2)=2
      paulival(3,1)=cone
      paulival(3,2)=-cone
    
      return
      end subroutine setPauliMatrices
      end module paulimodule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module determinants
      use pacc
      use numbers
      implicit none
      contains
      function det2(M)
      implicit none
      complex(prc) :: M(2,2)
      complex(prc) :: det2

      det2=M(1,1)*M(2,2)-M(1,2)*M(2,1)
      return
      end function det2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function det3(M)
      implicit none
      complex(prc) :: M(3,3)
      complex(prc) :: det3
      complex(prc) :: T(2,2)

      T=czero
      T=M(2:3,2:3)
      det3=M(1,1)*det2(T)
      T=czero
      T(:,1)=M(2:3,3)
      T(:,2)=M(2:3,1)
      det3=det3+M(1,2)*det2(T)
      T=czero
      T=M(2:3,1:2)
      det3=det3+M(1,3)*det2(T)

      return
      end function det3
      end module determinants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module gammas
      use numbers
      implicit none
      complex(prc) gamval(6,4)
      integer gamin(6,4)
      contains

      subroutine setGammas
      use arraysizes
      implicit none

      if (Ndc.eq.2) then
        call setGammas2d()
        return
      endif

c     gamma_1
      gamval(1,1)=-zi
      gamval(1,2)=-zi
      gamval(1,3)= zi
      gamval(1,4)= zi
c
      gamin(1,1)=4
      gamin(1,2)=3
      gamin(1,3)=2
      gamin(1,4)=1

c     gamma_2
      gamval(2,1)=-cone
      gamval(2,2)= cone
      gamval(2,3)= cone
      gamval(2,4)=-cone
c
      gamin(2,1)=4
      gamin(2,2)=3
      gamin(2,3)=2
      gamin(2,4)=1

c     gamma_3
      gamval(3,1)=-zi
      gamval(3,2)= zi
      gamval(3,3)= zi
      gamval(3,4)=-zi
c
      gamin(3,1)=3
      gamin(3,2)=4
      gamin(3,3)=1
      gamin(3,4)=2

c     gamma_4
      gamval(4,1)= one
      gamval(4,2)= one
      gamval(4,3)= -one
      gamval(4,4)= -one
c
      gamin(4,1)=1
      gamin(4,2)=2
      gamin(4,3)=3
      gamin(4,4)=4

c     gamma_5 = gamma_1 * gamma_2 * gamma_3 * gamma_4
      gamval(5,1)=-one
      gamval(5,2)=-one
      gamval(5,3)=-one
      gamval(5,4)=-one
c
      gamin(5,1)=3
      gamin(5,2)=4
      gamin(5,3)=1
      gamin(5,4)=2

c     gamma_4 * gamma_5 (called gamma_3 gamma_5 in notes)
      gamval(6,1)=-one
      gamval(6,2)=-one
      gamval(6,3)= one
      gamval(6,4)= one
c
      gamin(6,1)=3
      gamin(6,2)=4
      gamin(6,3)=1
      gamin(6,4)=2
      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setGammas2d
      use arraysizes
      implicit none

c     gamma_1
      gamval(1,1)=cone
      gamval(1,2)=cone
c
      gamin(1,1)=2
      gamin(1,2)=1

c     gamma_2
      gamval(2,1)=-cone
      gamval(2,2)= cone
c
      gamin(2,1)=2
      gamin(2,2)=1

c     gamma_3=gamma_4=gamma_5
      gamval(3,1)= cone
      gamval(3,2)=-cone
c
      gamin(3,1)=1
      gamin(3,2)=2

      gamval(4,1)= cone
      gamval(4,2)=-cone
c
      gamin(4,1)=1
      gamin(4,2)=2

      gamval(5,1)= cone
      gamval(5,2)=-cone
c
      gamin(5,1)=1
      gamin(5,2)=2

      return
      end subroutine setGammas2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mGmu(phi,mu)
      use pacc
      use arraysizes
      implicit none
      complex(prc) phi(Nv,4)
      integer mu
      complex(prc) tmp(4)
      integer v,d,di

      do v=1,Nv
        do d=1,4
          di=gamin(mu,d)
          tmp(d)=gamval(mu,d)*phi(v,di)
        end do
        phi(v,:)=tmp
      end do
        
      return
      end subroutine mGmu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mGmu4(phi,mu)
      use pacc
      use arraysizes
      implicit none
      complex(prc) phi(4)
      integer mu
      complex(prc) tmp(4)
      integer d,di

      do d=1,4
        di=gamin(mu,d)
        tmp(d)=gamval(mu,d)*phi(di)
      end do
      phi=tmp
        
      return
      end subroutine mGmu4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine matrixGmu(M,mu)
      use pacc
      use arraysizes
      implicit none
      complex(prc) M(4,4)
      integer mu
      complex(prc) tmp(4,4)
      integer d,di,i

      do d=1,4
        di=gamin(mu,d)
        do i=1,4
          tmp(di,i)=gamval(mu,d)*M(i,di)
        end do
      end do

      return
      end subroutine matrixGmu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mGmu5(phi,mu)
      use pacc
      use arraysizes
      use options
      implicit none
      complex(prc) phi(Nv,4,Ls)
      integer mu
      complex(prc) tmp(Nv,4)
      integer v

      do v=1,Ls
        tmp=phi(:,:,v)
        call mGmu(tmp,mu)
        phi(:,:,v)=tmp
      end do
        
      return
      end subroutine mGmu5
      end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module indices
      use arraysizes
      implicit none
      integer id(Nv,NDT),iu(Nv,NDT)
      contains

#ifndef TWODIMENSIONS
      subroutine setIndices
      implicit none
      integer j1,j2,j3
      integer ic
c*******************************************************************
c
c     loads the addresses required during the update
c
c*******************************************************************
      do 30 j3=1,Nt
      do 30 j2=1,Ns
      do 30 j1=1,Ns
      ic=((j3-1)*Ns+(j2-1))*Ns+j1
      call ia(j1-1,j2,j3,id(ic,1))
      call ia(j1+1,j2,j3,iu(ic,1))
      call ia(j1,j2-1,j3,id(ic,2))
      call ia(j1,j2+1,j3,iu(ic,2))
      call ia(j1,j2,j3-1,id(ic,3))
      call ia(j1,j2,j3+1,iu(ic,3))
  30  continue
      return
      end
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ia(i1,i2,i3,nnn)
      implicit none
      integer i1,i2,i3
      integer nnn
      integer n1,n2,n3
c*******************************************************************
c
c     address calculator
c
c*******************************************************************
      n1=i1
      n2=i2
      n3=i3 
      if(n1) 2,2,3 ! if n1 <=0 goto 2, otherwise goto 3
   2  n1=n1+Ns
      go to 4
   3  if(n1-Ns) 4,4,5 ! if n1-Ns <= 0 goto 4, otherwise goto 5
   5  n1=n1-Ns
   4  if(n2) 6,6,7
   6  n2=n2+Ns
      go to 8
   7  if(n2-Ns) 8,8,9
   9  n2=n2-Ns
   8  if(n3) 10,10,11
  10  n3=n3+Nt 
      go to 12
  11  if(n3-Nt) 12,12,13
  13  n3=n3-Nt   
  12  nnn=((n3-1)*Ns+(n2-1))*Ns+n1
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getLatticeCoords(nnn,ix,iy,it)
      implicit none
      integer ix,iy,it
      integer nnn
      integer tmp
! coords are given by nnn=(it-1)*Ns*Ns+(iy-1)*Ns+ix
      tmp=mod(nnn,Ns*Ns)
c      print *,tmp
      ix=mod(tmp,Ns)
      if (ix.eq.0) then
        ix=Ns
      endif
c      print *,"ix:",ix
      iy=(tmp-ix)/Ns+1
      if (iy.eq.0) then 
        iy=Ns
      endif
      it=(nnn-ix-(iy-1)*Ns)/(Ns*Ns)+1
      if (it.eq.0) then
        it=Nt
      endif
      return
      end subroutine getLatticeCoords
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer function np1(idx,N)
      implicit none
      integer idx,N

      if (idx.eq.N) then
        np1=1
        return
      endif
      np1=idx+1
      
      return
      end function np1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer function nm1(idx,N)
      implicit none
      integer idx,N

      if (idx.eq.1) then
        nm1=N
        return
      endif
      nm1=idx-1
      
      return
      end function nm1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer function dist(v,i,j,k)
      use arraysizes
      implicit none
      integer v,i,j,k
      integer vt,r
      integer v2t
      integer j1,j2,j3
      integer ri,rj,rk

c      print *,"v:",v,"i:",i,"j:",j,"k:",k

      j1=mod(v-1,Ns)+1
      vt=v-(j1-1)
      v2t=mod(vt-1,Ns*Ns)+1
      j2=v2t/Ns+1
      vt=vt-(j2-1)*Ns
      j3=vt/(Ns*Ns)+1
c      print *,"j1:",j1,"j2:",j2,"j3:",j3
      
c      print *,v,((j3-1)*Ns+(j2-1))*Ns+j1

      ri=abs(j1-i)
      if (2*ri > Ns) ri=Ns-ri
      rj=abs(j2-j)
      if (2*rj > Ns) rj=Ns-rj      
      rk=abs(j3-k)
c      print *,"rk:",rk
      if (2*rk > Nt) rk=Nt-rk      
c      print *,"ri:",ri,"rj:",rj,"rk:",rk

      dist=ri+rj+rk
      return 
      end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setIndices2d
      implicit none
      integer ix,it
      integer itp1,itm1,ixp1,ixm1
      integer ic,ip

      do it=1,Nt
        itp1= np1(it,Nt)
        itm1= nm1(it,Nt)
        do ix=1,Ns
          ixp1=np1(ix,Ns)
          ixm1=nm1(ix,Ns)
          ic=(it-1)*Ns+ix
          ip=(it-1)*Ns+ixm1
          id(ic,1)=ip
          ip=(itm1-1)*Ns+ix
          iu(ic,2)=ip
          ip=(itp1-1)*Ns+ix
          id(ic,2)=ip
        end do
      end do
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module indices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module cgcountmodule
      implicit none
      integer ncgits,ncgmax
      real conv
      end module
