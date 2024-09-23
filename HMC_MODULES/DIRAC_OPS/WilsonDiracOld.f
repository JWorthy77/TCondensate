      module WilsonDirac
      use pacc
      use arraysizes
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DWilson(R,DR,u,DAGGER,mass)
      use gammas
      use indices
      use options
      implicit none
c     calculates DR = DWilson*R
      complex(prc),intent(in) :: R(Nv,Ndc)
      complex(prc),intent(out) :: DR(Nv,Ndc)
      complex(prc),intent(in) :: u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      real(prc) mult
      integer i,idirac
      integer mu,igork

      mult=one
      if (DAGGER) then
        mult=-one
      end if

!     mass term 
      do idirac=1,Ndc
        do i=1,Nv
          DR(i,idirac) = mass*R(i,idirac)
        enddo
      enddo

!     Dirac term (anti-Hermitian)
      do mu=1,NDT
        do idirac=1,Ndc
          igork=gamin(mu,idirac)
          do i=1,Nv
            DR(i,idirac)=DR(i,idirac)
     &        +mult*gamval(mu,idirac)*
     &         (u(i,mu)*R(iu(i,mu),igork)
     &         -conjg(u(id(i,mu),mu))*R(id(i,mu),igork))/two
          enddo
        enddo
      enddo

!     Wilson term (Hermitian)
      do mu=1,NDT
        do idirac=1,Ndc
          do i=1,Nv
            DR(i,idirac)=DR(i,idirac)
     &        -( u(i,mu)*R(iu(i,mu),idirac) - two*R(i,idirac)
     &         +conjg(u(id(i,mu),mu))*R(id(i,mu),idirac))/two
          enddo
        enddo
      enddo

      return
      end subroutine DWilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DWilsonJW(R,DR,u,DAG,mass)
      use gammas
      use indices
      use options
      implicit none
c     calculates DR = DWilson*R
      complex(prc),intent(in) :: R(Nv,Ndc)
      complex(prc),intent(out) :: DR(Nv,Ndc)
      complex(prc),intent(in) :: u(Nv,NDT)
      logical,intent(in) :: DAG
      real(prc),intent(in) :: mass
      real(prc) mult
      complex(prc),dimension(Ndc) :: gRp1,gRm1
      integer i,ipmu,immu,mu

      mult=one
      if (DAG) then
        mult=-one
      end if
!     mass term 
      DR = mass*R

!     Dirac term (anti-Hermitian)
      do mu=1,NDT
        do i=1,Nv
          ipmu=iu(i,mu)
          immu=id(i,mu)
          gRp1=R(ipmu,:)
          call mGmu4(gRp1,mu)
          gRm1=R(immu,:)
          call mGmu4(gRm1,mu)
          DR(i,:)=DR(i,:) + mult*half*
     &     (u(i,mu)*gRp1 - conjg(u(immu,mu))*gRm1)
        enddo
      enddo

!     Wilson term (Hermitian)
      do mu=1,3
        do i=1,Nv
          ipmu=iu(i,mu)
          immu=id(i,mu)
          DR(i,:)=DR(i,:) - half*
     &     ( u(i,mu)*R(ipmu,:) - two*R(i,:)
     &         +conjg(u(immu,mu))*R(immu,:) )
        enddo
      enddo

      return
      end subroutine DWilsonJW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine G5DW(R,DR,u,DAGGER,mass)
      use gammas
      implicit none
c     calculates DR = G5.DWilson*R
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TMP(Nv,4)

      if(.not.DAGGER) then
        call DWilson(R,DR,u,DAGGER,mass)
        call mGmu(DR,5)
      else
       TMP=R
        call mGmu(TMP,5)
        call DWilson(TMP,DR,u,DAGGER,mass)
      end if
      
      return
      end subroutine G5DW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DdagD(R,DR,u,DAGGER,mass)
      implicit none
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

      call DWilson(R,TMP,u,DAGGER,mass)
      call DWilson(TMP,DR,u,.not.DAGGER,mass)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DdagDpC(R,DR,u,DAGGER,mass,add) 
      implicit none
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

      call DWilson(R,TMP,u,DAGGER,mass)
      call DWilson(TMP,DR,u,.NOT.DAGGER,mass)
      DR=DR+add*R

      return
      end subroutine DdagDpC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDdagD(RR,DR,u,DAGGER,mass)
!     solve DdagD.DR = RR
!     note DAGGER is for M, not for MdagM
      implicit none
      integer,parameter :: kferm = Ndc*Nv
      integer,parameter :: niterc=Nv*10
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call DWilson(RR,x1,u,DAGGER,mass)
      call DWilson(x1,x2,u,.not.DAGGER,mass)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call DWilson(p,x1,u,DAGGER,mass)
        call DWilson(x1,x2,u,.not.DAGGER,mass)
        alphan=sum(conjg(r)*r)
        alphad=sum(conjg(p)*x2)
        alpha=alphan/alphad

        DR=DR+alpha*p
        r=r-alpha*x2
        betan=sum(conjg(r)*r)
        if (betan.lt.resid) goto 8
        beta=betan/alphan
        p=r+beta*p
      end do

8     continue
!      print *,itercg,niterc,betan
      return
      end subroutine IDdagD     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDW(RR,DR,u,DAGGER,mass,add)
!     solve Dw.DR = RR or Dw(dag).DR = RR
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,4)
      procedure(),pointer :: Mptr

      if (.not. DAGGER) then
        call DWilson(RR,TMP,u,.true.,mass)
        call IDdagD(TMP,DR,u,.false.,mass)
      elseif (DAGGER) then
        call DWilson(RR,TMP,u,.false.,mass)
        call IDdagD(TMP,DR,u,.true.,mass)
      end if

      return
      end subroutine IDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      subroutine IG5DW(RR,DR,u,DAGGER,mass,add)
!!     solve Dw(u) DR = RR or DwD(u) DR = RR
!      use axbmodule1
!      implicit none
!      complex(prc) RR(Nv,4),DR(Nv,4)
!      complex(prc) u(Nv,3)
!      logical DAGGER
!      real(prc) mass
!      complex(prc) add
!      complex(prc) TMP(Nv,4)
!      procedure(),pointer :: Mptr
!
!      Mptr=>G5DW
!      if (.not. DAGGER) then
!        call Mptr(RR,TMP,u,.true.,mass,add)
!        call IMdagM(TMP,DR,u,.false.,mass,add,Mptr)
!      elseif (DAGGER) then
!        call IMdagM(RR,TMP,u,.true.,mass,add,Mptr)
!        call Mptr(TMP,DR,u,.false.,mass,add)
!      end if
!
!      return
!      end subroutine IG5DW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WilsonDerivsSJH(dSdA,eta,nu,DAG) 
      use numbers
      use gammas
      use indices
      implicit none
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4),intent(in) ::  eta,nu
      logical :: DAG
      complex(prc),dimension(4) :: lhs1,lhs2,rhs1,rhs2
      integer mu,idirac,i
      integer pm1

      pm1=one
      if (DAG) then
        pm1=-one
      endif

      dSdA=0
      do mu=1,3
        do i=1,Nv
          lhs1=eta(i,:)
          rhs1=nu(iu(i,mu),:)
          lhs2=eta(iu(i,mu),:)
          rhs2=nu(i,:)
          dSdA(i,mu)=real( half*zi*(dot_product(lhs1,rhs1)
     &                             -dot_product(lhs2,rhs2) ) ) 
          call mGmu4(rhs1,mu)
          call mGmu4(rhs2,mu)
          dSdA(i,mu)=dSdA(i,mu)-pm1*real( half*zi*(
     &                            dot_product(lhs1,rhs1)
     &                           +dot_product(lhs2,rhs2) ) ) 
        end do
      end do
        
      return
      end subroutine WilsonDerivsSJH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WilsonDerivsJW(dSdA,eta,nu,DAG) ! dSdA = eta^dag.dDwdA.nu
      use numbers
      use gammas
      use indices
      implicit none
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4),intent(in) ::  eta,nu
      logical :: DAG
      complex(prc),dimension(4) :: lhs1,lhs2,rhs1,rhs2
      integer mu,idirac,i
      integer pm1

      pm1=one
      if (DAG) then
        pm1=-one
      endif

      dSdA=0
      do mu=1,3
        do i=1,Nv

          lhs1=eta(i,:)
          rhs1=nu(iu(i,mu),:)
          dSdA(i,mu)=dSdA(i,mu)+zi*half*dot_product(lhs1,rhs1)

          lhs2=eta(iu(i,mu),:)
          rhs2=nu(i,:)
          dSdA(i,mu)=dSdA(i,mu)-zi*half*dot_product(lhs2,rhs2)

          lhs1=eta(i,:)
          rhs1=nu(iu(i,mu),:)
          call mGmu4(rhs1,mu)
          dSdA(i,mu)=dSdA(i,mu)-pm1*zi*half*dot_product(lhs1,rhs1)

          lhs2=eta(iu(i,mu),:)
          rhs2=nu(i,:)
          call mGmu4(rhs2,mu)
          dSdA(i,mu)=dSdA(i,mu)-pm1*zi*half*dot_product(lhs2,rhs2)

        end do
      end do
        
      return
      end subroutine WilsonDerivsJW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WilsonDerivsComplexJW(dSdA,eta,nu,DAG) ! dSdA = eta^dag.dDwdA.nu
      use numbers
      use gammas
      use indices
      implicit none
      complex(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4),intent(in) ::  eta,nu
      logical :: DAG
      complex(prc),dimension(4) :: lhs1,lhs2,rhs1,rhs2
      integer mu,idirac,i
      integer pm1

      pm1=one
      if (DAG) then
        pm1=-one
      endif

      dSdA=0
      do mu=1,3
        do i=1,Nv

          lhs1=eta(i,:)
          rhs1=nu(iu(i,mu),:)
          dSdA(i,mu)=dSdA(i,mu)+zi*half*dot_product(lhs1,rhs1)

          lhs2=eta(iu(i,mu),:)
          rhs2=nu(i,:)
          dSdA(i,mu)=dSdA(i,mu)-zi*half*dot_product(lhs2,rhs2)

          lhs1=eta(i,:)
          rhs1=nu(iu(i,mu),:)
          call mGmu4(rhs1,mu)
          dSdA(i,mu)=dSdA(i,mu)-pm1*zi*half*dot_product(lhs1,rhs1)

          lhs2=eta(iu(i,mu),:)
          rhs2=nu(i,:)
          call mGmu4(rhs2,mu)
          dSdA(i,mu)=dSdA(i,mu)-pm1*zi*half*dot_product(lhs2,rhs2)

        end do
      end do
        
      return
      end subroutine WilsonDerivsComplexJW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WilsonDerivs(dSdA,eta,nu,DAG) ! dSdA = eta^dag.dDwdA.nu
      use numbers
      use gammas
      use indices
      implicit none
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4),intent(in) ::  eta,nu
      logical :: DAG
      complex(prc),dimension(4) :: lhs,rhs1,rhs2
      integer mu,idirac,i
      integer pm1

      pm1=one
      if (DAG) then
        pm1=-one
      endif

      dSdA=0
      do mu=1,3
        do i=1,Nv
          lhs=eta(i,:)
          rhs1=pm1*zi*half*nu(iu(i,mu),:)
          call mGmu4(rhs1,mu)
          rhs2=  -zi*half*nu(iu(i,mu),:)
          dSdA(i,mu)=dot_product(lhs,rhs1+rhs2)
          lhs=eta(iu(i,mu),:)
          rhs1=pm1*zi*half*nu(i,:)
          call mGmu4(rhs1,mu)
          rhs2=  zi*half*nu(i,:)
          dSdA(i,mu)=dSdA(i,mu)+dot_product(lhs,rhs1+rhs2)
        end do
      end do
        
      return
      end subroutine WilsonDerivs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WilsonDerivsComplex(dSdA,eta,nu,DAG)
      use numbers
      use gammas
      use indices
      implicit none
      complex(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4),intent(in) ::  eta,nu
      logical :: DAG
      complex(prc),dimension(4) :: lhs,rhs1,rhs2
      integer mu,idirac,i
      integer pm1

      pm1=one
      if (DAG) then
        pm1=-one
      endif

      dSdA=0
      do mu=1,3
        do i=1,Nv
          lhs=eta(i,:)
          rhs1=pm1*zi*half*nu(iu(i,mu),:)
          call mGmu4(rhs1,mu)
          rhs2=  -zi*half*nu(iu(i,mu),:)
          dSdA(i,mu)=dot_product(lhs,rhs1+rhs2)
          lhs=eta(iu(i,mu),:)
          rhs1=pm1*zi*half*nu(i,:)
          call mGmu4(rhs1,mu)
          rhs2=  zi*half*nu(i,:)
          dSdA(i,mu)=dSdA(i,mu)+dot_product(lhs,rhs1+rhs2)
        end do
      end do
        
      return
      end subroutine WilsonDerivsComplex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module WilsonDirac
