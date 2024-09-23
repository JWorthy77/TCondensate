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
      subroutine WilsonDerivs(dSdA,eta,nu,DAG) ! dSdA = eta^dag.dDwdA.nu
      use numbers
      use gammas
      use indices
      implicit none
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4),intent(in) ::  eta,nu
      logical :: DAG
      complex(prc) :: dSdAC(Nv,3)

      call WilsonDerivsComplex(dSdAC,eta,nu,DAG) 
      dSdA=dSdAC

      return
      end subroutine WilsonDerivs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WilsonDerivsComplex(dSdA,eta,nu,DAG) ! dSdA = eta^dag.dDwdA.nu
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
          dSdA(i,mu)=dSdA(i,mu)-zi*half*dot_product(lhs1,rhs1)

          lhs2=eta(iu(i,mu),:)
          rhs2=nu(i,:)
          dSdA(i,mu)=dSdA(i,mu)+zi*half*dot_product(lhs2,rhs2)

          lhs1=eta(i,:)
          rhs1=nu(iu(i,mu),:)
          call mGmu4(rhs1,mu)
          dSdA(i,mu)=dSdA(i,mu)+pm1*zi*half*dot_product(lhs1,rhs1)

          lhs2=eta(iu(i,mu),:)
          rhs2=nu(i,:)
          call mGmu4(rhs2,mu)
          dSdA(i,mu)=dSdA(i,mu)+pm1*zi*half*dot_product(lhs2,rhs2)

        end do
      end do
        
      return
      end subroutine WilsonDerivsComplex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module WilsonDirac
