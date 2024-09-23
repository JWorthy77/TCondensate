      module NaiveDirac
      use pacc
      use arraysizes
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DNF(R,DR,u,DAGGER,mass)
      use gammas
      use indices
      use options
      implicit none
c     calculates DR = DNaive*R
      complex(prc),intent(in) :: R(Nv,Ndc)
      complex(prc),intent(out) :: DR(Nv,Ndc)
      complex(prc),intent(in) :: u(Nv,NDT)
      logical,intent(in) :: DAGGER
      real(prc),intent(in) :: mass
      real(prc) mult
      integer i,idirac
      integer mu,igork

      mult=one
      if (DAGGER) then
        mult=-one
      end if

!     mass term 
      DR = mass*R

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

      return
      end subroutine DNF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DNFdagDNF(R,DR,u,DAGGER,mass)
      implicit none
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

      call DNF(R,TMP,u,DAGGER,mass)
      call DNF(TMP,DR,u,.not.DAGGER,mass)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDNFdagDNF(RR,DR,u,DAGGER,mass)
!     solve DNFdagDNF.DR = RR - note DAGGER is for DNF, not for DNFdagDNF, which would be
!     identical
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
      call DNF(RR,x1,u,DAGGER,mass)
      call DNF(x1,x2,u,.not.DAGGER,mass)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call DNF(p,x1,u,DAGGER,mass)
        call DNF(x1,x2,u,.not.DAGGER,mass)
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
      end subroutine IDNFdagDNF     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDNF(RR,DR,u,DAGGER,mass,add)
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
        call DNF(RR,TMP,u,.true.,mass)
        call IDNFdagDNF(TMP,DR,u,.false.,mass)
      elseif (DAGGER) then
        call DNF(RR,TMP,u,.false.,mass)
        call IDNFdagDNF(TMP,DR,u,.true.,mass)
      end if

      return
      end subroutine IDNF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DNFDerivs(dSdA,eta,nu,DAG)
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
!          rhs2=  -zi*half*nu(iu(i,mu),:)
!          rhs2=  -zi*nu(iu(i,mu),:)
          dSdA(i,mu)=dot_product(lhs,rhs1)
          lhs=eta(iu(i,mu),:)
          rhs1=pm1*zi*half*nu(i,:)
          call mGmu4(rhs1,mu)
!          rhs2=  zi*half*nu(i,:)
!          rhs2=  zi*nu(i,:)
          dSdA(i,mu)=dSdA(i,mu)+dot_product(lhs,rhs1)
        end do
      end do
        
      return
      end subroutine DNFDerivs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DNFDerivsSJH(dSdA,eta,nu,DAG)
      use gammas
      use indices
      use options
      use numbers
      implicit none
      complex(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4),intent(in) ::  eta,nu
      logical :: DAG
      complex(prc),dimension(4) :: lhs,rhs1,rhs2
      integer mu,idirac,i,igork
      real(prc) mult

      mult=one
      if (DAG) then
        mult=-one
      end if

!      dSdA=0
!      do mu=1,NDT
!        do idirac=1,Ndc
!          igork=gamin(mu,idirac)
!          do i=1,Nv
!            tot=tot+LHS(i,idirac)*(
!     &         mult*gamval(mu,idirac)*
!     &         (u(i,mu)*RHS(iu(i,mu),igork)
!     &         -conjg(u(id(i,mu),mu))*RHS(id(i,mu),igork))/two
!          enddo
!        enddo
!      enddo

      do mu=1,3
        do idirac=1,4
          igork=gamin(mu,idirac)
          do i=1,Nv
            dSdA(i,mu)=dSdA(i,mu) + eta(i,idirac)*
     &        mult*gamval(mu,idirac) * zi * nu(iu(i,mu),igork) * half

            dSdA(id(i,mu),mu)=dSdA(id(i,mu),mu) + eta(i,idirac) *
     &         mult*gamval(mu,idirac) * zi * nu(id(i,mu),igork) * half
          enddo
        enddo
      enddo

      return
      end subroutine DNFDerivsSJH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module NaiveDirac
