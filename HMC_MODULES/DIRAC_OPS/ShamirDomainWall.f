!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module ShamirDomWall
      use arraysizes
      use options
      use WilsonDirac
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_Shamir(R,DR,u,DAGGER,mass)
      use gammas
      use ratfuncs
      implicit none
c     calculates DR = DDW*R where DDW is the Shamir domain wall formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TR(Nv,4)
      integer s
      real(prc) as

c     diagonal blocks
      do s=1,Ls
        call DWilson(R(:,:,s),DR(:,:,s),u,DAGGER,-MDW)
      end do
      DR=DR+R

c     projection blocks, P_-.R  = R(:,3:4,s), P_+.R = R(:,1:2,s)
      if (.not.DAGGER) then
        do s=1,Ls-1
          DR(:,3:4,s)=DR(:,3:4,s)-R(:,3:4,s+1) 
          DR(:,1:2,s+1)=DR(:,1:2,s+1)-R(:,1:2,s)
        end do
      else
        do s=1,Ls-1
          DR(:,1:2,s)=DR(:,1:2,s)-R(:,1:2,s+1) 
          DR(:,3:4,s+1)=DR(:,3:4,s+1)-R(:,3:4,s)
        end do
      end if

c     mass terms
 
      if (MTYPE.eq.1) then
        if (.not.DAGGER) then
          DR(:,1:2,1)=DR(:,1:2,1)+mass*R(:,1:2,Ls)
          DR(:,3:4,Ls)=DR(:,3:4,Ls)+mass*R(:,3:4,1)
        else
          DR(:,3:4,1)=DR(:,3:4,1)+mass*R(:,3:4,Ls)
          DR(:,1:2,Ls)=DR(:,1:2,Ls)+mass*R(:,1:2,1)
        end if
      elseif (MTYPE.eq.4) then
        if (.not.DAGGER) then
          TR=zi*R(:,:,Ls)
          call mGmu(TR,4)
          DR(:,1:2,1)=DR(:,1:2,1)+mass*TR(:,1:2)
          TR=zi*R(:,:,1)
          call mGmu(TR,4)
          DR(:,3:4,Ls)=DR(:,3:4,Ls)+mass*TR(:,3:4)
        else
          TR=zi*R(:,:,Ls)
          call mGmu(TR,4)
          DR(:,3:4,1)=DR(:,3:4,1)-mass*TR(:,3:4)
          TR=zi*R(:,:,1)
          call mGmu(TR,4)
          DR(:,1:2,Ls)=DR(:,1:2,Ls)-mass*TR(:,1:2)
        endif
      elseif (MTYPE.eq.2) then
        if (.not.DAGGER) then
          DR(:,1:2,1)=DR(:,1:2,1)+zi*mass*R(:,1:2,Ls)
          DR(:,3:4,Ls)=DR(:,3:4,Ls)-zi*mass*R(:,3:4,1)
        else
          DR(:,3:4,1)=DR(:,3:4,1)+zi*mass*R(:,3:4,Ls)
          DR(:,1:2,Ls)=DR(:,1:2,Ls)-zi*mass*R(:,1:2,1)
        endif
      elseif (MTYPE.eq.3) then
        if (.not.DAGGER) then
          DR(:,1:2,1)=DR(:,1:2,1)-zi*mass*R(:,1:2,Ls)
          DR(:,3:4,Ls)=DR(:,3:4,Ls)+zi*mass*R(:,3:4,1)
        else
          DR(:,3:4,1)=DR(:,3:4,1)-zi*mass*R(:,3:4,Ls)
          DR(:,1:2,Ls)=DR(:,1:2,Ls)+zi*mass*R(:,1:2,1)
        endif
      endif

      return
      end subroutine DDW_Shamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ShamirDomainWallDerivs(dSdA,eta,nu,DAG)
      use numbers
!      use gammas
!      use indices
!      use WilsonDirac
      implicit none
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls),intent(in) ::  eta,nu
      logical,intent(in) :: DAG
      complex(prc) :: dSdAC(Nv,3)

      call ShamirDomainWallDerivsComplex(dSdAC,eta,nu,DAG)
      dSdA=dSdAC

      return
      end subroutine ShamirDomainWallDerivs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ShamirDomainWallDerivsComplex(dSdA,eta,nu,DAG)
      use numbers
!      use gammas
!      use indices
      use WilsonDirac
      implicit none
      complex(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls),intent(in) ::  eta,nu
      logical,intent(in) :: DAG
      complex(prc) :: tmp(Nv,3)
      complex(prc),dimension(Nv,4) ::  eta_l,nu_l
      integer l
      integer pm1

      dSdA=0
      do l=1,Ls  ! diagonal terms
        eta_l=eta(:,:,l)
        nu_l=nu(:,:,l)
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA+tmp
      end do

      return
      end subroutine ShamirDomainWallDerivsComplex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dslash(Phi,R,u,am,imass)
      use gammas
      use indices
      implicit none
      integer ksize,ksizet,kthird,kvol
      parameter(ksize=4,ksizet=4,kthird=8,kvol=ksizet*ksize*ksize)
      real(prc) akappa
      parameter(akappa=0.5)
      real(prc) beta,am3
      complex(prc) u(kvol,3)
      complex(prc) Phi(kvol,kthird,4),R(kvol,kthird,4)
      integer imass
      complex(prc) zkappa
      integer idirac,ithird,i,mu,igork
      real(prc) diag,am

      beta=gbeta
      am3=MDW
      am=baremass
c
c     diagonal term
      diag=(3.0-am3)+1.0
      do idirac=1,4
      do ithird=1,kthird
      do  i=1,kvol
      Phi(i,ithird,idirac)=diag*R(i,ithird,idirac)
      enddo
      enddo
      enddo
c
c     Wilson term
      do mu=1,3
      do idirac=1,4
      do ithird=1,kthird
      do i=1,kvol
      Phi(i,ithird,idirac)=Phi(i,ithird,idirac)
     &    -akappa*(      u(i,mu)*R(iu(i,mu),ithird,idirac)
     &         +conjg(u(id(i,mu),mu))*R(id(i,mu),ithird,idirac))
      enddo
      enddo
      enddo
      enddo
c
c     Dirac term
      do mu=1,3
      do idirac=1,4
      igork=gamin(mu,idirac)
      do ithird=1,kthird
      do i=1,kvol
      Phi(i,ithird,idirac)=Phi(i,ithird,idirac)
     &+half*gamval(mu,idirac)*
     &    (          u(i,mu)*R(iu(i,mu),ithird,igork)
     &         -conjg(u(id(i,mu),mu))*R(id(i,mu),ithird,igork))
      enddo
      enddo
      enddo
      enddo
c
c  s-like term exploiting projection
      do idirac=3,4
      do ithird=1,kthird-1
      do i=1,kvol
      Phi(i,ithird,idirac)=Phi(i,ithird,idirac)
     &   -R(i,ithird+1,idirac)
      enddo
      enddo
      enddo
      do idirac=1,2
      do ithird=1,kthird-1
      do i=1,kvol
      Phi(i,ithird+1,idirac)=Phi(i,ithird+1,idirac)
     &    -R(i,ithird,idirac)
      enddo
      enddo
      enddo
c
c  Mass term (couples the two walls unless imass=5)
      if(imass.eq.1)then
         zkappa=cmplx(am,0.0)
         do idirac=3,4
         do i=1,kvol
         Phi(i,kthird,idirac)=Phi(i,kthird,idirac)+zkappa*R(i,1,idirac)
         enddo
         enddo
         do idirac=1,2
         do i=1,kvol
         Phi(i,1,idirac)=Phi(i,1,idirac)+zkappa*R(i,kthird,idirac)
         enddo
         enddo
      elseif(imass.eq.3)then
         zkappa=cmplx(0.0,-am)
         do idirac=3,4
         do i=1,kvol
         Phi(i,kthird,idirac)=Phi(i,kthird,idirac)-zkappa*R(i,1,idirac)
         enddo
         enddo
         do idirac=1,2
         do i=1,kvol
         Phi(i,1,idirac)=Phi(i,1,idirac)+zkappa*R(i,kthird,idirac)
         enddo
         enddo
      elseif(imass.eq.5)then
         zkappa=cmplx(0.0,-am)
         do idirac=3,4
         do i=1,kvol
         Phi(i,kthird,idirac)=Phi(i,kthird,idirac)
     &                        -zkappa*R(i,kthird,idirac-2)
         enddo
         enddo
         do idirac=1,2
         do i=1,kvol
         Phi(i,1,idirac)=Phi(i,1,idirac)-zkappa*R(i,1,idirac+2)
         enddo
         enddo
      endif
c
      return
      end
c***********************************************************************
      subroutine dslashd(Phi,R,u,am,imass)
      use gammas
      use indices
      implicit none
      integer ksize,ksizet,kthird,kvol
      parameter(ksize=4,ksizet=4,kthird=8,kvol=ksizet*ksize*ksize)
      real(prc) akappa
      parameter(akappa=0.5)
      real(prc) beta,am3
      complex(prc) u(kvol,3)
      complex(prc) Phi(kvol,kthird,4),R(kvol,kthird,4)
      integer imass
      real(prc) am

      call dslashdA(Phi,R,u,am,imass)
      call dslashdB(Phi,R,u,am,imass)
 
      return
      end
c***********************************************************************
      subroutine dslashdA(Phi,R,u,am,imass)
      use gammas
      use indices
      implicit none
      integer ksize,ksizet,kthird,kvol
      parameter(ksize=4,ksizet=4,kthird=8,kvol=ksizet*ksize*ksize)
      real(prc) akappa
      parameter(akappa=0.5)
      real(prc) beta,am3
      complex(prc) u(kvol,3)
      complex(prc) Phi(kvol,kthird,4),R(kvol,kthird,4)
      integer imass
      complex(prc) zkappa
      integer idirac,ithird,i,mu,igork
      real(prc) diag,am

      beta=gbeta
      am3=MDW
      am=baremass
c
c     diagonal term (hermitian)
      diag=(3.0-am3)+1.0
      do idirac=1,4
      do ithird=1,kthird
      do  i=1,kvol
      Phi(i,ithird,idirac)=diag*R(i,ithird,idirac)
      enddo
      enddo
      enddo

c     Wilson term (hermitian)
      do mu=1,3
      do idirac=1,4
      do ithird=1,kthird
      do i=1,kvol
      Phi(i,ithird,idirac)=Phi(i,ithird,idirac)
     &    -akappa*(      u(i,mu)*R(iu(i,mu),ithird,idirac) 
     &         +conjg(u(id(i,mu),mu))*R(id(i,mu),ithird,idirac))
      enddo
      enddo
      enddo
      enddo
c
c     Dirac term (antihermitian)
      do mu=1,3
      do idirac=1,4
      igork=gamin(mu,idirac)
      do ithird=1,kthird
      do i=1,kvol
      Phi(i,ithird,idirac)=Phi(i,ithird,idirac)
     &-half*gamval(mu,idirac)*
     &    (          u(i,mu)*R(iu(i,mu),ithird,igork)
     &         -conjg(u(id(i,mu),mu))*R(id(i,mu),ithird,igork))
      enddo
      enddo
      enddo
      enddo

      return
      end
c***********************************************************************
      subroutine dslashdB(Phi,R,u,am,imass)
      use gammas
      use indices
      implicit none
      integer ksize,ksizet,kthird,kvol
      parameter(ksize=4,ksizet=4,kthird=8,kvol=ksizet*ksize*ksize)
      real(prc) akappa
      parameter(akappa=0.5)
      real(prc) beta,am3
      complex(prc) u(kvol,3)
      complex(prc) Phi(kvol,kthird,4),R(kvol,kthird,4)
      integer imass
      complex(prc) zkappa
      integer idirac,ithird,i,mu,igork
      real(prc) diag,am

      beta=gbeta
      am3=MDW
      am=baremass

c
c  s-like term exploiting projection
      do idirac=1,2
      do ithird=1,kthird-1
      do i=1,kvol
      Phi(i,ithird,idirac)=Phi(i,ithird,idirac)
     &   -R(i,ithird+1,idirac)
      enddo
      enddo
      enddo
      do idirac=3,4
      do ithird=1,kthird-1
      do i=1,kvol
      Phi(i,ithird+1,idirac)=Phi(i,ithird+1,idirac)
     &   -R(i,ithird,idirac)
      enddo
      enddo
      enddo
c
c  Mass term (couples the two walls unless imass=5) 
      if(imass.eq.1)then
         zkappa=cmplx(am,0.0)
         do idirac=1,2
         do i=1,kvol
         Phi(i,kthird,idirac)=Phi(i,kthird,idirac)+zkappa*R(i,1,idirac)
         enddo
         enddo
         do idirac=3,4
         do i=1,kvol
         Phi(i,1,idirac)=Phi(i,1,idirac)+zkappa*R(i,kthird,idirac)
         enddo
         enddo
      elseif(imass.eq.3)then
         zkappa=cmplx(0.0,am)
         do idirac=1,2
         do i=1,kvol
         Phi(i,kthird,idirac)=Phi(i,kthird,idirac)+zkappa*R(i,1,idirac)
         enddo
         enddo
         do idirac=3,4
         do i=1,kvol
         Phi(i,1,idirac)=Phi(i,1,idirac)-zkappa*R(i,kthird,idirac)
         enddo
         enddo
      elseif(imass.eq.5)then
         zkappa=cmplx(0.0,am)
         do idirac=1,2
         do i=1,kvol
         Phi(i,kthird,idirac)=Phi(i,kthird,idirac)
     &                        -zkappa*R(i,kthird,idirac+2)
         enddo
         enddo
         do idirac=3,4
         do i=1,kvol
         Phi(i,1,idirac)=Phi(i,1,idirac)-zkappa*R(i,1,idirac-2)
         enddo
         enddo
      endif
c
      return
      end
c***********************************************************************
      end module ShamirDomWall
