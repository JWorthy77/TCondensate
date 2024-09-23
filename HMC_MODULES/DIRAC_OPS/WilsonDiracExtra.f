      module WilsonExtraMod
      use pacc
      use arraysizes
      use numbers
      use WilsonDirac
      implicit none
      contains
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDW(RR,DR,u,DAGGER,mass)
!     solve Dw.DR = RR or Dw(dag).DR = RR
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DShamir(R,DR,u,DAGGER,mass)
      implicit none
c     calculates DR = DShamir*R = (2+Dw)^-1*Dw*R
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass,twoplusmass
      complex(prc) :: TMP(Nv,4)

      twoplusmass=two+mass
      call DWilson(R,TMP,u,DAGGER,mass)
      call IDW(TMP,DR,u,DAGGER,twoplusmass)

      return
      end subroutine DShamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDS(RR,DR,u,DAGGER,mass)
!     solve DS.DR = RR or DS(dag).DR = RR
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass,twoplusmass
      complex(prc) TMP(Nv,4)

      twoplusmass=two+mass
      call DWilson(RR,TMP,u,DAGGER,twoplusmass)
      call IDW(TMP,DR,u,DAGGER,mass)

      return
      end subroutine IDS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Dbasic(R,DR,u,DAGGER,mass)
      use options
      implicit none
c     chooses between Shamir and Wilson operator
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass

      if (dwkernel.eq.1) then
        call DShamir(R,DR,u,DAGGER,mass)
      elseif (dwkernel.eq.2) then
        call DWilson(R,DR,u,DAGGER,mass)
      end if
      return
      end subroutine Dbasic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDbasic(R,DR,u,DAGGER,mass)
      use options
      implicit none
c     chooses between Shamir and Wilson operator
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add

      if (dwkernel.eq.1) then
        call IDS(R,DR,u,DAGGER,mass)
      elseif (dwkernel.eq.2) then
        call IDW(R,DR,u,DAGGER,mass)
      end if
      return
      end subroutine IDbasic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Hkernel(R,DR,u,DAGGER,mass)
      use gammas
      implicit none
c     calculates Hw or Hs
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TMP(Nv,4)

      if(.not.DAGGER) then
        call Dbasic(R,DR,u,DAGGER,mass)
        call mGmu(DR,4)
      else
        TMP=R
        call mGmu(TMP,4)
        call Dbasic(TMP,DR,u,DAGGER,mass)
      endif
      return
      end subroutine Hkernel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IHkernel(R,DR,u,DAGGER,mass)
      use gammas
      implicit none
c     calculates Hw or Hs
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TMP(Nv,4)

      if(.not.DAGGER) then
        TMP=R
        call mGmu(TMP,4)
        call IDbasic(TMP,DR,u,DAGGER,mass)
      else
        call IDbasic(R,DR,u,DAGGER,mass)
        call mGmu(DR,4)
      endif
      return
      end subroutine IHkernel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module WilsonExtraMod
