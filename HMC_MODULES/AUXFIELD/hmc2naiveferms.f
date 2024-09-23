!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module hmc2naiveferms
      use pacc
      use arraysizes
      use numbers
      use indices
      use rvmodule
      use options
      implicit none
      logical,parameter :: VB_HNF=.true.
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine march2DNF(dH,thetat)
      ! march dAdt=P       (theta is A)
      !       dPdt=-dSdA
      use gaugefield
      implicit none
      real(prc),intent(out) :: dH
      real(prc),intent(inout) :: thetat(Nv,3)
      real(prc) F(Nv,3)
      real(prc) pp(Nv,3)
      complex(prc) ps(Nv,4),ut(Nv,3)
      real(prc) proby,ytest,avsteps,h0,h1
      integer mu,ts
      real(prc) :: etime,dt
      integer tsmax

      etime=HMC_etime ; dt=HMC_dt ; tsmax=HMC_tsmax

      avsteps=etime/dt
      proby=one/avsteps

      call coef(ut,thetat)
      call setGRVs(3*Nv,pp)  ! randomise starting momentum
      call setCGRVs(4*Nv,ps) ! randomise pseudo-fermion field

      h0=ham2DNF(thetat,ut,pp,ps)
      if (VB_HNF) then ; print *,"h0:",h0 ; end if
      call force2DNF(thetat,ut,ps,F)
      pp=pp-dt*half*F ! half time step before leap frog
      proby=1.0/avsteps
      do ts=1,tsmax ! time march
        thetat=thetat+dt*pp
        call coef(ut,thetat)
        call force2DNF(thetat,ut,ps,F)
        print *,"sumF:",sum(F)
        ytest=urv()
        if (ytest.lt.proby) then
          print *,"ts:",ts
          pp=pp-dt*half*F
          goto 501
        else
          pp=pp-dt*F
        endif
      end do

501   continue
      h1=ham2DNF(thetat,ut,pp,ps)
      dH=h0-h1
      return
      end subroutine march2DNF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine force2DNF(thetat,ut,ps,F)
      implicit none
      real(prc),intent(in) :: thetat(Nv,3)
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4)
      real(prc),intent(out) :: F(Nv,3)
      real(prc) :: dSdA(Nv,3)
      real(prc) :: dSdAComplex(Nv,3)

      call forceThirring3(thetat,F)
      if(.not.QUENCHED)then
        call fermionforceDNF(ut,ps,dSdAComplex) 
        dSdA=dSdAComplex
        F=F+dSdA
      endif
      return
      end subroutine force2DNF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine forceThirring3(thetat,F)
      implicit none
      real(prc) thetat(Nv,3),F(Nv,3)

      F=Nferms*gbeta*thetat

      return
      end subroutine forceThirring3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforceReal(ut,ps,dSdA) ! for Seff=1/2phi^dag (Dw^dag.Dw)^-1 phi
      use gammas                          ! dSdA = Re [ X^dag.Dw^dag.dDwdA.X ]
      use NaiveDirac                     ! X = (Dw^dag.Dw)^-1 phi
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4)
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc) :: dSdAC(Nv,3)
      complex(prc),dimension(Nv,4) :: eta,nu
      
      call IDNFdagDNF(ps,nu,ut,.false.,baremass)
      call DNF(nu,eta,ut,.false.,baremass)
      call DNFDerivs(dSdAC,eta,nu,.false.)
      dSdA=dSdAC
      return
      end subroutine fermionforceReal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforceDNF(ut,ps,dSdA) ! for Seff=1/2phi^dag (Dw^dag.Dw)^-1 phi
      use gammas                                 ! dSda = 1/2.X^dag.[dDw^dagdA.Dw + Dw^dag.dDwdA] X
      use NaiveDirac                            ! X = (Dw^dag.Dw)^-1 phi
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4)
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4) :: eta,nu
      complex(prc),dimension(Nv,3) :: dSdA1,dSdA2
      complex(prc),dimension(Nv,3) :: dSdA3
      
      call IDNFdagDNF(ps,nu,ut,.false.,baremass)
      call DNF(nu,eta,ut,.false.,baremass)
      call DNFDerivsSJH(dSdA1,eta,nu,.false.)
      eta=nu
      call DNF(eta,nu,ut,.false.,baremass)
      call DNFDerivsSJH(dSdA2,eta,nu,.true.)
      dSdA3=(dSdA1+dSdA2)/2.0
      dSdA=dSdA3

      return
      end subroutine fermionforceDNF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham2DNF(thetat,ut,pp,ps)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat,pp
      complex(prc),dimension(Nv,3),intent(in) :: ut
      complex(prc),dimension(Nv,4),intent(in) :: ps
      real(prc) hg,hp,hf

      hp=0.5*sum(pp*pp)
      hg=hamThirring(thetat)
      ham2DNF=(hg+hp)/Nv
      if (VB_HNF) print *,"hg:",hg/Nv
      if (VB_HNF) print *,"hp:",hp/Nv
      if (.not.QUENCHED) then
        hf=ham2NaiveFerms(ps,ut)
        ham2DNF=ham2DNF+hf/Nv
        if (VB_HNF) print *,"hf:",hf/Nv
      endif
      if (VB_HNF) print *,"h:",ham2DNF
      write(101,*) hg/Nv,hp/Nv,hf/Nv
      return
      end function ham2DNF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function hamThirring(thetat)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat

      hamThirring=Nferms*half*gbeta*sum(thetat*thetat)

      return
      end function hamThirring
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham2NaiveFerms(ps,ut)
      use NaiveDirac
      implicit none
      complex(prc),intent(in) :: ps(Nv*4)
      complex(prc),intent(in) :: ut(Nv*3)
      complex(prc) tmp(Nv*4)

      print *,'energy 2 Naive Fermions'
      call IDNFdagDNF(ps,tmp,ut,.false.,baremass)
      ham2NaiveFerms=half*dot_product(ps,tmp)

      return
      end function ham2NaiveFerms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module hmc2naiveferms
