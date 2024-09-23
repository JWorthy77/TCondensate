!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module hmc2wilsonferms
      use pacc
      use arraysizes
      use numbers
      use indices
      use rvmodule
      use options
      use WilsonExtraMod
      implicit none
      logical,parameter :: VB_H2=.true.
!      real(prc) :: etime,dt
!      integer tsmax
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine march2DW(dH,thetat)
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

      h0=ham2DW(thetat,ut,pp,ps)
      if (VB_H2) then ; print *,"h0:",h0 ; end if
      call force2DW(thetat,ut,ps,F)
      pp=pp-dt*half*F ! half time step before leap frog
      proby=1.0/avsteps
      do ts=1,tsmax ! time march
        thetat=thetat+dt*pp
        call coef(ut,thetat)
        call force2DW(thetat,ut,ps,F)
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
      h1=ham2DW(thetat,ut,pp,ps)
      dH=h0-h1
      return
      end subroutine march2DW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine force2DW(thetat,ut,ps,F)
      implicit none
      real(prc),intent(in) :: thetat(Nv,3)
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4)
      real(prc),intent(out) :: F(Nv,3)
      real(prc) :: dSdA(Nv,3)
      real(prc) :: dSdAComplex(Nv,3)

      call forceThirring3(thetat,F)
      return
      call fermionforce(ut,ps,dSdA) 
      if (.false.) then
        call fermionforceComplex(ut,ps,dSdAComplex) 
        print *,"diff:",maxval(abs(dSdAComplex-dSdA))
      endif

!      F=F+half*dSdA
      F=F+dSdA
      return
      end subroutine force2DW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine forceThirring3(thetat,dSdA)
      implicit none
      real(prc) thetat(Nv,3),dSdA(Nv,3)

      dSdA=Nferms*gbeta*thetat

      return
      end subroutine forceThirring3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforce(ut,ps,dSdA) ! for Seff=1/2phi^dag (Dw^dag.Dw)^-1 phi
      use gammas                          ! dSdA = Re [ X^dag.Dw^dag.dDwdA.X ]
      use WilsonDirac                     ! X = (Dw^dag.Dw)^-1 phi
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4)
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4) :: eta,nu,etaT
      
      call IDdagD(ps,nu,ut,.false.,baremass)
      call DWilson(nu,eta,ut,.false.,baremass)
      if (.false.) then
        print *,"test:"
        call DWilsonJW(nu,etaT,ut,.false.,baremass)
        print *,"diff:",maxval(abs(etaT-eta))
      end if
      call WilsonDerivs(dSdA,eta,nu,.false.)

      return
      end subroutine fermionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforceComplex(ut,ps,dSdA) ! for Seff=1/2phi^dag (Dw^dag.Dw)^-1 phi
      use gammas                                 ! dSda = 1/2.X^dag.[dDw^dagdA.Dw + Dw^dag.dDwdA] X
      use WilsonDirac                            ! X = (Dw^dag.Dw)^-1 phi
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4)
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4) :: eta,nu
      complex(prc),dimension(Nv,3) :: dSdA1,dSdA2
      complex(prc),dimension(Nv,3) :: dSdA3
      
      call IDdagD(ps,nu,ut,.false.,baremass)
      call DWilson(nu,eta,ut,.false.,baremass)
      call WilsonDerivsComplex(dSdA1,eta,nu,.false.)
      eta=nu
      call DWilson(eta,nu,ut,.false.,baremass)
      call WilsonDerivsComplex(dSdA2,eta,nu,.true.)
      dSdA3=(dSdA1+dSdA2)/2.0
      dSdA=dSdA3

      return
      end subroutine fermionforceComplex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham2DW(thetat,ut,pp,ps)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat,pp
      complex(prc),dimension(Nv,3),intent(in) :: ut
      complex(prc),dimension(Nv,4),intent(in) :: ps
      real(prc) hg,hp,hf

      hp=0.5*sum(pp*pp)
      hg=hamThirring(thetat)
      hf=ham2WilsonFerms(ps,ut)
      ham2DW=(hg+hp+hf)/Nv
      if (VB_H2) print *,"hg:",hg/Nv
      if (VB_H2) print *,"hp:",hp/Nv
      if (VB_H2) print *,"hf:",hf/Nv
      if (VB_H2) print *,"h:",ham2DW
      write(101,*) hg/Nv,hp/Nv,hf/Nv
      return
      end              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function hamThirring(thetat)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat

      hamThirring=Nferms*half*gbeta*sum(thetat*thetat)

      return
      end function hamThirring
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham2WilsonFerms(ps,ut)
      use WilsonDirac
      implicit none
      complex(prc),intent(in) :: ps(Nv*4)
      complex(prc),intent(in) :: ut(Nv*3)
      complex(prc) tmp(Nv*4)

      print *,'energy 2 Wilson Fermions'
      call IDdagD(ps,tmp,ut,.false.,baremass)
      ham2WilsonFerms=half*dot_product(ps,tmp)

      return
      end function ham2WilsonFerms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module hmc2wilsonferms
