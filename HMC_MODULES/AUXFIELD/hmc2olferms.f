!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module hmc2olferms
      use pacc
      use arraysizes
      use numbers
      use indices
      use rvmodule
      use options
      implicit none
      logical,parameter :: VB_OLF=.true.
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setPseudoFermField(ps,ut) 
      use domainwallmod       
      use gaugefield                  
      implicit none
      complex(prc),intent(out) :: ps(Nv,4,Ls)
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc) :: ATMP(Nv,4,Ls),BTMP(Nv,4,Ls)
      complex(prc) :: rvs(Nv,4)
      integer l,baseMTYPE

      baseMTYPE=MTYPE

      ATMP=0d0
      call setCGRVs(4*Nv,rvs) ! randomise pseudo-fermion field
      ATMP(:,:,1)=rvs
!      write(151,*) phi3

      ! ps = K^dag.rv = Pdag.Ddag(m)^-1.Ddag(1).P.rv
      call PermM(ATMP,BTMP,.false.,4)
      MTYPE=1
      call IDDW(BTMP,ATMP,ut,.true.,one)
      MTYPE=baseMTYPE
      call DDW(ATMP,BTMP,ut,.true.,baremass)
      call PermM(BTMP,ps,.true.,4)
!      write(152,*) ps
      return
      end subroutine setPseudoFermField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine march2OLFerms(dH,thetat) ! march dAdt=P (theta is A)
      use gaugefield                      !       dPdt=-dSdA
      implicit none
      real(prc),intent(out) :: dH
      real(prc),intent(inout) :: thetat(Nv,3)
      real(prc) F(Nv,3)
      real(prc) pp(Nv,3)
      complex(prc) ps(Nv,4,Ls),ut(Nv,3),tmp(Nv,4,Ls)
      complex(prc) phi3(Nv,4)
      real(prc) etime,proby,ytest,avsteps,dt,h0,h1
      integer mu,ts,tsmax,l

      if (VB_OLF) then ; print *,"march 2 OL fermions" ; endif

      etime=HMC_etime ; dt=HMC_dt ; tsmax=HMC_tsmax

      avsteps=etime/dt
      proby=one/avsteps

      call coef(ut,thetat)
      call setGRVs(3*Nv,pp)  ! randomise starting momentum
      call setPseudoFermField(ps,ut) ! randomise pseudo-fermion field

!      write(108,*) sum(abs(pp*pp))/(Nv*3),sum(abs(ps*ps))/(Nv*Ls*4)

      h0=ham2OLFerms(thetat,ut,pp,ps) ! initial hamiltonian energy
      if (VB_OLF) then ; print *,"h0:",h0 ; end if
      call force2OLFerms(thetat,ut,ps,F)
!      write(223,*) F
      if (VB_OLF) then ; print *,"sum(F):",sum(F),sum(F*F) ; endif
      pp=pp-dt*half*F ! half time step before leap frog
      do ts=1,tsmax  ! time march
        if (VB_OLF) then ; print *,ts ; endif
        thetat=thetat+dt*pp
!        if(VB_DWF) print *,"sum(thetat):",sum(thetat),sum(thetat*thetat)
        call coef(ut,thetat)
        call force2OLFerms(thetat,ut,ps,F)
        if (VB_OLF) print *,"sum(F):",sum(F),sum(F*F)
        ytest=urv()
!        if(VB_DWF) then ; print *,"ytest:",ytest,"prob:",proby ; endif
        if (ytest.lt.proby) then
          pp=pp-dt*half*F
          goto 501
        else
          pp=pp-dt*F
        endif
      end do
!      pp=pp+half*dt*F ! correction to half step if tsmax reached

501   continue
      h1=ham2OLFerms(thetat,ut,pp,ps) ! final hamiltonian energy
      if (VB_OLF) then ; print *,"h1:",h1 ; end if
!      stop
      dH=h0-h1
      return
      end subroutine march2OLFerms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine force2OLFerms(thetat,ut,ps,F)
      implicit none
      real(prc),intent(in) :: thetat(Nv,3)
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: F(Nv,3)
      real(prc) :: dSdA(Nv,3)
      real(prc) :: dSdAT(Nv,3)

      call forceThirring3(thetat,F)
      if (QUENCHED) then ; return ; endif
      call fermionforceOL(ut,ps,dSdA)
      F=F+dSdA

      return
      end subroutine force2OLFerms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine forceThirring3(thetat,dSdA)
      implicit none
      real(prc) thetat(Nv,3),dSdA(Nv,3)

      dSdA=Nferms*gbeta*thetat

      return
      end subroutine forceThirring3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforceOL(ut,ps,dSdA) 
      use gammas
      use WilsonDirac
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,3) :: dSdA1,dSdA2
      complex(prc),dimension(Nv,4,Ls) :: chi,eta
      complex(prc),dimension(Nv,4,Ls) :: TMP,TMP2
      complex(prc),dimension(Nv,4,Ls) :: rhs1,rhs2
      integer :: KTYPE
      integer baseMTYPE,i,ioffset,ind 

      KTYPE=kerneltype()
      baseMTYPE=MTYPE
      
      call PermM(ps,chi,.false.,4)
      call IDDW(chi,eta,ut,.true.,baremass) 

      MTYPE=1
      call DDW(eta,TMP,ut,.true.,one) 
      call DDW(TMP,TMP2,ut,.false.,one) 
      MTYPE=baseMTYPE
      call IDDW(TMP2,rhs1,ut,.false.,baremass) 
      call DomainWallDerivsComplex(dSdA1,eta,rhs1,.false., 
     &                                             KTYPE,baremass)

      MTYPE=1
      call DDW(eta,rhs2,ut,.true.,one) 
      call DomainWallDerivsComplex(dSdA2,eta,rhs2,.false.,KTYPE,one)

      dSdA=2*(-dSdA1+dSdA2)

!     anti-p.b.c. in timelike direction
      ioffset=(Nt-1)*Ns*Ns
      do i=1,Ns*Ns
        ind=ioffset+i
        dSdA(ind,3)=-dSdA(ind,3)
      enddo
      MTYPE=baseMTYPE
      return
      end subroutine fermionforceOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham2OLFerms(thetat,ut,pp,ps)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat,pp
      complex(prc),dimension(Nv,3),intent(in) :: ut
      complex(prc),dimension(Nv,4,Ls),intent(in) :: ps
      real(prc) hg,hp,hf

      hp=0.5*sum(pp*pp)
      hg=hamThirring(thetat)
      hf=0
      if (.not.QUENCHED) then
        hf=ham2Ferms(ps,ut)
      endif
      ham2OLFerms=(hg+hp+hf)
      if (VB_OLF) print *,"hg:",hg
      if (VB_OLF) print *,"hp:",hp
      if (VB_OLF) print *,"hf:",hf
      if (VB_OLF) print *,"h:",ham2OLFerms
!      write(101,'(4F12.4)') hg/Nv,hp/Nv,hf/Nv,sum(thetat*thetat)/Nv
!      if (FOUT) then
!        write(101,'(4F12.4)') hg/Nv,hp/Nv,hf/Nv,sum(thetat*thetat)/Nv
!      endif
      return
      end function ham2OLFerms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function hamThirring(thetat)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat

      hamThirring=Nferms*half*gbeta*sum(thetat*thetat)

      return
      end function hamThirring
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham2Ferms(ps,ut)
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ps(Nv*4*Ls)
      complex(prc),intent(in) :: ut(Nv*3)
      complex(prc),dimension(Nv*4*Ls) :: tmp,rhs

      if (VB_OLF) print *,'energy 2 OL Fermions'
      call IKDDW(ps,tmp,ut,.true.,baremass)
      call IKDDW(tmp,rhs,ut,.false.,baremass)
      ham2Ferms=dot_product(ps,rhs)

      return
      end function ham2Ferms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer function kerneltype()
      implicit none
      if (DWkernel.eq.1) then
        kerneltype=1
      elseif ((DWkernel.eq.2).or.(DWkernel.eq.3)) then
        kerneltype=2
      endif
      return
      end function kerneltype
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module hmc2olferms
