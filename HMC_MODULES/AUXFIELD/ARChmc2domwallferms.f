!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module hmc2domwallferms
      use pacc
      use arraysizes
      use numbers
      use indices
      use rvmodule
      use options
      implicit none
      logical,parameter :: VB_DWF=.true.
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setPseudoFermField(ps,ut) 
!     Pseudofermion fields: Phi =  Mdagger(1)^-1 * Mdagger(m) * R, where R is complex gaussian
!                               = M(1) * (Mdagger(1)M(1))^-1 * Mdagger(m) * R
      use domainwallmod       
      use gaugefield                  
      implicit none
      complex(prc),intent(out) :: ps(Nv,4,Ls)
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc) :: TMP(Nv,4,Ls),phi3(Nv,4)
      integer l

      TMP=0d0
!      call setCGRVs(4*Nv,phi3) ! randomise pseudo-fermion field
!      do l=1,Ls
      do l=1,1
        call setCGRVs(4*Nv,phi3) ! randomise pseudo-fermion field
        TMP(:,:,l)=phi3
      end do
      write(151,*) phi3

      ! ps = Q^dag.rv = Ddag(1)^-1.Ddag(m).rv
      call DDW(TMP,ps,ut,.true.,baremass)
      call IDDWdagDDW(ps,TMP,ut,one) 
      call DDW(TMP,ps,ut,.false.,one)  
      write(152,*) ps
!      stop

      return
      end subroutine setPseudoFermField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine march2DomWallFerms(dH,thetat) ! march dAdt=P (theta is A)
      use gaugefield                           !       dPdt=-dSdA
      implicit none
      real(prc),intent(out) :: dH
      real(prc),intent(inout) :: thetat(Nv,3)
      real(prc) F(Nv,3)
      real(prc) pp(Nv,3)
      complex(prc) ps(Nv,4,Ls),ut(Nv,3),tmp(Nv,4,Ls)
      complex(prc) phi3(Nv,4)
      real(prc) etime,proby,ytest,avsteps,dt,h0,h1
      integer mu,ts,tsmax,l

      etime=HMC_etime ; dt=HMC_dt ; tsmax=HMC_tsmax

      avsteps=etime/dt
      proby=one/avsteps

      call coef(ut,thetat)
      call setPseudoFermField(ps,ut) ! randomise pseudo-fermion field
      call setGRVs(3*Nv,pp)  ! randomise starting momentum
      write(153,*) pp

      h0=ham2DomWallFerms(thetat,ut,pp,ps) ! initial hamiltonian energy
      if (VB_DWF) then ; print *,"h0:",h0 ; end if
      call force2DomWallFerms(thetat,ut,ps,F)
      if (VB_DWF) then ; print *,"sum(F):",sum(F),sum(F*F) ; endif
      pp=pp-dt*half*F ! half time step before leap frog
      write(154,*) pp
      do ts=1,tsmax  ! time march
        if (VB_DWF) then ; print *,ts ; endif
        thetat=thetat+dt*pp
!        if(VB_DWF) print *,"sum(thetat):",sum(thetat),sum(thetat*thetat)
        call coef(ut,thetat)
        call force2DomWallFerms(thetat,ut,ps,F)
        if (VB_DWF) print *,"sum(F):",sum(F),sum(F*F)
        ytest=urv()
        if(VB_DWF) then ; print *,"ytest:",ytest,"prob:",proby ; endif
        if (ytest.lt.proby) then
          pp=pp-dt*half*F
          goto 501
        else
          pp=pp-dt*F
        endif
      end do
!      pp=pp+half*dt*F ! correction to half step if tsmax reached

501   continue
      h1=ham2DomWallFerms(thetat,ut,pp,ps) ! final hamiltonian energy
      if (VB_DWF) then ; print *,"h1:",h1 ; end if
      dH=h0-h1
      if (VB_DWF) then ; print *,"dH:",dH ; end if
      return
      end subroutine march2DomWallFerms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine force2DomWallFerms(thetat,ut,ps,F)
      implicit none
      real(prc),intent(in) :: thetat(Nv,3)
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: F(Nv,3)
      real(prc) :: dSdA(Nv,3)
      real(prc) :: dSdAT(Nv,3)

      call forceThirring3(thetat,F)
!      write(221,*) F
!      print *,"sum(forceAux):",sum(F),sum(F*F)
      if (QUENCHED) then
        return
      endif
      call fermionforceAlt(ut,ps,dSdA)
!      call fermionforceJW(ut,ps,dSdA)
!      call fermionforceJW2(ut,ps,dSdAT)
!      call fermionforceSJH(ut,ps,dSdA)

!      print *,"force err:",maxval(dSdAT-dSdA)
!      stop
!      print *,"sum(forceFerm):",sum(dSdA),sum(dSdA*dSdA)
!      write(222,*) dSdA
      F=F+dSdA
      return
      end subroutine force2DomWallFerms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine forceThirring3(thetat,dSdA)
      implicit none
      real(prc) thetat(Nv,3),dSdA(Nv,3)

      dSdA=Nferms*gbeta*thetat

      return
      end subroutine forceThirring3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforceAlt(ut,ps,dSdA) 
      use gammas
      use WilsonDirac
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      real(prc) :: dSdA1(Nv,3),dSdA2(Nv,3)
      complex(prc),dimension(Nv,4,Ls) :: TMP,R,QR
      integer KTYPE,baseMTYPE

      KTYPE=kerneltype()
      baseMTYPE=MTYPE
      dSdA=0

      MTYPE=1
      call DDW(ps,TMP,ut,.true.,one)
      MTYPE=baseMTYPE
      call IDDWdagDDW(TMP,R,ut,baremass) 
      call DDW(R,QR,ut,.false.,baremass)

!      eta=transToHere(R)
!      nu=transToHere(ps)
      call ShamirDomainWallDerivs(dSdA1,R,ps,.true.)

!      eta=transToHere(R)
!      nu=transToHere(QR)
      call ShamirDomainWallDerivs(dSdA2,R,QR,.true.)

      dSdA=2*(dSdA1-dSdA2)
      return
      end subroutine fermionforceAlt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforceJW(ut,ps,dSdA) 
      use gammas
      use WilsonDirac
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      real(prc) :: dSdA1(Nv,3),dSdA2(Nv,3)
      complex(prc),dimension(Nv,4,Ls) :: R,R1,X1,R2
      integer :: KTYPE
      integer baseMTYPE 
      integer ioffset,i,ind

      KTYPE=kerneltype()
      baseMTYPE=MTYPE

      MTYPE=1
      call DDW(ps,X1,ut,.true.,one) ! X1=Ddag(1).ps
      MTYPE=baseMTYPE
      call IDDWdagDDW(X1,R,ut,baremass) ! R=IDdagD(m).Ddag(1).ps

      call DDW(R,R2,ut,.false.,baremass) ! R2=D(m).R
      call DomainWallDerivs(dSdA2,R,R2,.true.,KTYPE,baremass)

      MTYPE=1
      call DomainWallDerivs(dSdA1,R,ps,.true.,KTYPE,one)

      dSdA=2*(dSdA1-dSdA2)

!     anti-p.b.c. in timelike direction
      ioffset=(Nt-1)*Ns*Ns
      do i=1,Ns*Ns
        ind=ioffset+i
        dSdA(ind,3)=-dSdA(ind,3)
      enddo
      MTYPE=baseMTYPE

      return
      end subroutine fermionforceJW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforceJW2(ut,ps,dSdA) 
      use gammas
      use WilsonDirac
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      real(prc),dimension(Nv,3) :: dSdA1,dSdA2,dSdA3,dSdA4
      complex(prc),dimension(Nv,4,Ls) :: R1,X1,R2
      integer :: KTYPE
      integer baseMTYPE 

      dSdA=0
      KTYPE=kerneltype()
      call DDW(ps,X1,ut,.true.,one) ! X1=Ddag(1).phi
      call IDDWdagDDW(X1,R1,ut,baremass) ! R1=IDdagD(m).Ddag(1).phi
      call DDW(R1,R2,ut,.false.,baremass) ! R2=D(m).R1

      call DomainWallDerivs(dSdA1,ps,R1,.false.,KTYPE,one)
      call DomainWallDerivs(dSdA2,R1,ps,.true.,KTYPE,one)

      call DomainWallDerivs(dSdA3,R1,R2,.true.,KTYPE,baremass)
      call DomainWallDerivs(dSdA4,R2,R1,.false.,KTYPE,baremass)

      dSdA=dSdA1+dSdA2-dSdA3-dSdA4
!      dSdA=-dSdA

      return
      end subroutine fermionforceJW2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforceSJH(ut,ps,dSdA) ! for Seff=1/2phi ...  phi
      use gammas                             ! for Shamir only
      use WilsonDirac
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls) :: R,X1,X2
      integer :: KTYPE
      integer baseMTYPE 

      dSdA=0
      call DDW(ps,R,ut,.true.,one) ! R=Ddag(1).Phi
      call IDDWdagDDW(R,X1,ut,baremass) ! X1=IDdagD(m).R
      call DDW(X1,X2,ut,.false.,baremass) ! X2=D(m).X1
      R=X2-ps  ! R=X2-Phi

c     dSdpi=dSdpi-Re(X1dagger *(d(Mdagger)dp)* R) 
      KTYPE=kerneltype()
      call DomainWallDerivs(dSdA,X1,R,.true.,KTYPE,zero)
      dSdA=2*dSdA

      return
      end subroutine fermionforceSJH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham2DomWallFerms(thetat,ut,pp,ps)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat,pp
      complex(prc),dimension(Nv,3),intent(in) :: ut
      complex(prc),dimension(Nv,4,Ls),intent(in) :: ps
!      logical FOUT
      real(prc) hg,hp,hf

      hp=0.5*sum(pp*pp)
      hg=hamThirring(thetat)
      hf=0
      if (.not.QUENCHED) then
        hf=ham2Ferms(ps,ut)
      endif
      ham2DomWallFerms=(hg+hp+hf)
      if (VB_DWF) print *,"hg:",hg
      if (VB_DWF) print *,"hp:",hp
      if (VB_DWF) print *,"hf:",hf
      if (VB_DWF) print *,"h:",ham2DomWallFerms
      write(101,'(4F12.4)') hg/Nv,hp/Nv,hf/Nv,sum(thetat*thetat)/Nv
!      if (FOUT) then
!        write(101,'(4F12.4)') hg/Nv,hp/Nv,hf/Nv,sum(thetat*thetat)/Nv
!      endif
      return
      end function ham2DomWallFerms
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
      complex(prc),dimension(Nv*4*Ls) :: lhs,rhs
      integer baseMTYPE 

      if (VB_DWF) print *,'energy 2 DomWall Fermions'
!     Seff = psdag.DDW(1).[DDWdag(m).DDW(m)]^-1.DDWdag(1).ps
      baseMTYPE=MTYPE
      MTYPE=1
      call DDW(ps,lhs,ut,.true.,one)
      MTYPE=baseMTYPE
      call IDDWdagDDW(lhs,rhs,ut,baremass)
      ham2Ferms=dot_product(lhs,rhs)

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
      end module hmc2domwallferms
