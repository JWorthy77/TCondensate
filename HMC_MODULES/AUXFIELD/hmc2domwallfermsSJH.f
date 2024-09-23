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
      function transToHere(There)
      implicit none
      complex(prc) :: transToHere(Nv,4,Ls)
      complex(prc) :: There(Nv,Ls,4)
      integer ii,id,il

      do ii=1,Nv
        do il=1,Ls
          do id=1,4
            transToHere(ii,id,il)=There(ii,il,id)
          end do
        end do
      end do

      return
      end function transToHere
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function transToThere(Here)
      implicit none
      complex(prc) :: Here(Nv,4,Ls)
      complex(prc) :: transToThere(Nv,Ls,4)
      integer ii,id,il

      do ii=1,Nv
        do il=1,Ls
          do id=1,4
            transToThere(ii,il,id)=Here(ii,id,il)
          end do
        end do
      end do

      return
      end function transToThere
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setPseudoFermField(ps,ut) 
!     Pseudofermion fields: Phi =  Mdagger(1)^-1 * Mdagger(m) * R, where R is complex gaussian
!                               = M(1) * (Mdagger(1)M(1))^-1 * Mdagger(m) * R
      use domainwallmod       
      use gaugefield                  
      implicit none
      complex(prc),intent(out) :: ps(Nv,4,Ls)
!      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc) :: ut(Nv,3)
      complex :: uTMP(Nv,3)
      real :: thetaTMP(Nv,3)
      real(prc) :: thetaTMP2(Nv,3)
      complex(prc) :: TMP(Nv,4,Ls),TMP2(Nv,4,Ls)
      complex(prc) :: rv(Nv,Ls,4),pff(Nv,4,Ls)
      integer ii,id,il

!      call setCGRVs(4*Nv*Ls,TMP) ! randomise pseudo-fermion field
      if (DUPLICATE) then
        TMP=ps
      else
        call setCGRVs(4*Nv*Ls,TMP) ! randomise pseudo-fermion field
      endif

!      open(unit=110,file='IN/fort.110',status='unknown')
!      read(110,*) rv
!      close(110)
!      TMP=transToHere(rv)
 
      call DDW(TMP,ps,ut,.true.,baremass)
!      write(121,*) transToThere(ps)
!      call dslashd(TMP2,rv,ut,baremass,3)
!      write(131,*) TMP2
      call IDDWdagDDW(ps,TMP,ut,one)
!      write(122,*) transToThere(TMP)
      call DDW(TMP,ps,ut,.false.,one)   ! SJH sets to 1?
!      write(123,*) transToThere(ps)

!     ps=DDW(1).DDW(1)^-1.DDWdag(1)^-1.DDWdag(m)

!      write(108,*) sum(abs(rv*rv))/(Nv*Ls*4),sum(abs(ps*ps))/(Nv*Ls*4)

!      call setCGRVs(4*Nv,phi3) ! uses this in RHMC
!      do l=1,Ls
!        ps(:,:,l)=phi3
!      end do


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
!      open(unit=112,file='IN/fort.112',status='unknown')
!      read(112,*) thetat,ut
!      close(112)

      call setGRVs(3*Nv,pp)  ! randomise starting momentum
!      open(unit=113,file='IN/fort.113',status='unknown')
!      read(113,*) pp
!      close(113)
      if (DUPLICATE) then
        read(12,*) l,tmp,pp
        ps=transToHere(tmp)
      endif
      call setPseudoFermField(ps,ut) ! randomise pseudo-fermion field

      write(108,*) sum(abs(pp*pp))/(Nv*3),sum(abs(ps*ps))/(Nv*Ls*4)

      h0=ham2DomWallFerms(thetat,ut,pp,ps) ! initial hamiltonian energy
      if (VB_DWF) then ; print *,"h0:",h0 ; end if
      call force2DomWallFerms(thetat,ut,ps,F)
!      write(223,*) F
      print *,"sum(F):",sum(F),sum(F*F)
      pp=pp-dt*half*F ! half time step before leap frog
      if (DUPLICATE) tsmax=10
      do ts=1,tsmax  ! time march
        print *,ts
        if (DUPLICATE) h1=ham2DomWallFerms(thetat,ut,pp,ps) ! initial hamiltonian energy
        thetat=thetat+dt*pp
        print *,"sum(thetat):",sum(thetat),sum(thetat*thetat)
        call coef(ut,thetat)
        call force2DomWallFerms(thetat,ut,ps,F)
        print *,"sum(F):",sum(F),sum(F*F)
        ytest=urv()
        if (DUPLICATE) ytest=2
        print *,"ytest:",ytest,"prob:",proby
        if (ytest.lt.proby) then
!          print *,"ts:",ts
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
!      stop
      dH=h0-h1
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
      write(221,*) F
      print *,"sum(forceAux):",sum(F),sum(F*F)
      if (QUENCHED) then
        return
      endif
!      call fermionforce(ut,ps,dSdA)
      call fermionforceSJH(ut,ps,dSdA)
      if (.false.) then 
        call fermionforceComplex(ut,ps,dSdAT)
        print *,"force diff:",maxval(abs(dSdAT-dSdA)) 
      endif
      print *,"sum(forceFerm):",sum(dSdA),sum(dSdA*dSdA)
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
      subroutine fermionforceSJH(ut,ps,dSdA) ! for Seff=1/2phi ...  phi
      use gammas
      use WilsonDirac
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      real(prc),dimension(Nv,3) :: dSdA1,dSdA2,dSdA3,dSdAsub
      real(prc),dimension(Nv,3) :: dSdA4,dSdA5,dSdA6
      real(prc) :: dTmpdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls) :: eta,nu,chi
      complex(prc),dimension(Nv,4,Ls) :: TMP,TMP1,R,X1,X2
      complex(prc),dimension(Nv,4) :: X1s,Rs
      integer :: KTYPE
      procedure(),pointer :: Dptr=>NULL()
      integer baseMTYPE 
      complex(prc) tzi
      integer igork1,mu,idirac,i,ithird
      complex(prc),dimension(4) :: lsp1,rsp1,lsp2,rsp2     


!      if (VB_DWF) print *,"fermion force"
      if (DWkernel.eq.1) then
        KTYPE=1
      elseif ((DWkernel.eq.2).or.(DWkernel.eq.3)) then
        KTYPE=2
      endif

      dSdA=0

!      R=Ddag(1).Phi
      call DDW(ps,R,ut,.true.,one)
!      X1=IDdagD(m).R
      call IDDWdagDDW(R,X1,ut,baremass)
!      X2=D(m).X1
      call DDW(X1,X2,ut,.false.,baremass)
!      R=X2-Phi
      R=X2-ps
!      write(224,*) transToThere(R)
!      print *,"sum(R):",sum(R),sum(R*R)
!      write(225,*) transToThere(X1)
!      print *,"sum(X1):",sum(X1),sum(X1*X1)

      dSdA1=0
c     dSdpi=dSdpi-Re(X1dagger *(d(Mdagger)dp)* R) 
      tzi=(0.0,2.0)
      do mu=1,3
        do idirac=1,4
          igork1=gamin(mu,idirac)
          do ithird=1,Ls
            do i=1,Nv
              dSdA1(i,mu)=dSdA1(i,mu)+half*real( tzi*
     &               (conjg(X1(i,idirac,ithird))*
     &                        R(iu(i,mu),idirac,ithird)
     &               -conjg(X1(iu(i,mu),idirac,ithird))*
     &                        R(i,idirac,ithird) )  )
              dSdA1(i,mu)=dSdA1(i,mu)+real(tzi*gamval(mu,idirac)*half*
     &               (conjg(X1(i,idirac,ithird))*
     &                        R(iu(i,mu),igork1,ithird)
     &               +conjg(X1(iu(i,mu),idirac,ithird))*
     &                        R(i,igork1,ithird) )  )
            enddo
          enddo
        enddo
      enddo

!      print *,"dSdA2"
      dSdA2=0
      tzi=(0.0,2.0)
      do mu=1,3
        do ithird=1,Ls
          do i=1,Nv
            lsp1=X1(i,:,ithird)
            rsp1=R(iu(i,mu),:,ithird)
            lsp2=X1(iu(i,mu),:,ithird)
            rsp2=R(i,:,ithird) 
            dSdA2(i,mu)=dSdA2(i,mu)+half*real( tzi*
     &             ( dot_product(lsp1,rsp1)
     &              -dot_product(lsp2,rsp2) ) )

            lsp1=X1(i,:,ithird)
            rsp1=R(iu(i,mu),:,ithird)
            call mGmu4(rsp1,mu)
            lsp2=X1(iu(i,mu),:,ithird)
            rsp2=R(i,:,ithird)
            call mGmu4(rsp2,mu)
            dSdA2(i,mu)=dSdA2(i,mu)+real(tzi*half*
     &             ( dot_product(lsp1,rsp1)
     &              +dot_product(lsp2,rsp2) ) )
          enddo
        enddo
      enddo

!      print *,"dSdA3"
      dSdA3=0
      tzi=(0.0,2.0)
      do ithird=1,Ls
        X1s=X1(:,:,ithird)
        Rs=R(:,:,ithird)
!        call WilsonDerivsSJH(dSdASub,X1s,Rs,.true.) 
        dSdAsub=2*dSdAsub
        dSdA3=dSdA3+dSdAsub
      end do

!      print *,"dSdA4"
      dSdA4=0
      tzi=(0.0,2.0)
      do ithird=1,Ls
        X1s=X1(:,:,ithird)
        Rs=R(:,:,ithird)
!        call WilsonDerivsJW(dSdASub,X1s,Rs,.true.) 
        dSdAsub=2*dSdAsub
        dSdA4=dSdA4+dSdAsub
      end do

!      print *,"dSdA5"
      dSdA5=0
      tzi=(0.0,2.0)
      do ithird=1,Ls
        X1s=X1(:,:,ithird)
        Rs=R(:,:,ithird)
        call WilsonDerivs(dSdASub,X1s,Rs,.true.) 
        dSdAsub=2*dSdAsub
        dSdA5=dSdA5+dSdAsub
      end do

      call DomainWallDerivs(dSdA6,X1,R,.true.,KTYPE,zero)
!      dSdA6=-2*dSdA6
      dSdA6=2*dSdA6
      print *,"sum(ferm dSdpi):",sum(dSdA1),sum(dSdA1*dSdA1)

      dSdA=-dSdA1
!      write(226,*) dSdA1
!      write(227,*) dSdA2
!      print *,"sum(ferm dSdA):",sum(dSdA1),sum(dSdA1*dSdA1)

!      print *,"diff 2-1:",maxval(abs(dSdA2-dSdA1))
!      print *,"diff 3-1:",maxval(abs(dSdA3-dSdA1))
!      print *,"diff 4-5:",maxval(abs(dSdA4+dSdA5))
!      print *,"diff 4-1:",maxval(abs(dSdA1-dSdA4))
!      print *,"diff 6-1:",maxval(abs(dSdA6-dSdA1))

!       stop
!      baseMTYPE=MTYPE
!      ! 2 Re [dOdagdA (Qdag.Q)^-1.O]
!      eta=ps ; 
!      MTYPE=1;
!      call DDW(ps,TMP,ut,.false.,one)
!      MTYPE=baseMTYPE;
!      call IDDWdagDDW(TMP,nu,ut,baremass)
!      call DomainWallDerivs(dSdA1,eta,nu,.true.,KTYPE)

!      ! 2.Re[Odag.(QdagQ)^-1.dQdagdA.Q.(QdagQ)^-1.O]
!      MTYPE=1;
!      call DDW(ps,TMP,ut,.false.,one)
!      MTYPE=baseMTYPE;
!      call IDDWdagDDW(TMP,chi,ut,baremass)
!      eta=chi
!      call DDW(chi,nu,ut,.false.,baremass)
!      call DomainWallDerivs(dSdA2,eta,nu,.true.,KTYPE)

!      dSdA=(dSdA1+dSdA2)/Ls ! *2*half
!      MTYPE=baseMTYPE

      return
      end subroutine fermionforceSJH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforceComplex(ut,ps,dSdA) ! for Seff=1/2phi ...  phi
      use gammas
      use WilsonDirac
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,3) :: dSdA1,dSdA2,dSdA3,dSdA4
      real(prc) :: dTmpdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls) :: eta,nu,chi
      complex(prc),dimension(Nv,4,Ls) :: TMP,TMP1
      integer :: KTYPE
      procedure(),pointer :: Dptr=>NULL()
      integer baseMTYPE 
     
      if (VB_DWF) print *,"fermion force Complex"
      if (DWkernel.eq.1) then
        KTYPE=1
      elseif ((DWkernel.eq.2).or.(DWkernel.eq.3)) then
        KTYPE=2
      endif

      baseMTYPE=MTYPE
      ! dOdagdA (Qdag.Q)^-1 O
      eta=ps ; 
      MTYPE=1;
      call DDW(ps,TMP,ut,.false.,one)
      MTYPE=baseMTYPE;
      call IDDWdagDDW(TMP,nu,ut,baremass)
      call DomainWallDerivsComplex(dSdA1,eta,nu,.true.,KTYPE,zero)
      ! Odag.(Q^dag.Q)^-1. dOdA
      nu=ps
      MTYPE=1;
      call DDW(ps,TMP,ut,.false.,one)
      MTYPE=baseMTYPE;
      call IDDWdagDDW(TMP,eta,ut,baremass)
      call DomainWallDerivsComplex(dSdA2,eta,nu,.false.,KTYPE,zero)

      ! Odag.(QdagQ)^-1.dQdagdA.Q.(QdagQ)^-1.O
      MTYPE=1;
      call DDW(ps,TMP,ut,.false.,one)
      MTYPE=baseMTYPE;
      call IDDWdagDDW(TMP,chi,ut,baremass)
      eta=chi
      call DDW(chi,nu,ut,.false.,baremass)
      call DomainWallDerivsComplex(dSdA3,eta,nu,.true.,KTYPE,zero)

      ! Odag.(QdagQ)^-1.Qdagd.dQdA.(QdagQ)^-1.O
      MTYPE=1;
      call DDW(ps,TMP,ut,.false.,one)
      MTYPE=baseMTYPE;
      call IDDWdagDDW(TMP,chi,ut,baremass)
      nu=chi
      call DDW(chi,eta,ut,.false.,baremass)
      call DomainWallDerivsComplex(dSdA4,eta,nu,.false.,KTYPE,zero)

!      print *,"cmplx 0:",maxval(abs(aimag(dSdA1+dSdA2)))
!      print *,"cmplx 0:",maxval(abs(aimag(dSdA1+dSdA2+dSdA3+dSdA4)))
      dSdA=half*(dSdA1+dSdA2+dSdA3+dSdA4)/Ls
      MTYPE=baseMTYPE

      if (VB_DWF) print *,"fermion force Complex done"

      return
      end subroutine fermionforceComplex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforce(ut,ps,dSdA) ! for Seff=1/2phi ...  phi
      use gammas
      use WilsonDirac
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      real(prc),dimension(Nv,3) :: dSdA1,dSdA2
      real(prc) :: dTmpdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls) :: eta,nu,chi
      complex(prc),dimension(Nv,4,Ls) :: TMP,TMP1
      integer :: KTYPE
      procedure(),pointer :: Dptr=>NULL()
      integer baseMTYPE 
     
      if (VB_DWF) print *,"fermion force"
      if (DWkernel.eq.1) then
        KTYPE=1
      elseif ((DWkernel.eq.2).or.(DWkernel.eq.3)) then
        KTYPE=2
      endif

      baseMTYPE=MTYPE
      ! 2 Re [dOdagdA (Qdag.Q)^-1.O]
      eta=ps ; 
!      MTYPE=1;
      call DDW(ps,TMP,ut,.false.,one)
      MTYPE=baseMTYPE;
      call IDDWdagDDW(TMP,nu,ut,baremass)
      call DomainWallDerivs(dSdA1,eta,nu,.true.,KTYPE,zero)

      ! 2.Re[Odag.(QdagQ)^-1.dQdagdA.Q.(QdagQ)^-1.O]
!      MTYPE=1;
      call DDW(ps,TMP,ut,.false.,one)
      MTYPE=baseMTYPE;
      call IDDWdagDDW(TMP,chi,ut,baremass)
      eta=chi
      call DDW(chi,nu,ut,.false.,baremass)
      call DomainWallDerivs(dSdA2,eta,nu,.true.,KTYPE,zero)

      dSdA=(dSdA1+dSdA2)/Ls ! *2*half
      MTYPE=baseMTYPE

      return
      end subroutine fermionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      real(prc) function ham2DomWallFerms(thetat,ut,pp,ps,FOUT)
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
      complex(prc),dimension(Nv*4*Ls) :: tmp,lhs,rhs
      procedure(),pointer :: Dptr=>NULL()

      Dptr => DDW_OWilson
      Dptr => DDW_Wilson
      Dptr => MDomWall
      if (VB_DWF) print *,'energy 2 DomWall Fermions'
!      call IMdagM_DWkernel(ps,tmp,ut,baremass,Dptr)
!      ham2Ferms=half*dot_product(ps,tmp)

!     SJH
      call DDW(ps,lhs,ut,.true.,one)
!      write(131,*) transToThere(lhs)
      call IDDWdagDDW(lhs,rhs,ut,baremass)
      ham2Ferms=dot_product(lhs,rhs)

      return
      end function ham2Ferms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module hmc2domwallferms
