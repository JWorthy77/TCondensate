!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module domainwallmod
      use arraysizes
      use options
      use WilsonDirac
      use ShamirDomWall
      use WilsonDomWall
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      subroutine DomWall(R,DR,u,DAG,mass)
!      implicit none
!      complex(prc),intent(in) :: R(Nv,4,Ls)
!      complex(prc),intent(out) :: DR(Nv,4,Ls)
!      complex(prc),intent(in) :: u(Nv,3)
!      logical,intent(in) :: DAG
!      real(prc),intent(in) :: mass
!      complex(prc) :: TMP(Nv,4,Ls)
!      integer baseMTYPE
!
!      baseMTYPE=MTYPE;
!      if (.not.DAG) then
!        call DDW(R,TMP,u,DAG,mass)
!        MTYPE=1
!        call IDDW(TMP,DR,u,DAG,one)
!        MTYPE=baseMTYPE
!      else 
!        MTYPE=1
!        call IDDW(R,TMP,u,DAG,one)
!        MTYPE=baseMTYPE
!        call DDW(TMP,DR,u,DAG,mass)
!      end if
!
!      return
!      end subroutine DomWall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine MDomWall(R,DR,u,DAG,mass)
      implicit none
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAG
      real(prc),intent(in) :: mass
      complex(prc) :: TMP(Nv,4,Ls)
      integer baseMTYPE

      baseMTYPE=MTYPE;
      if (.not.DAG) then
        call DDW(R,TMP,u,DAG,mass)
        MTYPE=1
        call IDDW(TMP,DR,u,DAG,one)
        MTYPE=baseMTYPE
      else 
        MTYPE=1
        call IDDW(R,TMP,u,DAG,one)
        MTYPE=baseMTYPE
        call DDW(TMP,DR,u,DAG,mass)
      end if

      return
      end subroutine MDomWall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMDomWall(R,DR,u,DAG,mass)
      implicit none
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAG
      real(prc),intent(in) :: mass
      complex(prc) :: TMP(Nv,4,Ls)
      integer baseMTYPE

      baseMTYPE=MTYPE;
      if (.not.DAG) then
        MTYPE=1
        call DDW(R,TMP,u,DAG,one)
        MTYPE=baseMTYPE
        call IDDW(TMP,DR,u,DAG,mass)
      else 
        call IDDW(R,TMP,u,DAG,mass)
        MTYPE=1
        call DDW(TMP,DR,u,DAG,one)
        MTYPE=baseMTYPE
      end if

      return
      end subroutine IMDomWall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testMDomWall(u)
      use rvmodule
      implicit none
      complex(prc),intent(in) :: u(Nv,3)
      complex(prc) :: R(Nv,4,Ls)
      complex(prc) :: DR(Nv,4,Ls)
      logical :: DAG
      real(prc) :: mass
      complex(prc) :: TMP(Nv,4,Ls)

      call setCGRVs(Nv*4*Ls,R);

      call IMDomWall(R,TMP,u,.false.,baremass)
      call MDomWall(TMP,DR,u,.false.,baremass)
      print *,"error:",maxval(abs(DR-R))
      call IMDomWall(R,TMP,u,.true.,baremass)
      call MDomWall(TMP,DR,u,.true.,baremass)
      print *,"error:",maxval(abs(DR-R))

      call MdagMDomWall(R,DR,u,baremass)
      call IMDomWall(DR,TMP,u,.true.,baremass)
      call IMDomWall(TMP,DR,u,.false.,baremass)
      print *,"error:",maxval(abs(DR-R))

      call IMdagMDomWall(R,TMP,u,baremass)
      call MdagMDomWall(TMP,DR,u,baremass)
      print *,"error:",maxval(abs(DR-R))

      return
      end subroutine testMDomWall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine MdagMDomWall(R,DR,u,mass)
      implicit none
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      real(prc),intent(in) :: mass
      complex(prc) :: TMP(Nv,4,Ls)

      call MDomWall(R,TMP,u,.false.,mass)
      call MDomWall(TMP,DR,u,.true.,mass)

      return
      end subroutine MdagMDomWall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMdagMDomWall(R,DR,u,mass)
      implicit none
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      real(prc),intent(in) :: mass
      complex(prc) :: TMP(Nv,4,Ls)

      call IMDomWall(R,TMP,u,.true.,mass)
      call IMDomWall(TMP,DR,u,.false.,mass)

      return
      end subroutine IMdagMDomWall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW(R,DR,u,DAGGER,mass)
      implicit none
c     calculates DR = DDW*R where DDW is the domain wall formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAGGER
      real(prc),intent(in) :: mass

!      print *,"dwkernel:",dwkernel

      if (DWkernel.eq.1) then
        call DDW_Shamir(R,DR,u,DAGGER,mass)
      elseif (DWkernel.eq.2) then
        call DDW_Wilson(R,DR,u,DAGGER,mass)
      elseif (DWkernel.eq.3) then
        call DDW_OWilson(R,DR,u,DAGGER,mass)
      else
        print *,"DWkernel not set properly"
        stop
      endif

      return
      end subroutine DDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDDW(RR,DR,u,DAGGER,mass)
!     solve DDW.DR = RR
      use options
      implicit none
      complex(prc) RR(Nv,4,Ls),DR(Nv,4,Ls)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) TMP(Nv,4,Ls)
      procedure(),pointer :: Dptr=>NULL()

      if (DWkernel.eq.1) then
        Dptr => DDW_Shamir
      elseif (DWkernel.eq.2) then
        Dptr => DDW_Wilson
      elseif (DWkernel.eq.3) then
        Dptr => DDW_OWilson
      else
        print *, "Domain Wall kernel not set properly"
        stop
      endif

      if (.not. DAGGER) then
        call Dptr(RR,TMP,u,.true.,mass)
        call IMdagM_DWkernel(TMP,DR,u,mass,Dptr)
        open(unit=105,file="fort.105",status='unknown')
        write(105,*) DR
        close(105)
      elseif (DAGGER) then
        call IMdagM_DWkernel(RR,TMP,u,mass,Dptr) 
        call Dptr(TMP,DR,u,.false.,mass)
      end if

      return
      end subroutine IDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDDWdagDDW(RR,DR,u,mass)
!     solve DDWdag.DDW.DR = RR
      use options
      implicit none
      complex(prc) RR(Nv,4,Ls),DR(Nv,4,Ls)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
!      complex(prc) TMP(Nv,4,Ls)
      procedure(),pointer :: Dptr=>NULL()

      if (DWkernel.eq.1) then
        Dptr => DDW_Shamir
      elseif (DWkernel.eq.2) then
        Dptr => DDW_Wilson
      elseif (DWkernel.eq.3) then
        Dptr => DDW_OWilson
      else
        print *, "Domain Wall kernel not set properly"
        stop
      endif

      call IMdagM_DWkernel(RR,DR,u,mass,Dptr)

      return
      end subroutine IDDWdagDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMdagM_DWkernel(RR,DR,u,mass,Dptr)
!     solve (MdagM).DR = RR for M=DDW
      use countmod
      implicit none
      integer,parameter :: kferm = 4*Nv*Ls
      integer,parameter :: niterc=Nv*Ls*100
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      real(prc) mass,shft
      procedure(),pointer :: Dptr

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call Dptr(RR,x1,u,.false.,mass)
      call Dptr(x1,x2,u,.true.,mass)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Dptr(p,x1,u,.false.,mass)
        call Dptr(x1,x2,u,.true.,mass)
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
      write(401,*) itercg,niterc,betan
      oc_idx=oc_idx+1
      outer_count=outer_count+itercg

      return
      end subroutine IMdagM_DWkernel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDDWdagDDWpC(RR,DR,u,mass,shft)
!     solve DDWdag.DDW.DR = RR
      use options
      implicit none
      complex(prc) RR(Nv,4,Ls),DR(Nv,4,Ls)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass,shft
      procedure(),pointer :: Dptr=>NULL()

      if (DWkernel.eq.1) then
        Dptr => DDW_Shamir
      elseif (DWkernel.eq.2) then
        Dptr => DDW_Wilson
      elseif (DWkernel.eq.3) then
        Dptr => DDW_OWilson
      else
        print *, "Domain Wall kernel not set properly"
        stop
      endif

      call IMdagMpC_DWkernel(RR,DR,u,mass,Dptr,shft)

      return
      end subroutine IDDWdagDDWpC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMdagMpC_DWkernel(RR,DR,u,mass,Dptr,shft)
!     solve (MdagM+C).DR = RR for M=DDW
      use countmod
      implicit none
      integer,parameter :: kferm = 4*Nv*Ls
      integer,parameter :: niterc=Nv*Ls
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      real(prc) mass,shft
      procedure(),pointer :: Dptr

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call Dptr(RR,x1,u,.false.,mass)
      call Dptr(x1,x2,u,.true.,mass)
      x2=x2+shft*RR
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Dptr(p,x1,u,.false.,mass)
        call Dptr(x1,x2,u,.true.,mass)
        x2=x2+shft*p
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
      write(401,*) itercg,niterc,betan
      oc_idx=oc_idx+1
      outer_count=outer_count+itercg

      return
      end subroutine IMdagMpC_DWkernel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_calcPhi(R,DR,u,DAGGER)
      use options
      implicit none
c     calculates DR = DDW*R where DDW is the domain wall formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      complex(prc) :: TMP1(Nv,4,Ls),TMP2(Nv,4,Ls)

      call DDW_Wilson(R,TMP1,u,.false.,baremass)
      DR(:,1:2)=TMP1(:,1:2,Ls)
      DR(:,3:4)=TMP1(:,3:4,Ls)

      end subroutine DDW_calcPhi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDDW_calcPhi(R,DR,u,DAGGER)
      use options
      implicit none
c     calculates DR = DDW*R where DDW is the domain wall formulation
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      complex(prc) :: TMP1(Nv,4,Ls),TMP2(Nv,4,Ls)

      call IDDW(R,TMP1,u,.false.,baremass)
      DR(:,1:2)=TMP1(:,1:2,Ls)
      DR(:,3:4)=TMP1(:,3:4,Ls)

      end subroutine IDDW_calcPhi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine PermM(R,DR,DAGGER,gi)
      use options
      implicit none
c     calculates DR(:,:,l) = Pminus.R(:,:,l)+Pplus.R(:,:,l+1)
c     if DAGGER=.true. then
c     calculates DR(:,:,l) = Pminus.R(:,:,l)+Pplus.R(:,:,l-1)
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      logical DAGGER
      integer gi
      complex(prc) :: TMP(Nv,4)
      integer l

      do l=1,Ls
        call Pminus(R(:,:,l),DR(:,:,l),gi)
      enddo
      do l=1,Ls
        call Pplus(R(:,:,l),TMP,gi)
        if (.not.DAGGER) then
          if (l.eq.1) then
            DR(:,:,Ls)=DR(:,:,Ls)+TMP
          else
            DR(:,:,l-1)=DR(:,:,l-1)+TMP
          endif
        else
          if (l.eq.Ls) then
            DR(:,:,1)=DR(:,:,1)+TMP
          else
            DR(:,:,l+1)=DR(:,:,l+1)+TMP
          endif
        endif
      enddo

      return
      end subroutine PermM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine KDDW(R,DR,u,DAGGER,mass)
      implicit none
c     calculates DR = Pdag.IDDW(1).DDW(m).P.R where P is the permutation matrix
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TMP(Nv,4,Ls)
      integer gi
      procedure(),pointer :: Dptr => NULL()
      integer MTMP

!      print *,"KDDW:",dwkernel

      if (DWkernel.eq.1) then
        Dptr=>DDW_Shamir
      elseif (DWkernel.eq.2) then
        Dptr=>DDW_Wilson
!        print *,"KDDW with Wilson"
      elseif (DWkernel.eq.3) then
        Dptr=>DDW_OWilson
!        print *,"KDDW with OWilson"
      else
        print *,"DW kernel not set properly"
        stop
      endif
!      print *,"MTYPE:",MTYPE

      MTMP=MTYPE
      gi=4
      if (.not.DAGGER) then
        open(unit=101,file="fort.101",status='unknown')
        write(101,*) R
        close(101)
        call PermM(R,TMP,.false.,gi)
        open(unit=102,file="fort.102",status='unknown')
        write(102,*) TMP
        close(102)
        call Dptr(TMP,DR,u,.false.,mass)
        open(unit=103,file="fort.103",status='unknown')
        write(103,*) DR
        close(103)
        MTYPE=1
        call IDDW(DR,TMP,u,.false.,one)
        open(unit=104,file="fort.104",status='unknown')
        write(104,*) TMP
        close(104)
        call PermM(TMP,DR,.true.,gi)
      elseif(DAGGER) then
        call PermM(R,TMP,.false.,gi)
        MTYPE=1
        call IDDW(TMP,DR,u,.true.,one)
        MTYPE=MTMP
        call Dptr(DR,TMP,u,.true.,mass)
        call PermM(TMP,DR,.true.,gi)
      endif
      MTYPE=MTMP

      return
      end subroutine KDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IKDDW(R,DR,u,DAGGER,mass)
      implicit none
c     calculates DR = Pdag.IDDW(m).DDW(1).P.R where P is the permutation matrix
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TMP(Nv,4,Ls)
      procedure(),pointer :: Dptr => NULL()
      integer gi
      integer MTMP

      if (DWkernel.eq.2) then
        Dptr=>DDW_Wilson
!        print *,"KDDW with Wilson"
      elseif (DWkernel.eq.3) then
        Dptr=>DDW_OWilson
!        print *,"KDDW with OWilson"
      else
        print *,"DW kernel not set properly"
        stop
      endif
      MTMP=MTYPE
      gi=4
      if (.not.DAGGER) then
        call PermM(R,TMP,DAGGER,gi)
        MTYPE=1
        call Dptr(TMP,DR,u,DAGGER,one)
        MTYPE=MTMP
        call IDDW(DR,TMP,u,DAGGER,mass)
        call PermM(TMP,DR,.not.DAGGER,gi)
      elseif (DAGGER) then
        call PermM(R,TMP,.not.DAGGER,gi)
        call IDDW(TMP,DR,u,DAGGER,mass)
        MTYPE=1
        call Dptr(DR,TMP,u,DAGGER,one)
        call PermM(TMP,DR,DAGGER,gi)
        MTYPE=MTMP
      endif

      return
      end subroutine IKDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine KDDW4(R,DR,u,DAGGER,mass)
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: R5(Nv,4,Ls)
      complex(prc) :: DR5(Nv,4,Ls)
      complex(prc) :: TMP(Nv,4,Ls)
      integer gi

      R5=czero
      R5(:,:,1)=R
      call KDDW(R5,DR5,u,DAGGER,mass)
      DR=DR5(:,:,1)
      
      return
      end subroutine KDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IKDDW4(R,DR,u,DAGGER,mass)
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: R5(Nv,4,Ls)
      complex(prc) :: DR5(Nv,4,Ls)
      complex(prc) :: TMP(Nv,4,Ls)
      integer gi

      R5=czero
      R5(:,:,1)=R
      call IKDDW(R5,DR5,u,DAGGER,mass)
      DR=DR5(:,:,1)
      
      return
      end subroutine IKDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DomainWallDerivs(dSdA,eta,nu,DAG,KTYPE,mass)
      implicit none
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls),intent(in) ::  eta,nu
      logical,intent(in) :: DAG
      integer,intent(in) :: KTYPE ! 1=Shamir,2=Wilson
      real(prc),intent(in) :: mass

      if (KTYPE.eq.1) then
        call ShamirDomainWallDerivs(dSdA,eta,nu,DAG)
      elseif (KTYPE.eq.2) then
        call WilsonDomainWallDerivs(dSdA,eta,nu,DAG,mass)
      else
        print *,"KTYPE not set properly"
        stop
      endif

      return
      end subroutine DomainWallDerivs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DomainWallDerivsComplex(dSdA,eta,nu,DAG,KTYPE,mass)
      implicit none
      complex(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls),intent(in) ::  eta,nu
      logical,intent(in) :: DAG
      integer,intent(in) :: KTYPE ! 1=Shamir,2=Wilson
      real(prc),intent(in) :: mass

      if (KTYPE.eq.1) then
        call ShamirDomainWallDerivsComplex(dSdA,eta,nu,DAG)
      elseif (KTYPE.eq.2) then
        call WilsonDomainWallDerivsComplex(dSdA,eta,nu,DAG,mass)
      else
        print *,"KTYPE not set properly"
        stop
      endif

      return
      end subroutine DomainWallDerivsComplex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module domainwallmod
