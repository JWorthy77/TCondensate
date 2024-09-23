      module overlap ! code used for main results in thesis
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      use WilsonDirac
      use WilsonExtraMod
      implicit none
      type ioptions
        type(sgnratfunc) :: SRF
        logical :: DAGGER
        real(prc) :: mass
        complex(prc) :: add
      end type
      logical,parameter :: VB_OL=.false.
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOLMW(R,DR,u,DAGGER,mass,SRF)
      use options
      use gammas
      implicit none
!     calculates DR = DOl*R  for Wilson kernel using multishift cg
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      complex(prc) :: TMP(Nv,4),VR(Nv,4),DL(Nv,4)

      if (VB_OL) then ; print *,"DOLMW" ; end if;

      if (MTYPE.eq.1) then ! DOL = (1+m) + (1-m).V

!        if (DAGGER .eq. .false.) then
          call VOLMWpf(R,DR,u,DAGGER,-MDW,SRF)
!        elseif (DAGGER .eq. .true.) then
!          call VOLMWpfDag(R,DR,u,DAGGER,-MDW,SRF)
!        endif
        
        DR=(one+mass)/two*R + (one-mass)/two*DR

      elseif (MTYPE.eq.3) then ! DOL = (1+i.m.g3) + (1-i.m.g3).V
!       this is not the form corresponding to the domain wall formulation MTYPE=3
        if (.not.DAGGER) then
          call VOLMWpf(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          DR=zi*mass*(R-VR)/two
          call mGmu(DR,4)
          DR=DL+DR
        elseif (DAGGER) then
          call VOLMWpf(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          TMP=R
          call mGmu(TMP,4)
          call VOLMWpf(TMP,VR,u,DAGGER,-MDW,SRF)
          DR=-zi*mass*(TMP-VR)/two
          DR=DL+DR
        endif

      elseif (MTYPE.eq.4) then 

        if (.not.DAGGER) then ! DOL=(1+im.g3) + V(1-im.g3)
          TMP=R
          call mGmu(TMP,4)
          DL=R+zi*mass*TMP
          VR=R-zi*mass*TMP
          call VOLMWpf(VR,DR,u,DAGGER,-MDW,SRF)
          DR=(DL+DR)/two
        elseif (DAGGER) then ! DOLdag=(1-im.g3) + (1+im.g3)Vdag
          TMP=R
          call mGmu(TMP,4)
          DL=R-zi*mass*TMP
          call VOLMWpf(R,VR,u,DAGGER,-MDW,SRF)
          TMP=VR
          call mGmu(TMP,4)
          DR=VR+zi*mass*TMP
          DR=(DL+DR)/two
        endif

      elseif (MTYPE.eq.5) then ! DOL = (1-i.m.g3) + (1+i.m.g3).V
!       this is not the form corresponding to the domain wall formulation MTYPE=3
        if (.not.DAGGER) then
          call VOLMWpf(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          DR=zi*mass*(-R+VR)/two
          call mGmu(DR,4)
          DR=DL+DR
        elseif (DAGGER) then
          call VOLMWpf(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          TMP=R
          call mGmu(TMP,4)
          call VOLMWpf(TMP,VR,u,DAGGER,-MDW,SRF)
          DR=-zi*mass*(-TMP+VR)/two
          DR=DL+DR
        endif

      endif

      return
      end subroutine DOLMW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLMWpf(RR,S,u,DAGGER,dwmass,SRF) 
      use options
!     approximate Wilson Voverlap using partial fraction rational functions with a multishift cg method
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      real(prc) mult
      if(VB_OL)then ; print *,"VOLMWpf" ; endif

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      if(VB_OL)then ; print *,"nn",nn,"nd",nd ; endif
      nf=nn-nd+1
      S=zero
      if (nf.ge.1) then ! it should only ever be 0 or 1
        S=front(1)*RR
      end if
      call MSCGW(RR,S,u,DAGGER,dwmass,SRF) ! S=S+...
      call DWilson(S,TMP1,u,DAGGER,dwmass)
      S=SRF%mult*TMP1

      end associate
      return
      end subroutine VOLMWpf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLMWpfDag(RR,S,u,DAGGER,dwmass,SRF) 
      use options
!     approximate Wilson Voverlap using partial fraction rational functions with a multishift cg method
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      real(prc) mult
      if(VB_OL)then ; print *,"VOLMWpf" ; endif

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      if(VB_OL)then ; print *,"nn",nn,"nd",nd ; endif
      nf=nn-nd+1
      call DWilson(RR,TMP1,u,DAGGER,dwmass)
      S=zero
      if (nf.ge.1) then ! it should only ever be 0 or 1
        S=front(1)*TMP1
      end if
      call MSCGW(TMP1,S,u,.false.,dwmass,SRF) ! S=S+...
      S=SRF%mult*S

      end associate
      return
      end subroutine VOLMWpfDag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine MSCGW(rhs,DR,u,DAGGER,dwmass,SRF)
      use options
      use countmod
!     a multishift cg method solves all (DdagD+sigma_i)Z_i=rhs
      implicit none
      complex(prc) rhs(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
c      procedure(),pointer :: Mptr
      complex(prc) Zm1ns(Nv,4)
      integer ns
      integer k,o,it,maxit,nsr
      complex(prc),dimension(:,:,:),allocatable :: Z,Zp1,P,Pp1
      real(prc),dimension(:),allocatable :: gamm1,gam,gamp1,sig,ocount
      complex(prc) R(Nv,4),Rp1(Nv,4),B(Nv,4)
      real(prc) alpham1,alpha,beta,betap1,dnm,conv,cns,drr,drp1rp1

      if(VB_OL)then ; print *,"MSCGW" ; endif
      maxit=10*Nv*4*10
      associate(sigma => SRF%pfrf%denom%zeros,
     &          mult => SRF%pfrf%pf)

      if(VB_OL)then ; print *,"sigma:",sigma ; endif
      if(VB_OL)then ; print *,"mult:",mult ; endif
      ns=size(sigma)
      nsr=ns
      allocate(Z(Nv,4,ns),P(Nv,4,ns))
      allocate(gamm1(ns),gam(ns),gamp1(ns),sig(ns),ocount(ns))
      Z=0
      sig=-sigma
      ocount=0
      do k=1,nsr
        P(:,:,k)=rhs
      end do
      R=rhs
      alpham1=one
      beta=zero
      gamm1=one
      gam=one
      gamp1=one
      dnm=dopr(rhs,rhs)
      if(VB_OL)then ; print *,"dnm:",dnm ; endif
      do it=1,maxit
        call DdagDpC(P(:,:,1),B,u,DAGGER,dwmass,czero)
        B=B+sig(1)*P(:,:,1)
        drr=dopr(R,R)
        alpha=drr/dopr(P(:,:,1),B)
        Rp1=R-alpha*B
        drp1rp1=dopr(Rp1,Rp1)
!        if(VB_OL)then ; print *,"drp1rp1:",drp1rp1 ; endif
        betap1=drp1rp1/drr
        Z(:,:,1)=Z(:,:,1)+alpha*P(:,:,1)
        P(:,:,1)=Rp1+betap1*P(:,:,1)
        Zm1ns=Z(:,:,nsr)
        do o=2,nsr
          gamp1(o)=gam(o)*gamm1(o)*alpham1/
     &      ( alpha*beta*(gamm1(o)-gam(o)) +
     &          gamm1(o)*alpham1*(one+alpha*(sig(o)-sig(1))) )
          Z(:,:,o)=Z(:,:,o)+alpha*gamp1(o)/gam(o)*P(:,:,o)
          P(:,:,o)=gamp1(o)*Rp1+betap1*(gamp1(o)/gam(o))**2*P(:,:,o)
        end do
        beta=betap1
        conv=sqrt(drp1rp1/dnm)
!        if(VB_OL)then ; print *,"beta:",beta,"conv:",conv ; endif
!        if(VB_OL) print *,maxval(abs(Z(:,:,nsr)-Zm1ns))
!        if (VB_OL) print *,maxval(abs(Z(:,:,nsr)))
        cns=maxval(abs(Z(:,:,nsr)-Zm1ns))/maxval(abs(Z(:,:,nsr)))
!        if (cns.lt.1d-30) then
!          if(VB_OL)then ; print *,"cns converged (1d30):",cns; endif
!          ocount(nsr)=it
!          nsr=nsr-1
!          if (nsr.eq.0) then
!            goto 501
!          endif
!        endif
        if (conv.lt.1d-12) then
          if(VB_OL)then ; print *,"conv converged (1d-12):",cns; endif
          if(VB_OL)then ; print *,"beta:",beta,"conv:",conv ; endif
          goto 501
        endif
        alpham1=alpha
        R=Rp1
        gamm1=gam
        gam=gamp1
      end do
      if(VB_OL)then ; print *,"not converged"; endif

501   continue
      do o=1,ns
        DR=DR+mult(o)*Z(:,:,o)
      end do

      ic_idx=ic_idx+1
      inner_count=inner_count+it

      deallocate(gamm1,gam,gamp1,sig,ocount)
      deallocate(Z,P)
      end associate
      return
      end subroutine MSCGW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function dopr(R1,R2) result(val)
      implicit none
      complex(prc) R1(Nv,4),R2(Nv,4)
      real(prc) val
      complex(prc) cval
      integer v,d
     
      cval=czero
      do d=1,4
        do v=1,Nv
          cval=cval+conjg(R1(v,d))*R2(v,d)
        end do
      end do
c      print *,"cval:",cval
      val=real(cval,prc)
      return 
      end function dopr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLMW(R,DR,u,DAGGER,mass,SRF)
      use options
      implicit none
!     solves DOL.DR = R 
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Mptr => NULL()
      type(ioptions)  iopts
      if(VB_OL)then ; print *,"IDOLMW" ; endif
      iopts%SRF=SRF
      iopts%mass=mass
      Mptr => DOLMW
      call IMnonsymOpts(R,DR,u,DAGGER,Mptr,iopts)

      return
      end subroutine IDOLMW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMnonsymOpts(RR,DR,u,DAGGER,Mptr,iopts) ! Mptr like SGNfactor
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      procedure(),pointer :: Mptr
      type(ioptions) iopts
      complex(prc) TMP(Nv,4)

      if (VB_OL) then ; print *,"IMnonsymOpts" ; end if ;
      if (.not. DAGGER) then
        call Mptr(RR,TMP,u,.not.DAGGER,iopts%mass,iopts%SRF)
        call IMdagMOpts(TMP,DR,u,DAGGER,Mptr,iopts)
      elseif (DAGGER) then
        call IMdagMOpts(RR,TMP,u,DAGGER,Mptr,iopts)
        call Mptr(TMP,DR,u,.not.DAGGER,iopts%mass,iopts%SRF)
      end if

      return
      end subroutine IMnonsymOpts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMdagMOpts(RR,DR,u,DAGGER,Mptr,iopts)
      use countmod
      implicit none
      integer,parameter :: kferm = 4*Nv
      integer,parameter :: niterc=Nv*10
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr
      type(ioptions) iopts

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

      if (VB_OL) then ; print *,"IMdagMOpts" ; end if ;
c     initialise
      DR=RR
!      call Mptr(RR,x1,u,DAGGER,iopts%mass,iopts%SRF)
!      call Mptr(x1,x2,u,.not.DAGGER,iopts%mass,iopts%SRF)
      call Mptr(RR,x1,u,.false.,iopts%mass,iopts%SRF)
      call Mptr(x1,x2,u,.true.,iopts%mass,iopts%SRF)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
!        call Mptr(p,x1,u,DAGGER,iopts%mass,iopts%SRF)
!        call Mptr(x1,x2,u,DAGGER,iopts%mass,iopts%SRF)
        call Mptr(p,x1,u,.false.,iopts%mass,iopts%SRF)
        call Mptr(x1,x2,u,.true.,iopts%mass,iopts%SRF)
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
      oc_idx=oc_idx+1
      outer_count=outer_count+i
!      print *,itercg,niterc,betan
      return
      end subroutine IMdagMOpts     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module overlap
