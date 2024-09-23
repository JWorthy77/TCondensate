!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module dwcoeffs
      use pacc
      use options
      implicit none
      real(zprc) omega(Ls)
      end module dwcoeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module WilsonDomWall
      use arraysizes
      use options
      use WilsonDirac
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_Wilson(R,DR,u,DAGGER,mass)
      use gammas
      use ShamirDomWall
      implicit none
c     calculates DR = DDW*R where DDW is the domain wall formulation
c     with Wilson kernel
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc),intent(in) :: mass
      complex(prc) :: TR(Nv,4),TR2(Nv,4)
      integer s,gi
      complex(prc) zkappa

c     diagonal blocks
      do s=1,Ls
        call DWilson(R(:,:,s),DR(:,:,s),u,DAGGER,-MDW)
        DR(:,:,s)=DR(:,:,s)+R(:,:,s)
      end do

      gi=4 ! use gamma3
      if (.not. DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call Pminus(R(:,:,s+1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
          DR(:,:,s)=DR(:,:,s)+TR2
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call Pplus(R(:,:,s),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
          DR(:,:,s+1)=DR(:,:,s+1)+TR2
        end do
      elseif (DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s+1),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,s+1)
          call Pplus(TR,TR2,gi)
          DR(:,:,s)=DR(:,:,s)+TR2
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,s)
          call Pminus(TR,TR2,gi)
          DR(:,:,s+1)=DR(:,:,s+1)+TR2
        end do
      end if

c     mass terms
      if (MTYPE.eq.1) then

        if (.not.DAGGER) then
          call Pplus(R(:,:,Ls),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
        elseif (DAGGER) then
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,Ls)
          call Pminus(TR,TR2,gi)
        endif
        DR(:,:,1)=DR(:,:,1)-mass*TR2

        if (.not.DAGGER) then
          call Pminus(R(:,:,1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
        elseif (DAGGER) then
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,1)
          call Pplus(TR,TR2,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)-mass*TR2

      elseif ((MTYPE.eq.2).or.(MTYPE.eq.3)) then

        if (.not.DAGGER) then
          zkappa=cmplx(0,-mass)
          call Pminus(R(:,:,1),TR,gi)
          call mGmu(TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
        elseif (DAGGER) then
          zkappa=cmplx(0,mass)
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,1)
          call mGmu(TR,gi)
          call Pplus(TR,TR2,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)+zkappa*TR2

        if (.not.DAGGER) then
          zkappa=cmplx(0,-mass)
          call Pplus(R(:,:,Ls),TR,gi)
          call mGmu(TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR2=TR2-TR
        elseif(DAGGER) then
          zkappa=cmplx(0,mass)
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW)
          TR=TR-R(:,:,Ls)
          call mGmu(TR,gi)
          call Pminus(TR,TR2,gi)
        endif
        DR(:,:,1)=DR(:,:,1)+zkappa*TR2

      endif

      return
      end subroutine DDW_Wilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDW_OWilson(R,DR,u,DAGGER,mass)
      use dwcoeffs
      implicit none
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) :: TR(Nv,4),TR2(Nv,4),TR3(Nv,4)
      integer s,gi
      complex(prc) zkappa

c     diagonal blocks
      do s=1,Ls
        call DWilson(R(:,:,s),TR,u,DAGGER,-MDW)
        DR(:,:,s)=omega(s)*TR+R(:,:,s)
      end do

      gi=4
      if (.not. DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call Pminus(R(:,:,s+1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          DR(:,:,s)=DR(:,:,s)+omega(s)*TR2-TR
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call Pplus(R(:,:,s),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          DR(:,:,s+1)=DR(:,:,s+1)+omega(s+1)*TR2-TR
        end do
      elseif (DAGGER) then
c     upper diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s+1),TR,u,DAGGER,-MDW)
          TR2=omega(s+1)*TR-R(:,:,s+1)
          call Pplus(TR2,TR3,gi)
          DR(:,:,s)=DR(:,:,s)+TR3
        end do
c     lower diagonal blocks
        do s=1,Ls-1
          call DWilson(R(:,:,s),TR,u,DAGGER,-MDW)
          TR2=omega(s)*TR-R(:,:,s)
          call Pminus(TR2,TR3,gi)
          DR(:,:,s+1)=DR(:,:,s+1)+TR3
        end do
      end if

c     mass terms
      if (MTYPE.eq.1) then

        if (.not.DAGGER) then
          call Pminus(R(:,:,1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR3=omega(Ls)*TR2-TR
        elseif (DAGGER) then
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW)
          TR2=omega(1)*TR-R(:,:,1)
          call Pplus(TR2,TR3,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)-mass*TR3

        if (.not.DAGGER) then
          call Pplus(R(:,:,Ls),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR3=omega(1)*TR2-TR
        elseif (DAGGER) then
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW)
          TR2=omega(Ls)*TR-R(:,:,Ls)
          call Pminus(TR2,TR3,gi)
        endif
        DR(:,:,1)=DR(:,:,1)-mass*TR3

      elseif (MTYPE.eq.3) then

        if (.not.DAGGER) then
          zkappa=cmplx(0,-mass)
          call Pminus(R(:,:,1),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR3=omega(Ls)*TR2-TR
        elseif (DAGGER) then
          zkappa=cmplx(0,-mass)
          call DWilson(R(:,:,1),TR,u,DAGGER,-MDW)
          TR2=omega(1)*TR-R(:,:,1)
          call Pplus(TR2,TR3,gi)
        endif
        DR(:,:,Ls)=DR(:,:,Ls)+zkappa*TR3

        if (.not.DAGGER) then
          zkappa=cmplx(0,mass)
          call Pplus(R(:,:,Ls),TR,gi)
          call DWilson(TR,TR2,u,DAGGER,-MDW)
          TR3=omega(1)*TR2-TR
        elseif(DAGGER) then
          zkappa=cmplx(0,mass)
          call DWilson(R(:,:,Ls),TR,u,DAGGER,-MDW)
          TR2=omega(Ls)*TR-R(:,:,Ls)
          call Pminus(TR2,TR3,gi)
        endif
        DR(:,:,1)=DR(:,:,1)+zkappa*TR3

      endif

      return
      end subroutine DDW_OWilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Pplus(R,DR,gi)
      use gammas
      implicit none
c     calculates DR = (1+gamma)/2 R (gi should be 4 typically)
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      integer gi
      integer v,d,di

      do v=1,Nv
        do d=1,4
          di=gamin(gi,d)
          DR(v,d)=(R(v,d)+gamval(gi,d)*R(v,di))/two
        enddo
      enddo
      
      return
      end subroutine Pplus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Pminus(R,DR,gi)
      use gammas
      implicit none
c     calculates DR = (1-gamma)/2 R (gi should be 4 typically)
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      integer gi
      integer v,d,di

      do v=1,Nv
        do d=1,4
          di=gamin(gi,d)
          DR(v,d)=(R(v,d)-gamval(gi,d)*R(v,di))/two
        enddo
      enddo
      
      return
      end subroutine Pminus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WilsonDomainWallDerivs(dSdA,eta,nu,DAG,mass)
      use numbers
      use gammas
      use WilsonDirac
      use ShamirDomWall
      implicit none
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls),intent(in) ::  eta,nu
      logical,intent(in) :: DAG
      real(prc),intent(in) :: mass
      complex(prc) :: dSdAC(Nv,3)

      call WilsonDomainWallDerivsComplex(dSdAC,eta,nu,DAG,mass)
      dSdA=dSdAC

      return
      end subroutine WilsonDomainWallDerivs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WilsonDomainWallDerivsComplex(dSdA,eta,nu,DAG,mass)
      use numbers
      use gammas
      use WilsonDirac
      use dwcoeffs
      implicit none
      complex(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls),intent(in) ::  eta,nu
      logical,intent(in) :: DAG
      real(prc),intent(in) :: mass
      complex(prc) :: tmp(Nv,3)
      complex(prc),dimension(Nv,4) ::  eta_l,nu_l,tmp_l
      integer l
      integer pm1

      dSdA=0
      do l=1,Ls  ! diagonal terms
        eta_l=eta(:,:,l)
        nu_l=nu(:,:,l)
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA+omega(l)*tmp
      end do

      do l=1,Ls-1  ! upper diagonal
        if (.not.DAG) then
          eta_l=omega(l)*eta(:,:,l)
          call Pminus(nu(:,:,l+1),nu_l,4)
        else
          call Pplus(eta(:,:,l),eta_l,4)
          nu_l=omega(l+1)*nu(:,:,l+1)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA+tmp
      end do

      do l=2,Ls  ! lower diagonal
        if (.not.DAG) then
          eta_l=omega(l)*eta(:,:,l)
          call Pplus(nu(:,:,l-1),nu_l,4)
        else
          call Pminus(eta(:,:,l),eta_l,4)
          nu_l=omega(l-1)*nu(:,:,l-1)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA+tmp
      end do

!     mass components
      if (.not.DAG) then

        eta_l=omega(1)*eta(:,:,1)
!        eta_l=eta(:,:,1) ! DB
        if (MTYPE.eq.1) then
          call Pplus(nu(:,:,Ls),nu_l,4)
        else if (MTYPE.eq.3) then
          nu_l=-nu(:,:,Ls)
          tmp_l=zi*nu_l
          call mGmu(tmp_l,4)
          call Pplus(tmp_l,nu_l,4)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA-mass*tmp

        eta_l=omega(Ls)*eta(:,:,Ls)
!        eta_l=eta(:,:,Ls) ! DB
        if (MTYPE.eq.1) then
          call Pminus(nu(:,:,1),nu_l,4)
        else if (MTYPE.eq.3) then
          nu_l=-nu(:,:,1)
          tmp_l=zi*nu_l
          call mGmu(tmp_l,4)
          call Pminus(tmp_l,nu_l,4)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA-mass*tmp

      else if(DAG) then 

        nu_l=omega(Ls)*nu(:,:,Ls)
!        nu_l=nu(:,:,Ls) ! DB
        if (MTYPE.eq.1) then
          call Pminus(eta(:,:,1),eta_l,4)
        else if (MTYPE.eq.3) then
          eta_l=-eta(:,:,1)
          tmp_l=zi*eta_l
          call mGmu(tmp_l,4)
          call Pminus(tmp_l,eta_l,4)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA-mass*tmp

        nu_l=omega(1)*nu(:,:,1)
!        nu_l=nu(:,:,1) ! DB
        if (MTYPE.eq.1) then
          call Pplus(eta(:,:,Ls),eta_l,4)
        else if (MTYPE.eq.3) then
          eta_l=-eta(:,:,Ls)
          tmp_l=zi*eta_l
          call mGmu(tmp_l,4)
          call Pplus(tmp_l,eta_l,4)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA-mass*tmp
      endif

      return
      end subroutine WilsonDomainWallDerivsComplex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WilsonDomainWallDerivsComplexNZ(dSdA,eta,nu,DAG,mass)
      use numbers
      use gammas
      use WilsonDirac
      use dwcoeffs
      implicit none
      complex(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls),intent(in) ::  eta,nu
      logical,intent(in) :: DAG
      real(prc),intent(in) :: mass
      complex(prc) :: tmp(Nv,3)
      complex(prc),dimension(Nv,4) ::  eta_l,nu_l,tmp_l
      integer l
      integer pm1

      dSdA=0
      do l=1,Ls  ! diagonal terms
        eta_l=eta(:,:,l)
        nu_l=nu(:,:,l)
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA+omega(l)*tmp
      end do

      do l=1,Ls-1  ! upper diagonal
        if (.not.DAG) then
          eta_l=omega(l)*eta(:,:,l)
          call Pminus(nu(:,:,l+1),nu_l,4)
        else
          call Pplus(eta(:,:,l),eta_l,4)
          nu_l=omega(l+1)*nu(:,:,l+1)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA+tmp
      end do

      do l=2,Ls  ! lower diagonal
        if (.not.DAG) then
          eta_l=omega(l)*eta(:,:,l)
          call Pplus(nu(:,:,l-1),nu_l,4)
        else
          call Pminus(eta(:,:,l),eta_l,4)
          nu_l=omega(l-1)*nu(:,:,l-1)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA+tmp
      end do

!     mass components
      if (.not.DAG) then

        eta_l=omega(1)*eta(:,:,1)
        if (MTYPE.eq.1) then
          call Pplus(nu(:,:,Ls),nu_l,4)
        else if (MTYPE.eq.3) then
          nu_l=nu(:,:,Ls)
          tmp_l=zi*nu_l
          call mGmu(tmp_l,4)
          call Pplus(tmp_l,nu_l,4)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA-mass*tmp

        eta_l=omega(Ls)*eta(:,:,Ls)
        if (MTYPE.eq.1) then
          call Pminus(nu(:,:,1),nu_l,4)
        else if (MTYPE.eq.3) then
          nu_l=nu(:,:,1)
          tmp_l=zi*nu_l
          call mGmu(tmp_l,4)
          call Pminus(tmp_l,nu_l,4)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA-mass*tmp

      else if(DAG) then 

        nu_l=omega(Ls)*nu(:,:,Ls)
        if (MTYPE.eq.1) then
          call Pminus(eta(:,:,1),eta_l,4)
        else if (MTYPE.eq.3) then
          eta_l=eta(:,:,1)
          tmp_l=zi*eta_l
          call mGmu(tmp_l,4)
          call Pminus(tmp_l,eta_l,4)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA-mass*tmp

        nu_l=omega(1)*nu(:,:,1)
        if (MTYPE.eq.1) then
          call Pplus(eta(:,:,Ls),eta_l,4)
        else if (MTYPE.eq.3) then
          eta_l=eta(:,:,Ls)
          tmp_l=zi*eta_l
          call mGmu(tmp_l,4)
          call Pplus(tmp_l,eta_l,4)
        endif
        call WilsonDerivsComplex(tmp,eta_l,nu_l,DAG)
        dSdA=dSdA-mass*tmp
      endif

      return
      end subroutine WilsonDomainWallDerivsComplexNZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module WilsonDomWall
