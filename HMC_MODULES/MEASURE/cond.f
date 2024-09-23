!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module condmod
      use pacc
      use arraysizes
      use numbers
      use options
      use IOmodule
      use gaugefield
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine measureCondensateInstance()
      implicit none
      complex(prc) pbp,pbptot
      integer n,Nnoise
   
      pbptot=0 ; Nnoise=10
      do n=1,Nnoise
        if (dwkernel.eq.1) then
          call evalCondNoisy_DomWall_Shamir(u,pbp)
        elseif ((dwkernel.eq.2).or.(dwkernel.eq.3))then
          call evalCondNoisy_DomWall_Wilson(u,pbp)
        endif
        pbptot=pbptot+pbp
      end do
      pbptot=pbptot/Nnoise
      write(200,*) real(pbptot),real(pbptot)/baremass

      return
      end subroutine measureCondensateInstance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalCondNoisy_DomWall_Shamir(u,pbp)
c      use rvmodule
      use domainwallmod
      implicit none
      complex(prc) u(Nv,3),pbp
      integer idx
      complex(prc),dimension(Nv,4,Ls) :: eta,IDR
      complex(prc),dimension(Nv) :: rvs,XDR
      complex(prc) trcomp,denom
      real(prc) m

      m=one
      if (MTYPE.eq.3) then
        m=-one
      endif
      pbp=czero
      do idx=1,2
        eta=czero
        call setCGRVs(Nv,rvs)
        eta(:,idx,1)=rvs
        call IDDW(eta,IDR,u,.false.,baremass)
        XDR=conjg(rvs)*IDR(:,idx,Ls)
        trcomp=sum(XDR)
!        print *,"trcomp:",trcomp
        pbp=pbp+trcomp
      end do        
      do idx=3,4
        eta=czero
        call setCGRVs(Nv,rvs)
        eta(:,idx,Ls)=rvs
        call IDDW(eta,IDR,u,.false.,baremass)
        XDR=conjg(rvs)*IDR(:,idx,1)
        trcomp=sum(XDR)
!        print *,"trcomp:",trcomp
        pbp=pbp+m*trcomp
      end do        
      pbp=pbp/Nv
      if (MTYPE.eq.3)then
        pbp=zi*pbp
      endif
      print *,pbp

      return
      end subroutine evalCondNoisy_DomWall_Shamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalCondNoisy_DomWall_Wilson(u,pbp)
c      use rvmodule
      use domainwallmod
      implicit none
      complex(prc) u(Nv,3),pbp
      integer idx
      complex(prc),dimension(Nv,4,Ls) :: eta,IDR
      complex(prc),dimension(Nv) :: rvs,XDR
      complex(prc),dimension(Nv,4) :: TMP,DR
      complex(prc) trcomp,denom
      real(prc) m

      m=one
      if (MTYPE.eq.3) then
        m=-one
      endif
      pbp=czero
!      do idx=1,2
      do idx=1,4
        eta=czero
        call setCGRVs(Nv,rvs)
        eta(:,idx,1)=rvs
        call IDDW(eta,IDR,u,.false.,baremass)
        TMP=czero
        TMP(:,1:2)=IDR(:,1:2,Ls)
!        call DWilson(TMP,DR,u,.false.,-MDW,-cone)
        call DWilson(TMP,DR,u,.false.,-MDW)
        DR=DR-TMP
        XDR=conjg(rvs)*DR(:,idx)
        trcomp=sum(XDR)
!        print *,"trcomp:",trcomp
        pbp=pbp+trcomp
      end do        
!      do idx=3,4
      do idx=1,4
        eta=czero
        call setCGRVs(Nv,rvs)
        eta(:,idx,Ls)=rvs
        call IDDW(eta,IDR,u,.false.,baremass)
        TMP=czero
        TMP(:,3:4)=IDR(:,3:4,1)
!        call DWilson(TMP,DR,u,.false.,-MDW,-cone)
        call DWilson(TMP,DR,u,.false.,-MDW)
        DR=DR-TMP
        XDR=conjg(rvs)*DR(:,idx)
        trcomp=sum(XDR)
!        print *,"trcomp:",trcomp
        pbp=pbp+m*trcomp
      end do        
      pbp=pbp/Nv
      if (MTYPE.eq.3)then
        pbp=zi*pbp
      endif
      print *,pbp

      return
      end subroutine evalCondNoisy_DomWall_Wilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module condmod
