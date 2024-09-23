!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module condmod
      use pacc
      use arraysizes
      use numbers
      use options
      use IOmodule
      use gaugefield
      use ratfuncs
      implicit none
      type(sgnratfunc) :: SRF
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
!          call evalCondNoisy_DomWall_Wilson(u,pbp)
          call evalCondNoisy_OL_Wilson(u,pbp,SRF)
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
      subroutine evalCondNoisy_OL_Wilson(u,pbp,SRF)
      use rvmodule
      use ratfuncs
      use overlap
      implicit none
      complex(prc) u(Nv,3)
      complex(prc) pbp
      integer j
      type(sgnratfunc) :: SRF
      complex(prc),dimension(Nv,4) :: eta,IDR
      integer idx
      complex(prc) trcomp,denom

      pbp=czero
      do idx=1,4
        eta=czero
        call setCGRVs(Nv,eta(:,idx))
!        call IDOLop3(eta,IDR,u,.false.,baremass,SRF) ! separated multishift version
        call IDOLMW(eta,IDR,u,.false.,baremass,SRF)
!        IDR=conjg(eta)*IDR
        IDR=conjg(eta)*(IDR-eta) ! note this crashes for free field - check this
        if (MTYPE.eq.1) then
!          trcomp=sum(IDR(:,idx)-one)/(one-baremass)
          trcomp=sum(IDR(:,idx))/(one-baremass)
        elseif (MTYPE.eq.3) then ! assumes g3=diag(1 1 -1 -1)
          if (idx.eq.1 .or. (idx.eq.2)) then
            denom=cmplx(-baremass,one)
          else
            denom=cmplx(-baremass,-one)
          endif
!          trcomp=sum(IDR(:,idx)-one)/denom
          trcomp=sum(IDR(:,idx))/denom
        elseif (MTYPE.eq.4) then ! assumes g3=diag(1 1 -1 -1)
          if (idx.eq.1 .or. (idx.eq.2)) then
            denom=-cmplx(baremass,one)
          else
            denom=-cmplx(baremass,-one)
          endif
!          trcomp=sum(IDR(:,idx)-one)/denom
          trcomp=sum(IDR(:,idx))/denom
        else
          print *,"MASS OPTION NOT AVAILABLE IN evalCondNoisy_OL",MTYPE
          stop
        endif
!        print *,trcomp
        pbp=pbp+trcomp
      end do
      pbp=pbp/Nv
      print *,"pbp: ",pbp 

      return
      end subroutine evalCondNoisy_OL_Wilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module condmod
