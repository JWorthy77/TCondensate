      module condnoisytoolsmod
      use pacc
      use arraysizes
      use numbers
      use rvmodule
      use options
      implicit none
      logical,parameter :: VB_CN=.false.
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalCondNoisy_OL(u,pbp,SRF)
      use rvmodule
      use ratfuncs
!      use overlapmoduledev
      use overlapmoduleprod
      implicit none
      complex(prc) u(Nv,3)
      complex(prc) pbp
      integer j
      type(sgnratfunc) :: SRF
      complex(prc),dimension(Nv,4) :: eta,IDR
      integer idx
      complex(prc) trcomp,denom

      if(VB_CN)then ; print *,"evalCondNoisy_OL" ; endif
      pbp=czero
      do idx=1,4
        if(VB_CN)then ; print *,"idx: ",idx ; endif
        print *,"idx: ",idx 
        eta=czero
        call setCGRVs(Nv,eta(:,idx))
        if(VB_CN)then ; print *,"eta: ",eta ; endif
!        call IDOLop(eta,IDR,u,.false.,baremass,SRF)
!        call IDOLop2(eta,IDR,u,.false.,baremass,SRF) ! multishift version
        if(VB_CN)then ; print *,"invert" ; endif
        call IDOLop3(eta,IDR,u,.false.,baremass,SRF) ! separated multishift version
        if(VB_CN)then ; print *,"inverted" ; endif
!        IDR=conjg(eta)*IDR
        IDR=conjg(eta)*(IDR-eta) ! note this crashes for free field
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
      if(VB_CN)then ; print *,"pbp: ",pbp ; endif
      open(unit=10,file="condOL.dat",status='unknown',access='append')
      print *,"pbp:",pbp
      write(10,*) "pbp:",pbp
      close(10)

      return
      end subroutine evalCondNoisy_OL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalCondNoisyB_OL(u,pbp,SRF)
      use rvmodule
      use ratfuncs
!      use overlapmoduledev
      use overlapmoduleprod
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
        call IDOLop3(eta,IDR,u,.false.,baremass,SRF) ! separated multishift version
        IDR=conjg(eta)*IDR
c        IDR=conjg(eta)*(IDR-eta) ! note this crashes for free field
        if (MTYPE.eq.1) then
          trcomp=sum(IDR(:,idx)-one)/(one-baremass)
c          trcomp=sum(IDR(:,idx))/(one-baremass)
        elseif (MTYPE.eq.3) then ! assumes g3=diag(1 1 -1 -1)
          if (idx.eq.1 .or. (idx.eq.2)) then
            denom=cmplx(-baremass,one)
          else
            denom=cmplx(-baremass,-one)
          endif
          trcomp=sum(IDR(:,idx)-one)/denom
c          trcomp=sum(IDR(:,idx))/denom
        elseif (MTYPE.eq.4) then ! assumes g3=diag(1 1 -1 -1)
          if (idx.eq.1 .or. (idx.eq.2)) then
            denom=-cmplx(baremass,one)
          else
            denom=-cmplx(baremass,-one)
          endif
          trcomp=sum(IDR(:,idx)-one)/denom
c          trcomp=sum(IDR(:,idx))/denom
        else
          print *,"MASS OPTION NOT AVAILABLE IN evalCondNoisy_OL",MTYPE
          stop
        endif
!        print *,trcomp
        pbp=pbp+trcomp
      end do
      pbp=pbp/Nv
!      print *,pbp

      return
      end subroutine evalCondNoisyB_OL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalSubCondNoisy_OL(u,pbp,SRF,idx)
      use rvmodule
      use ratfuncs
!      use overlapmoduledev
      use overlapmoduleprod
      implicit none
      complex(prc) u(Nv,3)
      complex(prc) pbp
      integer j
      type(sgnratfunc) :: SRF
      integer idx
      complex(prc),dimension(Nv,4) :: eta,IDR
      complex(prc) trcomp,denom

      pbp=czero
        eta=czero
        call setCGRVs(Nv,eta(:,idx))
c        call IDOLop(eta,IDR,u,.false.,baremass,SRF)
!        call IDOLop2(eta,IDR,u,.false.,baremass,SRF) ! multishift version
        call IDOLop3(eta,IDR,u,.false.,baremass,SRF) ! multishift version
c        IDR=conjg(eta)*IDR
        IDR=conjg(eta)*(IDR-eta)
        if (MTYPE.eq.1) then
c          trcomp=sum(IDR(:,idx)-one)/(one-baremass)
          trcomp=sum(IDR(:,idx))/(one-baremass)
        elseif (MTYPE.eq.3) then ! assumes g3=diag(1 1 -1 -1)
          if (idx.eq.1 .or. (idx.eq.2)) then
            denom=cmplx(-baremass,one)
          else
            denom=cmplx(-baremass,-one)
          endif
c          trcomp=sum(IDR(:,idx)-one)/denom
          trcomp=sum(IDR(:,idx))/denom
        elseif (MTYPE.eq.4) then ! assumes g3=diag(1 1 -1 -1)
          if (idx.eq.1 .or. (idx.eq.2)) then
            denom=-cmplx(baremass,one)
          else
            denom=-cmplx(baremass,-one)
          endif
c          trcomp=sum(IDR(:,idx)-one)/denom
          trcomp=sum(IDR(:,idx))/denom
        else
          print *,"MASS OPTION NOT AVAILABLE IN evalCondNoisy_OL"
          stop
        endif
!        print *,trcomp
        pbp=pbp+trcomp
      pbp=pbp/Nv
      print *,pbp

      return
      end subroutine evalSubCondNoisy_OL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalCondNoisy_KDDW4(u,pbp)
      use domainwallmod
      implicit none
      complex(prc) u(Nv,3)
      complex(prc) pbp
      integer j
      complex(prc),dimension(Nv,4) :: eta,IDR
      integer idx
      complex(prc) trcomp,denom

      if(VB_CN)then ; print *,"evalCondNoisy_KDDW4" ; endif
      pbp=czero
      do idx=1,4
        eta=czero
        call setCGRVs(Nv,eta(:,idx))
        call IKDDW4(eta,IDR,u,.false.,baremass)
!        IDR=conjg(eta)*IDR
        IDR=conjg(eta)*(IDR-eta) ! note this crashes for free field, um maybe?
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
          print *,"MASS OPTION NOT AVAILABLE IN evalCondNoisy_OL"
          stop
        endif
        print *,trcomp
        pbp=pbp+trcomp
      end do
      pbp=pbp/Nv
!      open(unit=10,file="condInstance.dat",status='unknown',
!     &                                             access='append')
!      print *,"pbp:",pbp
!      write(10,*) "pbp:",pbp
!      close(10)

      return
      end subroutine evalCondNoisy_KDDW4
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        call DWilson(TMP,DR,u,.false.,-MDW,-cone)
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
        call DWilson(TMP,DR,u,.false.,-MDW,-cone)
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
      end module condnoisytoolsmod
