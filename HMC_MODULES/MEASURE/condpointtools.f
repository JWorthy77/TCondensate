      module condpointtoolsmod
      use pacc
      use arraysizes
      use numbers
      use rvmodule
      use options
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalCondPoint_OL(u,pbp,j,SRF)
      use ratfuncs
!      use overlapmoduledev
      use overlapmoduleprod
      use statsmod
      implicit none
      complex(prc) u(Nv,3)
      complex(prc) pbp
      integer j
      type(sgnratfunc) :: SRF
      complex(prc),dimension(Nv,4) :: delta,IDR
      integer idx
      complex(prc) trcomp,denom

      pbp=czero
      do idx=1,4
        delta=czero
        delta(j,idx)=one
!        call IDOLop(delta,IDR,u,.false.,baremass,SRF)
        call IDOLop3(delta,IDR,u,.false.,baremass,SRF)
        if (MTYPE.eq.1) then
          trcomp=(IDR(j,idx)-one)/(one-baremass)
        elseif (MTYPE.eq.3) then ! assumes g3=diag(1 1 -1 -1)
          if (idx.eq.1 .or. (idx.eq.2)) then
            denom=cmplx(-baremass,one)
          else
            denom=cmplx(-baremass,-one)
          endif
          trcomp=(IDR(j,idx)-one)/denom
        elseif (MTYPE.eq.4) then ! assumes g3=diag(1 1 -1 -1)
          if (idx.eq.1 .or. (idx.eq.2)) then
            denom=-cmplx(baremass,one)
          else
            denom=-cmplx(baremass,-one)
          endif
          trcomp=(IDR(j,idx)-one)/denom
        else
          print *,"MASS OPTION NOT AVAILABLE IN evalCondPoint_OL"
          stop
        endif
!        print *,trcomp
        pbp=pbp+trcomp
      end do
      print *,"pbp(evalCondPoint_OL):",pbp

      return
      end subroutine evalCondPoint_OL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalIDDWcol(jV,jD,jL,u,mass,DAGGER,col)
!     calculate column of IDDW
      use options
      use ratfuncs
      use gammas
      use domainwallmod
      implicit none
      integer jV,jD,jL
      complex(prc) u(Nv,3)
      real(prc) mass
      logical DAGGER
      complex(prc),dimension(Nv,4,Ls) :: col,delta

      delta=czero
      delta(jV,jD,jL)=cone
      call IDDW(delta,col,u,DAGGER,mass)

      return
      end subroutine evalIDDWcol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalIDDWdiag(jV,jD,jL,iL,u,mass,DAGGER,trcomp)
!     calculate column of IDDW
      use options
      use ratfuncs
      use gammas
      use domainwallmod
      implicit none
      integer jV,jD,jL,iL
      complex(prc) u(Nv,3)
      real(prc) mass
      logical DAGGER
      complex(prc) trcomp
      complex(prc),dimension(Nv,4,Ls) :: col

      call evalIDDWcol(jV,jD,iL,u,mass,DAGGER,col)
      trcomp=col(jV,jD,jL)

      return
      end subroutine evalIDDWdiag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalCondPoint_DDW_Wilson(u,pbp,j)
!     calculate condensate for standard domain wall formulation
      use options
      use ratfuncs
      use gammas
      use domainwallmod
      implicit none
      ! pbp M1 S =  [ (Dw-1).P^+.IM_N1 + P^- IM_1N ] / V
      ! pbp M3 S = i[ P^+ IM_N1 - P^- IM_1N ] / V
      complex(prc) u(Nv,3)
      complex(prc) pbp
      integer j
      complex(prc) trcomp
      complex(prc) :: col(Nv,4,Ls),rhs(Nv,4),slc(Nv,4)
      integer jL,iL
      real(prc) mp
      logical,parameter :: DAGGER=.false.

      pbp=cmplx(zero,zero,prc)
      mp=one
      if (MTYPE.eq.3) then 
        mp=-one
      endif

      jL=Ls ; iL=1 ! M_N1
      call evalIDDWcol(j,1,iL,u,baremass,DAGGER,col)
!      print *,maxval(abs(col)),minval(abs(col))
!      stop
      call Pplus(col(:,:,jL),slc,4)
      call DWilson(slc,rhs,u,DAGGER,-MDW,-cone)
      trcomp=rhs(j,1)
      pbp=pbp+trcomp
!      print *,"trcomp:",trcomp
      call evalIDDWcol(j,2,iL,u,baremass,DAGGER,col)
      call Pplus(col(:,:,jL),slc,4)
      call DWilson(slc,rhs,u,DAGGER,-MDW,-cone)
      trcomp=rhs(j,2)
      pbp=pbp+trcomp
!      print *,"trcomp:",trcomp
      call evalIDDWcol(j,3,iL,u,baremass,DAGGER,col)
      call Pplus(col(:,:,jL),slc,4)
      call DWilson(slc,rhs,u,DAGGER,-MDW,-cone)
      trcomp=rhs(j,3)
      pbp=pbp+trcomp
!      print *,"trcomp:",trcomp
      call evalIDDWcol(j,4,iL,u,baremass,DAGGER,col)
      call Pplus(col(:,:,jL),slc,4)
      call DWilson(slc,rhs,u,DAGGER,-MDW,-cone)
      trcomp=rhs(j,4)
      pbp=pbp+trcomp
!      print *,"trcomp:",trcomp

      jL=1 ; iL=Ls 
      call evalIDDWcol(j,1,iL,u,baremass,DAGGER,col)
      call Pminus(col(:,:,jL),slc,4)
      call DWilson(slc,rhs,u,DAGGER,-MDW,-cone)
      trcomp=rhs(j,1)
      pbp=pbp+mp*trcomp
!      print *,"trcomp:",trcomp
      call evalIDDWcol(j,2,iL,u,baremass,DAGGER,col)
      call Pminus(col(:,:,jL),slc,4)
      call DWilson(slc,rhs,u,DAGGER,-MDW,-cone)
      trcomp=rhs(j,2)
      pbp=pbp+mp*trcomp
!      print *,"trcomp:",trcomp
      call evalIDDWcol(j,3,iL,u,baremass,DAGGER,col)
      call Pminus(col(:,:,jL),slc,4)
      call DWilson(slc,rhs,u,DAGGER,-MDW,-cone)
      trcomp=rhs(j,3)
      pbp=pbp+mp*trcomp
!      print *,"trcomp:",trcomp
      call evalIDDWcol(j,4,iL,u,baremass,DAGGER,col)
      call Pminus(col(:,:,jL),slc,4)
      call DWilson(slc,rhs,u,DAGGER,-MDW,-cone)
      trcomp=rhs(j,4)
      pbp=pbp+mp*trcomp
!      print *,"trcomp:",trcomp
       
      pbp=pbp
      if (MTYPE.eq.3) then
        pbp=zi*pbp
      endif

      if ((MTYPE.ne.1).and.(MTYPE.ne.3)) then
        print *,"MTYPE:",MTYPE,"not supported"
        stop
      endif

      print *,"MTYPE:",MTYPE,"pbp(evalCondPoint_DDW_Wilson):",pbp
      return
      end subroutine evalCondPoint_DDW_Wilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalCondPoint_DDW_Shamir(u,pbp,j)
!     calculate condensate for standard domain wall formulation
      use options
      use ratfuncs
      use gammas
      use domainwallmod
      implicit none
      ! pbp M1 S =  [ P^+ IM_N1 + P^- IM_1N ] / V
      ! pbp M3 S = i[ P^+ IM_N1 - P^- IM_1N ] / V
      complex(prc) u(Nv,3)
      complex(prc) pbp
      integer j
      complex(prc) trcomp
      integer jL,iL
      real(prc) mp

      pbp=cmplx(zero,zero,prc)
      mp=one
      if (MTYPE.eq.3) then 
        mp=-one
      endif

      jL=Ls ; iL=1 ! M_N1
      call evalIDDWdiag(j,1,jL,iL,u,baremass,.false.,trcomp)
      pbp=pbp+trcomp
      print *,"trcomp:",trcomp
      call evalIDDWdiag(j,2,jL,iL,u,baremass,.false.,trcomp)
      pbp=pbp+trcomp
      print *,"trcomp:",trcomp
      jL=1 ; iL=Ls 
      call evalIDDWdiag(j,3,jL,iL,u,baremass,.false.,trcomp)
      pbp=pbp+mp*trcomp
      print *,"trcomp:",trcomp
      call evalIDDWdiag(j,4,jL,iL,u,baremass,.false.,trcomp)
      pbp=pbp+mp*trcomp
      print *,"trcomp:",trcomp
       
      pbp=pbp
      if (MTYPE.eq.3) then
        pbp=zi*pbp
      endif

      if ((MTYPE.ne.1).and.(MTYPE.ne.3)) then
        print *,"MTYPE:",MTYPE,"not supported"
        stop
      endif

      print *,"MTYPE:",MTYPE,"pbp:",pbp
      return
      end subroutine evalCondPoint_DDW_Shamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module condpointtoolsmod
