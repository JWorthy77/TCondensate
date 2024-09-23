      module condensatetoolsmod
      use pacc
      use arraysizes
      use numbers
      use rvmodule
      use options
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalCondensateDW(u,mass,SWITCH,pbp,nnoise)
!     calculate condensate for Wilson Dirac operator
      use basicdiracopsmod
c      use rvmodule
      implicit none
      complex(prc) u(Nv,3)
      real(prc) mass
      integer SWITCH
      complex(prc) pbp
      integer nnoise
      complex(prc) pbpN,pbpP,pbpNP,pbpM,pbpNM
      complex(prc),dimension(Nv,4) :: delta,R,DR
      real(prc),dimension(Nv,2) :: eta
      integer j,id,n
      procedure(),pointer :: Mptr => NULL()
      procedure(),pointer :: IMptr => NULL()

      Mptr => DWilson
      IMptr => IDW
      pbp=cmplx(zero,zero,prc)
      if (SWITCH.eq.1) then ! point method
        print *,"Point Method"
        do id=1,4
          do j=1,Nv
c            print *,id,j,Nv
            delta=zero
            delta(j,id)=cone
            call IMptr(delta,DR,u,.false.,mass,czero)
            pbp=pbp+DR(j,id)
          end do
        end do
        pbp=pbp/Nv
      elseif (SWITCH.eq.2) then ! noisy estimator
        pbpNP=cmplx(zero,zero,prc)
        pbpNM=cmplx(zero,zero,prc)
        pbpN=cmplx(zero,zero,prc)
        do n=1,nnoise
c          print *,n,Nnoise
          do id=1,4
c            print *,id,4
            R=czero
            call gauss0(eta)
            R(:,id)=cmplx(eta(:,1),eta(:,2),prc)
            call IMptr(R,DR,u,.false.,mass,czero)
            pbpN=pbpN+dot_product(R(:,id),DR(:,id))/Nv
            pbpNP=pbpNP+dot_product(R(:,id),DR(:,id))/Nv
            pbpNM=pbpNM+dot_product(R(:,id),DR(:,id))/Nv
          end do
c          pbpN=pbp+pbpN/Nv
c          pbpP=pbpP+pbpNP/Nv
c          pbpM=pbpM+pbpNM/Nv
        end do
        pbp=pbpN/Nnoise
        pbpP=pbpP/Nnoise
        pbpM=pbpM/Nnoise
      endif
      print *,"condensate DW:",pbp
      return
      end subroutine evalCondensateDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalSingleCondensateOL(u,mass,SWITCH,pbp,SRF,IMptr,
     &              xskip,yskip,zskip,pbpVec,Npbp)
!     calculate overlap condensate with direct (partial fraction) 
!     or domain wall formulation with either Wilson or Shamir kernel
      use rvmodule
      use gammas
      use options
      use ratfuncs
      use overlapmoduledev
      use domainwallmod
      use statsmod
      use indices
      use IOmodule
      implicit none
      complex(prc) u(Nv,3)
      real(prc) mass
      integer SWITCH
      complex(prc) pbp
      type(sgnratfunc) :: SRF
      procedure(),pointer :: IMptr
      integer,intent(in) :: xskip,yskip,zskip
      complex(prc),dimension(4*Nv),intent(out) :: pbpVec
      integer,intent(out) :: Npbp
      complex(prc) pbpN,pbpP,pbpNP,pbpM,pbpNM
      complex(prc) pbpComp
      complex(prc),dimension(Nv,4) :: delta,R,DR,DG3R
      complex(prc),dimension(4) :: T4
      real(prc),dimension(Nv,2) :: eta
      integer j,idr,n,Np,ix,iy,it
      complex(prc) :: mean,stdev
      character*128 fname

      print *,"Eval single condensate"

      pbp=cmplx(zero,zero,prc)
      if (SWITCH.eq.1) then ! point method
      Np=0
      do ix=1,Ns,xskip
      do iy=1,Ns,yskip
      do it=1,Nt,zskip
        call ia(ix,iy,it,j)
        do idr=1,4
            Np=Np+1
            print *,Np," of ",4*Nv
            delta=zero
            delta(j,idr)=cone
            if(OLTYPE.eq.1)then
              call IMptr(delta,DR,u,.false.,mass,SRF)
            elseif(OLTYPE.eq.2)then
              call IMptr(delta,DR,u,.false.,mass)
            endif
            T4=DR(j,:)-one
            if (MTYPE.eq.1)then
              pbpComp=T4(idr)/(one-mass)
            elseif (MTYPE.eq.2) then
              pbpComp=T4(idr)/(-zi*gamval(4,idr)-mass)
            elseif (MTYPE.eq.3) then
              pbpComp=T4(idr)/(zi*gamval(4,idr)-mass)
            endif
            pbp=pbp+pbpComp
            pbpVec(Np)=pbpComp
            print *,idr,j,Nv,pbpComp
        end do
      end do
      end do
      end do
      Npbp=Np
      pbp=4*pbp/Np

      elseif (SWITCH.eq.2) then ! noisy estimator
        print *,"Noisy estimator not implemented for evalSingleCondensat
     &eOL"
      endif
      print *,pbp
      print *,"Eval condensate DONE"
      return
      end subroutine evalSingleCondensateOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalCondensateDDW2(u,pbp,unt)
!     calculate condensate for standard domain wall formulation
      use options
      use ratfuncs
      use gammas
      use domainwallmod
      implicit none
      complex(prc) u(Nv,3)
      complex(prc) pbp
      integer unt
      complex(prc) pbpcomp
      complex(prc),dimension(Nv,4,Ls) :: R,DR,delta
      integer id,j
      integer i1,i2
      integer skip,Nvals

      skip=Nv/20
      Nvals=0
      pbp=cmplx(zero,zero,prc)
      print *,"MTYPE:",MTYPE
      do j=1,Nv,skip
        Nvals=Nvals+1
        pbpcomp=czero
        do id=1,2
          print *,id,j,Nv
          delta=zero
          delta(j,id,1)=one
          DR=zero
          call IDDW(delta,DR,u,.false.,baremass)
          if (MTYPE.eq.3) then
            call mGmu5(DR,4)
          endif
          pbp=pbp+DR(j,id,Ls)
          pbpcomp=pbpcomp+DR(j,id,Ls)
          print *,DR(j,id,Ls)
          delta=zero
          delta(j,id+2,Ls)=one
          DR=zero
          call IDDW(delta,DR,u,.false.,baremass)
          if (MTYPE.eq.3) then
            call mGmu5(DR,4)
          endif
          pbp=pbp+DR(j,id+2,1)
          pbpcomp=pbpcomp+DR(j,id+2,1)
          print *,DR(j,id+2,1)
        end do
        if (MTYPE.eq.1) then
          print *,"pbpcomp",-pbpcomp
        elseif (MTYPE.eq.3) then
          print *,"pbpcomp",-zi*pbpcomp
        endif
        write(unt,*) pbpcomp
        flush(unt)
      end do
      if (MTYPE.eq.1) then
        pbp=-pbp/Nvals
      elseif (MTYPE.eq.3) then
        pbp=-zi*pbp/Nvals
      else
        print *,"evalCondensateDDW not implemented for MTYPE",MTYPE
      endif

      print *,"MTYPE:",MTYPE,"pbp:",pbp
      return
      end subroutine evalCondensateDDW2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module condensatetoolsmod
