!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module rhmc1domwallferm
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
      complex(prc) :: There(Ls,Nv,4)
      integer ii,id,il

      do ii=1,Nv
        do il=1,Ls
          do id=1,4
            transToHere(ii,id,il)=There(il,ii,id)
          end do
        end do
      end do

      return
      end function transToHere
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function transToThere(Here)
      implicit none
      complex(prc) :: Here(Nv,4,Ls)
      complex(prc) :: transToThere(Ls,Nv,4)
      integer ii,id,il

      do ii=1,Nv
        do il=1,Ls
          do id=1,4
            transToThere(il,ii,id)=Here(ii,id,il)
          end do
        end do
      end do

      return
      end function transToThere
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine initRHMC()
      use domainwallpowers
      use gaugefield                  
      use rvmodule
      implicit none 
      complex(prc) ut(Nv,3)
      real(prc) :: thetat(Nv,3)

      call setHalfPOWcoeffs()
      call setQuarterPOWcoeffs()
!      call testPowers()
!      call setGRVs(3*Nv,thetat)  ! randomise starting momentum
!      call coef(ut,thetat)
!      call testDWPowers(ut)

      return
      end subroutine initRHMC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setPseudoFermField(ps,ut) 
!     Pseudofermion fields: ps =  MdagM(1)^-1/4.MdagM(m)^1/4 * R, where R is complex gaussian
      use domainwallpowers       
      implicit none
      complex(prc),intent(out) :: ps(Nv,4,Ls)
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc) :: TMP(Nv,4,Ls)
      complex(prc) :: pstmp(Nv,4)
      integer baseMTYPE
      integer,parameter :: Npow=12
      real(prc) const,num(Npow),denom(Npow)
      integer l
     
      baseMTYPE=MTYPE

!      call setCGRVs(4*Nv*Ls,ps) ! randomise pseudo-fermion field
      call setCGRVs(4*Nv,pstmp) ! randomise pseudo-fermion field
!      ps=zero
!      ps(:,:,1)=pstmp ! DB
      do l=1,Ls
        ps(:,:,l)=pstmp
      end do
!      open(unit=19,file='fort.121',status='unknown')
!      read(19,*) TMP
!      close(19)
!      print *,"rv:",TMP
!      ps = transToHere(TMP)
!      print *,"ps:",transToThere(ps)
!      stop

      const=quarter_const
      num=quarter_num
      denom=quarter_denom
      call DDWdagDDWpow(ps,TMP,ut,baremass,Npow,const,num,denom)
!!      call DDWdagDDWpow(ps,TMP,ut,one,Npow,const,num,denom) ! debug

!      print *,"TMP:",transToThere(TMP)
      
c      print *,""
c      print *,""
c      print *,""
c      print *,"TMP:",transToThere(TMP)
c      print *,""
c      print *,""
c      print *,""

      const=mquarter_const
      num=mquarter_num
      denom=mquarter_denom
      MTYPE=1
      call DDWdagDDWpow(TMP,ps,ut,one,Npow,const,num,denom)
!!      call DDWdagDDWpow(TMP,ps,ut,one,Npow,const,num,denom) ! debug

!      print *,"ps:",transToThere(ps)
!      stop

      MTYPE=baseMTYPE
      return
      end subroutine setPseudoFermField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine march1DomWallFerm(dH,thetat) ! march dAdt=P (theta is A)
      use gaugefield                          !       dPdt=-dSdA
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
      call setGRVs(3*Nv,pp)  ! randomise starting momentum
!      open(unit=130,file='fort.130',status='unknown')
!      read(130,*) pp
!      close(13)

      call setPseudoFermField(ps,ut) ! randomise pseudo-fermion field

!      write(108,*) sum(abs(pp*pp))/(Nv*3),sum(abs(ps*ps))/(Nv*Ls*4)

      h0=ham1DomWallFerm(thetat,ut,pp,ps) ! initial hamiltonian energy
      if (VB_DWF) then ; print *,"h0:",h0 ; end if
      call force1DomWallFerm(thetat,ut,ps,F)
!      write(223,*) F
!      print *,"dSdp:",F
!      stop
      pp=pp-dt*half*F ! half time step before leap frog
      do ts=1,tsmax  ! time march
        if (VB_DWF) then; print *,ts-1,"sum(F):",sum(F),sum(F*F); endif
        thetat=thetat+dt*pp
!        if(VB_DWF) print *,"sum(thetat):",sum(thetat),sum(thetat*thetat)
        call coef(ut,thetat)
        call force1DomWallFerm(thetat,ut,ps,F)
!        if (VB_DWF) print *,"sum(F):",sum(F),sum(F*F)
        ytest=urv()
!        if(VB_DWF) then ; print *,"ytest:",ytest,"prob:",proby ; endif
        if (ytest.lt.proby) then
          pp=pp-dt*half*F
          goto 501
        else
          pp=pp-dt*F
        endif
      end do
!      pp=pp+half*dt*F ! correction to half step if tsmax reached

501   continue
      if (VB_DWF) then ; print *,ts,"sum(F):",sum(F),sum(F*F); endif
      h1=ham1DomWallFerm(thetat,ut,pp,ps) ! final hamiltonian energy
      if (VB_DWF) then ; print *,"h1:",h1 ; end if
      if (VB_DWF) then ; print *,"ts:",ts ; end if
!      stop
      dH=h0-h1
      if (VB_DWF) then ; print *,"dH:",dH ; end if
      return
      end subroutine march1DomWallFerm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine force1DomWallFerm(thetat,ut,ps,F)
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
      call fermion1force(ut,ps,dSdA)

      F=F+dSdA
      return
      end subroutine force1DomWallFerm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine forceThirring3(thetat,dSdA)
      implicit none
      real(prc) thetat(Nv,3),dSdA(Nv,3)

      dSdA=Nferms*gbeta*thetat

      return
      end subroutine forceThirring3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermion1force(ut,ps,dSdA) 
      use domainwallpowers
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,3) :: dSdA1,dSdA2
      complex(prc),dimension(Nv,3) :: dSdA1tot,dSdA2tot
      complex(prc),dimension(Nv,4,Ls) :: psihat,psibar
      complex(prc),dimension(Nv,4,Ls) :: lhs1,lhs2,rhs1,rhs2
      integer :: KTYPE
      integer baseMTYPE 
      integer,parameter :: Npow=12
      real(prc) const,num(Npow),denom(Npow)
      integer j,i,ind,ioffset

      dSdA=0
      KTYPE=kerneltype()
!      print *,"KTYPE:",KTYPE
      baseMTYPE=MTYPE

      MTYPE=1
      const=quarter_const ; num=quarter_num ; denom=quarter_denom
      call DDWdagDDWpow(ps,psihat,ut,one,Npow,const,num,denom)
!      call DDWdagDDWpow(ps,psihat,ut,one,Npow,const,num,denom) ! debug
!      print *,"psihat:",transToThere(psihat)
!      stop

      MTYPE=baseMTYPE
      const=mhalf_const ; num=mhalf_num ; denom=mhalf_denom
      call DDWdagDDWpow(psihat,psibar,ut,baremass,Npow,const,num,denom)
!      call DDWdagDDWpow(psihat,psibar,ut,baremass,Npow,const,num,denom) ! debug
!      print *,"psibar:",transToThere(psibar)
!      stop

      dSdA=czero
      dSdA1tot=zero
      dSdA2tot=zero
      do j=1,Npow
        MTYPE=1
        call IDDWdagDDWpC(ps,lhs1,ut,one,quarter_denom(j))
        call IDDWdagDDWpC(psibar,rhs1,ut,one,quarter_denom(j))
        call DDWdagDDWDerivs(dSdA1,lhs1,rhs1,ut,KTYPE,one)
!        call IDDWdagDDWpC(ps,lhs1,ut,one,quarter_denom(j)) ! debug
!        call IDDWdagDDWpC(psibar,rhs1,ut,one,quarter_denom(j)) ! debug
!        call DDWdagDDWDerivs(dSdA1,lhs1,rhs1,ut,KTYPE,one) ! debug

        MTYPE=baseMTYPE
        call IDDWdagDDWpC(psihat,rhs2,ut,baremass,mhalf_denom(j))
        call DDWdagDDWDerivs(dSdA2,rhs2,rhs2,ut,KTYPE,baremass)
!        call IDDWdagDDWpC(psihat,rhs2,ut,zero,mhalf_denom(j)) ! debug
!        call DDWdagDDWDerivs(dSdA2,rhs2,rhs2,ut,KTYPE,zero) ! debug

        dSdA=dSdA+2*quarter_num(j)*dSdA1+mhalf_num(j)*dSdA2
        dSdA1tot=dSdA1tot+quarter_num(j)*dSdA1
        dSdA2tot=dSdA2tot+mhalf_num(j)*dSdA2
      end do

      dSdA=-dSdA

!     anti-p.b.c. in timelike direction
      ioffset=(Nt-1)*Ns*Ns
      do i=1,Ns*Ns
        ind=ioffset+i
        dSdA(ind,3)=-dSdA(ind,3)
      enddo

      MTYPE=baseMTYPE
!      print *,"end of fforce"

!      print *,"dSdA2 ferm:",dSdA2tot
!      stop

!      print *,"dSdA ferm:",dSdA
!      stop

      return
      end subroutine fermion1force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham1DomWallFerm(thetat,ut,pp,ps)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat,pp
      complex(prc),dimension(Nv,3),intent(in) :: ut
      complex(prc),dimension(Nv,4,Ls),intent(in) :: ps
      real(prc) hg,hp,hf

      hp=0.5*sum(pp*pp)
      hg=hamThirring(thetat)
      hf=0
      if (.not.QUENCHED) then
        hf=ham1Ferm(ps,ut)
      endif
      ham1DomWallFerm=(hg+hp+hf)
      if (VB_DWF) print *,"hg:",hg
      if (VB_DWF) print *,"hp:",hp
      if (VB_DWF) print *,"hf:",hf
      if (VB_DWF) print *,"h:",ham1DomWallFerm
      write(101,'(4F12.4)') hg/Nv,hp/Nv,hf/Nv,sum(thetat*thetat)/Nv
      return
      end function ham1DomWallFerm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function hamThirring(thetat)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat

      hamThirring=Nferms*half*gbeta*sum(thetat*thetat)

      return
      end function hamThirring
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham1Ferm(ps,ut) ! psi^dag.(MdagM(1))^1/4.(MdagM(m))^-1/2.(MdagM(1))^1/4.psi
      use domainwallpowers
      implicit none
      complex(prc),intent(in) :: ps(Nv*4*Ls)
      complex(prc),intent(in) :: ut(Nv*3)
      complex(prc),dimension(Nv*4*Ls) :: lhs,rhs
      integer baseMTYPE 
      integer,parameter :: N=12
      real(prc) const,num(12),denom(12)

      if (VB_DWF) print *,'energy 1 DomWall Fermion'
      baseMTYPE=MTYPE

      MTYPE=1
      const=quarter_const ; num=quarter_num ; denom=quarter_denom
      call DDWdagDDWpow(ps,lhs,ut,one,N,const,num,denom)
!      call DDWdagDDWpow(ps,lhs,ut,one,N,const,num,denom) ! debug
!      print *,"lhs:",transToThere(lhs)
!      stop

      MTYPE=baseMTYPE
      const=mhalf_const ; num=mhalf_num ; denom=mhalf_denom
      call DDWdagDDWpow(lhs,rhs,ut,baremass,N,const,num,denom)
!      call DDWdagDDWpow(lhs,rhs,ut,baremass,N,const,num,denom) ! debug
!      print *,"rhs:",transToThere(rhs)
!      stop

      ham1Ferm=dot_product(lhs,rhs)

      return
      end function ham1Ferm
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
      end module rhmc1domwallferm
