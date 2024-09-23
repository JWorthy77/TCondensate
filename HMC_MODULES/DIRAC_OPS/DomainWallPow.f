!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module domainwallpowers
      use arraysizes
      use options
!      use WilsonDirac
!      use ShamirDomWall
!      use WilsonDomWall
      use domainwallmod
      implicit none

      real(prc) :: half_const
      real(prc) :: half_num(12)
      real(prc) :: half_denom(12)

      real(prc) :: mhalf_const
      real(prc) :: mhalf_num(12)
      real(prc) :: mhalf_denom(12)

      real(prc) :: quarter_const
      real(prc) :: quarter_num(12)
      real(prc) :: quarter_denom(12)

      real(prc) :: mquarter_const
      real(prc) :: mquarter_num(12)
      real(prc) :: mquarter_denom(12)

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDWdagDDWpow(R,DR,u,mass,N,const,num,denom)
      ! calculates DR = (DDWdagDDW)^pow.R
      implicit none
      complex(prc),intent(in) :: R(Nv,4,Ls)
      complex(prc),intent(out) :: DR(Nv,4,Ls)
      complex(prc),intent(in) :: u(Nv,3)
      real(prc),intent(in) :: mass
      integer N
      real(prc) const,num(N),denom(N)
      complex(prc) :: TMP(Nv,4,Ls)
      integer l

      DR=const*R
      do l=1,N
        call IDDWdagDDWpC(R,TMP,u,mass,denom(l))
        DR=DR+num(l)*TMP
      end do

      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DDWdagDDWDerivs(dSdA,lhs,rhs,ut,KTYPE,mass)
!     dMdagMdA=Mdag.dMdA+dMdagdA.M
      implicit none
      complex(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls),intent(in) :: lhs,rhs
      complex(prc),dimension(Nv,3),intent(in) :: ut
      integer :: KTYPE
      real(prc) :: mass
      complex(prc),dimension(Nv,4,Ls) :: lhs1,rhs1
      complex(prc),dimension(Nv,3) :: dSdA1,dSdA2
 
      call DDW(rhs,rhs1,ut,.false.,mass)
      call DomainWallDerivsComplex(dSdA1,lhs,rhs1,.true.,KTYPE,mass)
      call DDW(lhs,lhs1,ut,.false.,mass)
      call DomainWallDerivsComplex(dSdA2,lhs1,rhs,.false.,KTYPE,mass)

      dSdA=dSdA1+dSdA2

!      call DDW(lhs,lhs1,ut,.false.,mass)
!      print *,"sum(lhs1*rhs):",sum(conjg(lhs1)*rhs)
!      call DDW(rhs,rhs1,ut,.true.,mass)
!      print *,"sum(lhs*rhs1):",sum(conjg(lhs)*rhs1)
!      stop
      return
      end subroutine DDWdagDDWDerivs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setHalfPOWcoeffs()
      implicit none

      half_const = 3.4939402210660624e+01
      half_num(1) = -2.2595570037573332e-07 
      half_num(2) = -2.1018567078755487e-06 
      half_num(3) =-1.5112731454261532e-05
      half_num(4) =-1.0349689952587547e-04
      half_num(5) =-7.0043086782883921e-04
      half_num(6) =-4.7280433573365244e-03
      half_num(7) =-3.1962336099782780e-02
      half_num(8) =-2.1785728710770608e-01
      half_num(9) =-1.5306511467067245e+00
      half_num(10) =-1.2005797909776804e+01
      half_num(11) =-1.4305376570066056e+02
      half_num(12) =-1.3769554972593129e+04

      half_denom(1) = 4.6173174534437053e-05
      half_denom(2) = 2.6997165953746358e-04
      half_denom(3) = 1.0847426483492867e-03
      half_denom(4) = 3.9954125694289262e-03
      half_denom(5) = 1.4379118206005792e-02
      half_denom(6) = 5.1430954343411557e-02
      half_denom(7) = 1.8380164743730432e-01
      half_denom(8) = 6.5876654137294566e-01
      half_denom(9) = 2.3897526772383735e+00
      half_denom(10) = 9.0644346317094904e+00
      half_denom(11) = 4.0950110843585144e+01
      half_denom(12) = 4.7842097030875186e+02

      mhalf_const  = 2.8620981949567598e-02
      mhalf_num(1) = 4.2538727207099072e-03
      mhalf_num(2) = 6.0321777917329436e-03
      mhalf_num(3) = 1.0332242440400194e-02
      mhalf_num(4) = 1.8952011433362850e-02
      mhalf_num(5) = 3.5497194217913614e-02
      mhalf_num(6) = 6.6899784187265776e-02
      mhalf_num(7) = 1.2639133075134279e-01
      mhalf_num(8) = 2.3954426151265940e-01
      mhalf_num(9) = 4.5844694055900020e-01
      mhalf_num(10) = 9.0818544302681492e-01
      mhalf_num(11) = 2.0391621195489309e+00
      mhalf_num(12) = 7.4942613654523527e+00

      mhalf_denom(1) = 1.0451046902842116e-05
      mhalf_denom(2) = 1.2209979160003308e-04
      mhalf_denom(3) = 5.5160638287454170e-04
      mhalf_denom(4) = 2.0922667218341859e-03
      mhalf_denom(5) = 7.5899422420261714e-03
      mhalf_denom(6) = 2.7203238217467705e-02
      mhalf_denom(7) = 9.7217717692234806e-02
      mhalf_denom(8) = 3.4772646892294334e-01
      mhalf_denom(9) = 1.2514352180442438e+00
      mhalf_denom(10) = 4.6093882337979242e+00
      mhalf_denom(11) = 1.8520462512866683e+01
      mhalf_denom(12) = 1.0828798432022215e+02

      return
      end subroutine setHalfPOWcoeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setQuarterPOWcoeffs()
      implicit none

      quarter_const = 5.6498385400480871e+00
      quarter_num(1) = -1.5048280900635180e-06
      quarter_num(2) =-9.4742158074138799e-06
      quarter_num(3) = -4.8381491936974022e-05
      quarter_num(4) = -2.3914882558183337e-04
      quarter_num(5) = -1.1736612253436826e-03
      quarter_num(6) = -5.7518553830005889e-03
      quarter_num(7) = -2.8224764965678467e-02
      quarter_num(8) = -1.3937616807744102e-01
      quarter_num(9) = -7.0425716209846612e-01
      quarter_num(10) = -3.8663459742930737e+00
      quarter_num(11) = -2.8914380530992560e+01
      quarter_num(12) = -8.6183803495075233e+02

      quarter_denom(1) = 3.5185625599364660e-05
      quarter_denom(2) = 2.2542179339218258e-04
      quarter_denom(3) = 9.2254633714563801e-04
      quarter_denom(4) = 3.4113731261431065e-03
      quarter_denom(5) = 1.2279519127529441e-02
      quarter_denom(6) = 4.3882426976824672e-02
      quarter_denom(7) = 1.5661846452865200e-01
      quarter_denom(8) = 5.6025625272535784e-01
      quarter_denom(9) = 2.0245150863303447e+00
      quarter_denom(10) = 7.5937597706543736e+00
      quarter_denom(11) = 3.2877905381537445e+01
      quarter_denom(12) = 2.8696077476772416e+02


      mquarter_const = 1.7699620846005429e-01
      mquarter_num(1) = 1.9679229359758019e-04
      mquarter_num(2) = 5.0295961278630020e-04
      mquarter_num(3) = 1.2607080000709109e-03
      mquarter_num(4) = 3.2308495942848028e-03
      mquarter_num(5) = 8.3491564484584387e-03
      mquarter_num(6) = 2.1635747962193641e-02
      mquarter_num(7) = 5.6163502518072006e-02
      mquarter_num(8) = 1.4635499198455759e-01
      mquarter_num(9) = 3.8640016923667797e-01
      mquarter_num(10) = 1.0688835220783639e+00
      mquarter_num(11) = 3.5057341170825365e+00
      mquarter_num(12) = 2.2855137495681181e+01

      mquarter_denom(1) = 1.7423984180580676e-05
      mquarter_denom(2) = 1.5207781462890107e-04
      mquarter_denom(3) = 6.5843536680238422e-04
      mquarter_denom(4) = 2.4697272120915875e-03
      mquarter_denom(5) = 8.9244876352161696e-03
      mquarter_denom(6) = 3.1924715997233483e-02
      mquarter_denom(7) = 1.1394082653269423e-01
      mquarter_denom(8) = 4.0718206861948747e-01
      mquarter_denom(9) = 1.4656854630419729e+00
      mquarter_denom(10) = 5.4197819650664059e+00
      mquarter_denom(11) = 2.2180641564239263e+01
      mquarter_denom(12) = 1.4210348444366682e+02

      return
      end subroutine setQuarterPOWcoeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function evalRealParFrac(x,N,const,num,denom)
      implicit none
      real(prc) x
      integer N
      real(prc) const,num(N),denom(N)
      real(prc) val
      integer l

      val=const
      do l=1,N
        val=val+num(l)/(x+denom(l))
      end do
      evalRealParFrac=val
      return
      end function evalRealParFrac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testPowers()
      implicit none
      real(prc) x,val,pfval
      integer,parameter :: N=12
      integer nx
      real(prc) const,num(N),denom(N)
      real(prc) xL,xR,dx
      integer j

      xL=0.01 ; xR=10.0
      nx=20
      dx=(xR-xL)/(nx-1)

      const=half_const
      num=half_num
      denom=half_denom

      x=xL ; 
      do j=1,nx
        pfval=evalRealParFrac(x,N,const,num,denom)
        write(*,'(4f16.8)') x,x**0.5,pfval,abs(pfval-x**0.5)
        x=x+dx
      end do

      const=mhalf_const
      num=mhalf_num
      denom=mhalf_denom

      x=xL ; 
      do j=1,nx
        pfval=evalRealParFrac(x,N,const,num,denom)
        write(*,'(4f16.8)') x,x**(-0.5),pfval,abs(pfval-x**(-0.5))
        x=x+dx
      end do

      const=quarter_const
      num=quarter_num
      denom=quarter_denom

      x=xL ; 
      do j=1,nx
        pfval=evalRealParFrac(x,N,const,num,denom)
        write(*,'(4f16.8)') x,x**0.25,pfval,abs(pfval-x**0.25)
        x=x+dx
      end do

      const=mquarter_const
      num=mquarter_num
      denom=mquarter_denom

      x=xL ; 
      do j=1,nx
        pfval=evalRealParFrac(x,N,const,num,denom)
        write(*,'(4f16.8)') x,x**(-0.25),pfval,abs(pfval-x**(-0.25))
        x=x+dx
      end do

      return
      end subroutine testPowers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDWPowers(ut)
      use rvmodule
!      use gaugefield
      implicit none
      real(prc) x,val,pfval
      integer,parameter :: N=12
      integer nx
      real(prc) const,num(N),denom(N)
      complex(prc),dimension(Nv,4,Ls) :: R,DR,TMP
      complex(prc) ut(Nv,3)

      call setCGRVs(4*Nv*Ls,R) ! randomise pseudo-fermion field

      const=quarter_const
      num=quarter_num
      denom=quarter_denom
      call DDWdagDDWpow(R,TMP,ut,baremass,N,const,num,denom)
      const=mquarter_const
      num=mquarter_num
      denom=mquarter_denom
      call DDWdagDDWpow(TMP,DR,ut,baremass,N,const,num,denom)

      print *,"testDWPowers"
      print *,"err 1/4:",maxval(abs(DR-R))

      const=half_const
      num=half_num
      denom=half_denom
      call DDWdagDDWpow(R,TMP,ut,baremass,N,const,num,denom)
      const=mhalf_const
      num=mhalf_num
      denom=mhalf_denom
      call DDWdagDDWpow(TMP,DR,ut,baremass,N,const,num,denom)

      print *,"testDWPowers"
      print *,"err 1/2:",maxval(abs(DR-R))

      return
      end subroutine testDWPowers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module domainwallpowers
