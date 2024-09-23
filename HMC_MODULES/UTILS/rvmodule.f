      module rvmodule
      use arraysizes
      use numbers
      implicit none
      integer seedsize
      integer,allocatable,dimension(:) :: seed
      integer,allocatable,dimension(:) :: reseed

      contains
      subroutine initRVs(REPEATABLE,PROC_DISTINCT,aseed)
      implicit none
      logical REPEATABLE,PROC_DISTINCT
      integer aseed

!     requires fortran2018
c      call RANDOM_INIT(REPEATABLE,PROC_DISTINCT)
!     in meantime should use seed
c      call RANDOM_SEED()

      call RANDOM_SEED(size=seedsize)
      allocate(seed(seedsize))
      call RANDOM_SEED(get=seed)
      if (aseed.gt.0) then
        print *,"seed size:",seedsize
        print *,"initial seed:",seed
        print *,"seed(1):",seed(1)
        print *,"seed(2):",seed(2)
        allocate(reseed(seedsize))
        reseed=seed
        reseed(1)=aseed
        call RANDOM_SEED(put=reseed)
        call RANDOM_SEED(get=seed)
        print *,"new seed:",seed
      endif

      return
      end subroutine initRVs
c**********************************************************************
c calculate vector of gaussian random numbers with unit variance
c to refresh momenta
c   Numerical Recipes pp.203
c**********************************************************************
      subroutine gaussp(ps)
      implicit none
      real(prc) ps(Nv,2)
      integer il
      real(prc) theta
      do il=1,Nv
        ps(il,2)=sqrt(-2.0*log(urv()))
      end do
      do il=1,Nv
        theta=2*pi*urv()
        ps(il,1)=ps(il,2)*sin(theta)
        ps(il,2)=ps(il,2)*cos(theta)
      end do
      return
      end      
c**********************************************************************
c calculate vector of complex gaussian random numbers with unit variance
c to generate pseudofermion fields R (ie two independent fields with 0.5
c variance
c   Numerical Recipes pp.203
c**********************************************************************
      subroutine gauss0(ps)
      implicit none
      real(prc) ps(Nv,2)
      integer il
      real(prc) theta
      do  il=1,Nv
        ps(il,2)=sqrt(-log(urv()))
      end do
      do il=1,Nv
        theta=2*pi*urv()
        ps(il,1)=ps(il,2)*sin(theta)
        ps(il,2)=ps(il,2)*cos(theta)
      end do
      return
      end   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      real(prc) function urv()
      implicit none  

      call random_number(urv)
      
      return 
      end function urv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      real(prc) function grv()
      implicit none  

      real(prc) urv1,urv2,z0
      real(prc) mu,sigma

      call random_number(urv1)
      call random_number(urv2)
      
      grv=sqrt(-2.0*log(urv1))*cos(2.0*pi*urv2)

      return 
      end function grv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      real(prc) function hgrv()
!     gaussian rv with variance 0.5
      implicit none  

      real(prc) urv1,urv2,z0
      real(prc) mu,sigma

      call random_number(urv1)
      call random_number(urv2)
      
      hgrv=sqrt(-log(urv1))*cos(2.0*pi*urv2)

      return 
      end function hgrv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine setRVs(N,M)
      implicit none
      integer N
      complex(prc) :: M(N)
      integer j
      
      do j=1,N
        M(j)=cmplx(urv(),urv())
      end do

      return
      end subroutine setRVs   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine setGRVs(N,M)
      implicit none
      integer N
      real(prc) :: M(N)
      integer j
      
      do j=1,N
        M(j)=grv()
      end do

      return
      end subroutine setGRVs   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine setCGRVs(N,M)
      implicit none
      integer N
      complex(prc) :: M(N)
      integer j
      
      do j=1,N
        M(j)=cmplx(hgrv(),hgrv())
      end do

      return
      end subroutine setCGRVs   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine setRealRVs(N,M)
      implicit none
      integer N
      real(prc) :: M(N)
      integer j
      
      do j=1,N
        M(j)=urv()
      end do

      return
      end subroutine setRealRVs   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine testrvs()
      implicit none  

      integer j
      real(prc) x(1000),xbar,xvar,x2bar

      xbar=0
      x2bar=0
      xvar=0
      do j=1,1000
        x(j)=urv()
        xbar=xbar+x(j) 
        x2bar=x2bar+x(j)*x(j) 
      end do
      xbar=xbar/1000
      x2bar=x2bar/1000
      do j=1,1000
        xvar=xvar+(xbar-x(j))*(xbar-x(j))
      end do
      xvar=xvar/999
      print *,"urv mean:",xbar,0.5
      print *,"urv var:",xvar,1.0/12.0
      print *,"urv mean x^2:",x2bar,1.0/3.0

      return 
      end subroutine testrvs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine testgrvs()
      use statsmod
      implicit none  

      integer,parameter :: N=10000
      real(prc) r(N),rbar,rsd
      real(prc) h(N),hbar,hsd
      complex(prc) c(N),cbar,csd
      integer j

      do j=1,N
        r(j)=grv()
        h(j)=hgrv()
        c(j)=cmplx(hgrv(),hgrv())
      end do

      call calcVarReal(N,r,rbar,rsd)
      print *,"grv mean:",rbar,0.0
      print *,"grv sd:",rsd,1.0
      call calcVarReal(N,h,hbar,hsd)
      print *,"hgrv mean:",hbar,0.0
      print *,"hgrv sd:",hsd,sqrt(0.5)
      call calcVar(N,c,cbar,csd)
      print *,"cgrv mean:",cbar,0.0
      print *,"cgrv sd:",csd,1.0

      print *,sum(conjg(c)*c)/N

      return 
      end subroutine testgrvs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine make2piedgedist(Nedge,dist)
! make the quenched compact cosine cumulative distribution
      implicit none
      integer Nedge
      real(prc),dimension(Nedge) :: dist
      integer i
      real(prc) dx

      dx=2*pi/(Nedge-1)
      dist(1)=-pi
      do i=2,Nedge
        dist(i)=dist(i-1)+dx
      end do
      return
      end subroutine make2piedgedist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine makecumcosdist(beta,Nedge,cdist)
! make the quenched compact cosine cumulative distribution
      implicit none
      real(prc) beta
      integer Nedge
      real(prc),dimension(Nedge) :: cdist
      integer Ncentre,i
      real(prc) dtheta,ethetaL,ethetaR,ctheta
      real(prc),dimension(Nedge-1) :: dist

      Ncentre=Nedge-1
      dtheta=2*pi/Ncentre
      ethetaL=-pi
      do i=2,Nedge
        ethetaR=ethetaL+dtheta
        ctheta=(ethetaL+ethetaR)/two
        dist(i-1)=cospdf(beta,ctheta)
        ethetaL=ethetaL+dtheta
      end do
      dist=dist/sum(dist)
      cdist(1)=zero
      do i=2,Nedge
        cdist(i)=cdist(i-1)+dist(i-1)
      end do
      print *,"cdist:",cdist

      return
      end subroutine makecumcosdist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      real(prc) function crv(Nedge,cdist)
! make the cosine random random variable X with pdf
! 1/N exp(beta.cos(X))
      implicit none
      integer Nedge
      real(prc),dimension(Nedge) :: cdist
      integer idx
      real(prc) rv,mult,dtheta

      dtheta=2*pi/(Nedge-1)
      rv=urv()
      idx=distindex(rv,Nedge,cdist)
      print *,idx
      mult=(rv-cdist(idx))/(cdist(idx+1)-cdist(idx))
      crv = -pi+(idx-1)*dtheta+mult*dtheta
          
      end function crv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function cospdf(beta,x)
      implicit none
      real(prc) beta,x

      cospdf = exp(beta*(cos(x)-one))
      return
      end function cospdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer function distindex(rv,Nd,dist)
      implicit none
      real(prc) rv
      integer Nd
      real(prc),dimension(Nd) :: dist
      integer i

      do i=1,Nd-1
        if ((rv.ge.dist(i)).and.rv.le.dist(i+1)) then
          distindex=i
          return
        endif
      enddo
      print *,"index not found"
      stop
      return
      end function distindex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function gausspdf(var,x)
      implicit none
      real(prc) var,x

      gausspdf = one/sqrt(2*pi*var)*exp(-1/(2*var)*x*x)
      return
      end function gausspdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module rvmodule
