      module statsmod
!     calcs mean,var,sd,etc...
      use arraysizes
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcVar(N,vec,mean,sd)
      implicit none
      integer N
      complex(prc) vec(N)
      complex(prc) mean,sd
      complex(prc) dv(N),var
      integer i

      mean=sum(vec)/N
      var=czero
      do i=1,N
        dv(i)=mean-vec(i)
        var=var+dv(i)*conjg(dv(i))
      end do
      var=var/(N-1)
      sd=var**half

      return
      end subroutine calcVar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcVarReal(N,vec,mean,sd)
      implicit none
      integer N
      real(prc) vec(N)
      real(prc) mean,sd
      real(prc) dv(N),var
      integer i

      mean=sum(vec)/N
      var=czero
      do i=1,N
        dv(i)=mean-vec(i)
        var=var+dv(i)*dv(i)
      end do
      var=var/(N-1)
      sd=var**half

      return
      end subroutine calcVarReal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module statsmod  
