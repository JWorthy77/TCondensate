      module utilsmod
      use pacc
      use arraysizes
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function mag(V)
      implicit none
      complex(prc) :: V(Nv*Ndc)
      real(prc) :: rV(Nv*Ndc),iV(Nv*Ndc)

      rV=real(V,prc)
      iV=dimag(V)
      mag=dot_product(rV,rV)+dot_product(iV,iV)
      mag=sqrt(mag)
      return
      end function mag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function mag2(V)
      implicit none
      complex(prc) :: V(Nv*2)
      real(prc) :: rV(Nv*2),iV(Nv*2)

      rV=real(V,prc)
      iV=dimag(V)
      mag2=dot_product(rV,rV)+dot_product(iV,iV)
      mag2=sqrt(mag2)
      return
      end function mag2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine normalise(V)
      implicit none
      complex(prc) :: V(Nv*Ndc)
      real(prc) :: nrm

      nrm=mag(V)
      V=V/nrm
      return
      end subroutine normalise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine normalise2(V)
      implicit none
      complex(prc) :: V(Nv*2)
      real(prc) :: nrm

      nrm=mag2(V)
      V=V/nrm
      return
      end subroutine normalise2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function magv(N,V)
      implicit none
      integer N
      complex(prc) :: V(N)
      real(prc) :: rV(N),iV(N)

      rV=real(V,prc)
      iV=dimag(V)
      magv=dot_product(rV,rV)+dot_product(iV,iV)
      magv=sqrt(magv)
      return
      end function magv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine normalisev(N,V)
      implicit none
      integer N
      complex(prc) :: V(N)
      real(prc) :: nrm

      nrm=magv(N,V)
      V=V/nrm
      return
      end subroutine normalisev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer function intbin(val,Nedges,edges)
      implicit none
      real(prc) val
      integer Nedges
      real(prc),dimension(Nedges) :: edges
      integer j

      do j=2,Nedges
        if (val .le. edges(j)) then
         intbin=j-1
         return
        endif
      enddo
      print *,"bin not found"
      stop
      return
      end function intbin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module utilsmod
