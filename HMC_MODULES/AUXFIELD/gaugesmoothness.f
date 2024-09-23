      module gaugesmoothness
      use pacc
      use gaugefield
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine extractGaugeLines(Ng)
      use IOmodule
      implicit none
      integer Ng
      character*1024 fname
      character*80 pname

      call readConvertedThetaFileName(100,pname)
      call readConvertedThirringGaugeField(pname,theta)
      call coef(u,theta)
      fname='line1.dat'
      call writeGaugeLines(fname,theta)
      call readConvertedThetaFileName(105,pname)
      call readConvertedThirringGaugeField(pname,theta)
      call coef(u,theta)
      fname='line2.dat'
      call writeGaugeLines(fname,theta)
      call readConvertedThetaFileName(110,pname)
      call readConvertedThirringGaugeField(pname,theta)
      call coef(u,theta)
      fname='line3.dat'
      call writeGaugeLines(fname,theta)
      call readConvertedThetaFileName(115,pname)
      call readConvertedThirringGaugeField(pname,theta)
      call coef(u,theta)
      fname='line4.dat'
      call writeGaugeLines(fname,theta)

      return
      end subroutine extractGaugeLines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testFields()
      use gaugefield
      use options
      use IOmodule
      implicit none
      integer,parameter :: Ncf=100 ! 20 for NC 100 for C
      character(len=80) fname 
      integer icf,j
      real(prc) thetaAv(3),totAv(3),totAvMag(3)
      integer istart,skip
      print *,"testFields"
      istart=3000 ! 500 for NC 3000 for C
      skip=50 

      totAv=0
      totAvMag=0
      do icf=1,Ncf
        call readThetaFileName(istart+skip*icf,fname)
        call readMPIConFile(fname,theta)
        call coef(u,theta)

        do j=1,3  
          thetaAv(j)=sum(theta(:,j))/Nv
        end do
        totAv=totAv+thetaAv
        totAvMag=totAvMag+abs(thetaAv)
      end do
      totAv=totAv/Ncf

      open(unit=11,file='auxData.dat',access='append',
     &           status='unknown',form='formatted')
      open(unit=12,file='linkData.dat',access='append',
     &           status='unknown',form='formatted')
      write(11,*) fname
      write(11,*) "Average:",totAv(1),totAvMag(1),
     &                       totAv(2),totAvMag(2),
     &                       totAv(3),totAvMag(3)
      close(11)
      close(12)

      return
      end subroutine testFields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module gaugesmoothness
