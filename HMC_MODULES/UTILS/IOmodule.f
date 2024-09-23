      module IOmodule
!     we have 3 file types to read
      use arraysizes
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine readThetaFromFile(idx,theta)
!     read file location and type from 'confiledir.txt'
      implicit none
      integer idx
      real(prc) theta(Ns,Ns,Nt,3)
      character(len=80) fdir 
      integer ftype
      character(len=80) fname
      character(len=5) cnum
      integer i,nc

      call readThetaPathAndType(fdir,ftype)

      if (ftype.eq.1) then ! use standard (MPI) con file
        fname=trim(fdir)//"/con."
      elseif (ftype.eq.2) then ! use converted con(theta) file
        fname=trim(fdir)//"/theta"
      endif

      cnum=itoa(idx)
      nc=len(trim(cnum))
      do i=nc+1,5
        cnum="0"//cnum
      enddo
      fname=trim(fname)//cnum
 
      if (ftype.eq.2) then 
        fname=trim(fdir)//".dat"
      endif
      print *,"file name:",fname
      print *,"file type:",ftype

      if (ftype.eq.1) then 
        call readMPIConFile(fname,theta) 
      elseif (ftype.eq.2) then 
        call readConvertedThirringGaugeField(fname,theta) 
      endif

      return
      end subroutine readThetaFromFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function itoa(itgr) result(res)
        character*10 :: res
        integer,intent(in) :: itgr
        character(range(itgr)+2) :: tmp
        write(tmp,'(i0)') itgr
        print *,tmp
        res=trim(tmp)
      end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine readThetaPathAndType(fdir,ftype)
!     read the conidx.dat directory and return full path-name
      implicit none
      character(len=80) fdir 
      integer ftype
      logical FEXISTS
      integer i,nc

      open(unit=31,file='confiledir.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) fdir
        read(31,*) ftype
        close(31)
        print *,"file directory:",fdir
        print *,"file type:",ftype
      else
        print *,"couldn't open con dir file"
        stop
      endif

      return
      end subroutine readThetaPathAndType
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine readThetaFileName(idx,fname)
!     read the conidx.dat directory and return full path-name
      implicit none
      integer idx
      character(len=80) fname 
      character(len=80) fdir 
      character(len=5) cnum
      logical FEXISTS
      integer i,nc

      open(unit=31,file='confiledir.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) fdir
        close(31)
        print *,fdir
      else
        print *,"couldn't open con dir file"
        stop
      endif

      fname=trim(fdir)//"/con."

      cnum=itoa(idx)
      nc=len(trim(cnum))
      do i=nc+1,5
        cnum="0"//cnum
      enddo
      fname=trim(fname)//cnum
      print *,fname

      return
      end subroutine readThetaFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine readConvertedThetaFileName(idx,fname)
!     read the thetaidx.dat directory and return full path-name
      use arraysizes
      implicit none
      integer idx
      character(len=80) fname 
      character(len=80) fdir 
      character(len=5) cnum
      logical FEXISTS
      integer i,nc

      open(unit=31,file='confiledir.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) fdir
        close(31)
        print *,fdir
      else
        print *,"couldn't open con dir file"
        stop
      endif

      fname=trim(fdir)//"/theta"
      cnum=itoa(idx)
      fname=trim(fname)//trim(cnum)//".dat"

!      fname=trim(fdir)//"/theta."
!      cnum=itoa(idx)
!      nc=len(trim(cnum))
!      do i=nc+1,5
!        cnum="0"//cnum
!      enddo
!      fname=trim(fname)//cnum
      print *,fname

      return
      end subroutine readConvertedThetaFileName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine readThirringConFile(idx,theta)
c      use arraysizes
c      implicit none
c      integer idx
c      real(prc) theta(Ns,Ns,Nt,3)
c      character(len=80) fdir 
c      character(len=5) cnum
c      character(len=80) fname 
c      integer i,nc
c
c!      fname="/home/jude/2021/Sunbird/COMPACT/SIMS/12x12/Ls20/Bp25/mp01/g
c!     &3/con."
c!      fname="/home/jude/2021/Sunbird/SIMS/12x12/Ls20/Bp25/mp01/g3/con."
c!      fname="/home/jude/2021/Sunbird/SIMS/12x12/Ls20/Bp50/mp01/g3/con."
c      fname="/home/jude/2021/Sunbird/SIMS/12x12/Ls20/B2p0/mp00/g3/con."
c!      fname="/home/jude/2021/Sunbird/SIMS/12x12/Ls40/B1p0/mp01/g3/con."
c      cnum=itoa(idx)
cc      print *,cnum
c      nc=len(trim(cnum))
c      do i=nc+1,5
c        cnum="0"//cnum
cc        print *,cnum
c      enddo
c      fname=trim(fname)//cnum
c      print *,fname
c
cc      fname='con.00500'
c      call readMPIConFile(fname,theta)
c
c      return
c      end subroutine readThirringConFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine readConvertedThirringGaugeField(fname,theta) ! ftype=2
c      use arraysizes
      implicit none
c      integer idx
      character(len=80) fname 
      real(prc) theta(Ns,Ns,Nt,3)
      integer :: ip_x, ip_y, ip_t, np_x, np_y, np_t
      integer :: FNx, FNy, FNt
      integer :: sx,sy,st
      integer nprocs, rank, i
c      character*10 num

c      num=itoa(idx)
c      print *,num  
c      ThirringFile=trim(ThirringFileDir)//'theta'//trim(num)//'.dat'

      open(unit=10,file=fname,status='old')
        read(10,*) rank,nprocs
      close(10)
      open(unit=10,file=fname,status='old')
      do i=1,nprocs
        read(10,*) rank,nprocs
        read(10,*) ip_x,ip_y,ip_t,np_x,np_y,np_t
        read(10,*) FNx,FNy,FNt
        if ((FNx.ne.Ns).or.(FNy.ne.Ns).or.(FNt.ne.Nt)) then
          print *,"Thirring Gauge File not correct size"
          stop
        endif
        sx=FNx/np_x
        sy=FNy/np_y
        st=FNt/np_t
        read(10,*) theta(ip_x*sx+1:(ip_x+1)*sx,ip_y*sy+1:(ip_y+1)*sy,
     &                   ip_t*st+1:(ip_t+1)*st,1:3)
      end do
      close(10)

      return
      end subroutine readConvertedThirringGaugeField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine readThirringGaugeField2(idx,theta)
c      use arraysizes
c      implicit none
c      integer idx
c      real(prc) theta(Ns,Ns,Nt,3)
c      integer :: ip_x, ip_y, ip_t, np_x, np_y, np_t
c      integer :: FNx, FNy, FNt
c      integer :: sx,sy,st
c      integer nprocs, rank, i
c    
cc      stub='/home/jude/2019/Thirring/Sunbird/24x24/Ls30/Bp3/mp01/theta'
cc      stub='/home/jude/2019/Thirring/Sunbird/20x20/Ls24/Bp3/mp01/theta'
cc      stub='/home/jude/2019/Thirring/Sunbird/16x16/Ls24/Bp5/mp01/g3/thet
cc     &a'
cc      stub='/home/jude/2019/Thirring/Sunbird/16x16/Ls24/Bp6/mp01/theta'
cc      stub='/home/jude/2019/Thirring/Sunbird/12x24/Ls30/Bp3/mp01/theta'
cc      stub='/home/jude/2019/Thirring/Sunbird/converted/b_0.60/theta'
c      stub='/home/jude/2020/GaugeFields/8x8/Ls40/Bp28/mp005/g3/theta'
c      sdir='8x8/Ls40/Bp28/mp005/g3/'
c
c      if (idx < 10) then
c        write(sname,"(A56,I1,A4)") stub,idx,'.dat'
c!        write(sname,"(A58,I1,A4)") stub,idx,'.dat'
c!        write(sname,"(A61,I1,A4)") stub,idx,'.dat'
c      else if (idx < 100) then
c        write(sname,"(A56,I2,A4)") stub,idx,'.dat'
c!        write(sname,"(A58,I2,A4)") stub,idx,'.dat'
c!        write(sname,"(A61,I2,A4)") stub,idx,'.dat'
c      else if (idx < 1000) then
c        write(sname,"(A56,I3,A4)") stub,idx,'.dat'
c!        write(sname,"(A58,I3,A4)") stub,idx,'.dat'
c!        write(sname,"(A61,I3,A4)") stub,idx,'.dat'
c      else if (idx < 10000) then
c        write(sname,"(A56,I4,A4)") stub,idx,'.dat'
c!        write(sname,"(A58,I4,A4)") stub,idx,'.dat'
c!        write(sname,"(A61,I4,A4)") stub,idx,'.dat'
c      endif
c
c      print *,sname
c
c      open(unit=10,file=sname,status='old')
c        read(10,*) rank,nprocs
c      close(10)
c      open(unit=10,file=sname,status='old')
c      do i=1,nprocs
c        read(10,*) rank,nprocs
c        read(10,*) ip_x,ip_y,ip_t,np_x,np_y,np_t
c        read(10,*) FNx,FNy,FNt
c        if ((FNx.ne.Ns).or.(FNy.ne.Ns).or.(FNt.ne.Nt)) then
c          print *,"Thirring Gauge File not correct size"
c          stop
c        endif
c        sx=FNx/np_x
c        sy=FNy/np_y
c        st=FNt/np_t
c        read(10,*) theta(ip_x*sx+1:(ip_x+1)*sx,ip_y*sy+1:(ip_y+1)*sy,
c     &                   ip_t*st+1:(ip_t+1)*st,1:3)
c      end do
c      close(10)
c
c      return
c      end subroutine readThirringGaugeField2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine readMScQEDGaugeField(idx,theta) ! ftype=3
!     reads output from serial codes
      use arraysizes
      implicit none
      integer idx
      real(prc) theta(Ns,Ns,Nt,3)
      integer :: ip_x, ip_y, ip_t, np_x, np_y, np_t
      integer :: FNx, FNy, FNt
      integer :: sx,sy,st
      integer nprocs, rank, i
      character*1024 sname,stub
    
      stub='/home/jude/TRANSFER/Runs/Bp5/12x12/GAUGE/T3_20/theta'

      if (idx < 10) then
        write(sname,"(A52,I1,A4)") stub,idx,'.dat'
      else if (idx < 100) then
        write(sname,"(A52,I2,A4)") stub,idx,'.dat'
      else if (idx < 1000) then
        write(sname,"(A52,I2,A4)") stub,idx,'.dat'
      endif

      print *,sname

      open(unit=10,file=sname,status='old')
      read(10,*) theta
      close(10)

c      sname='slice1.dat'
c      call writeGaugeSlice(sname,theta)

      return
      end subroutine readMScQEDGaugeField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine readConFile(theta)
!     this reads con files outputted by the adjusted production code
!     and does not require mpi to read
c      use arraysizes
c      implicit none
c      integer idx
c      real(prc) theta(Ns,Ns,Nt,3)
c      real*4 stheta(Ns,Ns,Nt,3),seed
c    
c      open(unit=11,file='con',status='unknown',form='unformatted')
c      read(11) stheta,seed
c      close(11)
c      theta=stheta
c      return
c      end subroutine readConFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine readMPIConFile(fname,theta) ! ftype=1
!     this reads con files as outputted by the production code
      use arraysizes
#ifdef PARALLEL
      use basicparallelmod
      use mpi
#endif
      implicit none
      character(len=80) :: fname
      real(prc) theta(Ns,Ns,Nt,3)
      real stheta(Ns,Ns,Nt,3),seed
      integer MPI_THETA_TYPE,fht,ierr
      logical FEXIST   
 
      INQUIRE(file=fname,exist=FEXIST)
      if (.not.FEXIST) then
        print *,"cannot find file ",fname
        stop
      endif
#ifdef PARALLEL
      ! MPI array storage datatype
      call MPI_Type_Create_Subarray(4, ! dimensionality
     &             (/Ns,Ns,Nt,3/), ! global volume
     &             (/Ns,Ns,Nt,3/), ! local volume
     &             (/0,0,0,0/), ! start location
     &             MPI_Order_Fortran, ! array ordering
     &             MPI_Real, ! datatype to store
     &             MPI_THETA_TYPE, ierr) ! type descriptor for this subarray type
      call MPI_Type_Commit(MPI_THETA_TYPE, ierr)
          
      call MPI_FILE_OPEN(MPI_COMM_SELF,fname,MPI_MODE_RDONLY, 
     &               MPI_INFO_NULL, fht, ierr)
      call MPI_FILE_SET_VIEW(fht, 0_8, MPI_REAL, MPI_THETA_TYPE,  
     &                      'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_ALL(fht,stheta,Ns*Ns*Nt*3,MPI_REAL, 
     &                                  MPI_STATUS_IGNORE,ierr)
      call MPI_FILE_CLOSE(fht,ierr)

      theta=stheta

#else

      print *,"cannot open parallel type con file without parallel code"
      stop

#endif
      return
      end subroutine readMPIConFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module IOmodule

