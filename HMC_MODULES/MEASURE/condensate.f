      module condensatemodule
      use pacc
      use arraysizes
      use numbers
      use options
      use IOmodule
      use gaugefield
      use condnoisytoolsmod
      use condensatetoolsmod
      use ratfuncs
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine measureCondensate()
!     all inputs should be taken from input file
      implicit none
      complex(prc) pbp
      type(sgnratfunc) :: SRF
      integer fstart,fstop,fskip,nnoise
      complex(prc) cpbp
      integer fidx,inoise,np
      character(len=80) fname 
      logical FEXISTS
      integer Nterms
      logical ZOLO
      real(prc) lmin,lmax

      open(unit=31,file='noisycondopts.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) fstart,fstop,fskip,nnoise
        close(31)
        print *,fstart,fstop,fskip,nnoise
      else
        print *,"Condensate options file not found"
        stop
      endif
      close(31)

      open(unit=31,file='olopts.txt',status='unknown')
      inquire(unit=31,opened=FEXISTS)
      if (FEXISTS) then
        read(31,*) dwkernel,MTYPE,baremass,Nterms,ZOLO,lmin,lmax
        close(31)
        print *,dwkernel,MTYPE,baremass,Nterms,ZOLO,lmin,lmax
      else
        print *,"Overlap options file not found"
        stop
      endif
      close(31)
      if (ZOLO) then
        call setZoloCoeffs(Nterms,SRF,lmin,lmax)
      else
        call setHTcoeffs(Nterms,SRF)
      endif

      call measureNoisyCondensate(pbp,SRF,fstart,fstop,fskip,nnoise)

      open(unit=24,file='out.dat',status='unknown')
      write(24,*) "pbp:",pbp
      write(24,*) "kernel:",dwkernel
      write(24,*) "MTYPE:",MTYPE
      close(24)     

      return
      end subroutine measureCondensate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine condensatetask(fid,pbp,nnoise,SRF)
      implicit none
      integer fid
      complex(prc) pbp
      integer nnoise
      type(sgnratfunc) :: SRF
      integer inoise
      complex(prc) cpbp
  
      call readThetaFromFile(fid,theta)
      call coef(u,theta)
      print *,"condensatetask, fid:",fid,"nnoise:",nnoise
      pbp=czero
      do inoise=1,nnoise
        if (dwkernel.eq.1) then
          call evalCondNoisy_KDDW4(u,cpbp)
        elseif (dwkernel.eq.2) then
          call evalCondNoisy_OL(u,cpbp,SRF)
        endif
        pbp=pbp+cpbp
!        print *,"cpbp:",cpbp
!        print *,"pbp:",pbp
        write(88,*) inoise,real(cpbp,prc),dimag(cpbp)
        flush(88)
      enddo
      pbp=pbp/nnoise

      return
      end subroutine condensatetask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine qcondensatetask(pbp,nnoise,strength,SRF)
      implicit none
      complex(prc) pbp
      integer nnoise
      real(prc) strength
      type(sgnratfunc) :: SRF
      integer inoise
      complex(prc) cpbp
  
      call setGRVs(3*Nv,theta)
      theta=theta*strength
      call coef(u,theta)
      print *,"qcondensatetask, nnoise:",nnoise,"strength:",strength
      pbp=czero
      do inoise=1,nnoise
        if (dwkernel.eq.1) then
          call evalCondNoisy_KDDW4(u,cpbp)
        elseif (dwkernel.eq.2) then
          call evalCondNoisy_OL(u,cpbp,SRF)
        endif
        pbp=pbp+cpbp
!        write(88,*) inoise,real(cpbp,prc),dimag(cpbp)
!        flush(88)
      enddo
      pbp=pbp/nnoise

      return
      end subroutine qcondensatetask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine condensatesubtask(fid,pbp,didx,SRF)
      use cgcountmodule
      implicit none
      integer fid
      complex(prc) pbp
      integer nnoise
      type(sgnratfunc) :: SRF
      integer didx
      complex(prc) cpbp
  
      call readThetaFromFile(fid,theta)
      call coef(u,theta)
      print *,"condensatetask, fid:",fid
      print *,"didx:",didx
      pbp=czero
      if (dwkernel.eq.1) then
c        call evalCondNoisy_KDDW4(u,cpbp)
        print *,"not implemented"
        stop
      elseif (dwkernel.eq.2) then
        print *,"rank didx:",didx
        call evalSubCondNoisy_OL(u,cpbp,SRF,didx)
!        print *,"check evalSubCondNoisy has suitable didx range now"
!        stop
      endif
      pbp=pbp+cpbp
      write(88,'(2I9,E15.5,I5,2E19.10)') ncgits,ncgmax,conv,didx,
     &                          real(cpbp,prc),dimag(cpbp)
      flush(88)

      return
      end subroutine condensatesubtask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine measureNoisyCondensate(pbp,SRF,fstart,fstop,fskip,
     &                                                          nnoise)
!     calculate converted condensate 
      implicit none
      complex(prc) pbp
      type(sgnratfunc) :: SRF
      integer fstart,fstop,fskip,nnoise
      complex(prc) cpbp,npbp
      integer idx,inoise,np
      character(len=80) fname 

      pbp=czero
      np=0
      do idx=fstart,fstop,fskip
        call readThetaFromFile(idx,theta)
        call coef(u,theta)
        npbp=czero
        do inoise=1,nnoise
          if (dwkernel.eq.1) then
            call evalCondNoisy_KDDW4(u,cpbp)
          elseif (dwkernel.eq.2) then
            call evalCondNoisy_OL(u,cpbp,SRF)
c          elseif (dwkernel.eq.3) then
c            call evalCondNoisy_DomWall_Shamir(u,cpbp)
          endif
          pbp=pbp+cpbp
          npbp=npbp+cpbp
          np=np+1
          open(unit=25,file='full.dat',status='unknown',access='append')
          write(25,*) np,cpbp
          close(25)
        end do
        npbp=npbp/nnoise
        open(unit=23,file='cond.dat',status='unknown',access='append')
        write(23,*) idx,npbp
        close(23)
      end do
      pbp=pbp/np

      return
      end subroutine measureNoisyCondensate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine measureCondensateDomWall()
!     calculate condensate using the standard Shamir domain wall formulation
!     using Thirring gauge fields
      use options
      use gaugefield
      use domainwallmod
      use IOmodule
      use statsmod
      implicit none
      integer,parameter :: Ngf=50
      complex(prc) pbp,pbpvec(50)
      integer i,idx
      complex(prc) mean,std
      character(len=80) fname 

      print *,'Measure Thirring Condensate'
      open(unit=11,file='condDomWall.dat',
     &                         status='unknown',form='formatted')
      open(unit=12,file='condDomWall_Full.dat',
     &                         status='unknown',form='formatted')
      idx=0
      do i=5,5*Ngf,5
c        call readThetaFileName(i,fname)
c        call readMPIConFile(fname,theta)
        call readThetaFromFile(idx,theta)
c        call readThirringConFile(fname,theta) 
        call coef(u,theta)
        idx=idx+1
        call evalCondensateDDW2(u,pbp,12) 
        pbpvec(idx)=pbp
        print *,"pbp",pbp
        write(11,*) i,real(pbp,prc),dimag(pbp)
        flush(11)
      end do
      call calcVar(Ngf,pbpvec,mean,std)
      close(11)
      close(12)

      return
      end subroutine measureCondensateDomWall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine measureCondensateGW()
!     calculate condensate 
      use options
      use ratfuncs
      use gaugefield
      use gaugemodule
      use overlapmoduledev
      use domainwallmod
      use IOmodule
      use statsmod
      implicit none
      type(sgnratfunc) :: SRF
      complex(prc) pbp
      integer i
      integer SWITCH ! 1=point method 2=noisy
      procedure(),pointer :: IMptr => NULL()
      integer xskip,yskip,zskip
      complex(prc) pbpVec(4*Nv)
      integer Npbp
      complex(prc) mean,stdev
      character(len=80) fname 

      OLTYPE=1 ! 1 =direct,2=indirect
      DWkernel=2 ! 1=Shamir,2=Wilson
      MTYPE=3 ! 1=standard,3=alt3, 2=transpose alt3
c      ThirringFileDir=
c     &           '/home/jude/2019/Thirring/Sunbird/converted/b_0.30/'
c      ThirringCase='8x8/Ls40/Bp30/mp01/g3/'
c      ThirringFileDir='/home/jude/2020/GaugeFields/'//ThirringCase

      SWITCH=1

      call readConvertedThetaFileName(1,fname)
      call readConvertedThirringGaugeField(fname,theta)
      call coef(u,theta)
      print *,"measure GW condensate"

      IMptr => IKDDW4
      if (OLTYPE.eq.1) then
        call setHTcoeffs(Ls,SRF)
      endif

c      do i=1,4
        xskip=Ns !1 !Ns/i
        yskip=Ns !1 !Ns/i
        zskip=Ns !1 !Nt/i
        call evalSingleCondensateOL(u,baremass,SWITCH,pbp,SRF,IMptr,
     &         xskip,yskip,zskip,pbpVec,Npbp)
        print *,Npbp
        call calcVar(Npbp,pbpVec(1:Npbp),mean,stdev)
        open(unit=18,file='GWcondG3Ls90WilsonAll.dat',status='unknown')
        write(18,*) pbpVec
        close(18)
        
        print *,"pbpVec",pbpVec(1:Npbp)
        print *,"mean",mean,"sd",stdev,"err",(stdev*stdev/Npbp)**half
        open(unit=17,file='WG3Ls90FullDirect.dat',status='unknown',
     &                                          access='append')
!        write(17,*) Npbp, 
!     &            realpart(mean),imagpart(mean),
!     &            realpart(stdev),imagpart(stdev), 
!     &            realpart((stdev*stdev/Npbp)**half),
!     &            imagpart((stdev*stdev/Npbp)**half)
        write(17,*) Npbp, 
     &            real(mean,prc),dimag(mean),
     &            real(stdev,prc),dimag(stdev), 
     &            real((stdev*stdev/Npbp)**half,prc),
     &            dimag((stdev*stdev/Npbp)**half)
      close(17)
c      end do

      return
      end subroutine measureCondensateGW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef STAGGERED
      subroutine measureCondensateGN(nbos)
      use staggeredmodule
      use scalargaugemodule
      use statsmod
      implicit none
      integer nbos
      real(prc) sigma(Nv)
      logical FEXIST
      real(prc) pbp,pbpvecE(nbos),pbpVecN(nbos)
      real(prc) mean,sd,meanE,sdE,meanN,sdN
      real(prc) sig,sigvec(nbos)
      integer k

      ! check if there is a sigma file to open
      INQUIRE(FILE="sigma.dat",EXIST=FEXIST)
      ! make initial sigma array      
      if (FEXIST) then
        open(unit=22,file='sigma.dat',status='unknown',form='formatted')
        read(22,*) sigma
        close(22)
      else
        sigma=one
        do k=1,30
          call makeScalarBosonField(.false.,sigma)
        end do
      end if
      ! calculate condensate
      do k=1,nbos
        call makeScalarBosonField(.false.,sigma)
c        call calcStagGNCondensate(sigma,pbp)
c        pbpvecE(k)=pbp
        call calcStagGNCondensateNoisy(sigma,pbp)
        pbpvecN(k)=pbp
        sig=sum(sigma)/Nv
        sigvec(k)=sig
      end do

c      call calcVarReal(nbos,pbpvecE,meanE,sdE)
c      print *,"pbp:",meanE
c      print *,"sd:",sdE
c      print *,"err:",sdE/sqrt(real(nbos))
      call calcVarReal(nbos,pbpvecN,meanN,sdN)
      print *,"pbp:",meanN
      print *,"sd:",sdN
      print *,"err:",sdN/sqrt(real(nbos))
      call calcVarReal(nbos,sigvec,mean,sd)
      print *,"sig:",mean
      print *,"sd:",sd
      print *,"err:",sd/sqrt(real(nbos))


      open(unit=23,file='pbpN.dat',status='unknown',form='formatted')
      write(23,*) 1,pbpvecN(1),sigvec(1),meanN,sdN,mean,sd
      do k=2,nbos
        write(23,*) k,pbpvecN(k),sigvec(k)
      end do
c      open(unit=23,file='pbpE.dat',status='unknown',form='formatted')
c      do k=1,nbos
c        write(23,*) k,pbpvecE(k)
c      end do
      close(unit=23)
      open(unit=22,file='sigma.dat',status='unknown',form='formatted')
      write(22,*) sigma
      close(22)

      return
      end subroutine measureCondensateGN
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module condensatemodule
