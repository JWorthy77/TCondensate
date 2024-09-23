      module extpolyapproxmod
!     extended polys: c.x^-ol + ... + c.x^ou
      use pacc
      implicit none

      contains

      subroutine getExtPolyCoeffs(func,ol,ou,r1,r2,cfs)
      implicit none
      integer func,ol,ou
      real(prc) r1,r2
      real(prc),allocatable,dimension(:) :: cfs

      if(allocated(cfs)) then
        deallocate(cfs)
      endif
      allocate(cfs(ou-ol+1))

      if (func.eq.1) then ! 1/(x+2) - for Wilson range [1,5]
      elseif (func.eq.2) then ! 1/x^-0.5 - for Wilson range [1,5] 
        if((ol.eq.-2).and.(ou.eq.4)) then ! weighted - e=2e-5
          cfs(7)=0.000111829337855363
          cfs(6)=-0.00234686847754203
          cfs(5)=0.0213268050112455
          cfs(4)=-0.11698705563828
          cfs(3)=0.615883871126285
          cfs(2)=0.545537863494807
          cfs(1)=-0.0635396096916379
        elseif((ol.eq.-2).and.(ou.eq.2)) then ! weighted - e=3e-4
          cfs(5)=0.00248352661782352
          cfs(4)=-0.0430467797453998
          cfs(3)=0.466874924647038
          cfs(2)=0.691556000056989
          cfs(1)=-0.118097732919246
        elseif((ol.eq.-3).and.(ou.eq.3)) then ! weighted - e=8e-6
          cfs(7)=-0.000213186911060525
          cfs(6)=0.00467103481186696
          cfs(5)=-0.0482985174078031
          cfs(4)=0.454748581178275
          cfs(3)=0.760212656147282
          cfs(2)=-0.213966092313384
          cfs(1)=0.0428530884466567
        elseif((ol.eq.-4).and.(ou.eq.4)) then ! weighted - e=3e-7
          cfs(9)=1.92268732272259e-05
          cfs(8)=-0.000516532342509109
          cfs(7)=0.00636979675564852
          cfs(6)=-0.0513959321595631
          cfs(5)=0.448126618292868
          cfs(4)=0.800105218507219
          cfs(3)=-0.285681515694804
          cfs(2)=0.100788006372203
          cfs(1)=-0.0178151476163209
        elseif((ol.eq.-5).and.(ou.eq.5)) then ! weighted - e=1e-8

          cfs(11)=-1.78411923150305e-06
          cfs(10)=5.69938916767493e-05
          cfs(9)=-0.000836742068894304
          cfs(8)=0.00768652117673837
          cfs(7)=-0.0534406315754478
          cfs(6)=0.443960445024787
          cfs(5)=0.826220761761855
          cfs(4)=-0.340125115357108
          cfs(3)=0.160123744249576
          cfs(2)=-0.051521757153553
          cfs(1)=0.00787757343157861

        endif
      elseif (func.eq.3) then ! 1/x^0.5 - for Wilson range [1,25]
        if((ol.eq.-3).and.(ou.eq.8)) then ! weighted - e=2e-5
          cfs(12)=1.01780708427028e-11
          cfs(11)=-1.34883597541861e-09
          cfs(10)=7.76812683512139e-08
          cfs(9)=-2.56001526755114e-06
          cfs(8)=5.36382462990157e-05
          cfs(7)=-0.000754892738480525
          cfs(6)=0.00745112442430497
          cfs(5)=-0.0557621310192142
          cfs(4)=0.464548296391015
          cfs(3)=0.756119258583077
          cfs(2)=-0.216611644899582
          cfs(1)=0.044973813689521
        elseif((ol.eq.-2).and.(ou.eq.2)) then ! weighted - e=6e-3

          cfs(5)=0.000170992288225109
          cfs(4)=-0.0099728827925868
          cfs(3)=0.30441262228137
          cfs(2)=0.992843122977328
          cfs(1)=-0.29347831585989

        elseif((ol.eq.-3).and.(ou.eq.3)) then ! weighted - e=8e-4

          cfs(7)=-5.17408682771383e-06
          cfs(6)=0.000384698861267793
          cfs(5)=-0.0120729477454822
          cfs(4)=0.298689860625816
          cfs(3)=1.10796820938377
          cfs(2)=-0.585009617380238
          cfs(1)=0.190781112457482

        elseif((ol.eq.-4).and.(ou.eq.4)) then ! weighted - e=1e-4

          cfs(9)=1.64698386294396e-07
          cfs(8)=-1.50713894616698e-05
          cfs(7)=0.000583747302207984
          cfs(6)=-0.0133948437616912
          cfs(5)=0.295510186905655
          cfs(4)=1.17437572215953
          cfs(3)=-0.821354680458909
          cfs(2)=0.505515668080175
          cfs(1)=-0.141315174458789

        else
          print *,"bad choice of epoly orders"
          stop
        endif
      else
      print *,"bad choice of approximation function"
      stop
      endif

      return
      end subroutine getExtPolyCoeffs




      end module extpolyapproxmod










