      module remezapproxmod
!     rat polys: p(x)/q(x) ... p(x)=p0+p1.x+p2.x^2+..., q(x)=1+q1.x+q2.x^2+...
      use pacc
      implicit none

      contains

      subroutine getRemezCoeffs(func,op,oq,r1,r2,front,alpha,beta)
      implicit none
      integer func,op,oq
      real(prc) r1,r2
      real(prc),allocatable,dimension(:) :: p,q
      real(prc) front
      real(prc),allocatable,dimension(:) :: alpha,beta


      if(allocated(alpha)) then
        deallocate(alpha)
      endif
      if(allocated(beta)) then
        deallocate(beta)
      endif
      allocate(p(op+1),q(oq+1))
      allocate(alpha(op),beta(oq))

      if (func.eq.1) then ! 1/(x+2) - for Wilson range [1,5]
      elseif (func.eq.2) then ! 1/x^-0.5 - for Wilson range [1,5] 
      elseif (func.eq.3) then ! 1/x^0.5 - for Wilson range [1,25]

        if((op.eq.1).and.(oq.eq.1)) then ! weighted - e=3e-2

          p(1)=0.0894
          p(2)=1.5615
          q(1)=0.6986
          q(2)=1.0000

          front = 1.2808356311029631e-01
          alpha(1) = 2.0536046979873319e+00 
          beta(1) = 1.4314529429411673e+00


        elseif((op.eq.2).and.(oq.eq.2)) then ! weighted - e=1e-3

          p(1)=0.0179
          p(2)=1.1268
          p(3)=2.6040

          q(1)=0.2329
          q(2)=2.5197
          q(3)=1.0000

          front = 7.6804735054408699e-02
          alpha(1) = 9.2035205599151559e-01
          beta(1) = 4.1261136743482790e-01
          alpha(2) = 3.0868514372357883e+00
          beta(2) = 1.0405706860006751e+01

        elseif((op.eq.3).and.(oq.eq.3)) then ! weighted - e=4e-5

          front = 5.4860481441453827e-02
          alpha(1) = 6.0283006748355317e-01
          beta(1) = 1.9796761371608651e-01
          alpha(2) = 1.1564810350194275e+00
          beta(2) = 2.9717201118091126e+00
          alpha(3) = 4.0789792753288223e+00
          beta(3) = 2.6064690825592852e+01

        elseif((op.eq.4).and.(oq.eq.4)) then ! weighted - e=2e-6

          front = 4.2669263297759358e-02
          alpha(1) = 4.5196735109760777e-01
          beta(1) = 1.1676093783501723e-01
          alpha(2) = 6.8412993392797605e-01
          beta(2) = 1.4314529429404295e+00
          alpha(3) = 1.4169742296201886e+00
          beta(3) = 7.4865322768172540e+00
          alpha(4) = 5.0740529160227243e+00
          beta(4) = 4.7656792862619859e+01

        elseif((op.eq.5).and.(oq.eq.5)) then ! weighted - e=5e-8

          front = 3.4911215425473971e-02
          alpha(1) = 3.6289021021675444e-01
          beta(1) = 7.7164772476385549e-02
          alpha(2) = 4.8345810075591161e-01
          beta(2) = 8.5435964811180409e-01
          alpha(3) = 8.0036748313499995e-01
          beta(3) = 3.5951214784654200e+00
          alpha(4) = 1.6682583257060843e+00
          beta(4) = 1.3737039467286907e+01
          alpha(5) = 6.0817739785982941e+00
          beta(5) = 7.4926158879596414e+01

        elseif((op.eq.6).and.(oq.eq.6)) then ! weighted - e=2e-9

          front = 2.9540259207574544e-02
          alpha(1) = 3.0372719693732420e-01
          beta(1) = 5.4841271758336425e-02
          alpha(2) = 3.7466407478857572e-01
          beta(2) = 5.7260802506584174e-01
          alpha(3) = 5.4526527553917092e-01
          beta(3) = 2.1386238187959967e+00
          alpha(4) = 9.1620827089275680e-01
          beta(4) = 6.6084560779080048e+00
          alpha(5) = 1.9158563317687349e+00
          beta(5) = 2.1576983786437690e+01
          alpha(6) = 7.1003065519108430e+00
          beta(6) = 1.0777778135198650e+02

        elseif((op.eq.7).and.(oq.eq.7)) then ! weighted - e=7e-11
   
          front = 2.5601558004891313d-02
          alpha(1) = 2.6142544421313996d-01
          beta(1) = 4.1001224144249553d-02
          alpha(2) = 3.0678377396165396d-01
          beta(2) = 4.1261136573518747d-01
          alpha(3) = 4.1047795974121210d-01
          beta(3) = 1.4314529375511613d+00
          alpha(4) = 6.1126370296898858d-01
          beta(4) = 3.9267913598765207d+00
          alpha(5) = 1.0289496592779273d+00
          beta(5) = 1.0405706830410381d+01
          alpha(6) = 2.1633324415411654d+00
          beta(6) = 3.0921108869174930d+01
          alpha(7) = 8.1268678729269208d+00
          beta(7) = 1.4617121629726213d+02

        elseif((op.eq.8).and.(oq.eq.8)) then ! weighted - e=3e-12

          front = 2.2589610473336887d-02
          alpha(1) = 2.2960822735957476d-01
          beta(1) = 3.1822519142867566d-02
          alpha(2) = 2.6040271520000069d-01
          beta(2) = 3.1244033469565108d-01
          alpha(3) = 3.2858918209872684d-01
          beta(3) = 1.0319906118614628d+00
          alpha(4) = 4.5201276475200192d-01
          beta(4) = 2.6236251650802704d+00
          alpha(5) = 6.7643407233603114d-01
          beta(5) = 6.1877025890961548d+00
          alpha(6) = 1.1397643529852599e+00
          beta(6) = 1.4936230366691795d+01
          alpha(7) = 2.4120586892255740d+00
          beta(7) = 4.1720658787606439d+01
          alpha(8) = 9.1592877759912241d+00
          beta(8) = 1.9008714654958371d+02

        elseif((op.eq.9).and.(oq.eq.9)) then ! weighted - e=9e-14

          front = 2.0211765914238990d-02
          alpha(1) = 2.0477421138925295d-01
          beta(1) = 2.5420248256542183d-02
          alpha(2) = 2.2665249970546891d-01
          beta(2) = 2.4529733490998279d-01
          alpha(3) = 2.7406374141954942d-01
          beta(3) = 7.8281025264931647d-01
          alpha(4) = 3.5626303924177727d-01
          beta(4) = 1.8887073545308461d+00
          alpha(5) = 4.9427355087460145d-01 
          beta(5) = 4.1320957876700213d+00
          alpha(6) = 7.4022767062799999d-01
          beta(6) = 8.8934272406582853d+00
          alpha(7) = 1.2497724351720387d+00
          beta(7) = 2.0164108897484990d+01
          alpha(8) = 2.6624062437493525d+00
          beta(8) = 5.3947040405792862d+01
          alpha(9) = 1.0196026326446450d+01
          beta(9) = 2.3951526661449941d+02

        else
          print *,"bad choice of epoly orders"
          stop
        endif
      elseif (func.eq.4) then ! x^-0.5 over [0.2,8]


        if((op.eq.1).and.(oq.eq.1)) then ! weighted - e=4e-2

          front = 2.4375061489592584e-01
          alpha(1) = 1.0411532295794446e+00
          beta(1) = 3.4648028417661514e-01

        elseif((op.eq.2).and.(oq.eq.2)) then ! weighted - e=2e-3

          front = 1.4607044680863407e-01
          alpha(1) = 4.5316508496797248e-01
          beta(1) = 9.7163369581680176e-02
          alpha(2) = 1.6292339444185830e+00
          beta(2) = 2.7052598353989885e+00

        elseif((op.eq.3).and.(oq.eq.3)) then ! weighted - e=9e-5

          front = 1.0433575158630061e-01
          alpha(1) = 2.9372134061200883e-01
          beta(1) = 4.6236220271923048e-02
          alpha(2) = 5.9711555830425045e-01
          beta(2) = 7.3717257572479611e-01
          alpha(3) = 2.1613851122274363e+00
          beta(3) = 6.9677349589445070e+00

        elseif((op.eq.4).and.(oq.eq.4)) then ! weighted - e=4e-6

          front = 8.1150028530492524e-02
          alpha(1) = 2.1916984574796936e-01
          beta(1) = 2.7175504159093519e-02
          alpha(2) = 3.4662318419611321e-01
          beta(2) = 3.4648028417661514e-01
          alpha(3) = 7.4497297807967311e-01
          beta(3) = 1.9232042092093882e+00
          alpha(4) = 2.6859471434756075e+00
          beta(4) = 1.2902596034085194e+01

        else
          print *,"bad choice of epoly orders"
          stop
        endif

      else
        print *,"bad choice of approximation function"
        stop
      endif

      return
      end subroutine getRemezCoeffs


      end module remezapproxmod

