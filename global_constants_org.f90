!====================================================================!
      module  Com_var
!--------------------------------------------------------------------!
!     Definitions of parameters.                                     !
!--------------------------------------------------------------------!
      implicit none
      integer :: NVAR, T=1, nmax=5, Lmax=8
      integer, save :: FFN=12    ! The first file number 
      real*8, save  :: epsr = 1.0d-15
      real*8  ::epsa, epsr_manip, dxout, xmax, dE
      real*8  :: r0, Mn, hbar
      real*8 :: A
      real*8  :: V0, Vso, ALPHA, Z
      logical :: Show_Title         = .false., &
               & Show_Results       = .true., &
               & Make_Wavefunction  = .true., &
               & Make_Energy_Scheme = .true.
      character, save :: C(0:9)
      data C / 's','p','d','f','g','h','i','j','k','l' /
      parameter(NVAR=2, epsa=1.0d-300, epsr_manip=1.0d-10)
      parameter(xmax=16.0d0, dxout=0.05d0, dE=6.0d0, Z=82.0d0)
      parameter(r0=1.27d0, Mn=939.6d0, hbar=197.3d0)
      parameter(A=208.0d0, V0=44.0d0, Vso=19.4d0, ALPHA=0.67d0)
      end module 
!====================================================================!
      module Cash_Karp_Constants
!--------------------------------------------------------------------!
!     Definition of constants which appear in the Cash_Karp's        !
!     formula.                                                       !
!     Since substitution is not allowed in module subprogram and     !
!     operations such as division are not allowd in data sentence,   !
!     we have to define temporary variables to initialize arrays     !
!     used in other program units.                                   !
!--------------------------------------------------------------------!
      implicit none
      real*8, public  :: w(1:6), dw(1:6), alph(2:6), beta(2:6, 5)
      real*8, private ::  w1,  w2,  w3,  w4,  w5,  w6
      real*8, private :: dw1, dw2, dw3, dw4, dw5, dw6
      real*8, private :: alph2,  alph3,  alph4,  alph5,  alph6
      real*8, private :: beta21, beta22, beta23, beta24, beta25,     &
                       & beta31, beta32, beta33, beta34, beta35,     &
                       & beta41, beta42, beta43, beta44, beta45,     &
                       & beta51, beta52, beta53, beta54, beta55,     &
                       & beta61, beta62, beta63, beta64, beta65
      parameter(                                                     &
        & w1 =  37.0d0/378.0d0, w2 =   0.0d0,                        &
        & w3 = 250.0d0/621.0d0, w4 = 125.0d0/ 594.0d0,               &
        & w5 =   0.0d0,         w6 = 512.0d0/1771.0d0 )
      parameter(                                                     &
        & dw1 = - 277.0d0/ 64512.0d0, dw2 =      0.0d0,              &
        & dw3 =  6925.0d0/370944.0d0, dw4 = - 6925.0d0/202752.0d0,   &
        & dw5 = - 277.0d0/ 14336.0d0, dw6 =    277.0d0/  7084.0d0 )
      parameter(                                                     &
        & alph2 = 0.2d0,   alph3 = 0.3d0,                            &
        & alph4 = 0.6d0,   alph5 = 1.0d0,                            &
        & alph6 = 0.875d0 )
      parameter(                                                     &
        & beta21 =    0.2d0,           beta31 = 0.075d0,             &
        & beta32 =  0.225d0,           beta41 =   0.3d0,             &
        & beta42 = -  0.9d0,           beta43 =   1.2d0,             &
        & beta51 = - 11.0d0/   54.0d0, beta52 =   2.5d0,             &
        & beta53 = - 70.0d0/   27.0d0, beta54 =  35.0d0/27.0d0,      &
        & beta61 = 1631.0d0/55296.0d0, beta62 = 0.341796875d0,       &
        & beta63 =  575.0d0/13824.0d0, beta64 = 44275.0d0/110592.0d0,&
        & beta65 =  0.061767578125d0 )
      data w    / w1,  w2, w3,  w4,  w5,  w6/
      data dw   /dw1, dw2,dw3, dw4, dw5, dw6/
      data alph /alph2, alph3, alph4, alph5, alph6/
      data beta(2,:) /beta21, 4*0.0d0/
      data beta(3,:) /beta31, beta32, 3*0.0d0/
      data beta(4,:) /beta41, beta42, beta43, 2*0.0d0/
      data beta(5,:) /beta51, beta52, beta53, beta54, 0.0d0/
      data beta(6,:) /beta61, beta62, beta63, beta64, beta65/

      end module
