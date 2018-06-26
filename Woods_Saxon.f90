!====================================================================!
!     Title  : Woods_Saxon_LS.f90                                    !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-4-11-Fri ~ 2008-4-23-Thu                         !
!     Last modified : 2008-4-23-Thu                                  !
!                                                                    !
!     The shell model is a model which describes the structure of a  !
!     nucleus from the microscopic level. In this model, a nucleon   !
!     in the nucleus feels a potential which is created by other     !
!     nucleons including itself and they move independently in this  !
!     potential.                                                     !
!     The behavior of the nucleon is described by Schroedinger       !
!     equation with proper boundary conditions. By solving this      !
!     equation, we can obtain energy levels of the bound states and  !
!     wave functions of the nucleon.                                 !
!     The shape of the potential is looks like the Fermi-Dirac       !
!     distribution and is well described by Woods-Saxon potential.   !
!     For this potential, we cannot solve the Schroedinger equation  !
!     exactly and have to solve it numerically.                      !
!     In addition, we have to add the spin-orbit coupling term to    !
!     the potential to reproduce the magic numbers. This term depends!
!     on the relative orientation of spin and orbital angular        !
!     momentum and splits the energy levels into two except for L=0  !
!     state.                                                         !
!                                                                    !
!     In this program, we solve the Schroedinger equation with       !
!     proper boundary conditions and find energy levels and wave     !
!     functions. We compute the problem with and without LS coupling !
!     and see how the energy levels split to produce the magic       !
!     numbers.                                                       !
!     To find the energy levels, we solves the equation with energy  !
!     E and E+dE. If E is different from the one of the eigenvalue,  !
!     the solution diverges as x tends to infinity. The manner of    !
!     the divergence is related to the number of nodes of the wave   !
!     function. For example, if the true energy level E* is located  !
!     between E and E+dE, and E diverges to positive infinity at     !
!     large x, then the wave function with E+dE diverges toward      !
!     negative infinity. Therefore, we can see whether the E* is     !
!     situated between E and E+dE by observing the sign of the       !
!     product of the solutions with E and E+dE at large x. Once the  !
!     solution turned out to be in [E,E+dE], we can find it by using !
!     bisection method. After finding the energy level, we increment !
!     E and seek next energy level. We start from small E and        !
!     increment it up to E = 0.                                      !
!                                                                    !
!     *** Structure of this program ***                              !
!                                                                    !
!     program     main                                               !
!     subroutine  E_loop                                             !
!     module      sub_init                                           !
!         subroutine  Init_var     --> contained sub_init            !
!     module      Solve_Eq                                           !
!         subroutine  Solve        --> contained Solve_Eq            !
!     subroutine  Find_E                                             !
!     subroutine  Bisection                                          !
!     subroutine  Init_width                                         !
!     subroutine  Init_E                                             !
!     subroutine  Next_y                                             !
!                                                                    !
!     subroutine  func                                               !
!     function    f                                                  !
!     function    V                                                  !
!     function    Vpara                                              !
!     function    Vapara                                             !
!                                                                    !
!     subroutine  Standard_out                                       !
!     subroutine  line                                               !
!     subroutine  Show_Eq                                            !
!     subroutine  Title                                              !
!                                                                    !
!     subroutine  Waves                                              !
!     subroutine  WS_potential                                       !
!     subroutine  Gnuplot                                            !
!     subroutine  Gnuplot_head                                       !
!     subroutine  Gnuplot_label                                      !
!     subroutine  Gnuplot_arrow                                      !
!     subroutine  Gnuplot_plot                                       !
!     subroutine  Gnuplot_NoLS                                       !
!     subroutine  Gnuplot_LS                                         ! 
!====================================================================!
      program main
!--------------------------------------------------------------------!
!     Main program.                                                  !
!--------------------------------------------------------------------!
      use  Com_var
      implicit none
      integer :: L=0, n
      integer, allocatable :: P_num(:)
      real*8 , allocatable :: En(:,:,:)
      real*8  :: E

      if ( Show_Title ) then
          call Title
          call Show_Eq
      end if
      allocate( En(nmax, 0:Lmax, 3), P_num(0:Lmax) )
      En = 0.0d0
      !call WS_potential(11)

      do T=1, 3      
          call Init_E(0, E, En)
          do while(L <= Lmax)
              n = 0
              call E_loop(L, n, E, En)
              call Init_E(L, E, En)
              P_num(L) = n
              L = L + 1
          end do    
          L = 1
      end do    

      En(:,0,2) = En(:,0,1)

      if ( Show_Results ) call Standard_out(En)
      if ( Make_Energy_Scheme )  call Gnuplot(P_num, En, FFN, 'gnu')

      stop
      end program
!====================================================================!
      subroutine E_loop(L, n, E, En)
!--------------------------------------------------------------------!
!     This program computes energy levels for given 'L'.             !
!     After calling 'Find_E', we have obtained the energy level 'E'  !
!     and by calling 'Waves', we make '.dat' files in a certain      !
!     directory. Data of these files are the numerical values of the !
!     wave functions of each states.                                 !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(in) :: L
      integer, intent(inout) :: n
      real*8,  intent(inout) :: E, En(nmax, 0:Lmax, 3)
      real*8 :: E0

      do
          call Find_E(L, E)
          if (E >= 0.0d0) then
              E = E0
              return  
          end if
          n = n + 1
          if ( Make_Wavefunction ) call Waves(n, L, E)

          if (n > nmax) stop 'Error : n > nmax'

          En(n, L, T) = E
          E0 = E
          E  = E + 9.0d0
          if (E > 0.0d0) E = 0.0d0
      end do

      return
      end subroutine
!====================================================================!
      module sub_init
!--------------------------------------------------------------------!
!     A module subroutine which initializes variables.               !
!     Since the variable x is used in the denominator in some        !
!     program units, x must not be equal to 0.                       !
!     'yy' is optional.                                              !
!--------------------------------------------------------------------!
      implicit none
      contains
      subroutine  Init_var(L, x, xout, y, yy)
      use Com_var
      integer, intent(in) :: L
      real*8, intent(inout) :: x, y(NVAR), xout
      real*8, optional :: yy(NVAR)

      x = 1.0d-30
      xout = 0.0d0
      y(1) = x ** L
      y(2) = L * x ** (L + 1)
      if (present(yy)) yy = y

      end subroutine
      end module
!====================================================================!
      module Solve_Eq
!--------------------------------------------------------------------!
!     This module subroutine solves the differential equation for    !
!     given E and L. If 'yy' is 'present', we solves the equation    !
!     not only for E, but also for E+dE. If 'Fnum' is 'present', we  !
!     write down the numerical values of the wave function to the    !
!     file.                                                          !
!--------------------------------------------------------------------!
      contains
      subroutine  Solve(y, yy, L, E, Fnum, FNAME)
      use Com_var
      use sub_init
      implicit none
      integer, intent(in) :: L
      integer, intent(in), optional :: Fnum
      real*8,  intent(out) :: y(NVAR)
      real*8,  intent(out), optional :: yy(NVAR)
      real*8,  intent(in) :: E
      real*8 :: x2, h2, E2, h, x, xout
      character(len=*), intent(in), optional :: FNAME

      if ( present(Fnum) ) open(Fnum, file=FNAME)
      E2 = E + dE
      call Init_width(L, E, h)
      call Init_var(L, x, xout, y, yy)
      do while(xout < xmax)
          
          if (x >= xout) then
              if ( present(Fnum) ) write(Fnum,*) x, y(1)
              xout = xout + dxout
          end if
          
          x2 = x
          h2 = h
          call Next_y(x, y, h, L, E, xout)
          if ( present(yy) ) call Next_y(x2, yy, h2, L, E2, xout)

      end do
      if ( present(Fnum) ) close(Fnum)

      end subroutine
      end module
!====================================================================!
      subroutine  Find_E(L, E)
!--------------------------------------------------------------------!
!     This subroutine seeks the energy E such that the sign of the   !
!     product of wave functions with E and E+dE at xmax is negative. !
!     This means that the true energy level must be in [E, E+dE]     !
!     After that, we call the subroutine 'Bisection' to find the     !
!     energy level 'E'.                                              !
!--------------------------------------------------------------------!
      use Com_var
      use Solve_eq
      implicit none
      integer, intent(in) :: L
      real*8, intent(inout) :: E
      real*8 :: yy(NVAR), y(NVAR)

      epsr = 1.0d-1
      do
          call Solve(y=y, yy=yy, L=L, E=E)
            
          if (y(1) * yy(1) > 0.0d0) then
              E  = E + 0.5d0 * dE
              cycle
          else
              epsr = epsr_manip
              call Bisection(y, L, E)
              exit
          end if
      end do
      end subroutine
!====================================================================!
      subroutine Bisection(y, L, E)
!--------------------------------------------------------------------!
!     This subroutine seeks the energy level which is located        !
!     between E and E+dE by bisection method. We iterate the         !
!     procedure until E converges.                                   !
!     If E gives the positive y(1) at xmax before entering do loop,  !
!     E2 = E + dE gives negative y(1) at xmax. We define logical     !
!     variable 'num' to keep this situation.                         !
!--------------------------------------------------------------------!
      use Com_var
      use sub_init
      implicit none
      integer, intent(in) :: L
      real*8 :: x, xout
      real*8 :: half_E, yy(NVAR)
      real*8, intent(inout) :: E, y(NVAR)
      real*8 :: E2, eps, E0=0.0d0, h
      logical :: num

      num = (y(1) > 0.0d0)
      E2 = E + dE

 XX:  do 
          call Init_var(L=L, x=x, xout=xout, y=y)
          call Init_width(L, E, h)
          half_E = 0.5d0 * (E + E2)
          
          eps = epsa + epsr * (abs(E0) + abs(half_E))
          if (abs(E0 - half_E) < eps) then
              exit
          end if

          E0 = half_E
          do while (xout < xmax)
              if (x >= xout) xout = xout + dxout
              call Next_y(x, y, h, L, half_E, xout)
          end do

          if (y(1) > 0.0d0) then
              if (num) then
                  E = half_E
              else 
                  E2 = half_E
              end if
          else
              if (num) then
                  E2 = half_E
              else
                  E = half_E
              end if
          end if

      end do  XX
      E = half_E

      return
      end subroutine
!====================================================================!
      subroutine Init_width(L, E, h)
!--------------------------------------------------------------------!
!     Calculation of initial width 'h'.                              !
!     If this subrouinte is called for the first time, we have to    !
!     compute it by solving the equation. After this calculation, we !
!     save the value and use it in the following calling.            !
!--------------------------------------------------------------------!
      use  Com_var
      use  sub_init
      implicit none
      integer, intent(in) :: L
      integer, save :: L_temp=10
      real*8,  intent(in) :: E
      real*8 :: x, xout, y(NVAR)
      real*8, intent(out) :: h
      real*8, save :: h0, eps_temp
      logical, save :: yet=.true.

      yet = (L /= L_temp)
      L_temp = L
      if (yet) then
          eps_temp = epsr
          epsr = epsr_manip
          h = 1.0d0
          call Init_var(L=L, x=x, xout=xout, y=y)
          call Next_y(x, y, h, L, E, 0.1d0)
          epsr = eps_temp
          h0 = h
      else
          h = h0
      end if

      return
      end subroutine
!====================================================================!
      subroutine Init_E(L, E, En)
      use Com_var
!--------------------------------------------------------------------!
!     A module subroutine which initalizes the energy 'E'.           !
!     We initialize 'E' depending on when we initialize.             !
!     If this subroutine is called for the first time, then we set   !
!     'E' equal to -V0. If it is called with L = 0 and T=2 or 3(i.e. !
!     calculation with LS coupling), we set 'E' equal to the 1s      !
!     energy level computed without LS coupling. If it is called     !
!     other situation, we assign 'E' to the value that is the lowest !
!     one with L-1 and same T.                                       !
!--------------------------------------------------------------------!
      implicit none
      integer, intent(in)  :: L
      real*8,  intent(out) :: E
      real*8,  intent(in)  :: En(nmax,0:Lmax,3)
      logical, save :: yet = .true.

      if (yet) then
          E   = - V0
          yet = .false.
          return
      end if
      if (L == 0 .and. T /= 0) then
          E = En(1, 0, 1)
      else 
          E = En(1, L, T) - 0.1d0
      end if

      return
      end subroutine
!====================================================================!
      subroutine  Next_y(x, y, h, L, E, xout)
!--------------------------------------------------------------------!
!     We compute next step width 'h'.                                !
!     Since the step width is not fixed, we change 'h' to 'xout - x' !
!     if 'h' plus 'x' is greater than 'xout'. Then we can print out  !
!     the results at 'xout'.                                         !
!--------------------------------------------------------------------!
      use  Cash_Karp_Constants
      use  Com_var
      implicit none
      integer :: i, j, k
      integer, intent(in) :: L
      real*8, intent(in)    :: xout, E
      real*8, intent(inout) :: x, y(NVAR), h
      real*8 :: x0, y0(NVAR), kk(NVAR, 1:6)
      real*8 :: yerr, dy(NVAR), ynorm, eps

      if ( x + h > xout )  h = xout - x
      x0 = x
      y0 = y
      
      call func(x, y, L, E, kk(:,1))

  AA: do j=2, 6
          x = x0 + alph(j) * h
!$omp parallel do
   B:     do i=1, NVAR
              y(i) = Dot_product(beta(j,1:j-1), kk(i,1:j-1))
              y(i) = y0(i) + h * y(i)
          end do   B
!$omp end parallel do

          call func(x, y, L, E, kk(:,j))
      
      end do   AA

      x = x0 + h
      ynorm = 0.0d0
      yerr  = 0.0d0
  CC: do i=1, NVAR
          y(i)  = 0.0d0
          dy(i) = 0.0d0
   D:     do k=1, 6
              y(i)  = y(i)  +  w(k) * kk(i,k)
              dy(i) = dy(i) + dw(k) * kk(i,k)
          end do   D
          y(i)  = y0(i) + h * y(i)
          ynorm = max(ynorm, abs(y(i)))
          yerr  = max(yerr,  abs(dy(i)))
      end do   CC

      eps  = epsa + epsr * ynorm
      h    = 0.9d0 * h * (eps / (h * yerr) ) ** 0.2d0
      
      return
      end  subroutine
!===================================================================!
      subroutine func(x, y, L, E, f)
!--------------------------------------------------------------------!
!     Definition of the RHS of differential equations.               !
!     T = 0 corresponds to the potential without LS coupling.        !
!     T = 1 corresponds to the potential with LS coupling coupled    !
!     in parallel.                                                   !
!     T = 2 corresponds to the potential with LS coupling coupled    !
!     in antiparallel.                                               !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: i, L
      real*8, intent(inout) :: x, y(NVAR), f(NVAR)
      real*8, intent(in) :: E
      real*8 :: U
      real*8, external :: V, Vpara, Vapara

      select case(T)
          case(1); U = V(x)
          case(2); U = Vpara(L, x)
          case(3) 
              if (L == 0) then
                  U = V(x)
              else 
                  U = Vapara(L, x)
              end if
      end select

      if (x == 0) stop 'Division by 0.'

      f(1) = y(2) / (x * x)
      f(2) = L * (L + 1) * y(1) + 2.0d0 * Mn / (hbar * hbar)  &
           &   * x * x * (U - E) * y(1)

      return
      end subroutine
!====================================================================!
      function f(x)  result(FD)
!--------------------------------------------------------------------!
!     Definition of Fermi-Dirac function which appears in the Woods- !
!     Saxon potential.                                               !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      real*8, intent(in) :: x
      real*8 :: FD, R
      
      R = r0 * A ** (1.0d0/3.0d0)
      FD = 1.0d0 / (1.0d0 + exp((x - R) / ALPHA))

      return
      end function
!====================================================================!
      function V(x) result(WS)
!--------------------------------------------------------------------!
!     Definition of Woods-Saxon potential.                           !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      real*8, intent(in) :: x
      real*8, external :: f
      real*8 :: WS

      WS = - V0 * f(x)

      return
      end function
!====================================================================!
      function Vpara(L, x)  result(Vp)
!--------------------------------------------------------------------!
!     Definition of the potential when there is a LS coupling and    !
!     couples with 'parallel'.                                       !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(in) :: L
      real*8, intent(in) :: x
      real*8, external :: f, V
      real*8 :: Vp, R
      
      R  = r0 * A ** (1.0d0/3.0d0)
      Vp = V(x) - 0.5d0 * Vso * r0 * r0 * exp((x - R) / ALPHA) &
         &  * f(x) * f(x) * L / (ALPHA * x)

      return
      end function
!====================================================================!
      function Vapara(L, x)  result(Va)
!--------------------------------------------------------------------!
!     Definition of the potential when there is a LS coupling and    !
!     couples with 'antiparallel'.                                   !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(in) :: L
      real*8, intent(in) :: x
      real*8, external :: f, V
      real*8 :: Va, R
      
      R  = r0 * A ** (1.0d0/3.0d0)
      Va = V(x) + 0.5d0 * Vso * r0 * r0 * exp((x - R) / ALPHA) &
         &  * f(x) * f(x) * (L + 1) / (ALPHA * x)

      return
      end function
!====================================================================!
      subroutine Standard_out(En)
!--------------------------------------------------------------------!
!     Standard output.                                               !
!     We output the header and results of the calculation.           !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: L, i, J
      real*8, intent(in)   :: En(nmax,0:Lmax,3)
      real*8 :: En2(nmax,0:Lmax,3)
      integer :: Loc(3)
      character :: FM1*50, FM2*50
      parameter(FM1='(1x,a,i1,a,a,i1,a,i1,a,f8.4,a)')
      parameter(FM2='(1x,a,i1,a,a,i1,a,i1,a,i2,a,f8.4,a)')

      En2 = En
!
      write(6,*) '***********************************'
      write(6,*) ' Energy Levels without LS Coupling '
      write(6,*) '***********************************'
      call line(1)
      write(6,*) '| state | n | L | Energy level | '
      call line(1)
      do 
          Loc(1:2) = MINLOC(En2(:,:,1))
          Loc(2) = Loc(2) - 1
          if (En2(Loc(1), Loc(2), 1) == 0.0d0) exit
          write(6,FM1) '|  ',Loc(1),C(Loc(2)),'   | ',Loc(1),' | ',&
                     & Loc(2),' | ' , &
                     & En2(Loc(1),Loc(2),1),' Mev |'
          En2(Loc(1), Loc(2) ,1) = 0.0d0
      end do
      call line(1)
!
      write(6,*)
      write(6,*) '********************************'
      write(6,*) ' Energy Levels with LS Coupling '
      write(6,*) '********************************'
      call line(2)
      write(6,*) '| state | n | L |   j  | Energy level | '
      call line(2)
      do
          Loc = MINLOC(En2(:,:,2:3))
          Loc(2) = Loc(2) - 1
          Loc(3) = Loc(3) + 1
          select case(Loc(3))
              case(2); J = 2 * Loc(2) + 1
              case(3); J = 2 * Loc(2) - 1
          end select
          if ( En2(Loc(1), Loc(2), Loc(3)) == 0.0d0 ) exit
          write(6,FM2) '|  ',Loc(1),C(Loc(2)),'   | ', Loc(1),' | ',&
                     & Loc(2),  ' | ' ,J, '/2 | ', &
                     & En2(Loc(1),Loc(2),Loc(3)), ' Mev |'
          En2(Loc(1), Loc(2) ,Loc(3)) = 0.0d0
      end do
      call line(2)

      write(6,*)
      write(6,*) ' %> gnuplot < gnu'
      write(6,*) ' %> cd TEX_Files'
      write(6,*) ' %> make_pdf_file'
      return
      end subroutine
!====================================================================!
      subroutine  line(k)
      implicit none
      integer :: k
      
      select case(k)
          case(1)
              write(6,*) '--------------------------------'
          case(2)
              write(6,*) '---------------------------------------'
      end select
      
      return 
      end subroutine

!====================================================================!
      subroutine  Show_Eq
      use Com_var
      implicit none

      write(6,*)
      write(6,*) '|*********************************************************|'
      write(6,*) '| Schroedinger Equation :                                 |'
      write(6,*) '|                                                         |'
      write(6,*) '|  d^2R     2   dR     L(L+1)        2m                   |'
      write(6,*) '| ------ + --- ---- - ------- R - --------(V(r) - E)R = 0 |'
      write(6,*) '|  dx^2     r   dr      r^2        hbar^2                 |'
      write(6,*) '|                                                         |'
      write(6,*) '| where                                                   |'
      write(6,*) '|                        1   df  --> -->                  |'
      write(6,*) '|  V(r) = V0.f(r) + Vso --- ----  L . s                   |'
      write(6,*) '|                        r   df                           |'
      write(6,*) '| with                                                    |'
      write(6,*) '|                  1                                      |'
      write(6,*) '|  f(r) = ---------------------    .                      |'
      write(6,*) '|          1 + exp((r - R)/a)                             |'
      write(6,*) '|                                                         |'
      write(6,*) '| Boundary Conditions :                                   |'
      write(6,*) '|                                                         |'
      write(6,*) '|   R ~ r^L  near the center  r ~ 0.                      |'
      write(6,*) '|   R ~ exp(- kr)  as r tends to infinity.                |'
      write(6,*) '|                                                         |'
      write(6,*) '***********************************************************'
      write(6,*)

      end subroutine

!====================================================================!
      subroutine Title

      write(6,*)
      write(6,*)'******************************************************'
      write(6,*)' Numerical calculation of energy levels of an nucleon '
      write(6,*)' based on the shell model with Woods-Saxon potential. '
      write(6,*)'******************************************************'

      return 
      end subroutine
!====================================================================!
      subroutine  Waves(n, L, E)
!--------------------------------------------------------------------!
!     This subroutine makes the files which we write down the        !
!     numerical values of the wave functions by calling subroutine   !
!     'Solve'. This subroutine also makes the file name from 'n', L' !
!     and T by using intrinsic function 'CHAR'.                      !
!--------------------------------------------------------------------!
      use Com_var
      use Solve_eq
      implicit none
      integer, save :: m = 0
      integer :: i
      integer, intent(in) :: n, L
      real*8,  intent(in) :: E
      real*8 :: y(NVAR)
      character(len=40) :: FNAME
      character(len=14), save :: di(2)
      character, save :: Q(0:19)*2
      data di / 'Waves/No_LS/', 'Waves/With_LS/'/ 
      logical, save :: set = .true.

      if (set) then
          do i=0, 9
              Q(i)    = CHAR(i+48)
              Q(i+10) = CHAR(49) // CHAR(i+48)
          end do
          set = .false.
      end if

      if (T == 1) then
          m = m + 1
          FNAME = TRIM(di(1)) // TRIM(Q(n)) // C(L) // '.dat'
          call Solve(y=y, L=L, E=E, Fnum=m+FFN, Fname=FNAME)
          if (L == 0) then
              m = m + 1
              FNAME = di(2) //TRIM(Q(n)) // C(L) // '1_2.dat'
              call Solve(y=y, L=L, E=E, Fnum=m+FFN, Fname=FNAME)
          end if
      else if(T == 2) then
          m = m + 1
          FNAME = di(2)//TRIM(Q(n))//C(L)//TRIM(Q(2*L+1))//'_2.dat'
          call Solve(y=y, L=L, E=E, Fnum=m+FFN, Fname=FNAME)
      else if(T == 3) then
          m = m + 1
          FNAME = di(2)//TRIM(Q(n))//C(L)//TRIM(Q(2*L-1))//'_2.dat'
          call Solve(y=y, L=L, E=E, Fnum=m+FFN, Fname=FNAME)
      end if

      return
      end  subroutine
!====================================================================!
      subroutine WS_potential(Fnum)
      implicit none
      integer, intent(in) :: Fnum
      real*8 :: x = 0.0d0
      real*8, external :: V
      
      open(Fnum,file='WSpot.dat')
      do while(x < 12.0d0)
          write(Fnum,*) x, V(x)
          x = x + 0.1d0
      end do
      close(Fnum)

      return
      end subroutine
!====================================================================!
      subroutine Gnuplot(P_num, En, FNUM, FNAME)
!--------------------------------------------------------------------!
!     This subroutine makes a file whose name is 'FNAME'.            !
!     This is a file for gnuplot. After the execution of this        !
!     program, redirect the file 'FNAME' :                           !
!           %> gnuplot < FNAME                                       !
!     then eps file will be created.                                 !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(in) :: P_num(0:Lmax), FNUM
      real*8,  intent(in) :: En(nmax, 0:Lmax, 3)
      character(len=*), intent(in) :: FNAME

      open(FNUM, file=FNAME)
      call Gnuplot_Wave_NoLS(En, FNUM)
      call Gnuplot_Wave_LS(En, Fnum)
      call Gnuplot_head(FNUM)
      call Gnuplot_label(P_num, En, FNUM)
      call Gnuplot_arrow(P_num, En, FNUM)
      call Gnuplot_plot(P_num, En, FNUM)
      close(FNUM)

      return
      end
!====================================================================!
      subroutine Gnuplot_head(Fnum)
!--------------------------------------------------------------------!
!     Set options of gnuplot.                                        !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(in) :: Fnum
      integer :: i, j, lstyle(nmax*Lmax)
      integer :: tmp(7) = (/2,3,4,5,7,8,9/)

      j = 1
      do i=1, nmax * Lmax
          lstyle(i) = tmp(j)
          if (j == 7) j = 0
          j = j + 1
      end do

      write(Fnum,*) 'set terminal postscript  eps enhanced  color '
      write(Fnum,*) 'q(x)=1/(1+exp(10**10*x))'
      write(Fnum,*) 'f(x)=q(-x-5)*q(x+3)  '
      write(Fnum,*) 'g(x)=q(-x-8)*q(x+6)  '       
      write(Fnum,*) 'unset xtics          '    
      write(Fnum,*) 'unset key            '   
      write(Fnum,*) 'set ylabel "E(MeV)"  '    
      write(Fnum,*) 'set xrange [-10.8:2]  '   
      write(Fnum,*) 'set yrange [-45:5]  '   
      write(Fnum,*) 'set nolabel           ' 
      do i=1, nmax * Lmax
          write(Fnum,*) 'set style line',i,'lt',lstyle(i)
      end do
      write(Fnum,*) 'set output "EPS_Files/Energy_Levels.eps"  '

      return
      end subroutine
!====================================================================!
      subroutine Gnuplot_label(P_num, En, Fnum)
!--------------------------------------------------------------------!
!     Set lebel options of gnuplot.                                  !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(in) :: Fnum, P_num(0:Lmax)
      integer ::  i=0, L
      real*8,  intent(in) :: En(nmax,0:Lmax, 3)
      character :: F1*20
      parameter(F1='(a,i1,a,a,f10.5)')

      do L=0, Lmax
          do i=1, P_num(L)
              if (P_num(L) > nmax) stop 'Error3'
              write(Fnum,F1) 'set label "', i, C(L),'" at first &
                          & -8.6, first ',En(i, L, 1)
          end do
      end do

      return
      end subroutine
!====================================================================!
      subroutine  Gnuplot_arrow(P_num, En, Fnum)
!--------------------------------------------------------------------!
!     Set arrow options of gnuplot.                                  !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: i, L
      integer, intent(in) :: Fnum, P_num(0:Lmax)
      real*8,  intent(in) :: En(nmax, 0:Lmax, 3)

      write(Fnum,*) 'unset arrow'
      do L=0, Lmax
          do i=1, P_num(L)
               write(Fnum,*) 'set arrow from -6,', En(i, L, 1), &
                          & 'to -5,', En(i, L, 2)
               if(L == 0) cycle 
               write(Fnum,*) 'set arrow from -6,', En(i, L, 1), &
                          & 'to -5,', En(i, L, 3)
          end do
      end do

      return
      end subroutine
!====================================================================!
      subroutine  Gnuplot_plot(P_num, En, Fnum)
!--------------------------------------------------------------------!
!     Plot options.                                                  !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(in) :: Fnum, P_num(0:Lmax)
      real*8,  intent(in) :: En(nmax, 0:Lmax, 3)
      integer :: i, n=0, L

      write(Fnum,*) 'plot \'
      
      do L=0, Lmax
          do i=1, nmax
              if (En(i, L, 2) == 0.0d0) then
                  exit
              else if (En(i, L, 2) /= 0.0d0) then
                  write(Fnum,*) En(i, L, 2),'*f(x) ls ',n + i,',\'
              end if
              if (En(i, L, 1) /= 0.0d0) then
                  write(Fnum,*) En(i, L, 1),'*g(x) ls ',n + i,',\'
              end if
              if(L /= 0 .and. En(i, L, 3) /= 0.0d0) then
                  write(Fnum,*) En(i, L, 3),'*f(x) ls ',n + i,',\'
              end if
          end do
          n = n + P_num(L)
      end do

      write(Fnum,*) '0 ls 1'
      write(Fnum,*) 'set term x11'

      return
      end subroutine
!====================================================================!
      subroutine  Gnuplot_Wave_NoLS(En, Fnum)
      use Com_var
      implicit none
      integer :: n, L
      integer, intent(in) :: Fnum
      real*8,  intent(in) :: En(nmax, 0:Lmax, 3)
      character(len=2), save :: Q(0:19)
      character :: FNAME*20
      character(len=30) :: FM1
      character(len=30), save :: di(2)
      logical, save :: set = .true.
      data di /'Waves/No_LS/', 'EPS_Files/'/
      parameter(FM1='(a,a,i1,a,a)')

      write(Fnum,*) 'set terminal postscript eps enhanced  color '
      write(Fnum,*) 'unset key'
      write(Fnum,*) "set xlabel 'r [fm]'"
      write(Fnum,*) "set ylabel 'Radial Wave Functions'"
      write(Fnum,*) 'set size 0.5,0.5'
      do L=0, Lmax
          if (En(1, L, 1) == 0.0d0) then
               exit
          else
               write(Fnum,*) "set title '",C(L) ," states'"
               write(Fnum,*) "set output '",TRIM(di(2)),C(L),&
                             & "_states.eps'"
               write(Fnum,*) 'plot \'
          end if
          do n=1, nmax
              if (En(n, L, 1) == 0.0d0) then
                  write(Fnum,*) '0'
                  exit
              end if
              write(Fnum,FM1) "'",TRIM(di(1)),n,C(L), &
                              & ".dat'with line ,\"
          end do
      end do
      write(Fnum,*) 'unset xlabel'
      write(Fnum,*) 'unset ylabel'
      write(Fnum,*) 'unset title'
      write(Fnum,*) 'set size 1,1'
      write(Fnum,*) 'set term x11'

      return
      end subroutine
!====================================================================!
      subroutine  Gnuplot_Wave_LS(En, Fnum)
      use Com_var
      implicit none
      integer :: n, L, i
      integer, intent(in) :: Fnum
      real*8,  intent(in) :: En(nmax, 0:Lmax, 3)
      character(len=2), save :: Q(0:19)
      character :: FNAME*20
      character(len=30) :: FM1
      character(len=30), save :: di(2)
      logical, save :: set = .true.
      data di /'Waves/With_LS/', 'EPS_Files/'/
      parameter(FM1='(a,i1,a,a)')

      if (set) then
          do i=0, 9
              Q(i)    = CHAR(i+48)
              Q(i+10) = CHAR(49) // CHAR(i+48)
          end do
          set = .false.
      end if

      write(Fnum,*) 'set terminal postscript eps enhanced  color '
      write(Fnum,*) 'unset key'
      write(Fnum,*) "set xlabel 'r [fm]'"
      write(Fnum,*) "set ylabel 'Radial Wave Functions'"
      write(Fnum,*) 'set size 0.5,0.5'

      do L=0, Lmax
          if (En(1, L, 2) == 0.0d0 .and. En(1, L, 3) == 0.0d0) then
               exit
          else
               write(Fnum,*) "set title '",C(L) ," states'"
               write(Fnum,*) "set output '",TRIM(di(2)),C(L),&
                             & "_states_LS.eps'"
               write(Fnum,*) 'plot \'
          end if
          do n=1, nmax
              if (En(n, L, 2) == 0.0d0  .and. &
                 &En(n, L, 3) == 0.0d0) then
                  write(Fnum,*) '0'
                  exit
              else
                  if (En(n, L, 2) /= 0.0d0) then
                  write(Fnum,FM1) "'" // TRIM(di(1)),n,C(L), &
                       & TRIM(Q(2*L+1))// "_2.dat'with line ,\"
                  end if
                  if (En(n, L, 3) /= 0.0d0) then
                  write(Fnum,FM1) "'" // TRIM(di(1)),n,C(L), &
                       & TRIM(Q(2*L-1))// "_2.dat'with line ,\"
                  end if
              end if
          end do
      end do
      write(Fnum,*) 'unset xlabel'
      write(Fnum,*) 'unset ylabel'
      write(Fnum,*) 'unset title'
      write(Fnum,*) 'set size 1,1'
      write(Fnum,*) 'set term x11'

      return
      end subroutine
