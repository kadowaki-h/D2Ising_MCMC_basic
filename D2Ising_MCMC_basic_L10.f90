! D2Ising_MCMC_basic_L10.f90  2023/7/1
!
! gfortran -O -o D2Ising_MCMC_basic_L10.ex mt19937_64_OMP.f90 D2Ising_MCMC_basic_L10.f90
! ifort -O -o D2Ising_MCMC_basic_L10.ex mt19937_64_OMP.f90 D2Ising_MCMC_basic_L10.f90
! ./D2Ising_MCMC_basic_L10.ex
!
!**********************************************************************
! Copyright (c) 2023, Hiroaki Kadowaki [email: Kadowaki.Hiroaki.G at gmail.com (remove spaces)]
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, 
! this list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice, 
! this list of conditions and the following disclaimer in the documentation and/or 
! other materials provided with the distribution.
! 
! 3. Neither the name of the copyright holder nor the names of its contributors may 
! be used to endorse or promote products derived from this software without specific 
! prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
! NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
! POSSIBILITY OF SUCH DAMAGE.
! 
!    (See https://opensource.org/licenses/BSD-3-Clause )
!
!**********************************************************************
!
! A fortran90 program of a MCMC simulation 
! which calculates <E(T)>, <C(T)>, <|M|(T)> etc of a 2d Ising model.
!-----
!  Markov chain Monte Carlo method
!  Metropolis algorithm
!  single spin flip
!  Ising model
!  d=2 LxL square lattice
!  E = - Jnn * sum_{x = 1,L; y = 1,L} [ spin(x,y)*spin(x+1,y) + spin(x,y)*spin(x,y+1) ]
!     spin(x,y) = 1 or -1
!     spin(x, L+1) = spin(x, 1)  ;  spin(L+1, y) = spin(1, y)
!  M = sum_{x=1,L; y=1,L} spin(x,y)
!  temperature (k_B * T) = temperature_list( 1, 2, ... , N_temperature_list )
!-----
PROGRAM main
!
  use,intrinsic :: iso_fortran_env
  use mt19937_64 , only : init_genrand64, genrand64_real3
  IMPLICIT NONE
  integer, allocatable :: spin( :,: )
  integer :: L, N_spin, N_mcs_average, N_mcs_discard, nf
  integer :: iter, j_T, i_spin, x, y, N_temperature_list, N_observable
  real(REAL64) :: Jnn, dE, temperature, E, M, accept
  real(REAL64), allocatable :: temperature_list( : ), sum_observables( : )
  character(128) :: file_table
  integer(int64) :: seed_64bit
!
! initial settings
    seed_64bit = 21478288371993_int64 + 202306251056_int64
!!!!!    seed_64bit = seed_64bit + 202306291046_int64       !! any 64-bit integer of the seed
    call init_genrand64( seed_64bit )         !! initialize pseudo-random generator mt19937_64
    L = 10                                    !! LxL square lattice, L = 10  16  32  64  128 ...
    Jnn = 1.0d0                               !! exchange constant, Jnn > 0 ferromagnetic coupling
    file_table = "D2Ising_MCMC_basic_L10.tab" !! output file, "D2Ising_MCMC_basic_L16.tab" .....
    N_temperature_list = 31  ;  allocate( temperature_list( N_temperature_list ) )  !! number of temperatures
    temperature_list( 1:N_temperature_list ) = &
      & (/ 4.0d0 &
        &  , 3.9d0, 3.8d0, 3.7d0, 3.6d0, 3.5d0,    3.4d0, 3.3d0, 3.2d0, 3.1d0, 3.0d0 &
        &  , 2.9d0, 2.8d0, 2.7d0, 2.6d0, 2.5d0,    2.4d0, 2.3d0, 2.2d0, 2.1d0, 2.0d0 &
        &  , 1.9d0, 1.8d0, 1.7d0, 1.6d0, 1.5d0,    1.4d0, 1.3d0, 1.2d0, 1.1d0, 1.0d0    /)
    N_mcs_discard = 500000 !! 5000000         !! MC steps (per spin) to be discarded (for thermalization)
    N_mcs_average = 500000 !! 5000000         !! MC steps (per spin) for averaging
    allocate( spin( L,L ) )
    N_spin = L**2                             !! number of spins
    N_observable = 16  ;  allocate( sum_observables( N_observable ) )  !! number of observables
!
! open file_table and write several lines
    open( newunit = nf , file= TRIM(file_table) )
    write(nf,'(A,A)')     '  file_table         = ', TRIM(file_table)
    write(nf,'(A,5G15.7)')'  Jnn                = ', Jnn
    write(nf,'(A,5I10)')  '  L, N_spin          = ', L, N_spin
    write(nf,'(A,5I10)')  '  N_temperature_list = ', N_temperature_list
    write(nf,'(A)')       '  temperature_list   = '
    write(nf,'(10X,5G15.7)') temperature_list(1:N_temperature_list)
    write(nf,'(A,5I20)')  '  seed_64bit         = ', seed_64bit
    write(nf,'(A,5I10)')  '  N_mcs_discard, N_mcs_average = ', N_mcs_discard, N_mcs_average
    write(nf,'(A,A)') "  temperature  <E/N_spin>  <C/N_spin>  <chi/N_spin>  <abs_chi/N_spin>" &
    & , "  <M/N_spin>  <|M|/N_spin>  <(M/N_spin)**2>  Binder_<M**4>/<M**2>**2  acceptance"
!
  call init_spin
!
  do j_T = 1, N_temperature_list
    temperature = temperature_list( j_T )
!!!!!    call init_spin
!
    DO iter = 1, N_mcs_discard            !! MC steps to be discarded (for thermalization)
      DO i_spin = 1, N_spin               !! one Monte Carlo step per spin
        x = int( L*genrand64_real3()+1 )  !! (x,y) of a trial spin
        y = int( L*genrand64_real3()+1 )
        call calc_dE(x, y, dE)            !! dE = E_{after trial} - E_{before trial}
        IF( genrand64_real3() <= exp( -dE/temperature ) ) THEN  !! accept  ( Metropolis )
          spin(x,y) = - spin(x,y)
     !! ELSE                                                    !! reject
        END IF
      END DO
    END DO
!
    sum_observables( : ) = 0.0d0
    DO iter = 1, N_mcs_average            !! MC steps for averaging
      accept = 0.0d0
      DO i_spin = 1, N_spin               !! one Monte Carlo step per spin
        x = int( L*genrand64_real3()+1 )  !! (x,y) of a trial spin
        y = int( L*genrand64_real3()+1 )
        call calc_dE(x, y, dE)            !! dE = E_{after trial} - E_{before trial}
        IF( genrand64_real3() <= exp( - dE/temperature ) ) THEN  !! accept  ( Metropolis )
          spin(x,y) = - spin(x,y)
          accept = accept + 1.0d0
     !! ELSE                                                     !! reject
        END IF
      END DO
      call calc_E(E)
      M = DBLE( SUM( spin(1:L, 1:L) ) )  !! call calc_M(M)
      CALL sum_observable_data( E, M, accept )
    END DO
    CALL output_table( nf )
  END DO
!
  close(nf)
!
contains
!-----
  subroutine init_spin
    use,intrinsic :: iso_fortran_env
    use mt19937_64 , only : genrand64_real3
    IMPLICIT NONE
    integer :: x, y
    spin(:,:) = -1
    DO y = 1, L
      DO x = 1, L
        IF( genrand64_real3() < 0.5D0 ) spin(x, y) = 1  !!  random spin configuration
      END DO
    END DO
  end subroutine init_spin
!
!-----
  subroutine calc_E(E)
    use,intrinsic :: iso_fortran_env
    IMPLICIT NONE
    real(REAL64) :: E
    integer :: x, y, xp1, yp1
    E = 0.0d0
    DO y = 1, L
      yp1 = y + 1  ;  IF(y == L) yp1 = 1
      DO x = 1, L
        xp1 = x + 1  ;  IF(x == L) xp1 = 1
        E = E + DBLE( spin(x, y) * ( spin(x, yp1) + spin(xp1, y) ) )
      END DO
    END DO
    E = - Jnn * E
  end subroutine calc_E
!
!-----
  subroutine calc_dE(x, y, dE)
    use,intrinsic :: iso_fortran_env
    IMPLICIT NONE
    integer :: x, y
    real(REAL64) :: dE
    integer :: xp1, xm1, yp1, ym1
    !! calc dE = E_{after trial} - E_{before trial}
      xm1 = x - 1  ;  IF(x == 1) xm1 = L
      xp1 = x + 1  ;  IF(x == L) xp1 = 1
      ym1 = y - 1  ;  IF(y == 1) ym1 = L
      yp1 = y + 1  ;  IF(y == L) yp1 = 1
    dE = Jnn * DBLE( 2 * spin(x,y) * ( spin(xm1,y) + spin(xp1,y) + spin(x,ym1) + spin(x,yp1) ) )
  end subroutine calc_dE
!
!-----
SUBROUTINE sum_observable_data( E, M, accept )
  use,intrinsic :: iso_fortran_env
!  summation of observables
  IMPLICIT NONE
  real(REAL64), intent(in) :: E, M, accept
  real(REAL64) :: magnetization
  sum_observables(1) = sum_observables(1) + accept / DBLE(N_spin)  !! sum accept/N_spin
  sum_observables(2) = sum_observables(2) + E                      !! sum E
  sum_observables(3) = sum_observables(3) + E*E                    !! sum E**2
  magnetization = M / DBLE(N_spin)
  sum_observables(4) = sum_observables(4) + magnetization          !! sum M/N_spin
  sum_observables(5) = sum_observables(5) + abs(magnetization)     !! sum |M|/N_spin
  sum_observables(6) = sum_observables(6) + magnetization**2       !! sum (M/N_spin)**2
  sum_observables(7) = sum_observables(7) + magnetization**4       !! sum (M/N_spin)**4
END SUBROUTINE sum_observable_data
!
!-----
SUBROUTINE output_table( nf )
  use,intrinsic :: iso_fortran_env
  IMPLICIT NONE
  integer, intent(in) :: nf
  real(REAL64) :: acceptance, energy_av, energy_sq_av, specific_heat
  real(REAL64) :: m_av, m_abs_av, m_sq_av, m_sqsq_av, chi, abs_chi, binder
!
  acceptance    = sum_observables(1) / DBLE(N_mcs_average)
  energy_av     = sum_observables(2) / DBLE(N_mcs_average)  !! <E>
  energy_sq_av  = sum_observables(3) / DBLE(N_mcs_average)  !! <E**2>
  specific_heat = (energy_sq_av - energy_av * energy_av) / ((temperature**2) * DBLE(N_spin) )
!                                                           !! <C/N_spin> = (1/T)**2 * (<E**2>-<E>**2)/N_spin
  energy_av     = energy_av / DBLE(N_spin)                  !! <E/N_spin>
  m_av          = sum_observables(4) / DBLE(N_mcs_average)  !! <M/N_spin>
  m_abs_av      = sum_observables(5) / DBLE(N_mcs_average)  !! <|M|/N_spin>
  m_sq_av       = sum_observables(6) / DBLE(N_mcs_average)  !! <(M/N_spin)**2>
  m_sqsq_av     = sum_observables(7) / DBLE(N_mcs_average)  !! <(M/N_spin)**4>
  chi     = DBLE(N_spin) * ( m_sq_av - (m_av**2)     ) / temperature !! (1/T) * ( <M**2> - <M>**2 )/ Nspin
  abs_chi = DBLE(N_spin) * ( m_sq_av - (m_abs_av**2) ) / temperature !! (1/T) * ( <M**2> - <|M|>**2 )/ Nspin
  binder  = m_sqsq_av / ( m_sq_av**2 )                               !! <M**4>/<M**2>**2
  write(nf,'(20G15.7)') &
  & temperature, energy_av, specific_heat, chi, abs_chi, m_av, m_abs_av, m_sq_av, binder, acceptance
END SUBROUTINE output_table
!
END PROGRAM main
!
