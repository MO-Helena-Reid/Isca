!-----------------------------------------------------------------------------
! Crown Copyright (c) Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! This file was written by Helena Reid
! as part of the initial release of the 
! PAPILLON stochastic physics scheme,
! under the BSD 3-clause license,
! which should have been provided with this file.
!
! A description of the scheme may be found in preprint here:
! https://doi.org/10.5194/egusphere-2025-6312
!

! Created by Helena Reid on 17/11/2023.
module ml_constants_mod
    use, intrinsic :: iso_fortran_env, only : int8, int16, int32, int64, &
                                            real32, real64, real128
    implicit none
    private

    public :: r_ml, r_um, ml_nlev, ml_output_len, ml_input_len, &
              reference_pressure, poisson_constant, z_papillon

    integer, parameter :: r_ml = real64
    integer, parameter :: r_um = real64
    integer, parameter :: ml_nlev = 70
    integer, parameter :: ml_output_len = 70
    integer, parameter :: ml_input_len = 213
    real(kind=r_um), parameter :: reference_pressure = 100000
    real(kind=r_um), parameter :: poisson_constant = 0.2854
    real,parameter,dimension(ml_nlev) :: z_papillon = [&
    20,&
    53,&
    100,&
    160,&
    233,&
    320,&
    420,&
    533,&
    660,&
    800,&
    953,&
    1120,&
    1300,&
    1493,&
    1700,&
    1920,&
    2153,&
    2400,&
    2660,&
    2933,&
    3220,&
    3520,&
    3833,&
    4160,&
    4500,&
    4853,&
    5220,&
    5600,&
    5993,&
    6400,&
    6820,&
    7253,&
    7700,&
    8160,&
    8634,&
    9121,&
    9622,&
    10137,&
    10667,&
    11213,&
    11775,&
    12355,&
    12954,&
    13575,&
    14221,&
    14895,&
    15602,&
    16348,&
    17137,&
    17980,&
    18884,&
    19861,&
    20923,&
    22087,&
    23369,&
    24789,&
    26371,&
    28141,&
    30130,&
    32371,&
    34904,&
    37771,&
    41022,&
    44712,&
    48902,&
    53659,&
    59060,&
    65187,&
    72133,&
    80000]
end module ml_constants_mod
