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
module papillon_config_mod
use ml_constants_mod, only: ml_nlev

#ifdef INTERNAL_FILE_NML
  use mpp_mod,        only: input_nml_file
#else
  use fms_mod,        only: open_namelist_file, close_file
#endif
use fms_mod,          only: stdlog, NOTE, error_mesg, &
                            uppercase, check_nml_error, file_exist
use time_manager_mod, only: time_type
 
implicit none
private
public :: sampling_method, sampling_method_correlated, &
    sampling_method_uncorrelated, &
    sampling_method_development, sampling_method_constant, &
    constant_t_sd_profile, &
    planet_radius, time_scale_factor, radius_scale_factor,&
    lat_scale_factor,lon_scale_factor,height_scale_factor,&
    noise_scale_factor,noise_centre_t,noise_centre_x,&
    noise_centre_y,noise_centre_z,papillon_init
integer, parameter :: sampling_method_correlated = 4
integer, parameter :: sampling_method_uncorrelated = 3
integer, parameter :: sampling_method_development = 2
integer, parameter :: sampling_method_constant = 1
integer            :: sampling_method = sampling_method_constant
real,parameter     :: planet_radius = 6371229.0
real               :: time_scale_factor = 4.62962963e-5
real               :: radius_scale_factor = 1.56955588945e-7
real               :: lat_scale_factor = 1.0
real               :: lon_scale_factor = 1.0
real               :: height_scale_factor = 1e-4
real               :: noise_scale_factor = 0.2
real               :: noise_centre_t=0.0,&
        noise_centre_x=0.0,&
        noise_centre_y=0.0,&
        noise_centre_z=0.0
real, parameter, dimension(ml_nlev) :: constant_t_sd_profile = [&
                        7.474816271768951470e-01,&
                        7.220088264402453326e-01,&
                        6.913767786022448902e-01,&
                        6.604720008469491255e-01,&
                        6.274260727768707913e-01,&
                        5.935843298504931420e-01,&
                        5.581323506608593110e-01,&
                        5.269029965448590591e-01,&
                        5.031533210446547111e-01,&
                        4.843463299555669233e-01,&
                        4.654120837871334460e-01,&
                        4.457017772443185866e-01,&
                        4.293024659581302083e-01,&
                        4.056554760987349240e-01,&
                        3.790849683043594487e-01,&
                        3.557921600403540219e-01,&
                        3.356802608285804901e-01,&
                        3.172699391939735047e-01,&
                        2.994834551311509885e-01,&
                        2.848621886974697093e-01,&
                        2.709924779630836222e-01,&
                        2.572476792558234804e-01,&
                        2.440186341901312839e-01,&
                        2.313308519787665718e-01,&
                        2.192662743952543958e-01,&
                        2.058378770425789150e-01,&
                        1.924132092201106536e-01,&
                        1.787514127334494940e-01,&
                        1.639840535199276195e-01,&
                        1.492921843433178131e-01,&
                        1.342466989060429050e-01,&
                        1.191909874240073275e-01,&
                        1.059360932657187615e-01,&
                        9.405187588018477929e-02,&
                        8.625684619950001186e-02,&
                        8.144014220540768401e-02,&
                        7.810128342839928184e-02,&
                        7.707802331457197509e-02,&
                        7.592792672160778022e-02,&
                        7.384283649908479630e-02,&
                        7.322961327378130214e-02,&
                        7.284487029770946032e-02,&
                        7.199318196461582109e-02,&
                        7.005515989035870916e-02,&
                        7.144544068708627571e-02,&
                        7.864479822799268216e-02,&
                        9.508975457712660895e-02,&
                        1.198609616134928141e-01,&
                        1.497640657192059743e-01,&
                        1.673066648051236560e-01,&
                        1.618827229052252348e-01,&
                        1.252618299910949196e-01,&
                        9.223363294064997053e-02,&
                        8.075789250709988765e-02,&
                        7.602602076690394284e-02,&
                        7.804772237575226257e-02,&
                        8.323806104333643374e-02,&
                        9.705793310236544846e-02,&
                        1.093643411308970886e-01,&
                        1.291602296345874534e-01,&
                        1.561289107357137773e-01,&
                        1.889193525029131404e-01,&
                        2.743245694254396461e-01,&
                        2.903155901354112500e-01,&
                        2.994331656197566915e-01,&
                        3.943417601417406049e-01,&
                        4.572177456707915422e-01,&
                        3.187616441406314616e-01,&
                        1.366351249509512866e-01,&
                        5.126769960918496627e-02]
character(len=14), parameter :: mod_name_ppl = "papillon_conf"
namelist /papillon_nml/ &
  sampling_method, time_scale_factor,radius_scale_factor,&
  height_scale_factor,lat_scale_factor,lon_scale_factor,&
  noise_scale_factor,noise_centre_x,noise_centre_y,&
  noise_centre_z,noise_centre_t
contains
subroutine papillon_init(axes, Time)
  type(time_type), intent(in)       :: Time
  integer, intent(in), dimension(4) :: axes
  integer :: io, ierr, nml_unit, stdlog_unit

#ifdef INTERNAL_FILE_NML
  read(input_nml_file, nml=papillon_nml, iostat=io)
  ierr = check_nml_error(io, 'papillon_nml')
#else
  if (file_exist('input.nml')) then
    nml_unit = open_namelist_file()
    ierr = 1
    do while (ierr /= 0)
      read(nml_unit, nml=papillon_nml, iostat=io, end=10)
      ierr = check_nml_error(io, 'papillon_nml')
    enddo
  10    call close_file(nml_unit)
  endif
#endif
  stdlog_unit = stdlog()
  write(stdlog_unit, papillon_nml)

  call error_mesg(mod_name_ppl, 'Using PAPILLON stochastic physics scheme', NOTE)
end subroutine papillon_init
end module papillon_config_mod
