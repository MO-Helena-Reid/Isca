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
module papillon_alg_mod
use time_manager_mod, only: time_type, get_time, print_time
use simplex_noise_mod, only: snoise4d
use spline_interp_mod, only: spline_type, spline_set_coeffs, spline_evaluate, spline_evaluate_1d
use ml_constants_mod, only: ml_nlev, ml_input_len, z_papillon
use papillon_config_mod, only: planet_radius, time_scale_factor, radius_scale_factor, &
  lat_scale_factor, lon_scale_factor, height_scale_factor, noise_scale_factor, noise_centre_t, noise_centre_x,&
  noise_centre_y, noise_centre_z, sampling_method_constant, sampling_method_development, sampling_method,&
  constant_t_sd_profile,output_papillon_diags,id_papillon_noise,id_papillon_t_sd,id_papillon_t_pert
use thetasds00v001_mod, only: thetasds00v001
use normalisation_mod, only: normalise_inputs
use sample_multivariate_normal_mod, only: sample_from_snoise
implicit none
contains
subroutine papillon_alg(&
  tpert,&
  t,&
  pfull,&
  q,&
  z_full,&
  lat,&
  lon,&
  fracland,&
  orog,&
  sd_orog,&
  Time)
  use reverse_array_mod, only: reverse
  implicit none
  real, intent(out), dimension(:,:,:) :: tpert
  real, intent(in),  dimension(:,:,:) :: t,pfull,q,z_full
  real, intent(in),  dimension(:,:)   :: lat,lon,fracland,orog,sd_orog
  type(time_type), intent(in)         :: Time

  integer                             :: seconds
  real, dimension(4)                  :: noise_loc
  real, parameter                     :: half_pi=1.570796327
  real,dimension(size(t,1),size(t,2)) :: sin_latitudes
  real,dimension(size(t,1),size(t,2)) :: cos_latitudes
  real,dimension(size(t,1),size(t,2)) :: sin_longitudes
  real,dimension(size(t,1),size(t,2)) :: cos_longitudes
  real,dimension(size(t,1),size(t,2),size(t,3)) :: noise, t_sd
  real                                :: noise_radius
  integer                             :: i,j,k
  type(spline_type)                   :: spline
  real,dimension(ml_nlev)             :: t_papillon, q_papillon, &
                                         p_papillon, tpert_papillon,&
                                         noise_papillon, t_sd_papillon
  real,dimension(ml_input_len)        :: inputs, normalised_inputs
  real,dimension(size(t,3))           :: z_full_r,pfull_r,q_r,t_r,tpert_r
  
  !print*,GETPID(),"BEGIN SUBROUTINE compute_t_pert"
  !tpert(:,:,:) = 0.0
  !noise(:,:,:) = 0.0
  !noise_papillon(:) = 0.0
  call get_time(Time, seconds)
  noise_loc = [0.0,0.0,0.0,time_scale_factor*seconds-noise_centre_t]
  sin_latitudes = sin(lat+half_pi)
  cos_latitudes = cos(lat+half_pi)
  sin_longitudes = sin(lon)
  cos_longitudes = cos(lon)
  do i=1, size(t,1)
    do j=1, size(t,2)
      call reverse(z_full(i,j,:),z_full_r)
      ! generate noise on papillon zlevs
      do k=1, ml_nlev
        noise_radius = (radius_scale_factor*planet_radius + height_scale_factor*z_papillon(k))
        noise_loc(1) = noise_radius*sin_latitudes(i,j)*cos_longitudes(i,j) - noise_centre_x
        noise_loc(2) = noise_radius*sin_latitudes(i,j)*sin_longitudes(i,j) - noise_centre_y
        noise_loc(3) = noise_radius*cos_latitudes(i,j) - noise_centre_z
        ! noise_papillon(k) = noise_scale_factor*snoise4d(noise_loc)
        noise_papillon(k) = snoise4d(noise_loc)
      end do

      ! only interpolate noise back to ISCA grid if it is requested as diagnostic
      if (id_papillon_noise > 0) then
        call spline_set_coeffs(z_papillon,noise_papillon,ml_nlev,spline)
        noise(i,j,:) = spline_evaluate_1d(z_full_r,spline)
        call reverse(noise(i,j,:),noise(i,j,:))
      end if
      
      if (any(isnan(noise_papillon))) print*, GETPID(), "WARN: NaNs detected in noise_papillon"

      select case (sampling_method)
        case (sampling_method_constant)
          call sample_from_snoise(noise_scale_factor*noise_papillon, constant_t_sd_profile, tpert_papillon)
        case (sampling_method_development)
          call reverse(z_full(i,j,:),z_full_r)
          call reverse(pfull(i,j,:),pfull_r)
          call reverse(q(i,j,:),q_r)
          call reverse(t(i,j,:),t_r)
          call spline_set_coeffs(z_full_r,pfull_r,size(t,3),spline)
          p_papillon = spline_evaluate(z_papillon, spline)
          call spline_set_coeffs(z_full_r,q_r,size(t,3),spline)
          q_papillon = spline_evaluate(z_papillon, spline)
          call spline_set_coeffs(z_full_r,t_r,size(t,3),spline)
          t_papillon = spline_evaluate(z_papillon, spline)
          inputs(1:70)=q_papillon
          inputs(71:140)=p_papillon
          inputs(141:210)=t_papillon
          inputs(211) = fracland(i,j)
          inputs(212) = orog(i,j)
          ! TODO: actually pipe sd orog into this subroutine, setting to zero for now
          inputs(213) = 0.0
      !    print*,GETPID(),"set inputs"
          call normalise_inputs(inputs, normalised_inputs)
          call thetasds00v001(normalised_inputs, t_sd_papillon)
          ! TODO: does this ml model output sd of theta or t? check
          call sample_from_snoise(noise_scale_factor*noise_papillon, t_sd_papillon, tpert_papillon)

          if (id_papillon_t_sd > 0) then
            call spline_set_coeffs(z_papillon, t_sd_papillon, ml_nlev, spline)
            t_sd(i,j,:) = spline_evaluate(z_full_r,spline)
            call reverse(t_sd(i,j,:),t_sd(i,j,:))
          end if
      end select
      !print*,GETPID(),"shape z zr zp t tp tpp",size(z_full(i,j,:)),size(z_full_r),size(z_papillon),size(t_r),size(t_papillon),size(tpert_papillon),"z_full_r",z_full_r,"t_papillon:",t_papillon,"temperature:",t_r,"t+tpert_papillon:",tpert_papillon
    
      !  print*, GETPID(), "get_alternative_temperature done"
      !tpert_papillon = tpert_papillon - t_papillon
      if (any(isnan(tpert_papillon))) print*, GETPID(), "WARN: NaNs detected in tpert_papillon"
      
      ! Interpolate perturbation back to ISCA grid
      call spline_set_coeffs(z_papillon,tpert_papillon,ml_nlev,spline)
      tpert_r = spline_evaluate_1d(z_full_r,spline)
      if (any(isnan(tpert_r))) then
        print*, GETPID(), "WARN: NaNs detected in tpert_r"
        ! print *, GETPID(), "spline n x y bcd lookup_index lookup_inv_dx", spline%n, spline%x, spline%y, spline%bcd, spline%lookup_index, spline%lookup_inv_dx
      endif
      !print*, GETPID(), "t pert spline evaluated"
      call reverse(tpert_r,tpert(i,j,:))
      !do k=1, ml_nlev
      !  print*,GETPID(),k,z_papillon(k),"tp",tpert_papillon(k),"n",noise_papillon(k),"t",t_papillon(k)
      !end do
      !do k=1, size(t,3)
      !  print*,GETPID(),k,z_full(i,j,k),"tp",tpert(i,j,k),"n",noise(i,j,k),"t",t(i,j,k)
      !end do
      !error stop "print papillon diags"
      !print*,GETPID(),"tpert_papillon",tpert_papillon,"tpert_r:",tpert_r,"tpert:",tpert(i,j,:)
      !print*,GETPID(),"noise:",noise(i,j,:)
      !error stop "print papillon diags"
      !print *, GETPID(), "reversed"
    end do
  end do
  call output_papillon_diags(noise, tpert, t_sd, Time)
  ! call print_time(Time)
  ! print*, GETPID(), "tpert min max",MINVAL(tpert),MAXVAL(tpert)
  ! print*, GETPID(), "temp min max",MINVAL(t+tpert),MAXVAL(t+tpert)
end subroutine papillon_alg
end module papillon_alg_mod
