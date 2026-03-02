module papillon_alg_mod
use    time_manager_mod, only: time_type, get_time
use simplex_noise_mod, only: snoise4d
use spline_interp_mod, only: spline_type, spline_set_coeffs, spline_evaluate
use ml_constants_mod, only: ml_nlev, ml_input_len, z_papillon
use plausible_alternative_mod, only: get_alternative_temperature
implicit none
contains
subroutine compute_t_pert(&
    noise,&
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
    implicit none
    real, intent(out), dimension(:,:,:) :: noise,tpert
    real, intent(in),  dimension(:,:,:) :: t,pfull,q,z_full
    real, intent(in),  dimension(:,:)   :: lat,lon,fracland,orog,sd_orog
    type(time_type), intent(in)         :: Time

    integer                         :: seconds
    real, dimension(4) :: noise_loc
    real,parameter     :: planet_radius = 6371229.0
    real,parameter     :: time_scale_factor = 4.62962963e-5
    real,parameter     :: radius_scale_factor = 1.56955588945e-7
    real,parameter     :: lat_scale_factor = 1.0
    real,parameter     :: lon_scale_factor = 1.0
    real,parameter     :: height_scale_factor = 1e-4
    real,parameter     :: noise_scale_factor = 2.0
    real,parameter     :: noise_centre_t=0.0,&
            noise_centre_x=0.0,&
            noise_centre_y=0.0,&
            noise_centre_z=0.0,&
            half_pi=1.570796327
    real,dimension(size(t,1),size(t,2))::sin_latitudes
    real,dimension(size(t,1),size(t,2))::cos_latitudes
    real,dimension(size(t,1),size(t,2))::sin_longitudes
    real,dimension(size(t,1),size(t,2))::cos_longitudes
    real :: noise_radius
    integer :: i,j,k
    type(spline_type) :: spline
    real,dimension(ml_nlev) :: t_papillon, q_papillon, p_papillon, tpert_papillon
    real,dimension(ml_input_len) :: inputs

    
    tpert(:,:,:) = 0.0
    noise(:,:,:) = 0.0
    call get_time(Time, seconds)
    noise_loc = [0.0,0.0,0.0,time_scale_factor*seconds-noise_centre_t]
    sin_latitudes = sin(lat+half_pi)
    cos_latitudes = cos(lat+half_pi)
    sin_longitudes = sin(lon)
    cos_longitudes = cos(lon)
    do i=1, size(t,1)
        do j=1, size(t,2)
            call spline_set_coeffs(z_full(i,j,:),pfull(i,j,:),size(t,3),spline)
            p_papillon = spline_evaluate(z_papillon, spline)
            call spline_set_coeffs(z_full(i,j,:),q(i,j,:),size(t,3),spline)
            q_papillon = spline_evaluate(z_papillon, spline)
            call spline_set_coeffs(z_full(i,j,:),t(i,j,:),size(t,3),spline)
            t_papillon = spline_evaluate(z_papillon, spline)
            inputs(1:70)=q_papillon
            inputs(71:140)=p_papillon
            inputs(141:210)=t_papillon
            inputs(211) = fracland(i,j)
            inputs(212) = orog(i,j)
            ! TODO: actually pipe sd orog into this subroutine, setting to zero for now
            inputs(213) = 0.0
            do k=1, size(t,3)
              noise_radius = (radius_scale_factor*planet_radius + height_scale_factor*z_full(i,j,k))
              noise_loc(1) = noise_radius*sin_latitudes(i,j)*cos_longitudes(i,j) - noise_centre_x
              noise_loc(2) = noise_radius*sin_latitudes(i,j)*sin_longitudes(i,j) - noise_centre_y
              noise_loc(3) = noise_radius*cos_latitudes(i,j) - noise_centre_z
              noise(i,j,k) = noise_scale_factor*snoise4d(noise_loc)
            end do
            call get_alternative_temperature(inputs, noise, t_papillon, tpert_papillon)
            ! TODO: does this ml model output sd of theta or t? check
            tpert(i,j,:) = tpert_papillon - t_papillon
        end do
    end do
end subroutine compute_t_pert
end module papillon_alg_mod
