! Created by  on 08/11/2023.

module plausible_alternative_mod
use ml_constants_mod, only: r_um, r_ml, ml_nlev, ml_output_len, ml_input_len
use papillon_config_mod, only: sampling_method, sampling_method_correlated, &
        sampling_method_uncorrelated, sampling_method_multivariate, &
        sampling_method_development, sampling_method_constant, &
        constant_t_sd_profile
use log_mod, only: log_event, LOG_LEVEL_ERROR, log_scratch_space
implicit none
contains
subroutine get_alternative_temperature(inputs, noise, temperature, alt_temperature)
    use thetasds00v001_mod, only: thetasds00v001
    use normalisation_mod, only: normalise_inputs
    use sample_multivariate_normal_mod, only: sample_independent_normal, sample_uniform_normal, &
                                              sample_from_snoise
    use fill_covariance_matrix_mod, only: fill_covariance_matrix
    implicit none
    ! ML local variables
    real(kind=r_ml), parameter :: epsilon = 0.00001
!    If epsilon is too large the covariance matrix might be perturbed more than necessary
!    too small and the algorithm is unstable, for example calling tiny(1.) results in NaNs
!    in the cholesky decomp.

    real(kind=r_ml), intent(in) :: inputs(ml_input_len)
    real(kind=r_um), intent(in) :: noise(ml_nlev)
    real(kind=r_um), intent(in) :: temperature(ml_nlev)
    real(kind=r_um), intent(out) :: alt_temperature(ml_nlev)
    integer :: i
    real(kind=r_ml) :: normalised_inputs(ml_input_len)
    real(kind=r_ml) :: y_outputs(ml_output_len)
    real(kind=r_ml) :: t_std(ml_nlev)
    ! program
    !     initialisation
    alt_temperature = 0
    y_outputs = 0
    t_std = 0
    call normalise_inputs(inputs, normalised_inputs)
!    print*, 'normalised inputs', normalised_inputs
    call thetasds00v001(normalised_inputs, y_outputs)
!    do i=1, ml_nlev
!        print*, 'y_outputs', i, y_outputs(i)
!    end do
    y_outputs = y_outputs + epsilon
    t_std = y_outputs(:ml_nlev)
!    call fill_covariance_matrix(y_outputs, t_covariance_matrix)
!    call sample_multivariate_normal(real(temperature, r_ml), &
!            t_covariance_matrix, t_plausible_alternative)
    select case (sampling_method)
        case (sampling_method_constant)
            ! TODO: can we not cast the sd profile
            t_std = real(constant_t_sd_profile, r_ml)
            call sample_from_snoise(noise, temperature, t_std, alt_temperature)
        case (sampling_method_correlated)
            call sample_uniform_normal(temperature, t_std, alt_temperature)
        case (sampling_method_uncorrelated)
            call sample_independent_normal(temperature, t_std, alt_temperature)
        case (sampling_method_development)
            call sample_from_snoise(noise, temperature, t_std, alt_temperature)
        case (sampling_method_multivariate)
            write(log_scratch_space, &
                    '(A)') "Not implemented error, no multivariate normal sampling yet with PAPILLON."
            call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end select
end subroutine get_alternative_temperature
end module plausible_alternative_mod