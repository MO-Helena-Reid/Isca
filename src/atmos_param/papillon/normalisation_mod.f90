! Created by  on 04/12/2023.

module normalisation_mod
use ml_constants_mod, only: r_ml
implicit none
contains
subroutine normalise_inputs(inputs, outputs)
    implicit none
    real(kind=r_ml), dimension(213), intent(in) :: inputs
    real(kind=r_ml), dimension(213), intent(out) :: outputs
    integer :: i
    do i = 1, 70
        call normalise_humidity(inputs(i), outputs(i)) ! q
        call normalise_pressure(inputs(i + 70), outputs(i + 70)) ! P
        call normalise_temperature(inputs(i + 140), outputs(i + 140)) ! T
    end do
    outputs(211) = inputs(211) ! lsm
    call normalise_orography(inputs(212), outputs(212)) ! orog
    call normalise_orography(inputs(213), outputs(213)) ! std orog
end subroutine normalise_inputs

subroutine normalise_temperature(input, output)
    implicit none
    real(kind=r_ml), intent(in) :: input
    real(kind=r_ml), intent(out) :: output
    real(kind=r_ml), parameter :: min_val = 140.
    real(kind=r_ml), parameter :: max_val = 320.
    call normalise(input, output, min_val, max_val)
end subroutine normalise_temperature

subroutine normalise_pressure(input, output)
    implicit none
    real(kind=r_ml), intent(in) :: input
    real(kind=r_ml), intent(out) :: output
    real(kind=r_ml), parameter :: min_val = 0.
    real(kind=r_ml), parameter :: max_val = 106000.
    call normalise(input, output, min_val, max_val)
end subroutine normalise_pressure
    
subroutine normalise_humidity(input, output)
    implicit none
    real(kind=r_ml), intent(in) :: input
    real(kind=r_ml), intent(out) :: output
    real(kind=r_ml), parameter :: min_val = 0.
    real(kind=r_ml), parameter :: max_val = 0.025
    call normalise(input, output, min_val, max_val)
end subroutine normalise_humidity
    
subroutine normalise_orography(input, output)
    implicit none
    real(kind=r_ml), intent(in) :: input
    real(kind=r_ml), intent(out) :: output
    real(kind=r_ml), parameter :: min_val = 0.
    real(kind=r_ml), parameter :: max_val = 4000.
    call normalise(input, output, min_val, max_val)
end subroutine normalise_orography
    
subroutine normalise(input, output, minval, maxval)
    implicit none
    real(kind=r_ml), intent(in) :: input
    real(kind=r_ml), intent(in) :: minval, maxval
    real(kind=r_ml), intent(out) :: output
    output = (input - minval) / (maxval - minval)
end subroutine normalise
end module normalisation_mod