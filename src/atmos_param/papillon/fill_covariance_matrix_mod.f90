! Created by  on 08/11/2023.

module fill_covariance_matrix_mod
use ml_constants_mod, only: r_ml, ml_nlev, ml_output_len
implicit none
    contains
    subroutine fill_covariance_matrix(diag_70_69_68, covariance_matrix)
        implicit none
        real(kind=r_ml), dimension(ml_output_len), intent(in) :: diag_70_69_68
        real(kind=r_ml), dimension(ml_nlev, ml_nlev), intent(out) :: covariance_matrix
        integer i
        covariance_matrix = 0_r_ml
        do i=1,ml_nlev
            covariance_matrix(i,i) = diag_70_69_68(i)
        end do
    end subroutine fill_covariance_matrix
end module fill_covariance_matrix_mod