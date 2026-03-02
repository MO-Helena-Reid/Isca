! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module sample_multivariate_normal_mod
use ml_constants_mod, only: r_ml, r_um

implicit none
contains
subroutine sample_multivariate_normal(mean, covariance_matrix, output_vector)
    ! Makes output_vector a sample from a multivariate normal distribution
    ! of the specified mean and covariance matrix.
    implicit none
    integer, parameter :: precision = r_ml
    ! arguments
    real(kind=precision), dimension(:), intent(in) :: mean
    real(kind=precision), dimension(:, :), intent(in) :: covariance_matrix
    real(kind=precision), dimension(:), intent(out) :: output_vector
    ! locals
    real(kind=precision), dimension(:), allocatable :: z
    real(kind=precision), dimension(:, :), allocatable :: decomp
    integer :: len_mean, i

    ! Initialise z to be filled with random normals
    len_mean = size(mean)
    allocate(z(len_mean))
    do i = 1, len_mean
        z(i) = rnorm()
    end do
    allocate(decomp(len_mean, len_mean))
    call cholesky_banachiewicz(covariance_matrix, decomp)
    output_vector = mean + matmul(decomp, z)
    deallocate(z)
    deallocate(decomp)
end subroutine sample_multivariate_normal

subroutine cholesky_banachiewicz(positive_definite, decomp)
    ! Implements the Cholesky-Banachiewicz algorithm for performing Cholesky decomposition.
    ! This is a decomposition of a Hermitian, positive-definite matrix into
    ! the product of a lower triangular matrix and its conjugate transpose,
    ! which is useful for efficient numerical solutions.
    implicit none
    ! constants
    integer, parameter :: precision = r_ml
    ! arguments
    real(kind=precision), dimension(:,:), intent(in) :: positive_definite
    real(kind=precision), dimension(:, :), intent(out) :: decomp
    ! locals
    real(kind=precision) :: sum
    integer i, j, k
    complex(kind=precision) :: aijminussigma, logdecomp
    decomp = 0
    do i = 1, size(positive_definite,1)
        do j = 1, i
            sum = 0
            do k = 1, j-1
                sum = sum + decomp(i,k) * decomp(j,k)
            end do
            if (i == j) then
                decomp(i, j) = sqrt(positive_definite(i,i) - sum)
                if (isnan(decomp(i, j))) then
                    print*, 'decomp(', i, j, ') was nan. Dumping vars...'
                    print*, 'sum', sum
                    print*, 'positive_definite(i,i)', positive_definite(i,i)
                    print*, 'positive_definite(i,i)-sum', positive_definite(i,i) - sum
                    sum = 0
                    do k = 1, j-1
                        sum = sum + decomp(i,k) * decomp(j,k)
                        print*, 'k=', k, 'increment', decomp(i,k) * decomp(j,k), 'sum=', sum
                    end do
                    print*, 'exiting...'
                    call exit(1)
                endif
            else
                aijminussigma = cmplx(positive_definite(i, j) - sum)
                logdecomp = log(aijminussigma) - log(decomp(j, j))
                decomp(i, j) = real(exp(logdecomp))
                if (isnan(decomp(i, j))) then
                    print*, 'decomp(', i, j, ') was nan. Dumping vars...'
                    print*, 'sum', sum
                    print*, 'positive_definite(i, j)', positive_definite(i, j)
                    print*, 'decomp(j, j)', decomp(j, j)
                    print*, 'log(decomp(j, j))', log(decomp(j, j))
                    print*, 'aijminussigma', aijminussigma
                    print*, 'log(aijminussigma)', log(aijminussigma)
                    print*, 'logdecomp', logdecomp
                    print*, 'exiting...'
                    call exit(1)
                end if
            end if
        end do
    end do
end subroutine cholesky_banachiewicz

subroutine sample_independent_normal(mean, stds, output_vector)
    ! Makes output_vector a sample from a multivariate normal distribution
    ! of the specified mean and standard deviations, where all other
    ! covariances are set to zero, to give completely uncorrelated
    ! samples from level to level.
    implicit none
    ! arguments
    real(kind=r_um), dimension(:), intent(in) :: mean
    real(kind=r_ml), dimension(:), intent(in) :: stds
    real(kind=r_um), dimension(:), intent(out) :: output_vector
    ! locals
    integer :: i

    do i = 1, size(output_vector)
        output_vector(i) = real(stds(i)*rnorm(), r_um) + mean(i)
    end do
end subroutine sample_independent_normal

subroutine sample_from_snoise(noise, mean, stds, output_vector)
    ! Makes output_vector from noise with the specified mean and standard deviations.
    implicit none
    ! arguments
    real(kind=r_um), dimension(:), intent(in) :: noise
    real(kind=r_um), dimension(:), intent(in) :: mean
    real(kind=r_ml), dimension(:), intent(in) :: stds
    real(kind=r_um), dimension(:), intent(out) :: output_vector
    ! locals
    integer :: i
    real(kind=r_um), parameter :: one_over_sd_snoise = 1.0_r_um/0.30210421270734245_r_um
    do i = 1, size(output_vector)
        ! divide by the standard deviation of simplex noise to get something with sd of 1,
        ! then multiply by the desired standard deviation and add the mean.
        output_vector(i) = one_over_sd_snoise*noise(i)*real(stds(i), r_um) + mean(i)
    end do
end subroutine sample_from_snoise

subroutine sample_uniform_normal(mean, stds, output_vector)
    ! Makes output_vector a sample where the same quantile is chosen
    ! for all levels, but the means and standard deviations are different.
    implicit none
    ! arguments
    real(kind=r_um), dimension(:), intent(in) :: mean
    real(kind=r_ml), dimension(:), intent(in) :: stds
    real(kind=r_um), dimension(:), intent(out) :: output_vector
    real(kind=r_ml) :: quantile
    ! locals
    integer :: i

    quantile = rnorm()
    do i = 1, size(output_vector)
        output_vector(i) = real(stds(i)*quantile, r_um) + mean(i)
    end do
end subroutine sample_uniform_normal

function rnorm() result( fn_val )
    !   Adapted from code released into the public domain by Alan Miller
    !   https://jblevins.org/mirror/amiller/rnorm.f90
    !   This version doubles the computations required for many calls
    !   but is thread-safe.
    !   Generate a random normal deviate using the polar method.
    !   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
    !              normal variables', Siam Rev., vol.6, 260-264, 1964.

    implicit none
    real(r_ml)  :: fn_val

    ! Local variables
    real(r_ml)            :: u, v, sum, sln
    real(r_ml), parameter :: one = 1.0, vsmall = tiny( one )

    ! Generate a pair of random normals
    do
    call random_number( u )
    call random_number( v )
    u = scale( u, 1 ) - one
    v = scale( v, 1 ) - one
    sum = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
    if (sum < one) exit
    end do
    sln = sqrt(- scale( log(sum), 1 ) / sum)
    fn_val = u*sln
return
end function rnorm
end module sample_multivariate_normal_mod