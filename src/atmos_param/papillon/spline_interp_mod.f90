!> Module for cubic spline interpolation
!>
!> Released into the public domain by Jannis Teunissen
!>
!> https://github.com/jannisteunissen/spline_interpolation_fortran
module spline_interp_mod
  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type spline_type
     !> Number of tabulated points
     integer               :: n
     !> Tabulated x values
     real(dp), allocatable :: x(:)
     !> Tabulated y values
     real(dp), allocatable :: y(:)
     !> Spline coefficients such that s(x) = y(i) + b(i)*(x-x(i)) +
     !> c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
     real(dp), allocatable :: bcd(:, :)
     !> Array to quickly get index for a coordinate
     integer, allocatable  :: lookup_index(:)
     !> Inverse spacing in lookup array
     real(dp)              :: lookup_inv_dx
  end type spline_type

  public :: spline_type
  public :: spline_set_coeffs
  public :: spline_evaluate
  public :: spline_evaluate_1d

contains


  !> Calculate coefficients for cubic spline interpolation
  subroutine spline_set_coeffs(x, y, n, spl)
    !> Size of tabulated data
    integer, intent(in)         :: n
    !> Tabulated coordinates (in strictly increasing order)
    real(dp), intent(in)        :: x(n)
    !> Tabulated values
    real(dp), intent(in)        :: y(n)
    type(spline_type), intent(out) :: spl
    real(dp), allocatable       :: b(:), c(:), d(:)
    integer                     :: i, j, lookup_size
    real(dp)                    :: h, x_lookup, min_dx
    !print*,GETPID(),"BEGIN SUBROUTINE spline_set_coeffs"
    if (n < 2) &
         error stop "spline_set_coeffs requires n >= 2"
    if (any(x(2:n) <= x(1:n-1))) &
         error stop "spline_set_coeffs: x(:) not strictly increasing"

    allocate(spl%x(n), spl%y(n), spl%bcd(3, n), b(n), c(n), d(n))
    spl%n = n
    spl%x = x
    spl%y = y
    !print*,GETPID(),"allocated and set n x y"
    if (n < 3) then
     ! print*,GETPID(),"special case n < 3 early exit ahead"
       ! Handle special case by linear interpolation
       b(1) = (y(2)-y(1))/(x(2)-x(1))
       c(1) = 0
       d(1) = 0
       b(2) = b(1)
       c(2) = 0
       d(2) = 0
       return
    end if

    ! set up tridiagonal system
    ! b = diagonal, d = offdiagonal, c = right hand side.
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, n-1
       d(i) = x(i+1) - x(i)
       b(i) = 2*(d(i-1) + d(i))
       c(i+1) = (y(i+1) - y(i))/d(i)
       c(i) = c(i+1) - c(i)
    end do
    !print*,GETPID(),"setup tri diagonal system"
    ! end conditions.  third derivatives at  x(1)  and  x(n)
    ! obtained from divided differences
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0
    c(n) = 0
    if(n /= 3) then
       c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
       c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
       c(1) = c(1)*d(1)**2/(x(4)-x(1))
       c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if
    !print*,GETPID(),"obtained third derivatives"
    ! forward elimination
    do i = 2, n
       h = d(i-1)/b(i-1)
       b(i) = b(i) - h*d(i-1)
       c(i) = c(i) - h*c(i-1)
    end do
    !print*,GETPID(),"did forward elimination"
    ! back substitution
    c(n) = c(n)/b(n)
    do j = 1, n-1
       i = n-j
       c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do
    !print*,GETPID(),"did back substitution"
    ! compute spline coefficients
    b(n) = (y(n) - y(n-1))/d(n-1) + d(n-1)*(c(n-1) + 2*c(n))
    do i = 1, n-1
       b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2*c(i))
       d(i) = (c(i+1) - c(i))/d(i)
       c(i) = 3.*c(i)
    end do
    c(n) = 3*c(n)
    d(n) = d(n-1)
    !print*,GETPID(),"computed spline coefficients"
    spl%bcd(1, :) = b
    spl%bcd(2, :) = c
    spl%bcd(3, :) = d
    !print*,GETPID(),"set spline bcd"
    ! Create linear lookup table to find location between data points more
    ! quickly. First determine good size for the lookup table.
    min_dx      = minval(x(2:n) - x(1:n-1))
    lookup_size = min(4 * n, nint(1 + (x(n) - x(1))/min_dx))
    !print*,GETPID(),"create linear lookup table"
    ! The lookup table will have a regular (linear) spacing
    h = (x(n) - x(1)) / (lookup_size - 1)
    spl%lookup_inv_dx = 1/h
    allocate(spl%lookup_index(lookup_size))
    !print*,GETPID(),"allocate lookup table"
    ! At location z, the table index is i = ceiling((z - x(1)) * inv_dx)
    ! The tabulated points should then be x(spl%lookup_index(i))
    spl%lookup_index(1) = 1
    do i = 2, lookup_size
       x_lookup = x(1) + (i-1) * h
       if (size(spl%lookup_index)<i-1) then
         error stop "i-1 exceeded size of spl lookup index"
       end if
       spl%lookup_index(i) = spl%lookup_index(i-1)
       do while (x_lookup > x(spl%lookup_index(i)+1))
          spl%lookup_index(i) = spl%lookup_index(i) + 1
          if (size(x)<spl%lookup_index(i)+1) then
            !print*,GETPID(),"size(x):",size(x),"spl%lookup_index(i)+1",spl%lookup_index(i)+1, "spl lookup index exceeded size of x"
            exit
          end if
       end do
    end do
  ! print*,GETPID(),"END SUBROUTINE spline_set_coeffs"
  end subroutine spline_set_coeffs

  ! Evaluate the cubic spline interpolation at point u
  ! result = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
  ! where  x(i) <= u <= x(i+1)
  pure elemental function spline_evaluate(u, spl) result(spline_value)
    !> Evaluate at this coordinate
    real(dp), intent(in)       :: u
    !> Spline data
    type(spline_type), intent(in) :: spl
    integer                    :: i, j
    real(dp)                   :: spline_value, dx

    associate (n=>spl%n, x=>spl%x, y=>spl%y, bcd=>spl%bcd)
      if(u <= x(1)) then
         spline_value = y(1)
         return
      else if(u >= x(n)) then
         spline_value = y(n)
         return
      end if

      j = ceiling(spl%lookup_inv_dx * (u - x(1)))

      ! Even with the checks above, we probably have to check whether j < 1 or
      ! j > size(spl%lookup_index) due to numerical round-off error
      if (j < 1) then
         j = 1
      else if (j > size(spl%lookup_index)) then
         j = size(spl%lookup_index)
      end if

      ! TODO: Not sure whether in practical applications it could be useful to
      ! do a binary search here instead of a linear one

      ! Find index i so that x(i) <= u <= x(i+1)
      do i = spl%lookup_index(j), n-1
         if (u <= x(i+1)) exit
      end do

      ! evaluate spline interpolation
      dx = u - x(i)
      spline_value = y(i) + dx*(bcd(1, i) + dx*(bcd(2, i) + dx*bcd(3, i)))
    end associate
  end function spline_evaluate

 ! Evaluate the cubic spline interpolation at point u
  ! result = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
  ! where  x(i) <= u <= x(i+1)
  function spline_evaluate_1d(u, spl) result(spline_value)
    !> Evaluate at this coordinate
    real(dp), intent(in), dimension(:) :: u
    !> Spline data
    type(spline_type), intent(in)      :: spl
    integer                            :: i, j, k
    real(dp), dimension(size(u))             :: spline_value
    real(dp), dimension(size(u))       :: dx
    do k = 1, size(u)
      associate (n=>spl%n, x=>spl%x, y=>spl%y, bcd=>spl%bcd)
        if(u(k) <= x(1)) then
          spline_value = y(1)
          return
        else if(u(k) >= x(n)) then
          spline_value = y(n)
          return
        end if

        j = ceiling(spl%lookup_inv_dx * (u(k) - x(1)))

        ! Even with the checks above, we probably have to check whether j < 1 or
        ! j > size(spl%lookup_index) due to numerical round-off error
        if (j < 1) then
          j = 1
        else if (j > size(spl%lookup_index)) then
          j = size(spl%lookup_index)
        end if

        ! TODO: Not sure whether in practical applications it could be useful to
        ! do a binary search here instead of a linear one

        ! Find index i so that x(i) <= u <= x(i+1)
        do i = spl%lookup_index(j), n-1
          if (u(k) <= x(i+1)) exit
        end do

        ! evaluate spline interpolation
        dx(k) = u(k) - x(i)
        spline_value(k) = y(i) + dx(k)*(bcd(1, i) + dx(k)*(bcd(2, i) + dx(k)*bcd(3, i)))
        if (isnan(dx(k))) print*, GETPID(), "NaN detected in dx at", k
        if (isnan(y(i))) print*, GETPID(), "NaN detected in y at", i
        if (isnan(dx(k)*bcd(1, i))) print*, GETPID(), "NaN detected in d*bcd(1,i) at", k
        if (isnan(dx(k)*bcd(2, i))) print*, GETPID(), "NaN detected in d*bcd(2,i) at", k
        if (isnan(dx(k)*bcd(3, i))) print*, GETPID(), "NaN detected in d*bcd(3,i) at", k
        if (isnan(dx(k)*(bcd(2, i) + dx(k)*bcd(3, i)))) print*, GETPID(), "NaN detected in dx(k)*(bcd(2, i) + dx(k)*bcd(3, i))) at", k
        if (isnan(spline_value(k))) print*, GETPID(), "NaN detected in spline_value at", k
      end associate
    end do
    if (any(isnan(spline_value))) print*, GETPID(), "NaNs detected in spline_value"
  end function spline_evaluate_1d


end module spline_interp_mod
