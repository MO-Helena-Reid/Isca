! ---------------------------------------------------------------
!    Reverses the order of an input array. Based on https://staff.cs.manchester.ac.uk/~fumie/internal/fortranguide/chap08/reverse.html
! ---------------------------------------------------------------
! modified to output to a different array to the input
MODULE reverse_array_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE reverse(a,b)
    IMPLICIT NONE
    REAL, INTENT(IN), DIMENSION(:) :: a  ! input array
    REAL, INTENT(OUT), DIMENSION(:) :: b ! output array
    INTEGER            :: Head           ! pointer moving forward
    INTEGER            :: Tail           ! pointer moving backward
    REAL               :: Temp

    b(:) = a(:)
    Head = 1                            ! start with the beginning
    Tail = SIZE(a)                      ! start with the end
    DO                                  ! for each pair...
      IF (Head >= Tail)  EXIT           !    if Head crosses Tail, exit
      Temp    = b(Head)                 !    otherwise, swap them
      b(Head) = b(Tail)
      b(Tail) = Temp
      Head    = Head + 1                !    move forward
      Tail    = Tail - 1                !    move backward
    END DO                              ! loop back
END SUBROUTINE reverse
END MODULE reverse_array_mod
