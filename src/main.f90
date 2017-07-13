
PROGRAM TEST
IMPLICIT NONE

  WRITE (*,*) 'TESTING...'
  
END PROGRAM TEST



SUBROUTINE GAUSS_JORDAN_REDUCTION(A, B, NDIM, N, ERROR)
IMPLICIT NONE

  INTEGER, INTENT(IN) :: NDIM
  REAL, INTENT(INOUT), DIMENSION(NDIM, NDIM) :: A
  REAL, INTENT(INOUT), DIMENSION(NDIM) :: B
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(OUT) :: ERROR
  
  REAL, PARAMETER :: EPSILON = 1.0E-6
  
  REAL :: FACTOR
  
  INTEGER :: IROW
  INTEGER :: IPEAK
  INTEGER :: JROW
  INTEGER :: KCOL
  
  REAL :: TEMP
  
  ! Process N times to get all equations
  MAIN_LOOP:  DO IROW = 1, N
  
    ! Find peak pivot for column irow in rows irow to n
    IPEAK = IROW
    MAX_PIVOT: DO JROW = IROW + 1, N
      IF (ABS(A(JROW,IROW)) > ABS(A(IPEAK,IROW))) THEN
        IPEAK = JROW
      END IF
    END DO MAX_PIVOT
    
    ! Check for singular equations
  
  END DO MAIN_LOOP
  
END SUBROUTINE GAUSS_JORDAN_REDUCTION
