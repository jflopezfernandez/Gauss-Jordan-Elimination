
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
    SINGULAR: IF (ABS(A(IPEAK,IROW)) < EPSILON) THEN
      ERROR = 1
      RETURN
    END IF singular
      
    ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
    SWAP_EQUATION: IF (IPEAK /= IROW) THEN
      DO KCOL = 1,N
        TEMP = A(IPEAK, KCOL)
        A(IPEAK,KCOL) = A(IROW,KCOL)
        A(IROW,KCOL) = TEMP
      END DO
      
      TEMP = B(IPEAK)
      B(IPEAK) = B(IROW)
      B(IROW) = TEMP
    END IF SWAP_EQUATION
    
    ! Multiply equation irow by -a(jrow,irow)/a(irow,irow), and add
    ! it to equation jrow (for all eqns except irow itself)
    ELIMINATE: DO JROW = 1, N
      IF (JROW /= IROW) THEN
        FACTOR = -A(JROW,IROW)/A(IROW,IR0W)
        
        DO KCOL = 1, N
          A(JROW,KCOL) = A(IROW,KCOL) * FACTOR + A(JROW,KCOL)
        END DO
        
        B(JROW = B(IROW) * FACTOR + B(JROW)
      END IF
    END DO ELIMINATE
    
  END DO MAIN_LOOP
  
  
  ! End of mainloop over all equations. All off-diagonal terms are now zero. To
  ! the final answer, we must divide each equation by the coefficient of its on-diagonal
  ! term.
  
  DIVIDE: DO IROW = 1, N
    B(IROW) = B(IROW) / A(IROW,IROW)
    A(IROW,IROW) = 1.
  END DO DIVIDE
  
  ! Set error flag to zero and return
  ERROR = 0
  
END SUBROUTINE GAUSS_JORDAN_REDUCTION
