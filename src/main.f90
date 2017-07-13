! Author: Stephen J. Chapman
! Fortran for Scientists and Engineers

PROGRAM TEST
IMPLICIT NONE

  INTEGER, PARAMETER :: FILE_UNIT = 25

  INTEGER, PARAMETER :: MAX_SIZE = 10
  REAL, DIMENSION(MAX_SIZE,MAX_SIZE) :: A
  REAL,DIMENSION(MAX_SIZE) :: B
  INTEGER :: ERROR
  CHARACTER(LEN=20) :: FILE_NAME
  
  INTEGER :: I
  INTEGER :: J
  INTEGER :: N
  INTEGER :: ISTAT
  
  WRITE (*,"(' ENTER THE FILE NAME CONTAINING THE EQUATIONS: ')")
  READ (*,'(A20)') FILE_NAME
  
  OPEN (UNIT=FILE_UNIT, FILE=FILE_NAME, STATUS='OLD', ACTION='READ', IOSTAT=ISTAT)
  
  FILEOPEN: IF (ISTAT == 0) THEN
    READ (FILE_UNIT,*) N
    
    SIZE_OK: IF (N <= MAX_SIZE) THEN
      DO I = 1, N
        READ (FILE_UNIT,*) (A(I,J), J=1,N), B(I)
      END DO
      
      CALL GAUSS_JORDAN_REDUCTION(A, B, MAX_SIZE, N, ERROR)
      
      ERROR_CHECK: IF (ERROR /= 0) THEN
        WRITE (*,1010)
        1010 FORMAT (/1X,'ZERO PIVOT ENCOUNTERED!',//1X, 'THERE IS NO UNIQUE SOLUTION TO THIS SYSTEM.')
      ELSE ERROR_CHECK
        WRITE (*,"(/,1X,'THE SOLUTIONS ARE:')")
        
        DO I = 1, N
          WRITE (*,"(3X,'X(',I2,') = ',F16.6)") I, B(I)
        END DO
        
      END IF ERROR_CHECK
    END IF SIZE_OK
  ELSE FILEOPEN
    WRITE (*,1020) ISTAT
    1020 FORMAT (1X,'FILE OPEN FAILED--STATUS = ', I6)
  END IF FILEOPEN
  
  CLOSE (FILE_UNIT)
  
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
        FACTOR = -A(JROW,IROW)/A(IROW,IROW)
        
        DO KCOL = 1, N
          A(JROW,KCOL) = A(IROW,KCOL) * FACTOR + A(JROW,KCOL)
        END DO
        
        B(JROW) = B(IROW) * FACTOR + B(JROW)
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
