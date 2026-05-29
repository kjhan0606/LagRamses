SUBROUTINE quick_sort(list, order, n)
  
  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  
  use amr_parameters, ONLY: qdp
  IMPLICIT NONE
  INTEGER :: n
  REAL(qdp), DIMENSION (1:n), INTENT(INOUT)  :: list
  INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order
  
  ! Local variable
  INTEGER :: i
  
  DO i = 1, n
     order(i) = i
  END DO
  
  CALL quick_sort_1(1, n)
  
CONTAINS
  
  RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)
    
    use amr_parameters, ONLY: qdp
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    REAL(qdp)              :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6
    
    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort(left_end, right_end)
       
    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1
       
       DO
          ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (list(i) >= reference) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (list(j) <= reference) EXIT
          END DO
          
          
          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO
       
       IF (left_end < j) CALL quick_sort_1(left_end, j)
       IF (i < right_end) CALL quick_sort_1(i, right_end)
    END IF
    
  END SUBROUTINE quick_sort_1
  
  
  SUBROUTINE interchange_sort(left_end, right_end)
    
    use amr_parameters, ONLY: qdp
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    REAL(qdp)           :: temp
    
    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) > list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO
    
  END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort
!########################################################################
!########################################################################
!########################################################################
! OMP-parallel quick sort over qdp keys with integer companion array.
! Sorts list ascending and returns original-position indices in order.
! Strategy: serial Hoare quicksort with OMP TASKs spawned for the larger
! sub-range when its length exceeds task_threshold. Below seq_threshold
! the entire call falls back to plain quick_sort.
!########################################################################
SUBROUTINE quick_sort_omp(list, order, n)
  use amr_parameters, ONLY: qdp
  IMPLICIT NONE
  INTEGER :: n
  REAL(qdp), DIMENSION (1:n), INTENT(INOUT) :: list
  INTEGER,   DIMENSION (1:n), INTENT(OUT)   :: order

  INTEGER :: i
  INTEGER, PARAMETER :: seq_threshold  = 200000
  INTEGER, PARAMETER :: task_threshold = 50000

  IF (n <= 1) RETURN

  IF (n < seq_threshold) THEN
     CALL quick_sort(list, order, n)
     RETURN
  END IF

  !$OMP PARALLEL DO SCHEDULE(STATIC)
  DO i = 1, n
     order(i) = i
  END DO
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP SINGLE
  CALL quick_sort_task(list, order, 1, n, task_threshold)
  !$OMP END SINGLE
  !$OMP END PARALLEL

END SUBROUTINE quick_sort_omp

!------------------------------------------------------------------------
! Recursive task-spawning quicksort body. Operates in place on list/order
! over the [left_end, right_end] range. Spawns an OMP TASK for the larger
! partition when it exceeds task_threshold so the running thread keeps
! descending the smaller side (bounds task stack depth ~ log2(n)).
!------------------------------------------------------------------------
RECURSIVE SUBROUTINE quick_sort_task(list, order, left_end, right_end, task_threshold)
  use amr_parameters, ONLY: qdp
  IMPLICIT NONE
  ! Assumed-size dummies: caller (quick_sort_omp) is external with no explicit
  ! interface, so passing assumed-shape DIMENSION(:) here corrupts the array
  ! descriptor at the call site (segfault on first list((l+r)/2) read).
  REAL(qdp), DIMENSION(*), INTENT(INOUT) :: list
  INTEGER,   DIMENSION(*), INTENT(INOUT) :: order
  INTEGER, INTENT(IN) :: left_end, right_end, task_threshold

  INTEGER :: i, j, itemp
  REAL(qdp) :: reference, temp
  INTEGER, PARAMETER :: max_simple_sort_size = 6

  IF (right_end < left_end + max_simple_sort_size) THEN
     ! Interchange sort for small ranges
     DO i = left_end, right_end - 1
        DO j = i+1, right_end
           IF (list(i) > list(j)) THEN
              temp     = list(i);  list(i)  = list(j);  list(j)  = temp
              itemp    = order(i); order(i) = order(j); order(j) = itemp
           END IF
        END DO
     END DO
     RETURN
  END IF

  ! Hoare partition
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1
  DO
     DO; i = i + 1; IF (list(i) >= reference) EXIT; END DO
     DO; j = j - 1; IF (list(j) <= reference) EXIT; END DO
     IF (i < j) THEN
        temp  = list(i);  list(i)  = list(j);  list(j)  = temp
        itemp = order(i); order(i) = order(j); order(j) = itemp
     ELSE IF (i == j) THEN
        i = i + 1
        EXIT
     ELSE
        EXIT
     END IF
  END DO

  ! Spawn task for left half (if large), descend serially on right half.
  ! TASKWAIT before return so left-half writes complete before parent returns.
  IF (left_end < j) THEN
     IF (j - left_end > task_threshold) THEN
        !$OMP TASK SHARED(list, order) FIRSTPRIVATE(left_end, j, task_threshold)
        CALL quick_sort_task(list, order, left_end, j, task_threshold)
        !$OMP END TASK
     ELSE
        CALL quick_sort_task(list, order, left_end, j, task_threshold)
     END IF
  END IF
  IF (i < right_end) THEN
     IF (right_end - i > task_threshold) THEN
        !$OMP TASK SHARED(list, order) FIRSTPRIVATE(i, right_end, task_threshold)
        CALL quick_sort_task(list, order, i, right_end, task_threshold)
        !$OMP END TASK
     ELSE
        CALL quick_sort_task(list, order, i, right_end, task_threshold)
     END IF
  END IF
  !$OMP TASKWAIT
END SUBROUTINE quick_sort_task
!########################################################################
!########################################################################
!########################################################################
!########################################################################
!########################################################################
SUBROUTINE quick_sort_dp(list, order, n)
  use amr_parameters, ONLY: dp
  IMPLICIT NONE
  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.


  INTEGER :: n
  REAL(dp), DIMENSION (1:n), INTENT(INOUT)  :: list
  INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order

  ! Local variable
  INTEGER :: i

  DO i = 1, n
     order(i) = i
  END DO

  CALL quick_sort_1_dp(1, n)

CONTAINS

  RECURSIVE SUBROUTINE quick_sort_1_dp(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    !     Local variables
    INTEGER             :: i, j, itemp
    REAL(kind=8)              :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6

    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort_dp(left_end, right_end)

    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1

      DO
          ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (list(i) >= reference) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (list(j) <= reference) EXIT
          END DO


          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO
       IF (left_end < j) CALL quick_sort_1_dp(left_end, j)
       IF (i < right_end) CALL quick_sort_1_dp(i, right_end)
    END IF

  END SUBROUTINE quick_sort_1_dp


  SUBROUTINE interchange_sort_dp(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    !     Local variables                                                                                                                                                                           
    INTEGER             :: i, j, itemp
    REAL(kind=8)           :: temp

    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) > list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO

  END SUBROUTINE interchange_sort_dp

END SUBROUTINE quick_sort_dp




