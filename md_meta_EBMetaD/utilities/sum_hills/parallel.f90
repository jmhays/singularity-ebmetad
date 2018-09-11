SUBROUTINE PARALLEL_INIT(rank,npe)
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: rank,npe
  INTEGER              :: ierr
  INCLUDE "mpif.h"
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierr)
  
END SUBROUTINE PARALLEL_INIT

SUBROUTINE PARALLEL_SUM(input,output,n)
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: n
  DOUBLE PRECISION, INTENT(IN)  :: input(*)
  DOUBLE PRECISION, INTENT(OUT) :: output(*)
  INTEGER                       :: ierr
  INCLUDE "mpif.h"
  CALL MPI_ALLREDUCE(input,output,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
END SUBROUTINE

SUBROUTINE PARALLEL_FINALIZE()
  IMPLICIT NONE
  INTEGER :: ierr
  INCLUDE "mpif.h"
  CALL MPI_FINALIZE(ierr)
END SUBROUTINE PARALLEL_FINALIZE
