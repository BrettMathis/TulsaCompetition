EXECS=mpi_mpfr_det
MPICC?=mpicc

all: ${EXECS}

mpi_mpfr_det: mpi_mpfr_det.c
	${MPICC} -o mpi_mpfr_det mpi_mpfr_det.c -lmpfr -lgmp -lm
clean:
	rm ${EXECS}
