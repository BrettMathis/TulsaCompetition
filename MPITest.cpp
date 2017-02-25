// MPITest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Include\mpi.h"

int main(int argc, char** argv)
{
	MPI_Init(NULL, NULL);

	int num_processes;
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name;
	MPI_Get_processor_name(processor_name, &name);

	printf("Hello world from processor %s, rank %d" " out of %d processors\n", processor_name, rank, num_processes);

	MPI_Finalize();

    return 0;
}

