#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <time.h>

#define RND MPFR_RNDN
#define precision 70

void send_mpfr_row(mpfr_t *row, int count, int destination, int tag, MPI_Comm communicator){
	int numLimbs = ceil((float)precision / mp_bits_per_limb);
	int bytesPerLimb = mp_bits_per_limb / 8;
	int bytesPerVar = 1 + sizeof(mpfr_exp_t) + numLimbs * bytesPerLimb;

	uint8_t *toSend = (uint8_t *)malloc(sizeof(uint8_t) * bytesPerVar * count);

	int locationToStore = 0;

	//Add each variable to the byte array, in order.
	for(int currentNum = 0; currentNum < count; currentNum++){
		toSend[locationToStore] = mpfr_signbit(row[currentNum]) == 0 ? 0 : 1;
		locationToStore++;

		mpfr_exp_t exp = mpfr_get_exp(row[currentNum]);
		uint8_t *expBytes = (uint8_t *) &exp;
		for(int i = 0; i < sizeof(mpfr_exp_t); i++){
			toSend[locationToStore] = expBytes[i];
			locationToStore++;
		}


		for(int limb = 0; limb < numLimbs; limb++){
			uint8_t *limbBytes = (uint8_t *) &(row[currentNum][0]._mpfr_d[limb]);

			for(int byte = 0; byte < bytesPerLimb; byte++){
				toSend[locationToStore] = limbBytes[byte];
				locationToStore++;
			}
		}
	}

	//toSend is now fully made. Lets send it!

	MPI_Send(toSend, bytesPerVar * count, MPI_BYTE, destination, tag, MPI_COMM_WORLD);

	free(toSend);
}

//Sends the row to all nodes larger than this one.
void send_mpfr_row_above(mpfr_t *row, int count, int my_rank, int world_size, int tag, MPI_Comm communicator){
	int numLimbs = ceil((float)precision / mp_bits_per_limb);
	int bytesPerLimb = mp_bits_per_limb / 8;
	int bytesPerVar = 1 + sizeof(mpfr_exp_t) + numLimbs * bytesPerLimb;

	uint8_t *toSend = (uint8_t *)malloc(sizeof(uint8_t) * bytesPerVar * count);

	int locationToStore = 0;

	//mpfr_printf("Sending: %.12Rf\n", row[0]);

	for(int currentNum = 0; currentNum < count; currentNum++){
		toSend[locationToStore] = mpfr_signbit(row[currentNum]) == 0 ? 0 : 1;
		locationToStore++;

		mpfr_exp_t exp = mpfr_get_exp(row[currentNum]);
		uint8_t *expBytes = (uint8_t *) &exp;
		for(int i = 0; i < sizeof(mpfr_exp_t); i++){
			toSend[locationToStore] = expBytes[i];
			locationToStore++;
		}


		for(int limb = 0; limb < numLimbs; limb++){
			uint8_t *limbBytes = (uint8_t *) &(row[currentNum][0]._mpfr_d[limb]);

			for(int byte = 0; byte < bytesPerLimb; byte++){
				toSend[locationToStore] = limbBytes[byte];
				locationToStore++;
			}
		}
		/*if(currentNum == 0) {
			//mpfr_printf("%Nx\n", row[currentNum][0]._mpfr_d, numLimbs);
			printf("Sent  exp:");
			for(int i = 0; i < sizeof(mpfr_exp_t); i++){
				printf("%d", expBytes[i]);
			}
			printf("\n");
		}*/
		//if(currentNum == 0) mpfr_printf("Sent : %Nx\n", row[currentNum][0]._mpfr_d, numLimbs);
	}

	for(int i = world_size - 1; i > my_rank; i--){
		MPI_Send(toSend, bytesPerVar * count, MPI_BYTE, i, tag, MPI_COMM_WORLD);
	}

	free(toSend);

}





//row must be filled with initialized mpfr_t variables of the same precision as the ones being received.
void recv_mpfr_row(mpfr_t *row, int count, int source, int tag, MPI_Comm communicator, MPI_Status* status){

	int numLimbs = ceil((float)precision / mp_bits_per_limb);
	int bytesPerLimb = mp_bits_per_limb / 8;
	int bytesPerVar = 1 + sizeof(mpfr_exp_t) + numLimbs * bytesPerLimb;

	uint8_t *recvd = (uint8_t *)malloc(sizeof(uint8_t) * bytesPerVar * count);

	MPI_Recv(recvd, bytesPerVar * count, MPI_BYTE, source, tag, communicator, status);

	int startIndex;
	for(int currentNum = 0; currentNum < count; currentNum++){
		startIndex = currentNum * bytesPerVar;

		//First lets set the exponent
		uint8_t *expBytes = (uint8_t *)malloc(sizeof(mpfr_exp_t));
		for(int i = 0; i < sizeof(mpfr_exp_t); i++){
			expBytes[i] = *(recvd + startIndex + 1 + i);
		}
		mpfr_exp_t *expPointer = (mpfr_exp_t *)expBytes;

		if(mpfr_set_exp(row[currentNum], *expPointer) != 0){ //Exponent is some special value.
			if(expBytes[0] == 1) mpfr_set_zero(row[currentNum], recvd[startIndex] == 0 ? 1 : -1);
			else mpfr_set_inf(row[currentNum], recvd[startIndex] == 0 ? 1 : -1);

			free(expBytes);
			continue;
		}
		/*if(currentNum == 0) {
			//mpfr_printf("%Nx\n", row[currentNum][0]._mpfr_d, numLimbs);
			printf("Recvd exp:");
			for(int i = 0; i < sizeof(mpfr_exp_t); i++){
				printf("%d", expBytes[i]);
			}
			printf("\n");
		}*/
		free(expBytes);


		//Now set the mantissa/limbs
		uint8_t *limbBytes = (uint8_t *)malloc(bytesPerLimb * numLimbs);
		int limbLocation = startIndex + 1 + sizeof(mpfr_exp_t);
		for(int i = 0; i < bytesPerLimb * numLimbs; i++){
			limbBytes[i] = recvd[limbLocation + i];
		}
		mp_limb_t *limbsPointer = (mp_limb_t *)limbBytes;
		row[currentNum][0]._mpfr_d = limbsPointer;
		//free(limbBytes); //Don't do this! limbBytes is a pointer to the values we're pointing to with the MPFR variable!

		//Now set the sign
		if((recvd[startIndex] == 0) != (mpfr_signbit(row[currentNum]) == 0))  mpfr_neg(row[currentNum], row[currentNum], RND);
	}
	free(recvd);
}





int main(int agrc, char** argv){

	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); //Get number of nodes being run on currently.

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); //Get this node's number.


	mpfr_set_default_prec(precision);

	//Pre-reference the array we'll use for doing math with.
	mpfr_t **myRows;
	int numRows;
	int minRow = 0;
	int maxRow;
	int arrayDimension = 0032;

	clock_t t1, t2;

	if(world_rank == 0){
		char f_name[50];
		double **rows = (double **)malloc(sizeof(double *) * arrayDimension);
		//Create filename
		sprintf(f_name,"m0032x0032.bin");
		printf("Reading array file %s of size %dx%d\n",f_name,arrayDimension,arrayDimension);
		//Open file
		FILE *datafile=fopen(f_name,"rb");
		//Read elelements
		for (int i=0; i < arrayDimension; i++){
			rows[i] = (double *)malloc(sizeof(double) * arrayDimension);
			for (int j=0; j < arrayDimension; j++){
				fread(&rows[i][j],sizeof(double),1,datafile);
			}
		}


		//Send each node a copy of its rows - excluding this one.
		for(int node = 1; node < world_size; node++){
			int nodeMin = arrayDimension / world_size * node;
			int nodeMax = (arrayDimension / world_size * (node + 1));

			if(node == world_size - 1) nodeMax = arrayDimension; //If the array dimension isn't divisible by 8, give the extras to the last node.

			//Tell the node what arrayDimension is, so it can figure out how many rows/columns to expect.
			//We'll use 0 as the tag even though it breaks the general pattern for the rest of the program.
			//We'd use -1, but it errors on negative tags.
			MPI_Send(&arrayDimension, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
			printf("Arraydimension: %d", arrayDimension);

			for(int row = nodeMin; row < nodeMax; row++){
				MPI_Send(rows[row], arrayDimension, MPI_DOUBLE, node, row, MPI_COMM_WORLD);
			}
		}

		//Now, set up this node's mpfr_t array.
		numRows = arrayDimension / world_size;

		myRows = (mpfr_t **)malloc(sizeof(mpfr_t *) * numRows);

		for(int row = 0; row < numRows; row++){
			myRows[row] = (mpfr_t *)malloc(sizeof(mpfr_t) * arrayDimension);
			for(int column = 0; column < arrayDimension; column++){
				mpfr_init_set_d(myRows[row][column], rows[row][column], RND);
			}
		}

		//free the array of doubles.
		for(int i = 0; i < arrayDimension; i++){
			free(rows[i]);
		}
		free(rows);

	} else {
		arrayDimension;
		MPI_Recv(&arrayDimension, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//printf("Node %d received %d\n", world_rank, arrayDimension); //Check for the curious.

		numRows = arrayDimension / world_size;
		//Again, if this is the last node make sure to absorb any extra rows.
		if(world_rank == world_size - 1) numRows += arrayDimension % world_size;

		minRow = arrayDimension / world_size * world_rank;

		double **dRows = (double**)malloc(sizeof(double *) * numRows);
		for(int i = minRow; i < minRow + numRows; i++){
			dRows[i - minRow] = (double *)malloc(sizeof(double) * arrayDimension);
			MPI_Recv(dRows[i - minRow], arrayDimension, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		myRows = (mpfr_t **)malloc(sizeof(mpfr_t *) * numRows);
		for(int row = 0; row < numRows; row++){
			myRows[row] = (mpfr_t *)malloc(sizeof(mpfr_t) * arrayDimension);
			for(int column = 0; column < arrayDimension; column++){
				mpfr_init_set_d(myRows[row][column], dRows[row][column], RND);
			}
		}

		//free the array of doubles.
		for(int i = 0; i < numRows; i++){
			free(dRows[i]);
		}
		free(dRows);
	}

	//Now we start row reducing.
	maxRow = minRow + numRows;
	mpfr_t coefficient;
	mpfr_t temp;
	mpfr_init2(coefficient, precision);
	mpfr_init2(temp, precision);
	int startRow = 0;

	mpfr_t *reduRow = malloc(sizeof(mpfr_t) * arrayDimension);
	for(int i = 0; i < arrayDimension; i++) mpfr_init2(reduRow[i], precision);

	for(int column = 0; column < maxRow; column++){
		//if(world_rank == 7) printf("%d\n", column);
		//First, lets find the row we'll reduce with. This is the row# that matches the column#
		if(column >= minRow && column < maxRow){ //If this node has the row
			reduRow = myRows[column - minRow]; //My array has minRow starting at zero, account for this
			send_mpfr_row_above(reduRow, arrayDimension, world_rank, world_size, column, MPI_COMM_WORLD);
		} else {
			if(column == maxRow){
				reduRow = malloc(sizeof(mpfr_t) * arrayDimension);
				for(int i = 0; i < arrayDimension; i++) mpfr_init2(reduRow[i], precision);
			}
			recv_mpfr_row(reduRow, arrayDimension, MPI_ANY_SOURCE, column, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		if(column < minRow) startRow = 0;
		else startRow = column + 1 - minRow;

		for(int row = startRow; row < numRows; row++){
			mpfr_div(coefficient, myRows[row][column], reduRow[column], RND);

			for(int i = column; i < arrayDimension; i++){
				mpfr_mul(temp, coefficient, reduRow[i], RND);
				mpfr_sub(myRows[row][i], myRows[row][i], temp, RND);
			}
		}
	}

	//Now multiply the diagonals
	mpfr_t multiplied;
	mpfr_init_set(multiplied, myRows[0][minRow], RND);
	for(int i = 1; i < numRows; i++){
		mpfr_mul(multiplied, multiplied, myRows[i][minRow + i], RND);
	}

	//Now send the values to node 0
	if(world_rank == 0){
		for(int node = 1; node < world_size; node++){
			recv_mpfr_row(&temp, 1, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			mpfr_mul(multiplied, multiplied, temp, RND);
		}
		//Print the answer:
		mpfr_printf("Det: %.12Rf\n", multiplied);

		//Get log(abs(det))
		if(mpfr_signbit(multiplied) != 0) mpfr_neg(multiplied, multiplied, RND);
		mpfr_log10(multiplied, multiplied, RND);
		mpfr_printf("Log(|Det|): %.12Rf\n", multiplied);


		//Print the array

		/*for(int i = 0; i < 2; i++){
			for(int j = 0; j < arrayDimension; j++){
				mpfr_printf("%.5Rf, ", myRows[i][j]);
			}
			printf("\n");
		}

		mpfr_t *rowToPrint = (mpfr_t *)malloc(sizeof(mpfr_t) * arrayDimension);
		for(int i = 0; i < arrayDimension; i++) mpfr_init(rowToPrint[i]);

		for(int i = 2; i < arrayDimension; i++){
			recv_mpfr_row(rowToPrint, arrayDimension, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(int j = 0; j < arrayDimension; j++){
				mpfr_printf("%.5Rf, ", rowToPrint[j]);
			}
			printf("\n");
		}*/
	} else {
		send_mpfr_row(&multiplied, 1, 0, 0, MPI_COMM_WORLD);

		/*
		send_mpfr_row(myRows[0], arrayDimension, 0, minRow + 0, MPI_COMM_WORLD);
		send_mpfr_row(myRows[1], arrayDimension, 0, minRow + 1, MPI_COMM_WORLD);
		*/
	}


	// Finalize the MPI environment.
	MPI_Finalize();
}
