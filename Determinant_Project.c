#include <stdlib.h>
#include <stdio.h>

#define ARRAYSIZE 16 //needs to be changed with each matrix, 16 is default for the smalled .bin provided


//
//performs the (-1)^n operation without being obtuse by directly computing it
//

double sign(int input)
{
	if ((input % 2) == 0)
		return 1.0;
	else
		return -1.0;
}

//deletes a row from memspace and sends everything up a row
void removeRow(double** matrix, int row, int call)
{
	free(matrix[row]); //this works fine for whatever reason
	while (row < (ARRAYSIZE - call))
	{
		//move data up
		matrix[row] = matrix[row + 1];
		row++;
	}
}

//moves everything left and (tries) to delete the outermost column
void removeColumn(double** matrix, int col, int call)
{
	for (int i = 0;i < (ARRAYSIZE - call); i++)
	{
		while (col < (ARRAYSIZE - call))
		{
			//move data to the left
			matrix[i][col] = matrix[i][col + 1];
			col++;
		}
		//matrix[i] = realloc(matrix[i], sizeof(double*)*(ARRAYSIZE - call)); again, having issues with freeign memspace
	}
}

//
//recursive method
//

double cofactorMethod(double **matrix, int callNum)
{
	double cofactor = 0.0;
	double determinant = 0.0;

	//
	//printing the matrix to make sure it was accurate
	//

	//for (int i = 0; i< ARRAYSIZE; i++)
	//	for (int j = 0; j< ARRAYSIZE; j++)
	//	{
	//		printf("matrix[%d][%d]=%f\n", i, j, matrix[i][j]);
	//	}

	int colSize = (ARRAYSIZE - callNum);
	int rowSize = colSize;
	//printf("%d\n, %d\n", rowSize, colSize);


	//
	//base case
	//
	if (colSize <= 2)
	{
		determinant = (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
	}


	//
	//all other cases
	//
	else
	{
		for (int k = 0;k < (ARRAYSIZE - callNum);k++)
		{
			//
			//creating a minor matrix that's dynamically allocated (THIS was a pain in the butt)
			//
			double *minorMatrix[ARRAYSIZE];
			for (int i = 0; i< ARRAYSIZE; i++)
				minorMatrix[i] = (double *)malloc(ARRAYSIZE * sizeof(double));

			for (int l = 0;l < (ARRAYSIZE - callNum);l++)
			{
				for (int k = 0;k < (ARRAYSIZE - callNum);k++)
				{
					minorMatrix[k][l] = matrix[k][l];
				}
			}
			cofactor = matrix[k][0];
			removeRow(minorMatrix, k, callNum); //method above
			removeColumn(minorMatrix, k, callNum); //method above
			//
			//recursive call
			//
			determinant += cofactor * sign(k) * cofactorMethod(minorMatrix, (callNum + 1));
			printf("current det: %f\n", determinant);
			//free(minorMatrix); having issues with freeing up the memspace 
		}
	}
	return determinant;
}

int main(int argc, char *argv[])
{
	//
	//code given to read .bin file
	//
	int callCount = 0;
	char f_name[50];
	double a[ARRAYSIZE][ARRAYSIZE];
	int i, j;
	double det;
	//Create filename
	sprintf(f_name, "m0016x0016.bin");
	printf("Reading array file %s of size %dx%d\n", f_name, ARRAYSIZE, ARRAYSIZE);
	//Open file
	FILE *datafile = fopen(f_name, "rb");
	//Read elelements
	for (i = 0; i< ARRAYSIZE; i++)
		for (j = 0; j< ARRAYSIZE; j++)
		{
			fread(&a[i][j], sizeof(double), 1, datafile);
			//printf("a[%d][%d]=%f\n", i, j, a[i][j]);
		}
	printf("Matrix has been read.\n");

	//
	//making a dynamically matrix for C, otherwise shit goes DOWN when recursion hits (due to the
	//inital parameter that has to be passed in) in the main method
	//
	double *matrix[ARRAYSIZE];
	for (i = 0; i< ARRAYSIZE; i++)
		matrix[i] = (double *)malloc(ARRAYSIZE * sizeof(double));

	for (int l = 0;l < ARRAYSIZE;l++)
	{
		for (int k = 0;k < ARRAYSIZE;k++)
		{
			matrix[k][l] = a[k][l];
		}
	}

	//
	//calling the cofactor method
	//
	det = cofactorMethod(matrix, callCount);
	printf("The determinant is %f", det);

	return 0;

}


