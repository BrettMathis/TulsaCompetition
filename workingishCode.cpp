//Brett maed dis

#include <stdlib.h>
#include <stdio.h>

//max size allows with default stack limit in linux enviornment, hard limit needs to be drastically increased to allow for 5000x5000
//i.e. 'ulimit -s 16348' will double the stack size from a default user setting of 8192kB
#define ARRAYSIZE 1024

double rowReduce(double inputMatrix[ARRAYSIZE][ARRAYSIZE])
{
  //applies row reduction to all column's

  for(int column = 0;column < ARRAYSIZE;column++)
  {

      //////////////////////////////////////////////
      //////////////////////////////////////////////
    //row reduces a single column

      //sets first column entry of each row to 1.0
      for(int row = 0;row < ARRAYSIZE;row++)
      {
        double normalizer = 1.0;
        double diagonal = inputMatrix[row][column];
        double multiplyThrough = normalizer/diagonal;
        //inputMatrix[row][column] = 1.0;
        //multiplies through the row
        for(int j = 0;j < ARRAYSIZE;j++)
        {
            inputMatrix[row][j] = multiplyThrough * inputMatrix[row][j];
        }
      }

      //subtracts the first row from all other rows,
      // leavning the first column as zero
      for(int i = 0;i < ARRAYSIZE;i++)
      {
        for(int j = 0;j < ARRAYSIZE;j++)
        {
          if(i > column)
          {
            inputMatrix[i][j] = (inputMatrix[i][j] - inputMatrix[column][j]);
          }
        }
      }
      //////////////////////////////////////////////
      //////////////////////////////////////////////

      //printing shit, makes it run slower than my motivation

      for (int i=0; i< ARRAYSIZE; i++)
          for (int j=0; j< ARRAYSIZE; j++)
          {
            //if(i == j)
             //printf("inputMatrix[%d][%d]=%f\n",i,j,inputMatrix[i][j]);
          }
  }


//was originally in computeDet
  double topLeft = inputMatrix[0][0];

  for(int i = 0;i < ARRAYSIZE;i++)
  {
    if(i < (ARRAYSIZE-1))
    {
      printf("%f\n",topLeft);
      topLeft = (topLeft * inputMatrix[i+1][i+1]);
    }
  }

  topLeft = (topLeft / inputMatrix[0][0]);

  return topLeft;
}

double computeDet(double inputMatrix[ARRAYSIZE][ARRAYSIZE])
{
  double workingDet = rowReduce(inputMatrix);

//  printf("Im in the next method%f",workingDet);

  return workingDet;
}


int main(int argc, char *argv[])
{
    char f_name[50];
    double a[ARRAYSIZE][ARRAYSIZE];
    int i,j;
    double det;
    //Create filename
    sprintf(f_name,"m1024x1024.bin");
    printf("Reading array file %s of size %dx%d\n",f_name,ARRAYSIZE,ARRAYSIZE);
    //Open file
    FILE *datafile=fopen(f_name,"rb");
    //Read elelements
    for (i=0; i< ARRAYSIZE; i++)
        for (j=0; j< ARRAYSIZE; j++)
        {
            fread(&a[i][j],sizeof(double),1,datafile);
          //  printf("a[%d][%d]=%f\n",i,j,a[i][j]);
        }
        printf("Matrix has been read.\n");

    det = computeDet(a);
    printf("\nYou have %f determination.\n",det);

  return 0;
}
