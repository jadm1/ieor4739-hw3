#include <stdio.h>
#include <stdlib.h>
#include "utilities.h"
#include "power.h"

int cheap_rank1perturb(int n, double *scratch, double *matcopy, double *matrix, double scale);

int main(int argc, char *argv[])
{
  int code = 0, i, j, n;
  double *matrix = NULL, *vector = NULL, *newvector = NULL;
  double *matcopy = NULL, *scratch = NULL;
  double scale = 1.0, topeigenvalue;
  int quantity = 1;

  if(argc < 2){ 
    printf(" usage: power filename [-s scale] [-q quantity]\n");
    code = 1; goto BACK;
  }

  for(j = 2; j < argc; j++){
    if (0 == strcmp(argv[j],"-s")){
      j += 1;
      scale = atof(argv[j]);
    }
    else if (0 == strcmp(argv[j],"-q")){
      j += 1;
      quantity = atoi(argv[j]);
    }
    else{
      printf("bad option %s\n", argv[j]); code = 1; goto BACK;
    }
  }


  printf("will use scale %g and quantity %d\n", scale, quantity);

  if((code = PWRreadnload(argv[1], &n, &vector, &newvector, &matrix))) 
    goto BACK;

  /** now, allocate an extra matrix and a vector to use in 
      perturbation **/
  matcopy = (double *)calloc(n*n, sizeof(double));
  scratch = (double *)calloc(n, sizeof(double));
  if((!matcopy) || (!scratch)){
    printf("no memory\n"); code = NOMEMORY; goto BACK;
  }
  /** and copy matrix **/
  for(i = 0; i < n*n; i++)
    matcopy[i] = matrix[i];

  for(j = 0; j < quantity; j++){
    printf("\n************Experiment %d\n", j);
    if((code = cheap_rank1perturb(n, scratch, matcopy, matrix, scale)))
      goto BACK;


    PWRpoweralg(n, vector, newvector, matrix, &topeigenvalue);
    printf("top eigenvalue estimate: %.12e\n", topeigenvalue);
  }

  PWRfree(&vector);

	 
BACK:
  if(matcopy) free(matcopy);
  if(scratch) free(scratch);
  return code;
}



int cheap_rank1perturb(int n, double *scratch, double *matcopy, double *matrix, double scale)
{
  int retcode = 0, j, i;

  /** first, create a random vector **/
  for(j = 0; j < n; j++)
    scratch[j] = scale*((double) rand())/((double) RAND_MAX);

  printf("scale for random perturbation: %g\n", scale);

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      matrix[i*n + j] = scratch[i]*scratch[j] + matcopy[i*n + j];

  return retcode;
}























