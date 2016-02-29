#include <stdio.h>
#include <stdlib.h>
#include "utilities.h"
#include "power.h"

int cheap_rank1perturb(int n, double *scratch, double *matcopy, double *matrix, double scale);

int main(int argc, char *argv[])
{
  int code = 0, i, n;
  double *matrix = NULL, *vector = NULL, *newvector = NULL;
  double *matcopy = NULL, *scratch = NULL;
  double scale = 1.0;

  if(argc > 3){ 
    printf(" usage: power filename [scale]\n");
    code = 1; goto BACK;
  }
  if(argc == 3){
    scale = atof(argv[2]);
  }
  printf("will use scale %g\n", scale);

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

  if((code = cheap_rank1perturb(n, scratch, matcopy, matrix, scale)))
    goto BACK;


  PWRpoweralg(n, vector, newvector, matrix);

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























