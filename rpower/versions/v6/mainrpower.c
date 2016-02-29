#include <stdio.h>
#include <stdlib.h>
#include "utilities.h"
#include "power.h"

int cheap_rank1perturb(int n, double *scratch, double *matcopy, double *matrix, double scale);

void *PWR_wrapper(void *pvoidedbag);

int main(int argc, char *argv[])
{
  int code = 0, i, j, n;
  powerbag *pbag = NULL;
  double *matcopy = NULL, *scratch = NULL;
  double scale = 1.0;
  int quantity = 1;

  if(argc < 2){ 
    printf(" usage: rpower filename [-s scale] [-q quantity]\n");
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

  /*  if((code = PWRreadnload(argv[1], &n, &vector, &newvector, &matrix))) 
      goto BACK; */

  if((code = PWRreadnload_new(argv[1], 0, &pbag)))
    goto BACK;

  n = pbag->n;

  /** now, allocate an extra matrix and a vector to use in 
      perturbation **/
  matcopy = (double *)calloc(n*n, sizeof(double));
  scratch = (double *)calloc(n, sizeof(double));
  if((!matcopy) || (!scratch)){
    printf("no memory\n"); code = NOMEMORY; goto BACK;
  }
  /** and copy matrix **/
  for(i = 0; i < n*n; i++)
    matcopy[i] = pbag->matrix[i];

  for(j = 0; j < quantity; j++){
    printf("\n************Experiment %d\n", j);
    if((code = cheap_rank1perturb(n, scratch, matcopy, pbag->matrix, scale)))
      goto BACK;


    /*    PWRpoweralg_new(pbag);*/

    PWR_wrapper((void *) pbag);

    printf("top eigenvalue estimate: %.12e\n", pbag->topeigvalue);
  }
      
  PWRfreespace(&pbag);

	 
BACK:
  if(matcopy) free(matcopy);
  if(scratch) free(scratch);
  return code;
}



int cheap_rank1perturb(int n, double *scratch, double *matcopy, double *matrix, double scale)
{
  int retcode = 0, j, i;
  double sum2, invnorm;

  /** first, create a random vector **/
  for(j = 0; j < n; j++)
    scratch[j] = ((double) rand())/((double) RAND_MAX);

  /** next, convert to norm 1 **/
  sum2 = 0;
  for(j = 0; j < n; j++)
    sum2 += scratch[j]*scratch[j];

  invnorm = 1/sqrt(sum2);

  /** rescale **/
  for(j = 0; j < n; j++)
    scratch[j] *= scale*invnorm;


  printf("scale for random perturbation: %g\n", scale);

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      matrix[i*n + j] = scratch[i]*scratch[j] + matcopy[i*n + j];

  return retcode;
}

void *PWR_wrapper(void *pvoidedbag)
{
  powerbag *pbag = (powerbag *) pvoidedbag;

  PWRpoweralg_new(pbag);

  return (void *) &pbag->ID;
}





















