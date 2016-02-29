#include "utilities.h"
#include "power.h"



int PWRallocatespace(int n, double **pv, double **pnewvector, double **pmatrix)
{
  int retcode = 0;
  double *v = NULL;

  v = (double *) calloc(n + n + n*n, sizeof(double));
  if(!v){
    retcode = NOMEMORY; goto BACK;
  }
  printf("allocated vector at %p\n", (void *) v);
  *pv = v;
  *pnewvector = v + n;
  *pmatrix = v + 2*n;

 BACK:
  return retcode;
}

void PWRfreespace(powerbag **ppbag)
{
  powerbag *pbag = *ppbag;
  if(pbag == NULL) return;

  PWRfree(&pbag->vector);

  free(pbag);
  *ppbag = NULL;
}


void PWRfree(double **pvector)
{
  double *vector = *pvector;
  if(vector == NULL) return;
  printf("freeing vector at %p\n", (void *) vector);
  free(vector);
  *pvector = NULL; /** prevents double freeing **/
}

int PWRreadnload_new(char *filename, int ID, powerbag **ppbag)
{
  int code = 0, n;
  double *vector, *newvector, *matrix;
  powerbag *pbag = NULL;

  code = PWRreadnload(filename, &n, &vector, &newvector,  &matrix);
  if(code){
    goto BACK;
  }

 
  pbag = (powerbag *)calloc(1, sizeof(powerbag));
  if(!pbag){
    printf("cannot allocate bag for ID %d\n", ID); code = NOMEMORY; goto BACK;
  }
  *ppbag = pbag;

  pbag->n = n;
  pbag->ID = ID;
  pbag->vector = vector;
  pbag->newvector = newvector;
  pbag->matrix = matrix;


 BACK:
  return code;
}

int PWRreadnload(char *filename, int *pn, double **pvector, double **pnewvector, double **pmatrix)
{
  int retcode = 0, n, j;
  FILE *input = NULL;
  char buffer[100];

  input = fopen(filename, "r");
  if(!input){
    printf("cannot open file %s\n", filename); retcode = 1; goto BACK;
  }
  fscanf(input,"%s", buffer);
  fscanf(input,"%s", buffer);
  n = atoi(buffer);
  *pn = n;
  
  printf("n = %d\n", n);
  retcode = PWRallocatespace(n, pvector, pnewvector, pmatrix);
  if(retcode) goto BACK;


  fscanf(input,"%s", buffer);
  for(j = 0; j < n*n; j++){ 
    fscanf(input,"%s", buffer);
    (*pmatrix)[j] = atof(buffer);
  }

  /** ignore any vector and generate at random **/  
  for(j = 0; j < n; j++){ 
    (*pvector)[j] = rand()/((double) RAND_MAX);
  }

  fclose(input);

 BACK:
  printf("read and loaded data for n = %d with code %d\n", n, retcode);
  return retcode;
}




void PWRpoweriteration(int ID, int k, 
		    int n, double *vector, double *newvector, double *matrix,
		       double *peigvalue, double *perror)
{
  double norm2 = 0, mult, error;
  int i, j;
  for(i = 0; i <n; i++){
    newvector[i] = 0;
    for (j = 0; j < n; j++){
      newvector[i] += vector[j]*matrix[i*n + j];
    }
  }

  norm2 = 0;
  for(j = 0; j < n; j++) norm2 += newvector[j]*newvector[j];

  mult = 1/sqrt(norm2);

  for(j = 0; j < n; j++) newvector[j] = newvector[j]*mult;

  PWRcompute_error(n, &error, newvector, vector);

  if(0 == k%10){
    printf("ID %d at iteration %d, norm is %g, ", ID, k, 1/mult);
    printf("  L1(error) = %.9e\n", error);
  }

  *peigvalue = 1/mult;
  *perror = error;

  /** will need to map newvector into vector if not terminated **/

  for(j = 0; j < n; j++) vector[j] = newvector[j];

}

void PWRcompute_error(int n, double *perror, double *newvector, double *vector)
{
  int j;
  double error;

  error = 0;

  for(j = 0; j < n; j++){
    error += fabs(newvector[j] - vector[j]);
  }

  *perror = error;

}

void PWRpoweralg_new(powerbag *pbag)
{
  int n, ID;
  double *vector, *newvector, *matrix;
  int k;
  double error, tolerance;

  ID = pbag->ID;
  n = pbag->n;
  vector = pbag->vector;
  newvector = pbag->newvector;
  matrix = pbag->matrix;
  tolerance = .000001/n;
  if (tolerance < 1e-7) tolerance = 1e-7;

  for(k = 0; ; k++){
    PWRpoweriteration(ID, k, n, vector, newvector, matrix, &pbag->topeigvalue, &error);
    if(error < tolerance){
      printf(" ID %d converged to tolerance %g!\n", ID, tolerance); 
      printf(" ID %d top eigenvalue  %g!\n", ID, pbag->topeigvalue);

      break;
    }
  }

}

void PWRshowvector(int n, double *vector)
{
  int j;

  for (j = 0; j < n; j++){
    printf(" %g", vector[j]);
  }
  printf("\n");
}


