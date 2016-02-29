#include <pthread.h>
#include <unistd.h>
#include "utilities.h"
#include "power.h"

int cheap_rank1perturb(int n, double *scratch, double *matcopy, double *matrix, double scale);

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
  PWRfree(&pbag->matcopy);

  free(pbag);
  *ppbag = NULL;
}


void PWRfree(double **ppaddress)
{
  double *paddress = *ppaddress;
  if(paddress == NULL) return;
  printf("freeing double array at %p\n", (void *) paddress);
  free(paddress);
  *ppaddress = NULL; /** prevents double freeing **/
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
  pbag->status = WAITING;

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
		       double *peigvalue, double *perror,
		         pthread_mutex_t *poutputmutex)
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

  if(0 == k%1000){
    pthread_mutex_lock(poutputmutex);
    printf("ID %d at iteration %d, norm is %g, ", ID, k, 1/mult);
    printf("  L1(error) = %.9e\n", error);
    pthread_mutex_unlock(poutputmutex);
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
  int k, waitcount, retcode;
  double error, tolerance;
  char letsgo = 0, interrupting, forcedquit = 0;

  ID = pbag->ID;
  n = pbag->n;

  pthread_mutex_lock(pbag->poutputmutex);
  printf("ID %d starts\n", pbag->ID);
  pthread_mutex_unlock(pbag->poutputmutex);


  vector = pbag->vector;
  newvector = pbag->newvector;
  matrix = pbag->matrix;
  tolerance = 1e-20; /** excessive **/
  /*  if (tolerance < 1e-7) tolerance = 1e-7;*/

  for(;;){
    pthread_mutex_lock(pbag->poutputmutex);
    printf(" ID %d in big loop\n", pbag->ID);
    pthread_mutex_unlock(pbag->poutputmutex);

    letsgo = 0;
    waitcount = 0;
    while(letsgo == 0){
      /** wait until WORK signal **/
      usleep(10000);

      pthread_mutex_lock(pbag->psynchro);
      if(pbag->command == WORK){
	letsgo = 1;
      }
      else if(pbag->command == QUIT)
	letsgo = 2;
      pthread_mutex_unlock(pbag->psynchro);

      if (letsgo == 2) 
	goto DONE;

      if(0 == waitcount%20){
	pthread_mutex_lock(pbag->poutputmutex);
	printf("ID %d bag %p: wait %d for signal; right now have %d\n", pbag->ID, (void *) pbag, waitcount, pbag->command);
	pthread_mutex_unlock(pbag->poutputmutex);

      }
      ++waitcount;

    }

    pthread_mutex_lock(pbag->poutputmutex);
    printf("ID %d: got signal to start working\n", pbag->ID);
    pthread_mutex_unlock(pbag->poutputmutex);

    /** let's do the perturbation here **/
    if((retcode = cheap_rank1perturb(n, pbag->scratch, pbag->matcopy, pbag->matrix, 10.0)))
      goto DONE;

    /** initialize vector to random**/
    for(k = 0; k < n; k++){ 
      vector[k] = rand()/((double) RAND_MAX);
    }
    

    for(k = 0; ; k++){

      /*      PWRshowvector(n, vector);*/
      PWRpoweriteration(ID, k, n, vector, newvector, matrix, &pbag->topeigvalue, &error, pbag->poutputmutex);
      if(error < tolerance){
	pthread_mutex_lock(pbag->poutputmutex);
	printf(" ID %d converged to tolerance %g! on job %d at iteration %d\n", ID, tolerance,
	       pbag->jobnumber, k); 
	printf(" ID %d top eigenvalue  %g!\n", ID, pbag->topeigvalue);
	pthread_mutex_unlock(pbag->poutputmutex);

	break;
      }
      pbag->itercount = k;  /** well, in this case we don't really need k **/
      if(0 == k%1000){
	pthread_mutex_lock(pbag->psynchro);

	interrupting = 0;
	if(pbag->command == INTERRUPT || pbag->command == QUIT){
	  pthread_mutex_lock(pbag->poutputmutex);
	  printf(" ID %d interrupting after %d iterations\n", pbag->ID, k);
	  pthread_mutex_unlock(pbag->poutputmutex);

	  interrupting = 1;
	}

	pthread_mutex_unlock(pbag->psynchro);

	if(interrupting)
	  break; /** takes you outside of for loop **/
      }
    }

    /** first, let's check if we have been told to quit **/
    pthread_mutex_lock(pbag->psynchro);
    if(pbag->command == QUIT)
      forcedquit = 1;
    pthread_mutex_unlock(pbag->psynchro);

    if(forcedquit)
      break;

    pthread_mutex_lock(pbag->psynchro);
    pbag->status = DONEWITHWORK;
    pbag->command = STANDBY;
    pthread_mutex_unlock(pbag->psynchro);
  }

 DONE:
  pthread_mutex_lock(pbag->poutputmutex);
  printf(" ID %d quitting\n", pbag->ID);
  pthread_mutex_unlock(pbag->poutputmutex);
  
}




void PWRshowvector(int n, double *vector)
{
  int j;

  for (j = 0; j < n; j++){
    printf(" %g", vector[j]);
  }
  printf("\n");
}


