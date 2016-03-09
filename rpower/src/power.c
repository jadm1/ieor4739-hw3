#include <pthread.h>
#include <unistd.h>
#include "utilities.h"
#include "power.h"

int cheap_rank1perturb(int n, double *scratch, double *matcopy, double *qprime, unsigned int* pseed, double scale);


/** free memory of a bag**/
void PWRfreebag(powerbag **ppbag)
{
	powerbag *pbag = *ppbag;

	if (pbag == NULL) goto BACK;

	PWRfree((void**)&pbag->vector);
	PWRfree((void**)&pbag->qcopy);
	PWRfree((void**)&pbag);

	BACK:
	*ppbag = pbag;
}

/** free an address and set it to NULL to prevent double freeing**/
void PWRfree(void **paddress)
{
	void *address = *paddress;

	if (address == NULL) goto BACK;

	printf("freeing array at %p\n", address);
	free(address);
	address = NULL; /** prevents double freeing **/

	BACK:
	*paddress = address;
}

int PWRallocatebag(int ID, int n, int r, double *covmatrix, powerbag **ppbag, double scale, pthread_mutex_t *psyncmutex, pthread_mutex_t *poutputmutex)
{
	int retcode = 0;
	int i;
	powerbag *pbag = NULL;
	double *double_array = NULL;

	int status, command;
	double *vector, *newvector, *q, *qprime, *qcopy, *scratch, *eigenvalue;

	pbag = (powerbag *)calloc(1, sizeof(powerbag));
	if (pbag == NULL) {
		printf("cannot allocate bag for ID %d\n", ID);
		retcode = NOMEMORY; goto BACK;
	}

	status = PREANYTHING;
	command = STANDBY;

	double_array = calloc(n*r + n*r + n*n + n*n, sizeof(double));
	if (double_array == NULL) {
		retcode = NOMEMORY; goto BACK;
	}

	/** keep variables used for intensive computations close together in memory for more efficiency **/
	vector = &double_array[0];
	newvector = &double_array[n*r];
	qprime = &double_array[n*r + n*r];
	q = &double_array[n*r + n*r + n*n];

	/** ignore any vector and generate at random **/
	for(i = 0; i < r*n; i++) {
		vector[i] = rand()/((double) RAND_MAX);
	}

	/** now, allocate an extra matrix and a vector to use in perturbation **/
	/** should really do it in the power retcode since we are putting them in the bag **/
	qcopy = (double *)calloc(n*n, sizeof(double));
	scratch = (double *)calloc(n, sizeof(double));
	if ((qcopy == NULL) || (scratch == NULL)) {
		retcode = NOMEMORY; goto BACK;
	}
	/** and copy the covariance matrix **/
	for(i = 0; i < n*n; i++)
		qcopy[i] = covmatrix[i];

	eigenvalue = (double*)calloc(n, sizeof(double));

	BACK:
	if (pbag != NULL) {
		/** write bag contents **/
		pbag->ID = ID;
		pbag->n = n;
		pbag->r = r;
		pbag->command = command;
		pbag->status = status;
		pbag->q = q;
		pbag->qprime = qprime;
		pbag->psynchro = psyncmutex;
		pbag->poutputmutex = poutputmutex;
		pbag->qcopy = qcopy;
		pbag->scratch = scratch;
		pbag->scale = scale;
		pbag->eigenvalue = eigenvalue;
		pbag->vector = vector;
		pbag->newvector = newvector;
		pbag->rseed = ID; /** initialize a random seed for the thread using its ID **/
	}
	if (retcode != 0) {
		/** an error occured, cleanup **/
		PWRfreebag(&pbag);
	}
	*ppbag = pbag;

	return retcode;
}

/** Changed this function to only read the matrix and n from the file
 * This function returns the size of the cov matrix in *pn and the matrix (n*n)
 *
 * **/
int PWRreadnload(char *filename, int *pn, double **pmatrix)
{
	int retcode = 0, n, j;
	FILE *input = NULL;
	char buffer[100];
	double *matrix = NULL;

	input = fopen(filename, "r");
	if(!input){
		printf("cannot open file %s\n", filename); retcode = 1; goto BACK;
	}
	fscanf(input,"%s", buffer);
	fscanf(input,"%s", buffer);
	n = atoi(buffer);
	printf("n = %d\n", n);


	matrix = (double*)calloc(n*n, sizeof(double));
	if (matrix == NULL) {
		retcode = NOMEMORY;
		goto BACK;
	}

	fscanf(input, "%s", buffer);
	for(j = 0; j < n*n; j++){
		fscanf(input,"%s", buffer);
		matrix[j] = atof(buffer);
	}

	fclose(input);

	BACK:
	if (retcode == 0) {
		printf("read and loaded data for n = %d with code %d\n", n, retcode);
	}
	else {
		PWRfree((void**)&matrix);
		printf("failed to read/load data\n");
	}
	*pn = n;
	*pmatrix = matrix;
	return retcode;
}



/** Compute a power method iteration **/
void PWRpoweriteration(int ID, int k, 
		int n, int r, double *vector, double *newvector, double *q, double *qprime,
		double *eigenvalue, double *perror,
		pthread_mutex_t *poutputmutex)
{
	double norm2 = 0, mult, error;
	int i, j, f;

	/** The first time copy Q into Q' so that we only deal with Q' afterwards**/
	for (j = 0; j < n*n; j++)
		qprime[j] = q[j];

	/** loop for the first r == 2 eigen values**/
	for (f = 0; f < r; f++) {
		for(i = 0; i < n; i++){
			newvector[f*n + i] = 0;
			for (j = 0; j < n; j++) {
				newvector[f*n + i] += vector[f*n + j]*qprime[i*n + j];
			}
		}

		norm2 = 0;
		for(j = 0; j < n; j++)
			norm2 += newvector[f*n + j]*newvector[f*n + j];

		mult = 1/sqrt(norm2);

		for(j = 0; j < n; j++)
			newvector[f*n + j] = newvector[f*n + j]*mult;

		mult = 1/mult;

		eigenvalue[f] = mult;

		/** Set Q' = Q' - lambda w w^T **/
		for(i = 0; i < n; i++){
			for (j = 0; j < n; j++){
				qprime[i*n + j] -= mult*vector[f*n + i]*vector[f*n + j];
			}
		}
	}

	PWRcompute_error(n, r, &error, newvector, vector);

	if(0 == k%100){
		pthread_mutex_lock(poutputmutex);
		printf("ID %d at iteration %d, norm is %g, ", ID, k, eigenvalue[0]);
		printf("  L1(error) = %.9e\n", error);
		pthread_mutex_unlock(poutputmutex);
	}

	/** will need to map newvector into vector if not terminated **/
	for(j = 0; j < r*n; j++)
		vector[j] = newvector[j];

	*perror = error;
}


void PWRcompute_error(int n, int r, double *perror, double *newvector, double *vector)
{
	int j;
	int f;
	double error;

	error = 0;

	for (f = 0; f < r; f++) {
		for (j = 0; j < n; j++) {
			error += fabs(newvector[f*n + j] - vector[f*n + j]);
		}
	}
	error /= n;
	error /= r;

	*perror = error;

}

/** power method algorithm **/
void PWRpoweralg(powerbag *pbag)
{
	int n, r, ID;
	double *vector, *newvector;
	int k, waitcount, retcode;
	double error, tolerance;
	char letsgo = 0, interrupting, forcedquit = 0;

	ID = pbag->ID;
	n = pbag->n;
	r = pbag->r;

	pthread_mutex_lock(pbag->poutputmutex);
	printf("ID %d starts\n", pbag->ID);
	pthread_mutex_unlock(pbag->poutputmutex);


	vector = pbag->vector;
	newvector = pbag->newvector;

	tolerance = 1e-6;


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
		/** Q is initialized from qcopy at this line **/
		if((retcode = cheap_rank1perturb(n, pbag->scratch, pbag->qcopy, pbag->q, &pbag->rseed, pbag->scale)))
			goto DONE;

		/** initialize vector to random**/
		for(k = 0; k < n; k++){
			vector[k] = rand_r(&pbag->rseed)/((double) RAND_MAX);
		}


		for(k = 0; ; k++){

			/* PWRshowvector(n, vector);*/
			PWRpoweriteration(ID, k, n, r, vector, newvector, pbag->q, pbag->qprime, pbag->eigenvalue, &error, pbag->poutputmutex);
			if(error < tolerance){
				pthread_mutex_lock(pbag->poutputmutex);
				printf(" ID %d converged to tolerance %g! on job %d at iteration %d\n", ID, tolerance, pbag->jobnumber, k);
				/**printf(" ID %d top eigenvalue  %g!\n", ID, pbag->eigenvalue[0]); don't print eigen values twice**/
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

				if (interrupting)
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



/** print vector **/
void PWRshowvector(int n, double *vector)
{
	int j;

	for (j = 0; j < n; j++){
		printf(" %g", vector[j]);
	}
	printf("\n");
}


