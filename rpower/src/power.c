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

	PWRfree((void**)&pbag->vector0);
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

int PWRallocatebag(int ID, int n, int r, double *covmatrix, powerbag **ppbag, double scale, double tolerance, pthread_mutex_t *psyncmutex, pthread_mutex_t *poutputmutex)
{
	int retcode = 0;
	int i;
	powerbag *pbag = NULL;
	double *double_array = NULL;

	int status, command;
	double *vector, *vector0, *newvector, *q, *qprime, *qcopy, *scratch, *eigenvalue;

	pbag = (powerbag *)calloc(1, sizeof(powerbag));
	if (pbag == NULL) {
		printf("cannot allocate bag for ID %d\n", ID);
		retcode = NOMEMORY; goto BACK;
	}

	status = PREANYTHING;
	command = STANDBY;

	double_array = calloc(n*r + n*r + n*r + n*n + n*n, sizeof(double));
	if (double_array == NULL) {
		retcode = NOMEMORY; goto BACK;
	}

	/** keep variables used for intensive computations close together in memory for more efficiency **/
	vector0 = &double_array[0];
	vector = &double_array[n*r];
	newvector = &double_array[n*r + n*r];
	qprime = &double_array[n*r + n*r + n*r];
	q = &double_array[n*r + n*r + n*r + n*n];

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
		pbag->vector0 = vector0;
		pbag->newvector = newvector;
		pbag->rseed = ID; /** initialize a random seed for the thread using its ID **/
		pbag->tolerance = tolerance;
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
		int n, double *vector, double *newvector, double *q,
		double *peigenvalue, double *perror,
		pthread_mutex_t *poutputmutex)
{
	double norm2 = 0, mult, error;
	int i, j;

	/** w_k+1 = Q * w_k **/
	for(i = 0; i < n; i++){
		newvector[i] = 0;
		for (j = 0; j < n; j++) {
			newvector[i] += vector[j]*q[i*n + j];
		}
	}

	norm2 = 0;
	for(j = 0; j < n; j++)
		norm2 += newvector[j]*newvector[j];

	mult = 1.0/sqrt(norm2);

	for(j = 0; j < n; j++)
		newvector[j] = newvector[j]*mult;

	*peigenvalue = 1.0/mult;

	PWRcompute_error(n, &error, newvector, vector);

	if(0 == k%100){
		pthread_mutex_lock(poutputmutex);
		printf("ID %d at iteration %d, norm is %g, ", ID, k, *peigenvalue);
		printf("  L1(error) = %.9e\n", error);
		pthread_mutex_unlock(poutputmutex);
	}

	/** will need to map newvector into vector if not terminated **/
	for(j = 0; j < n; j++)
		vector[j] = newvector[j];

	*perror = error;
}


void PWRcompute_error(int n, double *perror, double *newvector, double *vector)
{
	int j;
	double error;

	error = 0;

	for (j = 0; j < n; j++) {
		error += fabs(newvector[j] - vector[j]);
	}
	error /= n;

	*perror = error;

}

/** power method algorithm **/
void PWRpoweralg(powerbag *pbag)
{
	int n, r, ID;
	int i, j, f;
	double *vector, *vector0, *newvector;
	int k, waitcount, retcode;
	double error, tolerance, sp;
	char letsgo = 0, interrupting, forcedquit = 0;

	ID = pbag->ID;
	n = pbag->n;
	r = pbag->r;
	interrupting = 0;

	pthread_mutex_lock(pbag->poutputmutex);
	printf("ID %d starts\n", pbag->ID);
	pthread_mutex_unlock(pbag->poutputmutex);


	vector = pbag->vector;
	vector0 = pbag->vector0;
	newvector = pbag->newvector;

	tolerance = pbag->tolerance;


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



		/** initialize first vector to random**/
		for(j = 0; j < n*1; j++){
			vector0[j] = rand_r(&pbag->rseed)/((double) RAND_MAX);
		}

		/** copy Q into Q'  so that we only deal with Q' and afterwards**/
		for (j = 0; j < n*n; j++)
			pbag->qprime[j] = pbag->q[j];

		for (f = 0; f < r; f++) {
			/** copy f-th column vector0 into vector **/
			for(j = 0; j < n; j++){
				vector[f*n + j] = vector0[f*n + j];
			}
			for(k = 0; ; k++) {

				/* PWRshowvector(n, vector);*/
				PWRpoweriteration(ID, k, n, &vector[f*n], &newvector[f*n], pbag->qprime, &pbag->eigenvalue[f], &error, pbag->poutputmutex);
				if(error < tolerance){
					/** finished to compute f-th eigen value **/


					/** Set Q' = Q' - lambda w w^T **/
					for(i = 0; i < n; i++){
						for (j = 0; j < n; j++){
							pbag->qprime[i*n + j] -= pbag->eigenvalue[f]*vector[f*n + i]*vector[f*n + j];
						}
					}

					/** Set w'_0 = w_0 - (w^T w_0) w **/
					if (f < r-1) {
						/** first compute sp = (w^T w_0)**/
						sp = 0.0;
						for (j = 0; j < n; j++) {
							sp += vector[f*n + j] * vector0[f*n + j];
						}
						/** Set w'_0 = w_0 - sp * w **/
						for (j = 0; j < n; j++) {
							vector0[(f+1)*n + j] = vector0[f*n + j] - sp * vector[f*n + j];
						}
					}


					pthread_mutex_lock(pbag->poutputmutex);
					printf(" ID %d converged to tolerance %g! on job %d at iteration %d\n", ID, tolerance, pbag->jobnumber, k);
					printf(" ID %d %d-th eigenvalue:  %g!\n", ID, f, pbag->eigenvalue[f]);
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
			if (interrupting)
				break; /** takes you outside of for loop **/
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


