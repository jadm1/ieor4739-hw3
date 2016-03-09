#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>
#include "utilities.h"
#include "power.h"

static powerbag **ppbagproxy = NULL;
static int numworkersproxy = 0;
static char *deadstatus = NULL;
static int activeworkers = 0;

int cheap_rank1perturb(int n, double *scratch, double *matcopy, double *matrix, unsigned int* pseed, double scale);

void *PWR_wrapper(void *pvoidedbag);
void (*sigset(int sig, void (*disp)(int)))(int);

void handlesigint(int i);


int main(int argc, char *argv[])
{
	int retcode = 0, j, n, initialruns, scheduledjobs;
	powerbag **ppbag = NULL, *pbag;
	double scale = 1.0;
	int quantity = 1, numworkers = 1, theworker;
	char gotone;
	double *covmatrix;
	int r;
	pthread_t *pthread;
	pthread_mutex_t outputmutex;
	pthread_mutex_t *psyncmutex;
	/**unsigned int rseed = 123;**/

	r = 2; /** default number of factors **/

	if(argc < 2){
		printf(" usage: rpower filename [-s scale] [-q quantity] [-w workers]\n");
		retcode = 1; goto BACK;
	}

	sigset(SIGINT, &handlesigint);



	for(j = 2; j < argc; j++){
		if (0 == strcmp(argv[j],"-s")){
			j += 1;
			scale = atof(argv[j]);
		}
		else if (0 == strcmp(argv[j],"-q")){
			j += 1;
			quantity = atoi(argv[j]);
		}
		else if (0 == strcmp(argv[j],"-w")){
			j += 1;
			numworkers = atoi(argv[j]);
		}
		else if (0 == strcmp(argv[j], "-r")){
			j += 1;
			r = atoi(argv[j]); /** number of eigen values we want to extract from the PCA on the cov matrix **/
		}
		else{
			printf("bad option %s\n", argv[j]); retcode = 1; goto BACK;
		}
	}

	printf("will use scale %g and quantity %d: %d workers\n", scale, quantity, numworkers);

	if ( numworkers > quantity ){
		numworkers = quantity;
		printf(" --> reset workers to %d\n", numworkers);
	}

	deadstatus = (char *) calloc(numworkers, sizeof(char));

	pthread_mutex_init(&outputmutex, NULL); /** common to everybody **/


	psyncmutex = (pthread_mutex_t *)calloc(numworkers, sizeof(pthread_mutex_t));
	if(!psyncmutex){
		printf("could not create mutex array\n"); retcode = NOMEMORY; goto BACK;
	}

	for(j = 0; j < numworkers; j++)
		pthread_mutex_init(&psyncmutex[j], NULL);

	ppbag = (powerbag **)calloc(numworkers, sizeof(powerbag *));
	if(!ppbag){
		printf("could not create bag array\n"); retcode = NOMEMORY; goto BACK;
	}
	numworkersproxy = numworkers;
	ppbagproxy = ppbag;

	pthread = (pthread_t *)calloc(numworkers, sizeof(pthread_t));
	if(!pthread){
		printf("could not create thread array\n"); retcode = NOMEMORY; goto BACK;
	}


	retcode = PWRreadnload(argv[1], &n, &covmatrix); /** read the data once **/
	if (retcode != 0)
		goto BACK;

	for(j = 0; j < numworkers; j++) {

		if((retcode = PWRallocatebag(j, n, r, covmatrix, &ppbag[j], scale, &psyncmutex[j], &outputmutex)))
			goto BACK;

		printf("about to launch thread for worker %d\n", j);

		pthread_create(&pthread[j], NULL, &PWR_wrapper, (void *) ppbag[j]);
	}

	initialruns = numworkers;
	if (initialruns > quantity) initialruns = quantity;

	for(theworker = 0; theworker < initialruns; theworker++){
		pbag = ppbag[theworker];
		/**retcode = cheap_rank1perturb(n, pbag->scratch, pbag->matcopy, pbag->matrix, &rseed, scale);
		if (retcode != 0)
			goto BACK;**/

		pthread_mutex_lock(&outputmutex);
		printf("*****master:  worker %d will run experiment %d\n", theworker, j);
		pthread_mutex_unlock(&outputmutex);

		/** tell the worker to work **/
		pthread_mutex_lock(&psyncmutex[theworker]);
		pbag->command = WORK;
		pbag->status = WORKING;
		pbag->jobnumber = theworker;
		pbag->itercount = 0;
		pthread_mutex_unlock(&psyncmutex[theworker]);

	}
	scheduledjobs = activeworkers = initialruns;

	while(activeworkers > 0) {
		/** check the workers' status **/
		gotone = 0;
		for(theworker = 0; theworker < numworkers; theworker++){

			pthread_mutex_lock(&psyncmutex[theworker]);
			pbag = ppbag[theworker];
			if(pbag->status == DONEWITHWORK){

				pthread_mutex_lock(&outputmutex);
				printf("master:  worker %d is done with job %d\n", pbag->ID, pbag->jobnumber);
				for (j = 0; j < r; j++) {
					printf("Job %d: Eigenvalue #%d estimate: %.12e\n", pbag->jobnumber, j+1, pbag->eigenvalue[j]);
				}
				/**for (j = 0; j < r; j++) {
					printf("Eigenvector #%d: ", j+1);
					PWRshowvector(n, &pbag->eigenvector[j*n]);
				} don't print the eigen vectors they take much room**/
				pthread_mutex_unlock(&outputmutex);

				if(scheduledjobs >= quantity){
					/** tell worker to quit **/
					pthread_mutex_lock(&outputmutex);
					printf("master: telling worker %d to quit\n", pbag->ID);
					pthread_mutex_unlock(&outputmutex);
					pbag->command = QUIT;
					pbag->status = QUIT;
					--activeworkers;
				}
				else {
					gotone = 1;
				}
			}
			else if(pbag->status == PREANYTHING) {
				pthread_mutex_lock(&outputmutex);
				printf("master:  worker %d is available\n", theworker);
				pthread_mutex_unlock(&outputmutex);
				gotone = 1;
			}
			else if( (pbag->status == WORKING) && (pbag->itercount > 100000)){
				pbag->command = INTERRUPT;
				pthread_mutex_lock(&outputmutex);
				printf("master: telling worker %d to interrupt\n", pbag->ID);
				pthread_mutex_unlock(&outputmutex);
			}
			else if(deadstatus[theworker]){
				pthread_mutex_lock(&outputmutex);
				printf("master: telling worker %d to quit\n", pbag->ID);
				pthread_mutex_unlock(&outputmutex);
				pbag->command = QUIT;
				pbag->status = QUIT;
				--activeworkers;
				/** and let's make sure we don't do it again **/
				deadstatus[theworker] = 0;
			}
			pthread_mutex_unlock(&psyncmutex[theworker]);
			if(gotone) break;
			usleep(10000);

		}
		/** at this point we have run through all workers **/

		if(gotone){
			/** if we are here, "theworker" can work **/
			pbag = ppbag[theworker];
			/**retcode = cheap_rank1perturb(n, pbag->scratch, pbag->qcopy, pbag->q, &rseed, scale);
			if (retcode != 0)
				goto BACK;**/

			pthread_mutex_lock(&outputmutex);
			printf("master:  worker %d will run experiment %d\n", theworker, scheduledjobs);
			pthread_mutex_unlock(&outputmutex);


			/** tell the worker to work **/
			pthread_mutex_lock(&psyncmutex[theworker]);
			pbag->command = WORK;
			pbag->status = WORKING;
			pbag->itercount = 0;
			pbag->jobnumber = scheduledjobs;
			pthread_mutex_unlock(&psyncmutex[theworker]);

			++scheduledjobs;
		}
	}



	/*  pthread_mutex_lock(&psynchro_array[theworker]);
  pbag->command = QUIT;
  pthread_mutex_unlock(&psynchro_array[theworker]);*/

	pthread_mutex_lock(&outputmutex);
	printf("master:  done with loop\n");
	pthread_mutex_unlock(&outputmutex);

	/** actually this is bad -- should wait for the threads to be done --
      but how **/
	for(j = 0; j < numworkers; j++){
		pthread_join(pthread[j], NULL);
		pthread_mutex_lock(&outputmutex);
		printf("master: joined with thread %d\n", j);
		pthread_mutex_unlock(&outputmutex);
		pbag = ppbag[j];
		PWRfreebag(&pbag);
	}
	free(ppbag);

	BACK:
	if (covmatrix != NULL) {
		free(covmatrix); covmatrix = NULL;
	}
	return retcode;
}



int cheap_rank1perturb(int n, double *scratch, double *matcopy, double *matrix, unsigned int* pseed, double scale)
{
	int retcode = 0, j, i;
	double sum2, invnorm;

	/** first, create a random vector **/
	for(j = 0; j < n; j++)
		scratch[j] = ((double) rand_r(pseed))/((double) RAND_MAX);

	/** next, convert to norm 1 **/
	sum2 = 0;
	for(j = 0; j < n; j++)
		sum2 += scratch[j]*scratch[j];

	invnorm = 1/sqrt(sum2);

	/** rescale **/
	invnorm *= scale;
	for(j = 0; j < n; j++)
		scratch[j] *= invnorm;


	printf("scale for random perturbation: %g\n", scale);

	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			matrix[i*n + j] = scratch[i]*scratch[j] + matcopy[i*n + j];

	return retcode;
}

void *PWR_wrapper(void *pvoidedbag)
{
	powerbag *pbag = (powerbag *) pvoidedbag;

	PWRpoweralg(pbag);

	return (void *) &pbag->ID;
}



void handlesigint(int signal)
{
	int j;
	printf("yo, what's happening\n");
	for(j = 0; j < numworkersproxy; j++){
		deadstatus[j] = 1;
	}
	/** brutal **/
}

















