#ifndef POWER
#define POWER


#define NOMEMORY 100


#define WAITING 100
#define WORKING 101
#define PREANYTHING 102
#define DONEWITHWORK 103


#define QUIT 200
#define WORK 201
#define STANDBY 202
#define INTERRUPT 203

typedef struct powerbag{
	int n;
	int r; /** number of eigen values and vectors we want to save in the pca (r == 2 for the homework)**/
	double *q; /** perturbed cov matrix (initially it was called q) used at the begining of a power method iteration**/
	double *qprime; /** Q' cov matrix used in the power method **/
	double *qcopy; /** initial covariance Q (initially it was called matcopy) **/
	double *scratch; /** vector used for the rank 1 perturbation **/
	double *eigenvalue; /** Array of eigen values sorted in decreasing order**/
	double *vector; /** Corresponding matrix of eigen vectors (r x n matrix) **/
	double *vector0; /** Corresponding matrix of eigen vectors at iteration 0 (r x n matrix) **/
	double *newvector; /** Matrix of new eigen vectors at the end of an iteration (r x n matrix)**/
	double scale; /** scale parameter for the rank 1 perturb **/
	double tolerance;

	int ID; /** worker thread ID **/
	int status; /** status code **/
	int command; /** command code **/
	int jobnumber;
	int itercount;
	pthread_mutex_t *psynchro; /** mutex pointer for communication with the master thread **/
	pthread_mutex_t *poutputmutex; /** mutex pointer for outputing text to the console**/
	unsigned int rseed; /** thread's random seed
	I used rand_r() inside threads because rand() is not thread safe and every time it is called, it updates
	an internal value so there is a non zero (though very low) risk of multiple core accessing the seed at the same
	time. That is why I prefer to give a specific seed to every thread. Moreover it allows result reproducibility.
	**/
}powerbag;


void PWRshowvector(int n, double *vector);
void PWRfree(void **paddress);
int PWRreadnload(char *filename, int *pn, double **pmatrix);
int PWRallocatebag(int ID, int n, int r, double *covmatrix, powerbag **ppbag, double scale, double tolerance, pthread_mutex_t *psyncmutex, pthread_mutex_t *poutputmutex);
void PWRfreebag(powerbag **ppbag);
void PWRpoweralg(powerbag *pbag);
void PWRpoweriteration(int ID, int k, 
		int n, double *vector, double *newvector, double *q,
		double *peigenvalue, double *perror,
		pthread_mutex_t *poutputmutex);
void PWRcompute_error(int n, double *perror, double *newvector, double *vector);

#endif

