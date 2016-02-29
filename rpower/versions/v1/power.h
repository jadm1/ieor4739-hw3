#ifndef POWER

#define POWER

void PWRfree(double **pvector);
int PWRallocatespace(int n, double **pvector, double **pnewvector, double **pmatrix);
int PWRreadnload(char *filename, int *pn, double **pvector, double **pnewvector, double **pmatrix);

void PWRpoweriteration(int k, 
		    int n, double *vector, double *newvector, double *matrix,
		    double *perror);

void PWRpoweralg(int n, double *vector, double *newvector, double *matrix);

void PWRshowvector(int n, double *vector);

void PWRcompute_error(int n, double *perror, double *newvector, double *vector);
#endif

/*class powerunit{
 public:
  powerunit(char *inputname);
  ~powerunit();
  int readnload(char *file);
  int allocatespace();
  void iterate();
  void showvector(){ for(int j = 0; j < n; j++) printf("%g ",vector[j]);
  printf("\n");}
  char *getname(){return name;}
 private:
  void releasespace();
  int n;
  double *matrix;
  double *vector;
  double *newvector;
  char *name;
};
*/
