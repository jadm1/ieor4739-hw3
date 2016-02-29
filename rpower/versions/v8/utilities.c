#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h> 
#include "utilities.h"



char does_it_exist(char *filename)
{
   struct stat buf;

/* the function stat returns 0 if the file exists*/

   if (0 == stat(filename, &buf)){
      return 1;
   }
   else return 0;
}

void gotosleep(int numseconds)
{
   sleep(numseconds);
}

void erasefile(char *filename)
{
  remove(filename);
}


double drawnormal(void)
{
  double U1, U2, drawn, pi;

  pi = 3.141592653589793;

  U1 = rand()/((double) RAND_MAX);
  U2 = rand()/((double) RAND_MAX);
   
  drawn = sqrt(-2*log(U1))*cos(2*pi*U2);

  return drawn;
}
