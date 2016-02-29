#include <stdio.h>
#include <stdlib.h>
#include "utilities.h"
#include "power.h"



int main(int argc, char *argv[])
{
  int code = 0, n;
  double *matrix = NULL, *vector = NULL, *newvector = NULL;

  if(argc != 2){ 
    printf(" usage: power filename\n");
    code = 1; goto BACK;
  }

  if((code = readnload(argv[1], &n, &vector, &newvector, &matrix))) goto BACK;

  poweralg(n, vector, newvector, matrix);

	 
BACK:
  return code;
}




























/*
	powerunit *punit = NULL;

	punit = new powerunit("two");
	punit->readnload(argv[1]);
	punit->iterate();

	delete punit;  */
 
