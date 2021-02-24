#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "c3drealdatawriter.h"

int c3drealdatawriter_(float *pdata, int *nx, int *ny, int *nz, char *filename, long int fnlen)
{
  char *prefix;

  prefix = malloc (6 + fnlen);
  strcpy (prefix, "cbin-");
  strncat (prefix, filename, fnlen);
  strcat (prefix, "\0");

  FILE *fp;
  fp = fopen(prefix, "wb");

  int npt;
  npt = (*nx) * (*ny) * (*nz); 

  if ( fwrite(pdata, sizeof(pdata[0]), npt, fp) != npt ) 
  {
    printf("Error when writing data in c3dbindatawriter\n");
    fclose(fp);
    return 1;
  }

  fclose(fp);

  free(prefix);
  return 0;
}
