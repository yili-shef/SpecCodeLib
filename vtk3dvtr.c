#include<stdio.h>
#include<stdlib.h>
#include"SwapEndian.h"
#include"vtkdata.h"


/* 
 * Swap the Edian of the numbers and pack them into velvtr
 */
int swapnpack(float *cmpnt, int nres3, int offset, float *velvtr)
{
    int i;
    float *ptr;

    for (i = 0; i < nres3; i++)
    { 
      ptr = (float *) SwapEndian ( cmpnt+i , 4 ); 
      *(velvtr + offset) = *ptr; 
      offset += 3 ;
    }
    return 0;
}


int vtk3dvtr_( float *ux, float *uy, float*uz, int *nres, char *outfilename, long int fnlen)
{
  int nres3;
  nres3 = (*nres) * (*nres) * (*nres);
  
  float *vtr; // To store the vector

  // Note the correct usage of malloc. Don't forget free.
  vtr = (float *) malloc( sizeof(float) * nres3 * 3 );

  swapnpack(ux, nres3, 0, vtr);
  swapnpack(uy, nres3, 1, vtr);
  swapnpack(uz, nres3, 2, vtr);

  // It is assumed that the string has reserved the last byte for \0.
  *(outfilename+fnlen-1) = '\0';

  FILE *fs;
  fs = fopen(outfilename, "wb");

  fprintf(fs, "# vtk DataFile Version 3.2\n");
  fprintf(fs, "3d vector field data\n");
  fprintf(fs, "BINARY\n");
  fprintf(fs, "DATASET STRUCTURED_POINTS\n");
  fprintf(fs, "DIMENSIONS %i %i %i\n", *nres, *nres, *nres);
  fprintf(fs, "ORIGIN %f %f %f\n", 0., 0., 0.);
  fprintf(fs, "SPACING %f %f %f\n", 1., 1., 1.);
  fprintf(fs, "POINT_DATA %i\n", nres3 );
  fprintf(fs, "VECTORS 3dvector float\n");
  /* for vector field no lookup table is needed */

  fwrite(vtr, sizeof( vtr[0] ), nres3*3, fs);

  fclose(fs);

  free(vtr);

  return 0;

}
