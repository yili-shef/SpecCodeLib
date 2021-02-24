#include<stdio.h>
#include<stdlib.h>
#include"SwapEndian.h"
#include"vtkdata.h"



int vtkscalar_( float *phi, int *nres, char *outfilename, long int fnlen)
{
  int nres3, i;
  nres3 = (*nres) * (*nres) * (*nres);

  /* Swapping the Endian */
  float *ptr;
  for (i = 0; i < nres3; i++)
  {
    ptr = (float *) SwapEndian ( phi+i, 4 ); 
    *(phi + i) = *ptr; 
  }

  /* It is assumed that the string has reserved the last byte for \0. */
  *(outfilename+fnlen-1) = '\0';

  FILE *fs;
  fs = fopen(outfilename, "wb");

  fprintf(fs, "# vtk DataFile Version 3.2\n");
  fprintf(fs, "Scalarfield data\n");
  fprintf(fs, "BINARY\n");
  fprintf(fs, "DATASET STRUCTURED_POINTS\n");
  fprintf(fs, "DIMENSIONS %i %i %i\n", *nres, *nres, *nres);
  fprintf(fs, "ORIGIN %f %f %f\n", 0., 0., 0.);
  fprintf(fs, "SPACING %f %f %f\n", 1., 1., 1.);
  fprintf(fs, "POINT_DATA %i\n", nres3 );
  fprintf(fs, "SCALARS phi float %i\n", 1);
  fprintf(fs, "LOOKUP_TABLE default\n");


  fwrite(phi, sizeof( phi[0] ), nres3, fs);

  fclose(fs);

  return 0;

}
