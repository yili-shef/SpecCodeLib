#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"SwapEndian.h"
#include"vtkdata.h"

#define mySINGLE 4
#define myDOUBLE 8

static FILE *fs;

/* Open vtk outout data file */
int vtkopenfile2d_(char *outfn, long int fnlen)
{
    *(outfn+fnlen-1) = '\0';
    fs = fopen(outfn, "wb");

    fprintf(fs, "# vtk DataFile Version 3.2\n");
    fprintf(fs, "2d data\n");
    fprintf(fs, "BINARY\n");

    return 0;
}

/* Write the header for Structured Points data files */
int vtk_sp_header2d_(int *nbytes, int *nx, int *ny, 
        void *ox, void *oy, void *spx, void *spy)
{
 
    fprintf(fs, "DATASET STRUCTURED_POINTS\n");
    fprintf(fs, "DIMENSIONS %i %i\n", *nx, *ny);
    if ( *nbytes == mySINGLE )
    { 
        fprintf(fs, "ORIGIN %f %f\n", *(float *)ox, *(float *)oy);
        fprintf(fs, "SPACING %f %f\n", *(float *)spx, *(float *)spy);
    }
    else if ( *nbytes == myDOUBLE )
    { 
        fprintf(fs, "ORIGIN %f %f\n", *(double *)ox, *(double *)oy);
        fprintf(fs, "SPACING %f %f\n", *(double *)spx, *(double *)spy);
    }
    else
    {
        printf("Wrong data type. Stopping.\n");
        exit(EXIT_FAILURE);
    }
    fprintf(fs, "POINT_DATA %i\n", (*nx) * (*ny) );
 
    return 0;
}

/* Write the header for Rectilinear Grid data files */
int vtk_rg_header2d_(void *xco, void *yco, int *nbytes, int *nx, int *ny)
{
    fprintf(fs, "DATASET RECTILINEAR_GRID\n");
    fprintf(fs, "DIMENSIONS %i %i\n", *nx, *ny);

    char *type;

    if ( *nbytes == mySINGLE )
    {
        type = "float";
    }
    else if ( *nbytes == myDOUBLE )
    {
        type = "double";
    }
    else
    {
        printf("Wrong data type. Stopping.\n");
        exit(EXIT_FAILURE);
    }
        

    char *ptrco, *pco, *ptmp;
    int i, j, nmax;

    i = *nx > *ny ? *nx : *ny;
    nmax = i > *nz ? i : *nz;

    ptrco = (char *) malloc( sizeof(char) * nmax * (*nbytes) );

    pco = ptrco;
    for (i = 0; i < *nx; i++)
    {
        ptmp = (char *) SwapEndian ( (char *)xco + i * (*nbytes), *nbytes ); 
        for (j=0; j < *nbytes; *(pco++) = *(ptmp + j), j++)
            ;
    }
    fprintf(fs, "\nX_COORDINATES %i %s\n", *nx, type);
    fwrite(ptrco, (size_t) *nbytes, *nx, fs); 

 
    pco = ptrco;
    for (i = 0; i < *ny; i++)
    {
        ptmp = (char *) SwapEndian ( (char *)yco + i * (*nbytes), *nbytes ); 
        for (j=0; j < *nbytes; *(pco++) = *(ptmp + j), j++)
            ;
    }
    fprintf(fs, "\nY_COORDINATES %i %s\n", *ny, type);
    fwrite(ptrco, (size_t) *nbytes, *ny, fs); 


    pco = ptrco;
    for (i = 0; i < *nz; i++)
    {
        ptmp = (char *) SwapEndian ( (char *)zco + i * (*nbytes), *nbytes ); 
        for (j=0; j < *nbytes; *(pco++) = *(ptmp + j), j++)
            ;
    }
    fprintf(fs, "\nZ_COORDINATES %i %s\n", *nz, type);
    fwrite(ptrco, (size_t) *nbytes, *nz, fs); 

    free(ptrco);

    fprintf(fs, "\nPOINT_DATA %i\n", (*nx)*(*ny)*(*nz) );
    return 0;
}

/* output vector data */
int vtk3dvtr_( void *ux, void *uy, void *uz, int *nbytes, 
        int *nx, int *ny, int *nz, char *name, long int namelen)
{
    int ntotal = (*nx) * (*ny) * (*nz);

    char *pv, *pvtr, *vtr;

    vtr = (char *) malloc( sizeof(char) * ntotal * 3 * (*nbytes) );
    pvtr = vtr;

    int i, j;
    for (i = 0; i < ntotal; i++)
    {
        pv = (char *) SwapEndian ( (char *) ux + i * (*nbytes), *nbytes ); 
        for (j = 0; j < *nbytes; *(pvtr++) = *(pv+j), j++) 
            ;

        pv = (char *) SwapEndian ( (char *) uy + i * (*nbytes), *nbytes ); 
        for (j = 0; j < *nbytes; *(pvtr++) = *(pv+j), j++) 
            ;

        pv = (char *) SwapEndian ( (char *) uz + i * (*nbytes), *nbytes ); 
        for (j = 0; j < *nbytes; *(pvtr++) = *(pv+j), j++) 
            ;
    }

    *(name + namelen-1) = '\0';
 
    if ( *nbytes == mySINGLE )
        fprintf(fs, "\nVECTORS %s float\n", name); 
    else if ( *nbytes == myDOUBLE )
        fprintf(fs, "\nVECTORS %s double\n", name);
    else
    {
        printf("Wrong data type. Stopping.\n");
        exit(1);
    }
 
    fwrite(vtr, (size_t) *nbytes, ntotal*3, fs);
 
    free(vtr);
 
    return 0;
}

/* Output scalar data. */
int vtkscalar_( void *phi, int *nbytes, int *nx, int *ny, int *nz, 
        char *name, long int namelen)
{
    int ntotal;
    ntotal = (*nx) * (*ny) * (*nz);
 
    char *ptr, *pv1, *pv2;
    ptr = (char *) malloc( sizeof(char) * ntotal * (*nbytes) );
    pv2 = ptr;

    // Swapping the Endian 
    int i, j;
    for (i = 0; i < ntotal; i++)
    {
        pv1 = (char *) SwapEndian ( (char *)phi + i * (*nbytes), *nbytes ); 
        for (j=0; j < *nbytes; *(pv2++) = *(pv1 + j), j++)
            ;
    }
 
    *(name + namelen-1) = '\0';
 
    if ( *nbytes == mySINGLE )
        fprintf(fs, "\nSCALARS %s float %i\n", name, 1);
    else if ( *nbytes == myDOUBLE )
        fprintf(fs, "\nSCALARS %s double %i\n", name, 1);
    else
    {
        printf("Wrong data type. Stopping.\n");
        exit(1);
    }

    fprintf(fs, "LOOKUP_TABLE default\n");
    fwrite(ptr, (size_t) *nbytes, ntotal, fs);
 
    free(ptr);
 
    return 0;
}

/* Close output data file */
int vtkclosefile_(void)
{
    fclose(fs);
    return 0;
}

