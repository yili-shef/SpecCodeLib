#ifndef _VTKDATA
#define _VTKDATA

extern int vtkopenfile_(char *, long int);
extern int vtk_sp_header_(int *, int *, int *, int *, 
        void *, void *, void *, void *, void *, void *);
extern int vtk_rg_header_(void *, void *, void *, int *, int *, int *, int *);
extern int vtk3dvtr_(void *, void *, void *, int *, int *, int *, int *, char *, long int);
extern int vtkscalar_(void *, int *, int *, int *, int *, char *, long int);
extern int vtkclosefile_(void);

#endif
