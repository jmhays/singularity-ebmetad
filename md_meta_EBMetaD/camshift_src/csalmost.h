#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "metadyn.h"
#include <string>
#include <vector>

#if defined (PLUMED_GROMACS)
#include "config.h"
#endif

#define cs_snew(ptr,nelem) (ptr)= (nelem==0 ? NULL : (typeof(ptr)) calloc(nelem,sizeof(*(ptr))))
#define cs_srenew(ptr,nelem) (ptr)= (typeof(ptr)) realloc(ptr,(nelem)*sizeof(*(ptr)))

double** PREFIX csd_2d_array_alloc(int ii,int jj);
int PREFIX csfree_2dr_array_alloc(double **xx,int ii);
void PREFIX camshift_sum(struct mtd_data_s *mtd_data,int nr,double r[]);
void PREFIX camshift_intersum(struct mtd_data_s *mtd_data,int nr, double r[]);
