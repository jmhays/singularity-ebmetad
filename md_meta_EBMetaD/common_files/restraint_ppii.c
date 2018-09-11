/*
*******************************************************************************
*                                                                             *
*                                PLUMED                                       *
*   A Portable Plugin for Free Energy Calculations with Molecular Dynamics    *
*                              VERSION 1.3                                    *
*                                                                             *
*******************************************************************************
*
*  
*  Copyright (c) 2010 The PLUMED team.
*  See http://merlino.mi.infn.it/plumed for more information. 
*
*  This file is part of PLUMED.
*
*  PLUMED is free software: you can redistribute it and/or modify
*  it under the terms of the GNU Lesser General Public License as 
*  published by the Free Software Foundation, either version 3 of 
*  the License, or (at your option) any later version.
*
*  PLUMED is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General
*  Public License along with PLUMED.  
*  If not, see <http://www.gnu.org/licenses/>.
*
*  For more info, see:  http://merlino.mi.infn.it/plumed
*  or subscribe to plumed-users@googlegroups.com
*
*/
#include "metadyn.h"

real check_ppii(int from, int num, real *phi, real *psi, real der[6]);

void PREFIX ppii_restraint(int i_c, struct mtd_data_s *mtd_data) 
{
  real ***d, *dih, *phi, *psi, der[6];
  real psisum, cos_psi;
  int i, j, k, ix, idtc, ndtc, count, num;
  int a0, a1, a2, a3;
  int t1, t2, t3;
  rvec vp1, vp2, vp3;
  real pi, nv1, nv2, n21, sc1, sc2;
  real twopi, in21, sign, tmp;
  rvec r01, r21, r23, m, n;

  twopi = M_2PI;
  pi = M_PI;

  ndtc = colvar.type[i_c];	// number of dih
  snew(d, ndtc);
  for(idtc=0;idtc<ndtc;idtc++) {
    snew(d[idtc], 4);
    for(i=0;i<4;i++) snew(d[idtc][i], 3);
  }
  snew(dih, ndtc);
  snew(phi, ndtc/2);
  snew(psi, ndtc/2);

  count=0;
  for(idtc=0;idtc<ndtc;idtc++){
    for(ix=0;ix<3;ix++) {
      colvar.myder[i_c][4*idtc+0][ix] = 0.;
      colvar.myder[i_c][4*idtc+1][ix] = 0.;
      colvar.myder[i_c][4*idtc+2][ix] = 0.;
      colvar.myder[i_c][4*idtc+3][ix] = 0.;
    }

    a0 = colvar.cvatoms[i_c][4*idtc+0];
    a1 = colvar.cvatoms[i_c][4*idtc+1];
    a2 = colvar.cvatoms[i_c][4*idtc+2];
    a3 = colvar.cvatoms[i_c][4*idtc+3];

#if defined (PLUMED_GROMACS)
#if defined (PLUMED_GROMACS45)
    dih[idtc] = dih_angle(mtd_data->pos[a0], mtd_data->pos[a1], mtd_data->pos[a2], mtd_data->pos[a3], &mtd_data->metapbc, 
                          r01, r21, r23, m, n, &sign, &t1, &t2, &t3);
#else
    dih[idtc] = dih_angle(mtd_data->pos[a0], mtd_data->pos[a1], mtd_data->pos[a2], mtd_data->pos[a3], &mtd_data->metapbc, 
                          r01, r21, r23, m, n, &cos_psi, &sign, &t1, &t2, &t3);
#endif
#else
    dih[idtc] = dih_angle(mtd_data->pos[a0], mtd_data->pos[a1], mtd_data->pos[a2], mtd_data->pos[a3],
                          r01, r21, r23, m, n, &cos_psi, &sign);
#endif

    if(dih[idtc]>=pi) dih[idtc] -= twopi;
    else if(dih[idtc] < -pi) dih[idtc] += twopi;

    if(!(idtc%2)) { phi[count] = dih[idtc]; }
    else { 
      psi[count] = dih[idtc];
      count++;
    }

    oprod(r01, r21, vp1);
    oprod(r21, r23, vp2);
    oprod(vp1, vp2, vp3);

    nv1 = norm(vp1);
    nv2 = norm(vp2);
    n21 = norm(r21);

    sc1 = iprod(r01, r21);
    sc2 = iprod(r23, r21);

    in21 = 1./(n21*n21);

    sc1 = sc1*in21;
    sc2 = sc2*in21;

    for(i=0;i<3;++i) {
      d[idtc][0][i] = -n21 * vp1[i] / (nv1*nv1);
      d[idtc][3][i] = n21 * vp2[i] / (nv2*nv2);
      d[idtc][1][i] = ((sc1-1.0)*d[idtc][0][i] - sc2*d[idtc][3][i]);
      d[idtc][2][i] = ((sc2-1.0)*d[idtc][3][i] - sc1*d[idtc][0][i]);
    }
  }

  psisum = 0.;

  for(i=0;i<ndtc/2;i=i+3) {
    num=min(ndtc/2-i,3);
    psisum += check_ppii(i,num,phi,psi,der);
    for(k=0;k<2*num;k++) {
      for(j=0;j<3;j++) {
       d[2*i+k][0][j] *= der[k];
       d[2*i+k][3][j] *= der[k];
       d[2*i+k][1][j] *= der[k];
       d[2*i+k][2][j] *= der[k];
      }
    }
  }

  colvar.ss0[i_c] = psisum;

  for(idtc=0;idtc<ndtc;idtc++) {
    for(ix=0;ix<3;ix++){
      colvar.myder[i_c][4*idtc+0][ix] += -d[idtc][0][ix];
      colvar.myder[i_c][4*idtc+1][ix] += -d[idtc][1][ix];
      colvar.myder[i_c][4*idtc+2][ix] += -d[idtc][2][ix];
      colvar.myder[i_c][4*idtc+3][ix] += -d[idtc][3][ix];
    }
  }

  for(idtc=0;idtc<ndtc;idtc++) {
    sfree(d[idtc][0]);
    sfree(d[idtc][1]);
    sfree(d[idtc][2]);
    sfree(d[idtc][3]);
    sfree(d[idtc]);
  }
  sfree(d);
  sfree(dih);
  sfree(phi);
  sfree(psi);
}

real check_ppii(int from, int num, real *phi, real *psi, real der[6])
{
  int  i;
  real isPP=0., t1=0., t2=0., t3=0., t4=0.;
  real mphi=0., mpsi=0.;

  for(i=from;i<from+num;i++)
  {
    mphi += phi[i];
    mpsi += psi[i];
  }
  mphi /= (real) num;
  mpsi /= (real) num;

  t1 = 1.+cos(mphi+1.3083); 
  t2 = 1.+cos(mpsi-2.5294);
  isPP = t1*t2*3./4.; 

  t3 = -0.25*t2*sin(1.3083+mphi);
  t4 =  0.25*t1*sin(2.5294-mpsi); 

  for(i=0;i<num;i++){
    der[2*i]   = t3;
    der[2*i+1] = t4;
  }

  return isPP;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int PREFIX read_ppii(char **word, int count, t_plumed_input *input,int *iline,FILE *fplog)
{
  int iat, atnr1, atnr2, atnr3, atnr4, iw;
  double delta = 0.0;
  char string[400];
  int help;

  help = 0;

  iw=seek_word(word,"NDIH");
  if(iw>=0) { 
      sscanf(word[iw+1],"%i", &colvar.type[count]);
      if(colvar.type[count]%2) help=1;
  } else {
    fprintf(fplog,"|- NEEDED NDIH KEYWORD FOR PPII\n");
    help=1;
  }
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

   if(help){
         fprintf(fplog,"|- PPII SYNTAX:\n");
         fprintf(fplog,"|- NDIH (even)            : number of dihedrals\n");
         fprintf(fplog,"|- SIGMA              : hills width for this cv\n");
         fprintf(fplog,"e.g. \n");
         fprintf(fplog,"PPII NDIH 2 SIGMA 0.1 \n");
         fprintf(fplog,"168 170 172 188 \n");
         fprintf(fplog,"170 172 188 190 \n");
         plumed_error("PluMeD dead with errors: check log file");         
  }

  colvar.type_s[count]   = 17;
  colvar.natoms[count]   = 4*colvar.type[count]; 

  fprintf(fplog, "\n%1i-PPII: %i DIHEDRAL; ", count+1, colvar.type[count]);
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 

  snew(colvar.cvatoms[count], colvar.natoms[count]);
  snew(colvar.myder[count], colvar.natoms[count]);

  for(iat=0;iat<colvar.type[count];iat++) {
    (*iline)++;
    sscanf(input->words[*iline][0],"%i",&atnr1);
    sscanf(input->words[*iline][1],"%i",&atnr2);
    sscanf(input->words[*iline][2],"%i",&atnr3);
    sscanf(input->words[*iline][3],"%i",&atnr4);
    atnr1--;
    atnr2--;
    atnr3--;
    atnr4--;
    colvar.cvatoms[count][4*iat]   = atnr1;
    colvar.cvatoms[count][4*iat+1] = atnr2;
    colvar.cvatoms[count][4*iat+2] = atnr3;
    colvar.cvatoms[count][4*iat+3] = atnr4;
    fprintf(fplog, "|--DIH %i, ATOMS: %i %i %i %i\n", iat+1, atnr1+1, atnr2+1, atnr3+1, atnr4+1);
  }

  fprintf(fplog,"\n");
  return colvar.type[count];
}
