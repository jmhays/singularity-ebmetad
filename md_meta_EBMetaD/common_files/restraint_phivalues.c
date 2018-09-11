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
*  Copyright (c) 2008-2011 The PLUMED team.
*  See http://www.plumed-code.org for more information. 
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
*  For more info, see:  http://www.plumed-code.org
*  or subscribe to plumed-users@googlegroups.com
*
*/

#include "metadyn.h"

void  PREFIX phivalues_restraint(int i_c, struct mtd_data_s *mtd_data) {

    int nres = colvar.type[i_c];
    int i, j, k, ix, aj, ak;
    real *phimat, tmp=0.,dfunc=0., tmp2=0., ss0=0.;
    real beta = colvar.beta[i_c];
    real mod_rij, totphi=0., itotphi=0.;
    rvec rij, *tmpder;

    snew(phimat, nres);
    snew(tmpder, colvar.natoms[i_c]);

    for(i=0;i<colvar.natoms[i_c];i++) for(ix=0;ix<3;ix++) colvar.myder[i_c][i][ix] = 0.;
    for(i=0; i<nres; i++) if(colvar.expphi[i_c][i]>-1.) totphi++;
    itotphi = 1./totphi;

    // this is a cycle of experimental phi-values
    for(i=0; i<nres; i++) {
      if(colvar.expphi[i_c][i]>-1.) {
        for(j=0;j<colvar.natoms[i_c]; j++) for(ix=0;ix<3;ix++) tmpder[j][ix] = 0.;
        // this is a cycle over residue[i].atoms
        // printf("RES %i PHI %lf NPH %lf\n", i, colvar.expphi[i_c][i], colvar.nphi[i_c][i]); fflush(stdout);
        for(j=0;j<colvar.natoms[i_c]; j++) {
          // this is a cycle over atoms in native contacts with residue[i]
          aj=colvar.cvatoms[i_c][j];
          if(colvar.ator[i_c][aj]==i) for(k=0;k<colvar.natoms[i_c]; k++) {
            ak=colvar.cvatoms[i_c][k];
            if(colvar.nmat[i_c][aj][ak]) {
              // fraction of native contacts formed
              if(colvar.cell_pbc[i_c]){
                minimal_image(mtd_data->pos[aj], mtd_data->pos[ak], &mod_rij, rij);
              } else {
                rij[0] = mtd_data->pos[aj][0]-mtd_data->pos[ak][0];
                rij[1] = mtd_data->pos[aj][1]-mtd_data->pos[ak][1];
                rij[2] = mtd_data->pos[aj][2]-mtd_data->pos[ak][2];
                mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
              }
              tmp = exp(beta*(mod_rij-colvar.r_0[i_c]));
              tmp2 = 1./(1.+tmp);
              phimat[i] += (tmp2/colvar.nphi[i_c][i]);
              dfunc = - 2.*(tmp2*tmp2)*beta*tmp/(colvar.nphi[i_c][i]*mod_rij);
              //printf("DFUNC %i %i %i %lf\n", i, j, k, dfunc); fflush(stdout);
              for(ix=0;ix<3;ix++) {
                tmpder[j][ix] += +dfunc*rij[ix];
                tmpder[k][ix] += -dfunc*rij[ix];
                //printf("DER A %i [%i] = %lf\n", j, ix, colvar.myder[i_c][j][ix]); fflush(stdout);
                //printf("DER A %i [%i] = %lf\n", k, ix, colvar.myder[i_c][k][ix]); fflush(stdout);
              }
            }
          }
        }

        ss0 += (phimat[i]-colvar.expphi[i_c][i])*(phimat[i]-colvar.expphi[i_c][i]);

        for(j=0;j<colvar.natoms[i_c]; j++) {
          for(ix=0;ix<3;ix++) {
            colvar.myder[i_c][j][ix] += tmpder[j][ix]*(phimat[i]-colvar.expphi[i_c][i])*itotphi;
          }
        }
      }
    }

    colvar.ss0[i_c] = ss0*itotphi;

    sfree(phimat);
    sfree(tmpder);
}

// ------------------------------------------------------------------------------------------------

int PREFIX read_phivalues(char **word, int count,t_plumed_input *input,FILE *fplog)
{

  int i,iw,j;
  int nres, resi, ai, resj, aj;
  FILE *nativa, *f_phi;
  double sigma = 0.0;
  char file_maps[129];
  char file_exp[129];
  int help;
  double r_0, phi;

  help=0; 
  colvar.cell_pbc[count]=1; // default is PBC

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &sigma); colvar.delta_r[count]  = (real) sigma; }
  iw = seek_word(word,"MAP");
  if(iw>=0) { sscanf(word[iw+1],"%s", file_maps);
  } else {
    fprintf(fplog,"|- NEEDED MAP KEYWORD FOR PHIVALUES\n");
    help=1;
  }
  iw = seek_word(word,"REFERENCE");
  if(iw>=0) { sscanf(word[iw+1],"%s", file_exp);
  } else {
    fprintf(fplog,"|- NEEDED REFERENCE KEYWORD FOR PHIVALUES\n");
    help=1;
  }
  iw = seek_word(word,"LIST");
  if(iw>=0){   
    j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
    colvar.natoms[count] = j;
  } else { fprintf(fplog,"|- NEEDED LIST KEYWORD FOR PHIVALUES\n"); help=1;}
 
  iw=seek_word(word,"NRES");
  if(iw>=0) {
       sscanf(word[iw+1],"%i", &nres);
       colvar.type[count]=nres;
  } else {
      fprintf(fplog,"|- NEEDED NRES KEYWORD FOR PHIVALUES\n");
      help=1;
  }
  iw=seek_word(word,"R_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &r_0); } else { fprintf(fplog,"|- NEEDED R_0 KEYWORD FOR PHIVALUES\n"); help=1;}
  colvar.r_0[count]      = (real) r_0;
  iw=seek_word(word,"LAMBDA");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &r_0); } else { fprintf(fplog,"|- NEEDED LAMBDA KEYWORD FOR PHIVALUES\n"); help=1;}
  colvar.beta[count]      = (real) r_0;
 
  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;}

  if(help){
         fprintf(fplog,"|- PHIVALUES SYNTAX:\n");
         fprintf(fplog,"|- MAP                : contact definition file \n");
         fprintf(fplog,"|- SIGMA              : hills width for this cv\n");
         fprintf(fplog,"|- e.g.\n");
         fprintf(fplog,"|- \n");
         fprintf(fplog,"|- PHIVALUES INDEX nat.ndx { NOPBC } SIGMA 1.0 \n");
         fprintf(fplog,"|- \n");
         plumed_error("PluMeD dead with errors: check log file");
  } 

  
  fprintf(fplog, "\n%1i-PHIVALUES: MAP file %s PHI file %s NATOMS %i NRES %i R_0 %lf", 
          count+1, file_maps, file_exp, colvar.natoms[count], colvar.type[count], colvar.r_0[count]);
  if(colvar.cell_pbc[count]) fprintf(fplog, " PBC ON ");
  else                       fprintf(fplog, " PBC OFF ");
  if (logical.do_hills){
        if (colvar.delta_r[count]>0){
                 fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
        }
  }
  else fprintf(fplog,"\n");

  snew(colvar.nphi[count], nres);
  snew(colvar.expphi[count], nres);
  for(i=0;i<nres;i++) colvar.expphi[count][i]=-1.;

  snew(colvar.ator[count], colvar.natoms[count]);
  snew(colvar.nmat[count], colvar.natoms[count]);
  for(i=0;i<colvar.natoms[count];i++) snew(colvar.nmat[count][i], colvar.natoms[count]);

  i=0;
  nativa=fopen(file_maps,"r");
  while(fscanf(nativa, "%i %i %i %i", &resi, &ai, &resj, &aj)!=EOF)
  {
    colvar.ator[count][ai]=resi;
    colvar.ator[count][aj]=resj;
    colvar.nmat[count][ai][aj]=1;
    i++;
  }
  fclose(nativa);

  fprintf(fplog, "|--%i NATIVE CONTACTS READED\n", i); 
 
  for(i=0;i<colvar.natoms[count];i++) 
    for(j=i;j<colvar.natoms[count];j++) {
      ai=colvar.cvatoms[count][i]; 
      aj=colvar.cvatoms[count][j]; 
      resi=colvar.ator[count][ai]; 
      resj=colvar.ator[count][aj]; 
      if(colvar.nmat[count][ai][aj]) {
        colvar.nphi[count][resi]+=1.; 
        colvar.nphi[count][resj]+=1.;
      }
    }

  i=0;
  f_phi=fopen(file_exp,"r");
  while(fscanf(f_phi, "%i %lf", &resi, &phi)!=EOF)
  {
    colvar.expphi[count][resi-1] = phi;
    i++;
  }
  fclose(f_phi);
  
  fprintf(fplog, "|--%i PHI-VALUES READED\n", i); 

  snew(colvar.myder[count], colvar.natoms[count]);

  fprintf(fplog,"\n");

  colvar.type_s[count]   = 65;

  return colvar.natoms[count];
}
