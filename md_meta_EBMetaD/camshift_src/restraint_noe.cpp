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

#ifdef HAVE_ALMOST

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <mdb.h>
#include <forcefield/const/meta.h>
#include <io/formostream.h>
#include <string>
#include <vector>
#include "csalmost.h"

using namespace std;
using namespace Almost;

extern "C" void PREFIX noe_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  static int readed[nconst_max];
  static vector<MetaD<double> > noe_list;
  char stringadb[512];
  char stringamdb[512];
  char stringapdb[512];
  char csfile[512];
  double energy = 0.0;
  bool printout = FALSE;
  Molecules molecules;

  if(firstTime&&logical.not_same_step) readed[i_c]=0;

  if(!readed[i_c]) {
    sprintf(stringadb, "%s/almnoe.dat", cam_shift.cam_data[i_c]);
    sprintf(stringamdb, "%s/%s", cam_shift.cam_data[i_c], cam_shift.cam_ff[i_c]);
    sprintf(stringapdb, "%s/template.pdb", cam_shift.cam_data[i_c]);
    MDB mdb(stringamdb);
    PDB pdb(stringapdb);
    //Molecules molecules;
    for(unsigned int i=0;i<pdb[0].size();i++){
      string str;
      str +='A'+i;
      Protein p(str);
      p.build_missing(pdb[0][i],mdb);
      if(cam_shift.disu[i_c]) p.auto_disu_bonds(3.5,mdb);
      molecules.add_protein(p);
    }

    MetaD<double> a=MetaD<double>(50,1000,10);
    a.read_upl(stringadb,molecules);
    noe_list.push_back(a);
    mol2pdb(molecules,"test-noe.pdb");
  }

  int N = colvar.natoms[i_c];
  Coor<double> coor(N); 
  Coor<double> forces(N);

  coor.clear();
  forces.clear();

  for (int i = 0; i < N; i++) {
     int ipos = 4 * i;
#ifdef PLUMED_GROMACS
     coor.coor[ipos  ] = (double) 10.*mtd_data->pos[i][0];
     coor.coor[ipos+1] = (double) 10.*mtd_data->pos[i][1];
     coor.coor[ipos+2] = (double) 10.*mtd_data->pos[i][2];
#else
     coor.coor[ipos  ] = (double) mtd_data->pos[i][0];
     coor.coor[ipos+1] = (double) mtd_data->pos[i][1];
     coor.coor[ipos+2] = (double) mtd_data->pos[i][2];
#endif
  }

  energy = noe_list[cam_shift.num[i_c]].Q_energy_force(coor, forces);

  for (int i = 0; i < N; i++)
  {
    int ipos = 4 * i;
#ifdef PLUMED_GROMACS
    colvar.myder[i_c][i][0] = (real) 41.86*forces.coor[ipos];  
    colvar.myder[i_c][i][1] = (real) 41.86*forces.coor[ipos+1];  
    colvar.myder[i_c][i][2] = (real) 41.86*forces.coor[ipos+2];
#else
    colvar.myder[i_c][i][0] = (real) forces.coor[ipos];  
    colvar.myder[i_c][i][1] = (real) forces.coor[ipos+1];  
    colvar.myder[i_c][i][2] = (real) forces.coor[ipos+2];
#endif
  }
  readed[i_c]++;
#ifdef PLUMED_GROMACS
  colvar.ss0[i_c] = (real) energy*4.186;
#else
  colvar.ss0[i_c] = (real) energy;
#endif
}

extern "C" int PREFIX read_noe(char **word, int count, t_plumed_input *input,FILE *fplog)
{
  int iw, j;
  int help;
  int static num = 0;
  double delta = 0.0;

  help = 0;
  cam_shift.mumo[count]=0;

  colvar.type_s[count] = 64;

  iw=seek_word(word,"DATA");
  if(iw>=0) {
       sscanf(word[iw+1],"%s", cam_shift.cam_data[count]);
  } else {
    fprintf(fplog,"|- NEEDED DATA KEYWORD FOR NOE\n");
    help=1;
  }

  iw=seek_word(word,"FF");
  if(iw>=0) {
       sscanf(word[iw+1],"%s", cam_shift.cam_ff[count]);
  } else {
    fprintf(fplog,"|- NEEDED FF KEYWORD FOR NOE SHIFTS\n");
    help=1;
  }

  iw=seek_word(word,"CYS-DISU");
  if(iw>=0) {
       cam_shift.disu[count]=1;
  } else {
       cam_shift.disu[count]=0;
  }
  iw = seek_word(word,"LIST");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR CAMSHIFT3-METHYL SHIFTS\n"); help=1;} 

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

   if(help){
         fprintf(fplog,"|- NOE SYNTAX:\n");
         fprintf(fplog,"|- DATA                   : parameters directory\n");
         fprintf(fplog,"|- LIST                   : group with all the atoms involved\n");
         fprintf(fplog,"|- SIGMA              : hills width for this cv\n");
         fprintf(fplog,"e.g. \n");
         fprintf(fplog,"NOE DATA data/ LIST <prot> SIGMA 0.1 \n");
         EXIT(); 
  }

  fprintf(fplog, "\n%1i-NOE: DATA %s, %i ATOMS ", count+1, cam_shift.cam_data[count], colvar.natoms[count]);
  if(logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n");

  cam_shift.num[count] = num; 
  cs_snew(colvar.myder[count], colvar.natoms[count]);

  num++; 

  return colvar.natoms[count];
}

#endif
