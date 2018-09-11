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

#if defined (PLUMED_GROMACS4) || defined (PLUMED_GROMACS45) || defined (NAMD)

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <mdb.h>
#include <forcefield/const/camshift2.h>
#include <io/formostream.h>
#include <string>
#include <vector>
#include "csalmost.h"

double** PREFIX csd_2d_array_alloc(int ii,int jj){
  double **xx;
  double *ptr;
  int i;
  ptr=(double *)calloc(ii*jj,sizeof(double));
  xx=(double **)calloc(ii,sizeof(double *));
  for (i=0;i<ii;i++)xx[i]=& ptr[jj*i];
  return xx;
};

int PREFIX csfree_2dr_array_alloc(double **xx,int ii){
  assert(xx[0]); free(xx[0]);
  assert(xx); free(xx);
  return 0;
};

using namespace std;
using namespace Almost;
//CamShift2 Static
string CamShift2::RingInfo::atomNames[4][6];
string CamShift2::RingInfo::types[4];
int    CamShift2::RingInfo::init_ = 0;

void PREFIX camshift_sum(struct mtd_data_s *mtd_data,int nr,double r[]){
  static double *buf=NULL;
  static int nalloc=0;
  int i;
#ifdef __PLUMED_MPI
  if (nr > nalloc) {
    nalloc = nr;
    cs_srenew(buf,nalloc);
  }

  MPI_Allreduce(r,buf,nr,MPI_DOUBLE,MPI_SUM,mtd_data->comm);

  for(i=0; i<nr; i++)
    r[i] = buf[i];
#endif
};

void PREFIX camshift_intersum(struct mtd_data_s *mtd_data,int nr, double r[]){
  static double *buf=NULL;
  static int nalloc=0;
  int i;
#ifdef __PLUMED_MPI
  if (nr > nalloc) {
    nalloc = nr;
    cs_srenew(buf,nalloc);
  }

  MPI_Allreduce(r,buf,nr,MPI_DOUBLE,MPI_SUM,mtd_data->intercomm);

  for(i=0; i<nr; i++)
    r[i] = buf[i];
#endif
};

extern "C" void PREFIX camshift_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  static int readed[nconst_max];
  static vector<CamShift2> cam_list;
  char stringadb[512];
  char stringamdb[512];
  char stringapdb[512];
  char csfile[512];
  double energy = 0.0;
  bool printout = FALSE;
  static Molecules molecules;

  if(firstTime&&logical.not_same_step) readed[i_c]=0;
  
  if(!readed[i_c]) {
    sprintf(stringadb, "%s/camshift.db", cam_shift.cam_data[i_c]);
    sprintf(stringamdb, "%s/%s", cam_shift.cam_data[i_c], cam_shift.cam_ff[i_c]);
    sprintf(stringapdb, "%s/template.pdb", cam_shift.cam_data[i_c]);
    MDB mdb(stringamdb);
    PDB pdb(stringapdb);
    for(unsigned int i=0;i<pdb[0].size();i++){
      string str;
      str +='A'+i;
      Protein p(str);
      if(cam_shift.autob[i_c]) p.build_missing(pdb[0][i],mdb);
      else p.build(pdb[0][i],mdb);
      if(cam_shift.disu[i_c]) p.auto_disu_bonds(3.5,mdb);
      molecules.add_protein(p);
    }
    //CamShift2 a = CamShift2(stringapdb, stringamdb, stringadb, cam_shift.disu[i_c]);
    CamShift2 a = CamShift2(molecules, stringadb);

    sprintf(stringadb, "%s/CAshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "CA");
    sprintf(stringadb, "%s/CBshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "CB");
    sprintf(stringadb, "%s/Cshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "C");
    sprintf(stringadb, "%s/HAshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "HA");
    sprintf(stringadb, "%s/Hshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "H");
    sprintf(stringadb, "%s/Nshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "N");
    a.set_flat_bottom_const(cam_shift.grains[i_c]);
    a.set_lambda(1);
    cam_list.push_back(a);
    mol2pdb(molecules,"test-camsh.pdb");
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

  energy = cam_list[cam_shift.num[i_c]].energy_force(coor, forces);
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

  printout = (colvar.logic[i_c]&&((logical.not_same_step)&&(!(colvar.it%colvar.nt_print))));
  if(printout) {
    sprintf(csfile, "cs%i.dat", colvar.it);
    cam_list[cam_shift.num[i_c]].printout_chemical_shifts(csfile);
  }

  readed[i_c]++;
#ifdef PLUMED_GROMACS
  colvar.ss0[i_c] = (real) energy*4.186;
#else
  colvar.ss0[i_c] = (real) energy;
#endif
}

extern "C" void PREFIX camshiftens_restraint(int i_c, struct mtd_data_s *mtd_data)
{
#ifdef __PLUMED_MPI
  static int readed[nconst_max]; 
  static int repliche;
  static int replica;
  static vector<CamShift2> cam_list;
  static double fact=1.0;
  char stringadb[512];
  char stringamdb[512];
  char stringapdb[512];
  char csfile[512];
  double energy = 0.0; 
  double **sh;
  double **tca, **tcb, **tco, **tha, **thn, **tnh;
  bool printout = FALSE;
  static Molecules molecules;

  if(firstTime&&logical.not_same_step) {
    if(logical.remd) { 
      if(mtd_data->repl!=-1) { 
        repliche=mtd_data->nrepl; 
        fact = 1.0/mtd_data->nrepl;
        replica = mtd_data->repl;
      } else if(mtd_data->mcr->ms->nsim>1) {
        repliche = mtd_data->mcr->ms->nsim; 
        fact = 1.0/mtd_data->mcr->ms->nsim;
        replica = mtd_data->mcr->ms->sim;
      } else {
        plumed_error("IN ORDER TO USE CAMSHIFT ENSEMBLE YOU MUST RUN GROMACS WITH MULTIPLE REPLICAS AND BIASXMD");
      }
    } else {
      plumed_error("IN ORDER TO USE CAMSHIFT ENSEMBLE YOU MUST RUN GROMACS WITH MULTIPLE REPLICAS AND BIASXMD");
    }
    readed[i_c]=0;
    if(!cam_shift.mumo[i_c]) fprintf(mtd_data->fplog, "|- CAMSHIFT AVERAGED OVER %i REPL WITH FACT %lf\n", repliche, fact);
    else fprintf(mtd_data->fplog, "|- CAMSHIFT AVERAGED OVER 2 REPL WITH FACT 0.5\n");
    fflush(mtd_data->fplog);
  }

  if(!readed[i_c]) {
    sprintf(stringadb, "%s/camshift.db", cam_shift.cam_data[i_c]);
    sprintf(stringamdb, "%s/%s", cam_shift.cam_data[i_c], cam_shift.cam_ff[i_c]);
    sprintf(stringapdb, "%s/template.pdb", cam_shift.cam_data[i_c]);
    MDB mdb(stringamdb);
    PDB pdb(stringapdb);
    for(unsigned int i=0;i<pdb[0].size();i++){
      string str;
      str +='A'+i;
      Protein p(str);
      if(cam_shift.autob[i_c]) p.build_missing(pdb[0][i],mdb);
      else p.build(pdb[0][i],mdb);
      if(cam_shift.disu[i_c]) p.auto_disu_bonds(3.5,mdb);
      molecules.add_protein(p);
    }
    //CamShift2 a = CamShift2(stringapdb, stringamdb, stringadb, cam_shift.disu[i_c]);
    CamShift2 a = CamShift2(molecules, stringadb);

    sprintf(stringadb, "%s/CAshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "CA");
    sprintf(stringadb, "%s/CBshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "CB");
    sprintf(stringadb, "%s/Cshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "C");
    sprintf(stringadb, "%s/HAshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "HA");
    sprintf(stringadb, "%s/Hshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "H");
    sprintf(stringadb, "%s/Nshifts.dat", cam_shift.cam_data[i_c]);
    a.read_cs(stringadb, "N");
    a.set_flat_bottom_const(cam_shift.grains[i_c]);
    a.set_lambda(1);
    cam_list.push_back(a);
    mol2pdb(molecules,"test-camsh.pdb");
  }

  int N = colvar.natoms[i_c];
  int numResidues = colvar.type[i_c];

  Coor<double> coor(N); 
  Coor<double> forces(N);

  coor.clear();
  forces.clear();

  sh = csd_2d_array_alloc(numResidues, 6);
  /* temp arrays for subreplica averaging */
  if(cam_shift.mumo[i_c]) {
    tca = csd_2d_array_alloc(repliche, numResidues);
    tcb = csd_2d_array_alloc(repliche, numResidues);
    tco = csd_2d_array_alloc(repliche, numResidues);
    tha = csd_2d_array_alloc(repliche, numResidues);
    thn = csd_2d_array_alloc(repliche, numResidues);
    tnh = csd_2d_array_alloc(repliche, numResidues);
  }
  /* end */

  for(int i = 0; i < N; i++) {
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

  // calculate all the chemical shifts
  cam_list[cam_shift.num[i_c]].ens_return_shifts(coor, sh);

  printout = (colvar.logic[i_c]&&((logical.not_same_step)&&(!(colvar.it%colvar.nt_print))));
  if(printout) {
    sprintf(csfile, "cs%i-%i.dat", replica, colvar.it);
    cam_list[cam_shift.num[i_c]].printout_chemical_shifts(csfile);
  }

  /* old way */
  if(!cam_shift.mumo[i_c]) { 
    // communication
    if(MASTER(mtd_data->mcr)) {
      // among replicas 
      camshift_intersum(mtd_data, numResidues*6, & sh[0][0]); 
      for(int i=0;i<6;i++) for(int j=0;j<numResidues;j++) sh[j][i] *= fact; 
    } else for(int i=0;i<6;i++) for(int j=0;j<numResidues;j++) sh[j][i] = 0.;
    // inside each replica
    camshift_sum(mtd_data, numResidues*6, & sh[0][0]);
  /* end */
  } else {
  /* new way */
    if(MASTER(mtd_data->mcr)) {
      // among replicas
      for(int j=0;j<numResidues;j++) {
        tca[replica][j] = sh[j][0]; 
        tcb[replica][j] = sh[j][1]; 
        tco[replica][j] = sh[j][2]; 
        tha[replica][j] = sh[j][3]; 
        thn[replica][j] = sh[j][4]; 
        tnh[replica][j] = sh[j][5];
      } 
      camshift_intersum(mtd_data, repliche*numResidues, & tca[0][0]); 
      camshift_intersum(mtd_data, repliche*numResidues, & tcb[0][0]); 
      camshift_intersum(mtd_data, repliche*numResidues, & tco[0][0]); 
      camshift_intersum(mtd_data, repliche*numResidues, & tha[0][0]); 
      camshift_intersum(mtd_data, repliche*numResidues, & thn[0][0]); 
      camshift_intersum(mtd_data, repliche*numResidues, & tnh[0][0]);
      for(int i=0; i<repliche; i=i+2) {
        for(int j=0;j<numResidues;j++) {
          tca[i][j] = (tca[i][j]+tca[i+1][j])*0.5;
          tcb[i][j] = (tcb[i][j]+tcb[i+1][j])*0.5;
          tco[i][j] = (tco[i][j]+tco[i+1][j])*0.5;
          tha[i][j] = (tha[i][j]+tha[i+1][j])*0.5;
          thn[i][j] = (thn[i][j]+thn[i+1][j])*0.5;
          tnh[i][j] = (tnh[i][j]+tnh[i+1][j])*0.5;
          tca[i+1][j] = tca[i][j];
          tcb[i+1][j] = tcb[i][j];
          tco[i+1][j] = tco[i][j];
          tha[i+1][j] = tha[i][j];
          thn[i+1][j] = thn[i][j];
          tnh[i+1][j] = tnh[i][j];
        }
      }
    } 
    for(int j=0;j<numResidues;j++) {
      sh[j][0] = tca[replica][j];
      sh[j][1] = tcb[replica][j];
      sh[j][2] = tco[replica][j];
      sh[j][3] = tha[replica][j];
      sh[j][4] = thn[replica][j];
      sh[j][5] = tnh[replica][j];
    } 
    camshift_sum(mtd_data, numResidues*6, & sh[0][0]);
    /* end */
  }
  // calculate all the forces
  energy = cam_list[cam_shift.num[i_c]].ens_energy_force(coor, forces, sh);

  for(int i = 0; i < N; i++) {
    int ipos = 4 * i;
#ifdef PLUMED_GROMACS
    colvar.myder[i_c][i][0] = (real) 41.86*fact*forces.coor[ipos];  
    colvar.myder[i_c][i][1] = (real) 41.86*fact*forces.coor[ipos+1];  
    colvar.myder[i_c][i][2] = (real) 41.86*fact*forces.coor[ipos+2];
#else
    colvar.myder[i_c][i][0] = (real) fact*forces.coor[ipos];  
    colvar.myder[i_c][i][1] = (real) fact*forces.coor[ipos+1];  
    colvar.myder[i_c][i][2] = (real) fact*forces.coor[ipos+2];
#endif
  }

  csfree_2dr_array_alloc(sh, numResidues);
  if(cam_shift.mumo[i_c]) {
    csfree_2dr_array_alloc(tca, repliche);
    csfree_2dr_array_alloc(tcb, repliche);
    csfree_2dr_array_alloc(tco, repliche);
    csfree_2dr_array_alloc(tha, repliche);
    csfree_2dr_array_alloc(thn, repliche);
    csfree_2dr_array_alloc(tnh, repliche);
  }

  readed[i_c]++;
#ifdef PLUMED_GROMACS
  colvar.ss0[i_c] = (real) energy*4.186;
#else
  colvar.ss0[i_c] = (real) energy;
#endif
#else
  plumed_error("IN ORDER TO USE CAMSHIFT ENSEMBLE YOU MUST COMPILE GROMACS WITH MPI!");
#endif
}

#endif

extern "C" int PREFIX read_camshift(char **word, int count, t_plumed_input *input,FILE *fplog)
{
  int iw, j;
  int help;
  int static num = 0;
  double delta = 0.0;

  help = 0;
  cam_shift.mumo[count]=0;

  iw=seek_word(word,"PRINT");
  if(iw>=0) colvar.logic[count] = 1;

  iw=seek_word(word,"ENSEMBLE");
  if(iw>=0) {
    colvar.type_s[count] = 61;
    iw=seek_word(word,"MUMO");
    if(iw>=0) cam_shift.mumo[count]=1;
  } else colvar.type_s[count] = 60;

  iw=seek_word(word,"DATA");
  if(iw>=0) {
       sscanf(word[iw+1],"%s", cam_shift.cam_data[count]);
  } else {
    fprintf(fplog,"|- NEEDED DATA KEYWORD FOR CAMSHIFT\n");
    help=1;
  }

  iw=seek_word(word,"FF");
  if(iw>=0) {
       sscanf(word[iw+1],"%s", cam_shift.cam_ff[count]);
  } else {
    fprintf(fplog,"|- NEEDED FF KEYWORD FOR CAMSHIFT\n");
    help=1;
  }

  iw=seek_word(word,"CYS-DISU");
  if(iw>=0) {
       cam_shift.disu[count]=1;
  } else {
       cam_shift.disu[count]=0;
  }
  iw=seek_word(word,"NOAUTO");
  if(iw>=0) {
       cam_shift.autob[count]=0;
  } else {
       cam_shift.autob[count]=1;
  }

  iw=seek_word(word,"NRES");
  if(iw>=0) {
       sscanf(word[iw+1],"%i", &colvar.type[count]);
  } else {
    if(colvar.type_s[count]==61) {
      fprintf(fplog,"|- NEEDED NRES KEYWORD FOR CAMSHIFT ENSEMBLE\n");
      help=1;
    }
  }

  iw=seek_word(word,"FLAT");
  if(iw>=0) {
       sscanf(word[iw+1],"%lf", &cam_shift.grains[count]);
  } else {
       cam_shift.grains[count] = 1.;
  }

  iw = seek_word(word,"LIST");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR CAMSHIFT\n"); help=1;} 

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

   if(help){
         fprintf(fplog,"|- CAMSHIFT SYNTAX:\n");
         fprintf(fplog,"|- DATA                   : parameters directory\n");
         fprintf(fplog,"|- LIST                   : group with all the atoms involved\n");
         fprintf(fplog,"|- SIGMA              : hills width for this cv\n");
         fprintf(fplog,"e.g. \n");
         fprintf(fplog,"CAMSHIFT DATA data/ LIST <prot> SIGMA 0.1 \n");
         EXIT(); 
  }

  fprintf(fplog, "\n%1i-CAMSHIFT: DATA %s, %i ATOMS ", count+1, cam_shift.cam_data[count], colvar.natoms[count]);
  if(colvar.type_s[count]==61) fprintf(fplog, "ENSEMBLE AVERAGING ON ");
  if(logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n");

  cam_shift.num[count] = num; 
  cs_snew(colvar.myder[count], colvar.natoms[count]);

  num++; 

  return colvar.natoms[count];
}

#endif
