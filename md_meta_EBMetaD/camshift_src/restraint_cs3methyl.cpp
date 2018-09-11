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
#include <camshift3/meth/methcs.h>
#include <io/formostream.h>
#include <string>
#include <vector>
#include "csalmost.h"

using namespace std;
using namespace Almost;

extern "C" void PREFIX cs3meth_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  static int readed[nconst_max];
  static vector<MethCS*> cam_list;
  char stringadb[512];
  char stringamdb[512];
  char stringapdb[512];
  char csfile[512];
  double energy = 0.0;
  bool printout = FALSE;
  static Molecules molecules;

  if(firstTime&&logical.not_same_step) readed[i_c]=0;

  if(!readed[i_c]) {
    sprintf(stringadb, "%s/methyls.cs", cam_shift.cam_data[i_c]);
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
    MethCS* a = new MethCS(molecules);
    a->read_cs(stringadb);
    a->set_w_cs(2);
    a->set_flat_bottom_const(cam_shift.grains[i_c]);
    cam_list.push_back(a);
    mol2pdb(molecules,"test-methyl.pdb");
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

  //energy = cam_list[cam_shift.num[i_c]]->calc_cs_force(coor, forces);
  energy = cam_list[cam_shift.num[i_c]]->calc_cs_and_force(coor, forces);

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
    sprintf(csfile, "csmethyl%i.dat", colvar.it);
    cam_list[cam_shift.num[i_c]]->write_cs(csfile);
  }

  readed[i_c]++;
#ifdef PLUMED_GROMACS
  colvar.ss0[i_c] = (real) energy*4.186;
#else
  colvar.ss0[i_c] = (real) energy;
#endif
}

extern "C" void PREFIX cs3methens_restraint(int i_c, struct mtd_data_s *mtd_data)
{
#ifdef __PLUMED_MPI
  static int readed[nconst_max]; 
  static int repliche;
  static int replica;
  static vector<MethCS*> cam_list;
  static double fact=1.0;
  char stringadb[512];
  char stringamdb[512];
  char stringapdb[512];
  char csfile[512];
  double energy = 0.0; 
  double **tala, **tile1, **tile2, **tleu1, **tleu2, **tthr, **tval1, **tval2;
  double **s1ala, **s1ile1, **s1ile2, **s1leu1, **s1leu2, **s1thr, **s1val1, **s1val2;
  double **s2ala, **s2ile1, **s2ile2, **s2leu1, **s2leu2, **s2thr, **s2val1, **s2val2;
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
        plumed_error("IN ORDER TO USE CAMSHIFT3-METHYL SHIFTS ENSEMBLE YOU MUST RUN GROMACS WITH MULTIPLE REPLICAS AND BIASXMD");
      }
    } else {
      plumed_error("IN ORDER TO USE CAMSHIFT3-METHYL SHIFTS ENSEMBLE YOU MUST RUN GROMACS WITH MULTIPLE REPLICAS AND BIASXMD");
    }
    readed[i_c]=0;
    if(!cam_shift.mumo[i_c]) fprintf(mtd_data->fplog, "|- CAMSHIFT3-METHYL SHIFTS AVERAGED OVER %i REPL WITH FACT %lf\n", repliche, fact);
    else fprintf(mtd_data->fplog, "|- CAMSHIFT3-METHYL SHIFTS AVERAGED OVER 2 REPL WITH FACT 0.5\n");
    fflush(mtd_data->fplog);
  }

  if(!readed[i_c]) {
    sprintf(stringadb, "%s/methyls.cs", cam_shift.cam_data[i_c]);
    sprintf(stringamdb, "%s/%s", cam_shift.cam_data[i_c], cam_shift.cam_ff[i_c]);
    sprintf(stringapdb, "%s/template.pdb", cam_shift.cam_data[i_c]);
    MDB mdb(stringamdb);
    PDB pdb(stringapdb);
    for(unsigned int i=0;i<pdb[0].size();i++){
      string str;
      str +='A'+i;
      Protein p(str);
      p.build_missing(pdb[0][i],mdb);
      if(cam_shift.disu[i_c]) p.auto_disu_bonds(3.5,mdb);
      molecules.add_protein(p);
    }
    MethCS* a = new MethCS(molecules);
    a->read_cs(stringadb);
    a->set_w_cs(2);
    a->set_flat_bottom_const(cam_shift.grains[i_c]);
    cam_list.push_back(a);
  }

  int N = colvar.natoms[i_c];
  int numResidues = colvar.type[i_c];

  Coor<double> coor(N); 
  Coor<double> forces(N);

  coor.clear();
  forces.clear();

  /* temp arrays for subreplica averaging */
  tala =  csd_2d_array_alloc(repliche, numResidues);
  tile1 = csd_2d_array_alloc(repliche, numResidues);
  tile2 = csd_2d_array_alloc(repliche, numResidues);
  tleu1 = csd_2d_array_alloc(repliche, numResidues);
  tleu2 = csd_2d_array_alloc(repliche, numResidues);
  tthr =  csd_2d_array_alloc(repliche, numResidues);
  tval1 = csd_2d_array_alloc(repliche, numResidues);
  tval2 = csd_2d_array_alloc(repliche, numResidues);
  /* end */

  if(cam_shift.mumo[i_c]==1) {
    s1ala =  csd_2d_array_alloc(repliche, numResidues);
    s1ile1 = csd_2d_array_alloc(repliche, numResidues);
    s1ile2 = csd_2d_array_alloc(repliche, numResidues);
    s1leu1 = csd_2d_array_alloc(repliche, numResidues);
    s1leu2 = csd_2d_array_alloc(repliche, numResidues);
    s1thr =  csd_2d_array_alloc(repliche, numResidues);
    s1val1 = csd_2d_array_alloc(repliche, numResidues);
    s1val2 = csd_2d_array_alloc(repliche, numResidues);
    s2ala =  csd_2d_array_alloc(repliche, numResidues);
    s2ile1 = csd_2d_array_alloc(repliche, numResidues);
    s2ile2 = csd_2d_array_alloc(repliche, numResidues);
    s2leu1 = csd_2d_array_alloc(repliche, numResidues);
    s2leu2 = csd_2d_array_alloc(repliche, numResidues);
    s2thr =  csd_2d_array_alloc(repliche, numResidues);
    s2val1 = csd_2d_array_alloc(repliche, numResidues);
    s2val2 = csd_2d_array_alloc(repliche, numResidues);
  } 

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
  cam_list[cam_shift.num[i_c]]->calc_cs(coor);

  printout = (colvar.logic[i_c]&&((logical.not_same_step)&&(!(colvar.it%colvar.nt_print))));
  if(printout) {
    sprintf(csfile, "csmethyl%i-%i.dat", replica, colvar.it);
    cam_list[cam_shift.num[i_c]]->write_cs(csfile);
  }

  if(MASTER(mtd_data->mcr)) {
    // among replicas
    int size = cam_list[cam_shift.num[i_c]]->ala_calc_hb.size();
    for(int j=0;j<size;j++) tala[replica][j] = cam_list[cam_shift.num[i_c]]->ala_calc_hb[j];
    size = cam_list[cam_shift.num[i_c]]->ile_calc_hd.size();
    for(int j=0;j<size;j++) tile1[replica][j] = cam_list[cam_shift.num[i_c]]->ile_calc_hd[j];
    size = cam_list[cam_shift.num[i_c]]->ile_calc_hg2.size();
    for(int j=0;j<size;j++) tile2[replica][j] = cam_list[cam_shift.num[i_c]]->ile_calc_hg2[j];
    size = cam_list[cam_shift.num[i_c]]->leu_calc_hd1.size();
    for(int j=0;j<size;j++) tleu1[replica][j] = cam_list[cam_shift.num[i_c]]->leu_calc_hd1[j];
    size = cam_list[cam_shift.num[i_c]]->leu_calc_hd2.size();
    for(int j=0;j<size;j++) tleu2[replica][j] = cam_list[cam_shift.num[i_c]]->leu_calc_hd2[j];
    size = cam_list[cam_shift.num[i_c]]->thr_calc_hg2.size();
    for(int j=0;j<size;j++) tthr[replica][j] = cam_list[cam_shift.num[i_c]]->thr_calc_hg2[j];
    size = cam_list[cam_shift.num[i_c]]->val_calc_hg1.size();
    for(int j=0;j<size;j++) tval1[replica][j] = cam_list[cam_shift.num[i_c]]->val_calc_hg1[j];
    size = cam_list[cam_shift.num[i_c]]->val_calc_hg2.size();
    for(int j=0;j<size;j++) tval2[replica][j] = cam_list[cam_shift.num[i_c]]->val_calc_hg2[j];

    camshift_intersum(mtd_data, repliche*numResidues, & tala[0][0]); 
    camshift_intersum(mtd_data, repliche*numResidues, & tile1[0][0]); 
    camshift_intersum(mtd_data, repliche*numResidues, & tile2[0][0]); 
    camshift_intersum(mtd_data, repliche*numResidues, & tleu1[0][0]); 
    camshift_intersum(mtd_data, repliche*numResidues, & tleu2[0][0]); 
    camshift_intersum(mtd_data, repliche*numResidues, & tthr[0][0]);
    camshift_intersum(mtd_data, repliche*numResidues, & tval1[0][0]); 
    camshift_intersum(mtd_data, repliche*numResidues, & tval2[0][0]);
 

    if(cam_shift.mumo[i_c]==1) {
      for(int i=0; i<repliche; i=i+2) {
        for(int j=0;j<numResidues;j++) {
          int k;
          if(i<repliche-1) k=i+1; else k=0;
          s1ala[i][j] = (tala[i][j]+tala[k][j])*0.5;
          s1ile1[i][j] = (tile1[i][j]+tile1[k][j])*0.5;
          s1ile2[i][j] = (tile2[i][j]+tile2[k][j])*0.5;
          s1leu1[i][j] = (tleu1[i][j]+tleu1[k][j])*0.5;
          s1leu2[i][j] = (tleu2[i][j]+tleu2[k][j])*0.5;
          s1thr[i][j] = (tthr[i][j]+tthr[k][j])*0.5;
          s1val1[i][j] = (tval1[i][j]+tval1[k][j])*0.5;
          s1val2[i][j] = (tval2[i][j]+tval2[k][j])*0.5;

          if(i>0) k=i-1; else k=repliche-1;
          s2ala[i][j] = (tala[i][j]+tala[k][j])*0.5;
          s2ile1[i][j] = (tile1[i][j]+tile1[k][j])*0.5;
          s2ile2[i][j] = (tile2[i][j]+tile2[k][j])*0.5;
          s2leu1[i][j] = (tleu1[i][j]+tleu1[k][j])*0.5;
          s2leu2[i][j] = (tleu2[i][j]+tleu2[k][j])*0.5;
          s2thr[i][j] = (tthr[i][j]+tthr[k][j])*0.5;
          s2val1[i][j] = (tval1[i][j]+tval1[k][j])*0.5;
          s2val2[i][j] = (tval2[i][j]+tval2[k][j])*0.5;
        }
      }
    } else if(cam_shift.mumo[i_c]==2) { 
      for(int i=0; i<repliche; i=i+2) {
        for(int j=0;j<numResidues;j++) {
          tala[i][j] = (tala[i][j]+tala[i+1][j])*0.5;
          tile1[i][j] = (tile1[i][j]+tile1[i+1][j])*0.5;
          tile2[i][j] = (tile2[i][j]+tile2[i+1][j])*0.5;
          tleu1[i][j] = (tleu1[i][j]+tleu1[i+1][j])*0.5;
          tleu2[i][j] = (tleu2[i][j]+tleu2[i+1][j])*0.5;
          tthr[i][j] = (tthr[i][j]+tthr[i+1][j])*0.5;
          tval1[i][j] = (tval1[i][j]+tval1[i+1][j])*0.5;
          tval2[i][j] = (tval2[i][j]+tval2[i+1][j])*0.5;
          tala[i+1][j] = tala[i][j];
          tile1[i+1][j] = tile1[i][j];
          tile2[i+1][j] = tile2[i][j];
          tleu1[i+1][j] = tleu1[i][j];
          tleu2[i+1][j] = tleu2[i][j];
          tthr[i+1][j] = tthr[i][j];
          tval1[i+1][j] = tval1[i][j];
          tval2[i+1][j] = tval2[i][j];
        }
      }
    } else {
      for(int i=1; i<repliche; i++) {
        for(int j=0;j<numResidues;j++) {
          tala[0][j]  += tala[i][j]; 
          tile1[0][j] += tile1[i][j];
          tile2[0][j] += tile2[i][j];
          tleu1[0][j] += tleu1[i][j];
          tleu2[0][j] += tleu2[i][j];
          tthr[0][j]  += tthr[i][j]; 
          tval1[0][j] += tval1[i][j];
          tval2[0][j] += tval2[i][j];
        }
      }
      for(int i=repliche-1; i>=0; i--) {
        for(int j=0;j<numResidues;j++) {
          tala[i][j]  = tala [0][j]*fact;
          tile1[i][j] = tile1[0][j]*fact;
          tile2[i][j] = tile2[0][j]*fact;
          tleu1[i][j] = tleu1[0][j]*fact;
          tleu2[i][j] = tleu2[0][j]*fact;
          tthr[i][j]  = tthr [0][j]*fact;
          tval1[i][j] = tval1[0][j]*fact;
          tval2[i][j] = tval2[0][j]*fact;
        }
      }
    }
  }
  // intra replica sharing 
  if(cam_shift.mumo[i_c]!=1) {
    camshift_sum(mtd_data, repliche*numResidues, & tala[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & tile1[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & tile2[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & tleu1[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & tleu2[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & tthr[0][0]);
    camshift_sum(mtd_data, repliche*numResidues, & tval1[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & tval2[0][0]);

    int size = cam_list[cam_shift.num[i_c]]->ala_calc_hb.size();
    for(int j=0;j<size;j++)  cam_list[cam_shift.num[i_c]]->ala_calc_hb[j]  = tala[replica][j] ;
    size = cam_list[cam_shift.num[i_c]]->ile_calc_hd.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->ile_calc_hd[j]  = tile1[replica][j];
    size = cam_list[cam_shift.num[i_c]]->ile_calc_hg2.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->ile_calc_hg2[j] = tile2[replica][j];
    size = cam_list[cam_shift.num[i_c]]->leu_calc_hd1.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->leu_calc_hd1[j] = tleu1[replica][j];
    size = cam_list[cam_shift.num[i_c]]->leu_calc_hd2.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->leu_calc_hd2[j] = tleu2[replica][j];
    size = cam_list[cam_shift.num[i_c]]->thr_calc_hg2.size();
    for(int j=0;j<size;j++)  cam_list[cam_shift.num[i_c]]->thr_calc_hg2[j] = tthr[replica][j] ;
    size = cam_list[cam_shift.num[i_c]]->val_calc_hg1.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->val_calc_hg1[j] = tval1[replica][j];
    size = cam_list[cam_shift.num[i_c]]->val_calc_hg2.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->val_calc_hg2[j] = tval2[replica][j];

    // calculate all the forces
    energy = cam_list[cam_shift.num[i_c]]->ens_calc_cs_force(coor, forces);

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

  } else {
    camshift_sum(mtd_data, repliche*numResidues, & s1ala[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s1ile1[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s1ile2[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s1leu1[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s1leu2[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s1thr[0][0]);
    camshift_sum(mtd_data, repliche*numResidues, & s1val1[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s1val2[0][0]);
    camshift_sum(mtd_data, repliche*numResidues, & s2ala[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s2ile1[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s2ile2[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s2leu1[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s2leu2[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s2thr[0][0]);
    camshift_sum(mtd_data, repliche*numResidues, & s2val1[0][0]); 
    camshift_sum(mtd_data, repliche*numResidues, & s2val2[0][0]);

    int size = cam_list[cam_shift.num[i_c]]->ala_calc_hb.size();
    for(int j=0;j<size;j++)  cam_list[cam_shift.num[i_c]]->ala_calc_hb[j]  = s1ala[replica][j] ;
    size = cam_list[cam_shift.num[i_c]]->ile_calc_hd.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->ile_calc_hd[j]  = s1ile1[replica][j];
    size = cam_list[cam_shift.num[i_c]]->ile_calc_hg2.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->ile_calc_hg2[j] = s1ile2[replica][j];
    size = cam_list[cam_shift.num[i_c]]->leu_calc_hd1.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->leu_calc_hd1[j] = s1leu1[replica][j];
    size = cam_list[cam_shift.num[i_c]]->leu_calc_hd2.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->leu_calc_hd2[j] = s1leu2[replica][j];
    size = cam_list[cam_shift.num[i_c]]->thr_calc_hg2.size();
    for(int j=0;j<size;j++)  cam_list[cam_shift.num[i_c]]->thr_calc_hg2[j] = s1thr[replica][j] ;
    size = cam_list[cam_shift.num[i_c]]->val_calc_hg1.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->val_calc_hg1[j] = s1val1[replica][j];
    size = cam_list[cam_shift.num[i_c]]->val_calc_hg2.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->val_calc_hg2[j] = s1val2[replica][j];

    // calculate all the forces from the first set
    energy = cam_list[cam_shift.num[i_c]]->ens_calc_cs_force(coor, forces);

    size = cam_list[cam_shift.num[i_c]]->ala_calc_hb.size();
    for(int j=0;j<size;j++)  cam_list[cam_shift.num[i_c]]->ala_calc_hb[j]  = s2ala[replica][j] ;
    size = cam_list[cam_shift.num[i_c]]->ile_calc_hd.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->ile_calc_hd[j]  = s2ile1[replica][j];
    size = cam_list[cam_shift.num[i_c]]->ile_calc_hg2.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->ile_calc_hg2[j] = s2ile2[replica][j];
    size = cam_list[cam_shift.num[i_c]]->leu_calc_hd1.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->leu_calc_hd1[j] = s2leu1[replica][j];
    size = cam_list[cam_shift.num[i_c]]->leu_calc_hd2.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->leu_calc_hd2[j] = s2leu2[replica][j];
    size = cam_list[cam_shift.num[i_c]]->thr_calc_hg2.size();
    for(int j=0;j<size;j++)  cam_list[cam_shift.num[i_c]]->thr_calc_hg2[j] = s2thr[replica][j] ;
    size = cam_list[cam_shift.num[i_c]]->val_calc_hg1.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->val_calc_hg1[j] = s2val1[replica][j];
    size = cam_list[cam_shift.num[i_c]]->val_calc_hg2.size();
    for(int j=0;j<size;j++) cam_list[cam_shift.num[i_c]]->val_calc_hg2[j] = s2val2[replica][j];

    // calculate all the forces from the second set
    energy = cam_list[cam_shift.num[i_c]]->ens_calc_cs_force(coor, forces);

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

  }

  csfree_2dr_array_alloc(tala, repliche);
  csfree_2dr_array_alloc(tile1, repliche);
  csfree_2dr_array_alloc(tile2, repliche);
  csfree_2dr_array_alloc(tleu1, repliche);
  csfree_2dr_array_alloc(tleu2, repliche);
  csfree_2dr_array_alloc(tthr, repliche);
  csfree_2dr_array_alloc(tval1, repliche);
  csfree_2dr_array_alloc(tval2, repliche);

  if(cam_shift.mumo[i_c]==1) {
    csfree_2dr_array_alloc(s1ala, repliche);
    csfree_2dr_array_alloc(s1ile1, repliche);
    csfree_2dr_array_alloc(s1ile2, repliche);
    csfree_2dr_array_alloc(s1leu1, repliche);
    csfree_2dr_array_alloc(s1leu2, repliche);
    csfree_2dr_array_alloc(s1thr, repliche);
    csfree_2dr_array_alloc(s1val1, repliche);
    csfree_2dr_array_alloc(s1val2, repliche);
    csfree_2dr_array_alloc(s2ala, repliche);
    csfree_2dr_array_alloc(s2ile1, repliche);
    csfree_2dr_array_alloc(s2ile2, repliche);
    csfree_2dr_array_alloc(s2leu1, repliche);
    csfree_2dr_array_alloc(s2leu2, repliche);
    csfree_2dr_array_alloc(s2thr, repliche);
    csfree_2dr_array_alloc(s2val1, repliche);
    csfree_2dr_array_alloc(s2val2, repliche);
  }

  readed[i_c]++;
#ifdef PLUMED_GROMACS
  colvar.ss0[i_c] = (real) energy*4.186;
#else
  colvar.ss0[i_c] = (real) energy;
#endif
#else
  plumed_error("IN ORDER TO USE CAMSHIFT3-METHYL SHIFTS ENSEMBLE YOU MUST COMPILE GROMACS WITH MPI!");
#endif
}

extern "C" int PREFIX read_cs3meth(char **word, int count, t_plumed_input *input,FILE *fplog)
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
    colvar.type_s[count] = 63;
    iw=seek_word(word,"WTE");
    if(iw>=0) cam_shift.mumo[count]=2;
    iw=seek_word(word,"MUMO");
    if(iw>=0) cam_shift.mumo[count]=1;
  } else colvar.type_s[count] = 62;

  iw=seek_word(word,"DATA");
  if(iw>=0) {
       sscanf(word[iw+1],"%s", cam_shift.cam_data[count]);
  } else {
    fprintf(fplog,"|- NEEDED DATA KEYWORD FOR CAMSHIFT3-METHYL SHIFTS\n");
    help=1;
  }

  iw=seek_word(word,"FF");
  if(iw>=0) {
       sscanf(word[iw+1],"%s", cam_shift.cam_ff[count]);
  } else {
    fprintf(fplog,"|- NEEDED FF KEYWORD FOR CAMSHIFT3-METHYL SHIFTS\n");
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
    if(colvar.type_s[count]==63) {
      fprintf(fplog,"|- NEEDED NRES KEYWORD FOR CAMSHIFT3-METHYL SHIFTS ENSEMBLE\n");
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
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR CAMSHIFT3-METHYL SHIFTS\n"); help=1;} 

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

   if(help){
         fprintf(fplog,"|- CAMSHIFT3-METHYL SHIFTS SYNTAX:\n");
         fprintf(fplog,"|- DATA                   : parameters directory\n");
         fprintf(fplog,"|- LIST                   : group with all the atoms involved\n");
         fprintf(fplog,"|- SIGMA              : hills width for this cv\n");
         fprintf(fplog,"e.g. \n");
         fprintf(fplog,"CAMSHIFT3-METHYL SHIFTS DATA data/ LIST <prot> SIGMA 0.1 \n");
         EXIT(); 
  }

  fprintf(fplog, "\n%1i-CAMSHIFT3-METHYL SHIFTS: DATA %s, %i ATOMS ", count+1, cam_shift.cam_data[count], colvar.natoms[count]);
  if(colvar.type_s[count]==63) fprintf(fplog, "ENSEMBLE AVERAGING ON ");
  if(cam_shift.mumo[count]==1) fprintf(fplog, "OVER A SUBSET OF REPLICAS (MUMO STYLE) ");
  if(cam_shift.mumo[count]==2) fprintf(fplog, "OVER A SUBSET OF REPLICAS (WTE STYLE) ");
  if(logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n");

  cam_shift.num[count] = num; 
  cs_snew(colvar.myder[count], colvar.natoms[count]);

  num++; 

  return colvar.natoms[count];
}

#endif
