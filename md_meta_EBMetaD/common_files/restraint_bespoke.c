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
#ifdef RECONMETAD  
#ifdef CVS

#include "metadyn.h"

void PREFIX bespoke_restraint(int i_c, struct mtd_data_s *mtd_data) {
  int i, j, k, m, n, nat, i_cv, ncolvar; 
  // We calculate all the bespoke collective coordinates in the first call
  ncolvar=colvar.nconst-colvar.nbespoke;
  if( i_c>ncolvar ){ return; }

  double cv_in[colvar.bespoke_ncv];
  double cv_out[colvar.nbespoke];
  double derivatives[colvar.bespoke_ncv*colvar.nbespoke]; 

  // Transfer the values of all the colvars to a local array
  for(i=0;i<colvar.bespoke_ncv;i++){ cv_in[i]=colvar.ss0[ colvar.bespoke_cvlist[i] ]; }

  // This is the wrapper to the c++
  calculate_bespoke_cvs( colvar.bespoke_ncv, colvar.nbespoke, cv_in, cv_out, derivatives, mybespokeObj );  

  k=0; m=0; 
  for(i=ncolvar;i<colvar.nconst;i++){ 
     colvar.ss0[i]=cv_out[k]; k++; nat=0; 
     for(j=0;j<colvar.bespoke_ncv;j++){
        i_cv=colvar.bespoke_cvlist[j]; 
        for(n=0;n<colvar.natoms[i_cv];n++){
           colvar.myder[i][nat][0]=derivatives[m]*colvar.myder[i_cv][n][0]; 
           colvar.myder[i][nat][1]=derivatives[m]*colvar.myder[i_cv][n][1];
           colvar.myder[i][nat][2]=derivatives[m]*colvar.myder[i_cv][n][2];      
           nat++;
        }
        m++; 
     }
  }
}

int PREFIX read_bespoke( char **word, int count, t_plumed_input *input, FILE *fplog ) {
  int i, j, k, iw, ncolvar, ngrid, nspline;  int* cvlist;
  double smooth, delta = 0.0; int help=0;

  ncolvar=0; cvlist=NULL;

  // This reads in all the CVs used for bespoke collective coordinates
  iw = seek_word(word,"CV_LIST"); 
  if(iw>=0){  ncolvar=plumed_get_group(word[iw+1],&cvlist,0,input,fplog); 
             if( ncolvar>count ){ plumed_error("NUMBER OF CVS USED TO CALCULATE BESPOKE CVS IS LARGER THAN NUMBER OF CVS IN INPUT FILE"); }
  } 
  else{ fprintf(fplog,"|- NEEDED CV_LIST KEYWORD FOR BESPOKE\n"); help=1;}

  if(colvar.nbespoke==0){
     colvar.bespoke_ncv=ncolvar; snew( colvar.bespoke_cvlist, colvar.bespoke_ncv );
     for(i=0;i<colvar.bespoke_ncv;++i){ colvar.bespoke_cvlist[i]=cvlist[i]; }
  // This is the sanity check for CV_LIST ( part 1 : check all the bespoke collective coordinates use the same CV_LIST )
  } else if( ncolvar!=colvar.bespoke_ncv ){
     fprintf(fplog,"NUMBER OF CVS USED TO CREATE EACH BESPOKE CV MUST BE THE SAME \n"); help=1;
  } else{
     for(i=0;i<colvar.bespoke_ncv;++i){
        if(cvlist[i]!=colvar.bespoke_cvlist[i]){ fprintf(fplog,"CVS USED MUST BE THE SAME FOR ALL BESPOKE CVS \n"); help=1; }
     }
  }
  // This is the sanity check for CV_LIST ( part 2 : check all the cvs used to construct bespoke cvs are declared before any BESPOKE cvs)
  for(i=0;i<colvar.bespoke_ncv;i++){
     if( colvar.bespoke_cvlist[i]>count ){fprintf(fplog,"BESPOKE COMMAND MUST COME AFTER ALL CVS USED TO CREATE BESPOKE CV \n"); help=1; }
     // Number of atoms involved in this CV (sum of number of atoms involved in all other cvs)
     colvar.natoms[count]+=colvar.natoms[ colvar.bespoke_cvlist[i] ];
  }

  // N.B. Please note that if some idiot tries to use the different filenames for different CVs this will be ignored.
  //      This is because this idiot doesn't know how to compare char* s
  iw=seek_word(word,"FILENAME");
  if(iw>=0){ sscanf(word[iw+1],"%s", bespoke_input.filename); }
  else{ fprintf(fplog,"| - NEEDED FILENAME KEYWORD FOR BESPOKE\n"); help=1; }

  // Sanity check for FILENAME missing because I don't know how to compare char*s

  iw = seek_word(word,"SMOOTH_PARAM");
  if(iw>=0) { sscanf(word[iw+1], "%lf", &smooth); 
             if(colvar.nbespoke>0 && smooth!=bespoke_input.smooth){fprintf(fplog,"PARAMETERS FOR ALL BESPOKE CVS SHOULD BE THE SAME \n"); help=1; }
             else{ bespoke_input.smooth = (real) smooth; }
  }
  else{ fprintf(fplog,"|- NEEDED SMOOTH_PARAM KEYWORD FOR BESPOKE\n"); help=1; }

  iw = seek_word(word,"NSPLINE");
  if(iw>=0) { sscanf(word[iw+1], "%i", &nspline); bespoke_input.nspline[colvar.nbespoke]=nspline; }
  else{ fprintf(fplog,"|- NEEDED NSPLINE KEYWORD FOR BESPOKE\n"); help=1; }
 
  iw = seek_word(word,"NGRID");
  if(iw>=0) { sscanf(word[iw+1], "%i", &ngrid); bespoke_input.ngrid[colvar.nbespoke]=ngrid; }
  else{ fprintf(fplog,"|- NEEDED NGRID KEYWORD FOR BESPOKE\n"); help=1; }

  // Setup array of atoms involved in this cv
  snew( colvar.cvatoms[count], colvar.natoms[count] ); k=0;
  for(i=0;i<colvar.bespoke_ncv;i++){
     for(j=0;j<colvar.natoms[ colvar.bespoke_cvlist[i] ];j++){
        colvar.cvatoms[count][k]=colvar.cvatoms[ colvar.bespoke_cvlist[i] ][j];
        k++;
     }
  }
  // Setup derivatives array
  snew( colvar.myder[count], colvar.natoms[count] );

  if(help){ 
     fprintf(fplog,"\n-BESPOKE CV: WRONG SYNTAX\n");
     fprintf(fplog,"e.g.:    \n");
     fprintf(fplog,"TORSION LIST <g1> <g2> <g3> <g4>    \n");
     fprintf(fplog,"TORSION LIST <g2> <g3> <g4> <g5>    \n");
     fprintf(fplog,"TORSION LIST <g3> <g4> <g5> <g6>    \n");
     fprintf(fplog,"TORSION LIST <g4> <g5> <g6> <g7>    \n");
     fprintf(fplog,"BESPOKE CV_LIST <cv_list> FILENAME bespoke.defs SMOOTH_PARAM 1E-3 NGRID 100 NSPLINE 10 \n");
     fprintf(fplog,"BESPOKE CV_LIST <cv_list> FILENAME bespoke.defs SMOOTH_PARAM 1E-3 NGRID 100 NSPLINE 10 \n");
     fprintf(fplog,"cvlist-> \n");
     fprintf(fplog,"1 2 3 4  \n");
     fprintf(fplog,"cvlist<- \n");
     fprintf(fplog,"         \n"); 
     plumed_error("PluMed dead with errors: check log file");
  }

  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

  fprintf(fplog, "%1i-BESPOKE COLLECTIVE COORDINATE CREATED FROM %i OTHER CVS \n", count+1, colvar.bespoke_ncv);
  fprintf(fplog, "|--PARAMETERS: NUM. POINTS FOR SPLINE = %i NUM. POINTS FOR INTEGRAL = %i \n", bespoke_input.nspline[colvar.nbespoke], bespoke_input.ngrid[colvar.nbespoke] );
  fprintf(fplog, "|--INPUT FILE: %s \n",bespoke_input.filename );
  fprintf(fplog,"|- COLVARS USED IN BESPOKE CV:");
  for(i=0;i<colvar.bespoke_ncv;i++){fprintf(fplog," %d ",colvar.bespoke_cvlist[i]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");

  return colvar.natoms[count];
}

// A separate test derivatives routine

void PREFIX bespoke_test_derivatives(struct mtd_data_s *mtd_data)
{ 
  real teststep = 0.00001, invstep = 100000. ; 
  double testforce[colvar.nbespoke], analder[colvar.nbespoke];
  double val1[colvar.nbespoke], val2[colvar.nbespoke];
  int ncolvar, i_c, i_oc, i_bc, ix, i, j, k, iat, it;

  ncolvar=colvar.nconst-colvar.nbespoke;

  for(i=0;i<colvar.natoms[ncolvar];i++) {   // Loop over all atoms in the CV 
      iat = colvar.cvatoms[ncolvar][i];
      for(ix=0;ix<3;ix++) {              // Loop over all coordinates 
          mtd_data->pos[iat][ix] += teststep;
          for(i_bc=0;i_bc<colvar.bespoke_ncv;i_bc++) { // Loop over all cvs used to construct bespoke cv
              i_c=colvar.bespoke_cvlist[i_bc];
              switch(colvar.type_s[i_c]){
                case 1: dist_restraint(i_c, mtd_data); break;
                case 2: mindist_restraint(i_c, mtd_data); break;
                case 3: coord_restraint(i_c, mtd_data); break;
                case 4: angle_restraint(i_c, mtd_data); break;
                case 5: torsion_restraint(i_c, mtd_data); break;
                case 6: alfabeta_restraint(i_c, mtd_data); break;
                case 7: hbonds_restraint(i_c, mtd_data); break;
                case 8: dipole_restraint(i_c, mtd_data); break;
                case 11: radgyr_restraint(i_c, mtd_data); break;
                case 16: dihcor_restraint(i_c, mtd_data); break;
                case 20: waterbridge_restraint(i_c, mtd_data); break;
                case 30: spath_restraint(i_c, mtd_data); break;
                case 31: zpath_restraint(i_c, mtd_data); break;
                case 32: position_restraint(i_c, mtd_data); break;
                case 33: elstpot_restraint(i_c, mtd_data); break;
                case 34: puckering_restraint(i_c, mtd_data); break;
                case 36: helix_restraint(i_c, mtd_data); break;
                case 37: alpharmsd_restraint(i_c, mtd_data); break;
                case 38: antibetarmsd_restraint(i_c, mtd_data); break;
                case 39: parabetarmsd_restraint(i_c, mtd_data); break;
                //case 40: camshift_restraint(i_c, mtd_data); break;
                case 45: cmap_restraint(i_c, mtd_data); break;
                case 47: rdf_restraint(i_c, 1, mtd_data); break;
                case 52: adf_restraint(i_c, 0, mtd_data); break;
              }
          }
          bespoke_restraint( ncolvar, mtd_data ); k=0;
          for(i_oc=ncolvar;i_oc<colvar.nconst;i_oc++){ val1[k] = colvar.ss0[ i_oc ]; k++; }

          mtd_data->pos[iat][ix] += -2.*teststep;
          for(i_bc=0;i_bc<colvar.bespoke_ncv;i_bc++) {
              i_c=colvar.bespoke_cvlist[i_bc];
              switch(colvar.type_s[i_c]){
               case 1: dist_restraint(i_c, mtd_data); break;
               case 2: mindist_restraint(i_c, mtd_data); break;
               case 3: coord_restraint(i_c, mtd_data); break;
               case 4: angle_restraint(i_c, mtd_data); break;
               case 5: torsion_restraint(i_c, mtd_data); break;
               case 6: alfabeta_restraint(i_c, mtd_data); break;
               case 7: hbonds_restraint(i_c, mtd_data); break;
               case 8: dipole_restraint(i_c, mtd_data); break;
               case 11: radgyr_restraint(i_c, mtd_data); break;
               case 16: dihcor_restraint(i_c, mtd_data); break;
               case 20: waterbridge_restraint(i_c, mtd_data); break;
               case 30: spath_restraint(i_c, mtd_data); break;
               case 31: zpath_restraint(i_c, mtd_data); break;
               case 32: position_restraint(i_c, mtd_data); break;
               case 33: elstpot_restraint(i_c, mtd_data); break;
               case 34: puckering_restraint(i_c, mtd_data); break;
               case 36: helix_restraint(i_c, mtd_data); break;
               case 37: alpharmsd_restraint(i_c, mtd_data); break;
               case 38: antibetarmsd_restraint(i_c, mtd_data); break;
               case 39: parabetarmsd_restraint(i_c, mtd_data); break;
               //case 40: camshift_restraint(i_c, mtd_data); break;
               case 45: cmap_restraint(i_c, mtd_data); break;
               case 47: rdf_restraint(i_c, 1, mtd_data); break;
               case 52: adf_restraint(i_c, 0, mtd_data); break;
              }
          }
          bespoke_restraint(ncolvar, mtd_data ); k=0;
          for(i_oc=ncolvar;i_oc<colvar.nconst;i_oc++){ val2[k] = colvar.ss0[ i_oc ]; k++; }

          for(k=0;k<colvar.nbespoke;k++){ testforce[k] = 0.5*((val1[k]*invstep)-(val2[k]*invstep)); }
          mtd_data->pos[iat][ix] += teststep;
          for(i_bc=0;i_bc<colvar.bespoke_ncv;i_bc++) {
              i_c=colvar.bespoke_cvlist[i_bc];
              switch(colvar.type_s[i_c]){
                case 1: dist_restraint(i_c, mtd_data); break;
                case 2: mindist_restraint(i_c, mtd_data); break;
                case 3: coord_restraint(i_c, mtd_data); break;
                case 4: angle_restraint(i_c, mtd_data); break;
                case 5: torsion_restraint(i_c, mtd_data); break;
                case 6: alfabeta_restraint(i_c, mtd_data); break;
                case 7: hbonds_restraint(i_c, mtd_data); break;
                case 8: dipole_restraint(i_c, mtd_data); break;
                case 11: radgyr_restraint(i_c, mtd_data); break;
                case 16: dihcor_restraint(i_c, mtd_data); break;
                case 20: waterbridge_restraint(i_c, mtd_data); break;
                case 30: spath_restraint(i_c, mtd_data); break;
                case 31: zpath_restraint(i_c, mtd_data); break;
                case 32: position_restraint(i_c, mtd_data); break;
                case 33: elstpot_restraint(i_c, mtd_data); break;
                case 34: puckering_restraint(i_c, mtd_data); break;
                case 36: helix_restraint(i_c, mtd_data); break;
                case 37: alpharmsd_restraint(i_c, mtd_data); break;
                case 38: antibetarmsd_restraint(i_c, mtd_data); break;
                case 39: parabetarmsd_restraint(i_c, mtd_data); break;
                //case 40: camshift_restraint(i_c, mtd_data); break;
                case 45: cmap_restraint(i_c, mtd_data); break;
                case 47: rdf_restraint(i_c, 1, mtd_data); break;
                case 52: adf_restraint(i_c, 0, mtd_data); break;
              }
          }
          bespoke_restraint( ncolvar, mtd_data );

          
          for(k=0;k<colvar.nbespoke;k++){ analder[k]=0; }
      
          for(j=0;j<colvar.natoms[ncolvar];j++) {
              if(colvar.cvatoms[ncolvar][j]==iat) {
                k=0; for(i_oc=ncolvar;i_oc<colvar.nconst;i_oc++){ analder[k]+=colvar.myder[i_oc][j][ix]; k++; }
              }
          }
 
          for(k=0;k<colvar.nbespoke;k++){
            printf("Bespoke CV : %i Force atom %5i[%i] ** analytic %15.10f ** numeric %15.10f *** DELTA %15.10f\n",
              k+colvar.bespoke_ncv+1, iat, ix, analder[k], testforce[k], analder[k]-testforce[k]); 
          }
      }
  }

}

#endif
#endif
