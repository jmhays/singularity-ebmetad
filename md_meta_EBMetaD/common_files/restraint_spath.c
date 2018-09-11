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

void  PREFIX spath_restraint(int i_c, struct mtd_data_s *mtd_data) {

        int iat,i,ii,j;
        real s,ci_vec,tmp1;
        struct coordinates_frameset *pmy_coord1;
        struct hybrid_frameset *hbd_pmy;
        struct sz_data *pmy_sz;
        struct cmap_inpack inpack;
        struct cmap_outpack work;
	real ds_temp_dr0[MAXFRAMES_PATH][3][MAXATOMS_PATH];
	real ds_dcm[MAXFRAMES_PATH][MAXDIM_CMAP]; 
	real ds_dr0[3][MAXATOMS_PATH]; 
        int start_avg = 0; 
        real ds_dr1[MAXFRAMES_PATH][3][MAXATOMS_PATH];
        real dmsd_dr1[3][MAXATOMS_PATH];
        int  tot_con;
        real *save_err; 
        int nneigh;

        pmy_sz=&my_sz_list[ic_to_sz[i_c]];

        // bernd alternative indexing:
        // a separate implementation is preferred so to avoid strange quirks
        if(pmy_sz->indexing_type==1){
           // call bernd's routine: pass the pointer to the structure  
           sbernd_restraint(i_c,mtd_data);
           return ;  
        }
 

// neigh list ?
        if((pmy_sz->neigh==1 && colvar.it%pmy_sz->neigh_time==0) || firstTime) {
            //fprintf(mtd_data->fplog,"|- CALCULATING NEIGHBOUR LIST AT STEP %d\n",colvar.it);
            for(i=0;i< pmy_sz->number;i++)pmy_sz->lneigh[i]=i;
               nneigh=pmy_sz->number;
               save_err=(real *)malloc(pmy_sz->number*sizeof(real));
        }else {
            nneigh=pmy_sz->nneigh;
        }

// delete vectors
        for (i=0;i<3;i++){
                for (j=0;j<colvar.natoms[i_c];j++) {
                        ds_dr0[i][j]=0.;
                }
        }
        ci_vec=0.;
        s=0.;


        if(pmy_sz->umb_on && (colvar.it-pmy_sz->umblagsteps>=0)) start_avg = 1;
                                                               
        if(strcmp(pmy_sz->path_type,"HYBRID") != 0){
               for (i=0;i<colvar.natoms[i_c];i++){
                     iat = colvar.cvatoms[i_c][i];
                     inpack.r0[i][0] = mtd_data->pos[iat][0];
                     inpack.r0[i][1] = mtd_data->pos[iat][1];
                     inpack.r0[i][2] = mtd_data->pos[iat][2];
              }
        }
        if(strcmp(pmy_sz->path_type,"CMAP") == 0){ 
         tot_con=pmy_sz->my_cmap_pack.number+pmy_sz->my_cmap_pack.gnumber;
         cmap_running(i_c, &inpack,&pmy_sz->my_cmap_pack);
        }
        
        if(strcmp(pmy_sz->path_type,"HYBRID") == 0){ 
            //retrive values of cv and derivatives (one set)
            // assuming they are all the same for all the framesets (generally it is) 
            // in case of distance where the comp weight is on the cv evaluation and not on the 
            // distance acquire the points in CV space 
            // eg: dist between two atoms-> collect the cvval (which is the coordinate) and the  
            // eg: rmsd between two struct-> collect the actual coordinates (which is the coordinate) 
            // note: if D is the distance function from the reference  
            // d D(r_0,r)/dx= dD(r_0,r)/dr dr/dx 
            //   
            // in case of rmsd there is no need for chain rule as the derivative is directly calculated with quaternion
            // d D(r_0,r)/dx via quaternion   (it's like dr/dx=1)
            //
            //
			
            hbd_collect_config(pmy_sz->hbd_running);  // this collect cv values and needed derivatives

			hbd_collect_jacobian(pmy_sz->hbd_running,pmy_sz->mathybrid,pmy_sz->myinverse,mtd_data->fplog,1,1);

			//calculates the projector
			
			//calc_projector_test  ( pmy_sz );
			
			// simple case of debugging of the metrics (need its own module)

			if(pmy_sz->debug_metrics)test_hbd_metrics_new(pmy_sz->hbd_running,pmy_sz->hbd_frameset[0],&work,pmy_sz->mathybrid,mtd_data->fplog);

			
        }
        for(ii=0;ii< nneigh;ii++){
                i=pmy_sz->lneigh[ii];
       
                if(strcmp(pmy_sz->path_type,"CMAP") == 0){
                 cmdist_eval(i_c, i,&inpack,&work,&pmy_sz->my_cmap_pack,start_avg); 
                } 
                if(strcmp(pmy_sz->path_type,"MSD") == 0){
                   pmy_coord1=pmy_sz->frameset[i];
                   msd_calculation(pmy_coord1,&inpack,&work,dmsd_dr1,pmy_sz->umb_on,pmy_sz->norot_on,pmy_sz->nocenter_on);
                }
                if(strcmp(pmy_sz->path_type,"DMSD") == 0){
                 pmy_coord1=pmy_sz->frameset[i];
                 dmsd_calculation(i_c,pmy_coord1,&inpack,&work,dmsd_dr1);
                }
                if(strcmp(pmy_sz->path_type,"HYBRID") == 0){
                    hbd_pmy=pmy_sz->hbd_frameset[i];
                    // calculate the distance between the two frames
					//legacy code:
                    //hbd_metrics(&pmy_sz->hbd_running,hbd_pmy,&work,pmy_sz->mathybrid);
					//fprintf(mtd_data->fplog,"FRAME I= %3d\n",i);
					//if(nneigh>5)check_hbd_vecmvec_ref(pmy_sz->hbd_running,pmy_sz->hbd_frameset[0],pmy_sz->hbd_frameset[1],pmy_sz->hbd_frameset[2],pmy_sz->hbd_frameset[3],pmy_sz->hbd_frameset[4],pmy_sz->mathybrid,mtd_data->fplog);
				    hbd_metrics_new(pmy_sz->hbd_running,pmy_sz->hbd_frameset[i],&work,pmy_sz->mathybrid,mtd_data->fplog);
					//EXIT();
                }
 
               // fprintf(mtd_data->fplog,"ERR %d %f \n",i,work.err);
               // fflush(mtd_data->fplog);
               if((pmy_sz->neigh==1 && colvar.it%pmy_sz->neigh_time==0) || firstTime ) save_err[i]=work.err;

             // sqrt option
                if(pmy_sz->sqrt_on){
                  // if(work.err<1.e-6){
                  //  char buf[1024];
                  //  sprintf(buf,"PATH. Too small error: %f",work.err);
                  //  plumed_error(buf);
                  // }
		   if(pmy_sz->targeted_on){
                       pmy_sz->lambda=1./sqrt(work.err);
                       if(work.err<1.e-15) pmy_sz->lambda=1.e-15;
                   }
                   if(work.err<1.e-15){
                        tmp1=0.;
                   }else{
                        tmp1=(0.5/sqrt(work.err))*exp(-pmy_sz->lambda*sqrt(work.err));
                   }
                   ci_vec+=exp(-pmy_sz->lambda*sqrt(work.err));
                   s+=(i+1.0)*exp(-pmy_sz->lambda*sqrt(work.err));
                }else{
                   //if(pmy_sz->lambda*outpack.err>1.e2)outpack.err=1.e2/pmy_sz->lambda; 
                   tmp1=exp(-pmy_sz->lambda*work.err);
                   ci_vec+=tmp1;
                   s+=(i+1.0)*tmp1;
                }
 
				for(j=0;j<colvar.natoms[i_c];j++){
					ds_temp_dr0[i][0][j]=(work.derr_dr0[0][j])*tmp1;
					ds_temp_dr0[i][1][j]=(work.derr_dr0[1][j])*tmp1;
					ds_temp_dr0[i][2][j]=(work.derr_dr0[2][j])*tmp1;
					if((strcmp(pmy_sz->path_type,"MSD") == 0 || strcmp(pmy_sz->path_type,"DMSD") == 0) && start_avg){
						ds_dr1[i][0][j]=(dmsd_dr1[0][j])*tmp1*pmy_sz->lambda;
						ds_dr1[i][1][j]=(dmsd_dr1[1][j])*tmp1*pmy_sz->lambda;
						ds_dr1[i][2][j]=(dmsd_dr1[2][j])*tmp1*pmy_sz->lambda;
					}
				}
			
                if(strcmp(pmy_sz->path_type,"CMAP") == 0 && start_avg){
                 for(j=0;j<tot_con;j++){
                         ds_dcm[i][j]=(work.derr_dcm[j])*tmp1*pmy_sz->lambda;
                        } 
                }
	}

        s=s/ci_vec;

	for(j=0;j<colvar.natoms[i_c];j++){
               for(ii=0;ii<nneigh;ii++){
                       i=pmy_sz->lneigh[ii];
        	       ds_dr0[0][j]+=(s-i-1.)*ds_temp_dr0[i][0][j]; 
        	       ds_dr0[1][j]+=(s-i-1.)*ds_temp_dr0[i][1][j]; 
        	       ds_dr0[2][j]+=(s-i-1.)*ds_temp_dr0[i][2][j]; 
                       if((strcmp(pmy_sz->path_type,"MSD") == 0 || strcmp(pmy_sz->path_type,"DMSD") == 0) && start_avg){ 
                        ds_dr1[i][0][j]=ds_dr1[i][0][j]*(s-i-1.)/ci_vec;
                        ds_dr1[i][1][j]=ds_dr1[i][1][j]*(s-i-1.)/ci_vec;
                        ds_dr1[i][2][j]=ds_dr1[i][2][j]*(s-i-1.)/ci_vec;

                       }
       	    	}
                ds_dr0[0][j]=((pmy_sz->lambda)/ci_vec)*ds_dr0[0][j]; 
       	    	ds_dr0[1][j]=((pmy_sz->lambda)/ci_vec)*ds_dr0[1][j]; 
       	     	ds_dr0[2][j]=((pmy_sz->lambda)/ci_vec)*ds_dr0[2][j]; 
	}

        if(strcmp(pmy_sz->path_type,"CMAP") == 0 && start_avg){
           for(i=0;i<pmy_sz->number;i++){           
               for(j=0;j<tot_con;j++){
          	        ds_dcm[i][j]=ds_dcm[i][j]*(s-i-1.)/ci_vec;
		       }
           }
        }

        for(i=0;i<colvar.natoms[i_c];i++) {
          colvar.myder[i_c][i][0] = ds_dr0[0][i];
          colvar.myder[i_c][i][1] = ds_dr0[1][i];
          colvar.myder[i_c][i][2] = ds_dr0[2][i];
			//fprintf(mtd_data->fplog,"V %12.6f %12.6f %12.6f\n",colvar.myder[i_c][i][0],colvar.myder[i_c][i][1],colvar.myder[i_c][i][2]); 

        }

        colvar.ss0[i_c]=s;

#ifdef PATHREF_FINDIFF
          for(j=0;j<colvar.natoms[i_c];j++){
                 for(ii=0;ii<pmy_sz->number;ii++){
                      pmy_sz->dpath_dr[0][j][ii]=0.; 
                      pmy_sz->dpath_dr[1][j][ii]=0.; 
                      pmy_sz->dpath_dr[2][j][ii]=0.; 
                 }
                 for(ii=0;ii<nneigh;ii++){
                      i=pmy_sz->lneigh[ii];
                      pmy_sz->dpath_dr[0][j][i]=ds_dr1[i][0][j]; 
                      pmy_sz->dpath_dr[1][j][i]=ds_dr1[i][1][j]; 
                      pmy_sz->dpath_dr[2][j][i]=ds_dr1[i][2][j]; 
  
                 } 
          } 
#endif  


// neigh list? do quicksort and deallocate save_err
        if((pmy_sz->neigh==1 && colvar.it%pmy_sz->neigh_time==0) || firstTime){
             //   for(i=0;i<nneigh;i++)printf("BEFORE SORTING %d %f\n",pmy_sz->lneigh[i],save_err[i]);
                realquicksort(save_err,pmy_sz->lneigh,0,nneigh-1);
             //   for(i=0;i<nneigh;i++)printf("AFTER SORTING %d %f\n",pmy_sz->lneigh[i],save_err[i]);
                free(save_err);
        }



        if(pmy_sz->umb_on==1){
          if(strcmp(pmy_sz->path_type,"CMAP") == 0) mean_map(pmy_sz,ds_dcm,i_c,mtd_data->fplog); 
          if(strcmp(pmy_sz->path_type,"MSD") == 0 || strcmp(pmy_sz->path_type,"DMSD") == 0) mean_rmsd(pmy_sz,ds_dr1,i_c,mtd_data->fplog); 
        }
        
        return;
}

// ------------------------------------------------------------------------------------------------

int PREFIX read_path(char **word, int count,t_plumed_input *input,FILE *fplog)
{

  int i,iw;
  double lambda, tol;
  double sigma = 0.0;
  char file_maps[129];
  char file_maps2[129];
  char file_group[129];
  char type[2];
  struct sz_data *my_sz;
  int path_help;
  int neigh;

  path_help=0; 
  neigh=0;

  my_sz=&(my_sz_list[nsz]);
  ic_to_sz[count]=nsz;

  my_sz->number = 0;
  my_sz->lambda = 0.;                
  my_sz->my_cmap_pack.logical_group = 0;
  my_sz->umb_on = 0;
  my_sz->grad_on = 0;
  my_sz->targeted_on = 0;
  my_sz->sqrt_on = 0;
  my_sz->norot_on = 0;
  my_sz->nocenter_on = 0;
  // indexing type 0: From Branduardi
  // indexing type 1: From Ensing
  my_sz->indexing_type = 0;
  my_sz->ievol = -1;
  my_sz->reset = 0;
  my_sz->fadefactor = 1.0;
  my_sz->debug_metrics = 0;
  my_sz->intraframe_dist=0;
  my_sz->intraframe_diff=0;

 
// targeted md ?
  iw=seek_word(word,"TARGETED");
  if(iw>=0) {
      fprintf(fplog,"|- TARGETED MD: only one frame needed \n"); 
      my_sz->targeted_on=1;
     // my_sz->sqrt_on=1; // targeted is always square rooted
  }
// specific sqrt keyword
  iw=seek_word(word,"SQRT");
  if(iw>=0) {
      fprintf(fplog,"|- SQRT enabled: the measure will be square rooted \n"); 
      my_sz->sqrt_on=1;
  }
// type
  iw=seek_word(word,"TYPE");
  if(iw>=0) {
      sscanf(word[iw+1],"%s",my_sz->path_type);
  } else {
      fprintf(fplog,"|- NEEDED \"TYPE\" KEWORD FOR PATH\n"); 
      path_help=1; 
  }
// common stuff
  // if targeted you dont need many frames 
  if(my_sz->targeted_on==0){
    
     iw=seek_word(word,"NFRAMES");
     if(iw>=0){
         sscanf(word[iw+1],"%i", &my_sz->number); 
     } else {
         fprintf(fplog,"|- NEEDED \"NFRAMES\" KEYWORD FOR PATH\n"); 
         path_help=1; 
     }
  } else { 
     my_sz->number=1;
  }
  iw=seek_word(word,"ENSING");
  if(iw>=0){
          my_sz->indexing_type=1;
  } 
  iw=seek_word(word,"I_EVOL");
  if(iw>=0){
          sscanf(word[iw+1],"%i", &my_sz->ievol); 
  } 
  iw=seek_word(word,"RESET");
  if(iw>=0){
          my_sz->reset=1; 
  } 
  iw=seek_word(word,"FADEFACTOR");
  if(iw>=0){
      //sscanf(word[iw+1],"%f", &my_sz->fadefactor);   
      my_sz->fadefactor=atof(word[iw+1]);
  }
  iw=seek_word(word,"DEBUGMETRICS");
  if(iw>=0){
		my_sz->debug_metrics=1;
  }

  // Bernd ensing path required to use a reset in the weights
  if(my_sz->reset==1 && my_sz->indexing_type!=1 ) plumed_error("ENSING KEYWORD NEEDED TO RESET WEIGHTS");
  // Bernd ensing path required to use a fadefactor
  if(my_sz->fadefactor!=1.0 && my_sz->indexing_type!=1 ) plumed_error("ENSING KEYWORD NEEDED TO USE FADEFACTOR");

  // if targeted  lambda can be simply 1.0 
  if(my_sz->targeted_on==0){
     iw=seek_word(word,"LAMBDA");
     if(iw>=0){ 
         sscanf(word[iw+1],"%lf", &lambda);
#ifdef STANDALONE
         lambda/=(mtd_data.ampli)*(mtd_data.ampli);
#endif
     } else {
         if( my_sz->indexing_type==0){
           fprintf(fplog,"|- NEEDED \"LAMBDA\" KEYWORD FOR PATH\n"); 
           path_help=1; 
         }
     }
  } else { 
     lambda=1.0;
  }
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &sigma);
             colvar.delta_r[count]  = (real) sigma; }
// rmsd parameters
  iw=seek_word(word,"FRAMESET");
  if(iw>=0) sscanf(word[iw+1],"%s", my_sz->names); 
// cmap parameters
  iw = seek_word(word,"INDEX");
  if(iw>=0) sscanf(word[iw+1],"%s", file_maps);
  iw = seek_word(word,"MAP");
  if(iw>=0) sscanf(word[iw+1],"%s", file_maps2);
  iw=seek_word(word,"GROUP");
  if(iw>=0) {
   sscanf(word[iw+1],"%s",file_group);
   my_sz->my_cmap_pack.logical_group = 1;
   }
// UMBRELLA stuff
  iw=seek_word(word,"UMB_LAG");
  if(iw>=0) {
    sscanf(word[iw+1],"%d",&(my_sz->umblagsteps));
    my_sz->umb_on=1; 
    my_sz->countperm=0;
  }
  iw=seek_word(word,"UMB_BLOCK");
  if(iw>=0) sscanf(word[iw+1],"%d",&(my_sz->umbblocksize)); 
  iw=seek_word(word,"UMB_STRIDE");
  if(iw>=0) sscanf(word[iw+1],"%d",&(my_sz->umbstride));
  iw=seek_word(word,"UMB_PERM");
  if(iw>=0) sscanf(word[iw+1],"%d",&(my_sz->umbpermanency));
  iw=seek_word(word,"UMB_TOL");
  if(iw>=0) sscanf(word[iw+1],"%lf",&tol);
  iw=seek_word(word,"NO_ROT");
  if(iw>=0) my_sz->norot_on=1;
  iw=seek_word(word,"NO_CENTER");
  if(iw>=0) my_sz->nocenter_on=1;



#ifdef PATHREF_FINDIFF
    my_sz->umblagsteps=0;
    my_sz->umbblocksize=10;
    my_sz->umb_on=1;
    my_sz->countperm=0;
    my_sz->umbstride=1;
    my_sz->umbpermanency=1000000;
    tol=0.001;
#endif
// NEIGHBOUR LIST STUFF
  my_sz->neigh=0;
  my_sz->lneigh=(int *)malloc( MAXFRAMES_PATH * sizeof(int)); 
  iw=seek_word(word,"NEIGHLIST");  
  if(iw>=0) {
       sscanf(word[iw+1],"%d",&(my_sz->neigh_time)); // time for list 
       sscanf(word[iw+2],"%d",&(my_sz->nneigh));     // number of neighbours 
       if(my_sz->nneigh>=my_sz->number){
              my_sz->nneigh=my_sz->number;
              my_sz->neigh=0;
              for(i=0;i<my_sz->nneigh;i++)my_sz->lneigh[i]=i;
              neigh=1;
       }else {   
              my_sz->neigh=1;
              for(i=0;i<my_sz->nneigh;i++)my_sz->lneigh[i]=i;
              neigh=2;
       } 
 
  } else {
       my_sz->neigh=0;
       my_sz->nneigh=my_sz->number;
       for(i=0;i<my_sz->nneigh;i++)my_sz->lneigh[i]=i; 
  } 


  my_sz->lambda  = (real) lambda; 
  my_sz->umbtolerance = (real) tol;

  if(colvar.type_s[count]==30) strcpy(type,"S");
  if(colvar.type_s[count]==31) strcpy(type,"Z");

  if(path_help){
         fprintf(fplog,"|- PATH/TARGETED SYNTAX:\n");
         fprintf(fplog,"|- TYPE               : can be MSD/CMAP/DMSD\n");
         fprintf(fplog,"|- NFRAMES            : the number of reference structures\n");
         fprintf(fplog,"|- LAMBDA             : the common prefactor in the exponential equation for the path\n");
         fprintf(fplog,"|- SIGMA              : hills width for this cv\n");
         fprintf(fplog,"|- NEIGHLIST (opt.)   : neighlist on the closest frames    NEIGHLIST (ntimesteps) (nframes) \n");
         fprintf(fplog,"|- ENSING   (opt.)    : Bernd Ensing's path \n");
         fprintf(fplog,"|- I_EVOL   (opt.)    : update Bernd Ensing's path each ievol timesteps (only if ENSING is defined )\n");
         fprintf(fplog,"|- RESET    (opt.)    : reset weights to zero every Bernd Ensing's path evolution \n");
         fprintf(fplog,"|-                      (only if ENSING is defined )\n");
         fprintf(fplog,"|- FADEFACTOR (opt.)  : fadefactor for weights in Bernd Ensing's path (only if ENSING is defined) \n");
         fprintf(fplog,"|-                      Should be 0<=fadefactor<=1. Default=1.    \n");
         fprintf(fplog,"|- .....many other keywords are specific to MSD/CMAP/DMSD paths... \n");
         fprintf(fplog,"|- (MSD)   FRAMESET :base for frameset name pippo_ will look for pippo_1.pdb pippo_2.pdb etc...  \n");
         fprintf(fplog,"|-     \n");
         fprintf(fplog,"|- e.g.\n");
         fprintf(fplog,"|- \n");
         fprintf(fplog,"|- Z_PATH   TYPE MSD FRAMESET frame_ NFRAMES 1 LAMBDA 1.0 SIGMA 1.0 NEIGHLIST 10 5 \n");
         fprintf(fplog,"|- \n");
         fprintf(fplog,"|- S_PATH   TYPE MSD FRAMESET frame_ NFRAMES 1 LAMBDA 1.0 SIGMA 1.0 NEIGHLIST 10 5 \n");
         fprintf(fplog,"|- \n");
         fprintf(fplog,"|- TARGETED TYPE MSD FRAMESET frame.pdb SIGMA 1.0 \n");
         fprintf(fplog,"|- \n");
         fprintf(fplog,"|- S_PATH   TYPE MSD FRAMESET frame_ NFRAMES 10 ENSING  RESET FADEFACTOR 0.9 I_EVOL 10000 \n");
         fprintf(fplog,"|- \n");
 
         plumed_error("PluMeD dead with errors: check log file");
  } 
  if(strcmp(my_sz->path_type,"CMAP") != 0 && strcmp(my_sz->path_type,"MSD") != 0
     && strcmp(my_sz->path_type,"DMSD") !=0 && strcmp(my_sz->path_type,"HYBRID") !=0  && strcmp(my_sz->path_type,"CHIRAL") !=0 ){
          char buf[1024];
          sprintf(buf,"%s Unknown type of path!!",my_sz->path_type);
          plumed_error(buf);
  }
  fprintf(fplog,"\n%1i-%s_PATH in %s space \n",count+1,type,my_sz->path_type);
  fprintf(fplog,"|--NFRAMES %i ",my_sz->number);
  if(my_sz->indexing_type==0)fprintf(fplog,"|--LAMBDA %f ",my_sz->lambda);
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  if(strcmp(my_sz->path_type,"CMAP") == 0){
   fprintf(fplog,"|--READING CONTACT MAPS INDEX FROM FILE %s AND VALUES FROM %s\n",file_maps,file_maps2); 
   if(my_sz->my_cmap_pack.logical_group) fprintf(fplog,"|--AND GROUP FROM FILE %s \n",file_group);
  }
  if(strcmp(my_sz->path_type,"MSD") == 0 || strcmp(my_sz->path_type,"DMSD") == 0  || strcmp(my_sz->path_type,"HYBRID") == 0 || strcmp(my_sz->path_type,"CHIRAL") == 0 ){
   fprintf(fplog,"|--BASENAME FOR FRAMES %s \n",my_sz->names);
  }
	if(my_sz->debug_metrics)fprintf(fplog,"|- DEBUGGING THE METRICS \n");

  if(neigh==0) fprintf(fplog,"|--NEIGHBOUR LIST OFF \n");  
  if(neigh==1) fprintf(fplog,"|--NO NEED FOR NEIGHBOUR LIST: THE LIST APPEARS TOO LARGE \n");  
  if(neigh==2) {fprintf(fplog,"|--NEIGHBOUR LIST ON : TIME FOR LIST: %d TSTEPS \n", my_sz->neigh_time);
                fprintf(fplog,"|--                  : LIST LENGTH  : %d FRAMES \n", my_sz->nneigh); } 

  
  if(my_sz->number>MAXFRAMES_PATH)
                plumed_error("Maximum number of frames in path CV exceeded. Increase MAXFRAMES_PATH and recompile");

   if(my_sz->number==0)   plumed_error("NUMBER KEYWORD NEEDED FOR PATH");  
   
   // lambda needed only for traditional indexing style
   if(my_sz->lambda==0. && my_sz->indexing_type==0 ) plumed_error("LAMBDA KEYWORD NEEDED FOR PATH");
  

   if(my_sz->umb_on){
    fprintf(fplog,"|--MEAN EVALUATION ON %s IS ACTIVE\n",type);
    fprintf(fplog,"|---LAGSTEPS %d - BLOCKSTEPS %d  - STRIDESTEPS %d  \n",my_sz->umblagsteps,my_sz->umbblocksize,my_sz->umbstride);
    fprintf(fplog,"|---PERMSTEPS %d - TOLERANCESTEPS %lf \n",my_sz->umbpermanency,my_sz->umbtolerance);
   }
   if(my_sz->norot_on){
    fprintf(fplog,"|--NO ROTATION OF REFERENCE/RUNNING FRAME\n");
    if(strcmp(my_sz->path_type,"MSD") != 0) plumed_error("NO_ROT CAN BE USED ONLY WITH TYPE MSD");
   }
   if(my_sz->nocenter_on){
    fprintf(fplog,"|--NO COM CENTERING OF REFERENCE/RUNNING FRAME\n");
    if(strcmp(my_sz->path_type,"MSD") != 0) plumed_error("NO_CENTER CAN BE USED ONLY WITH TYPE MSD");
    if(!my_sz->norot_on) plumed_error("NO_CENTER CAN BE USED ONLY IF NO_ROT IS ACTIVE");
   }
   switch (my_sz->indexing_type) {
      case 0:
        fprintf(fplog,"|--TRADITIONAL INDEXING STYLE FOR PATHWAY  IS ON\n");
        break; 
      case 1:
        fprintf(fplog,"|--BERND ENSING PATHWAY INDEXING STYLE IS ON\n");
        if(my_sz->ievol>0){ 
           fprintf(fplog,"|--EVOLVE PATHWAY EACH %d TIMESTEPS\n",my_sz->ievol);
        }else{
           fprintf(fplog,"|--NEVER EVOLVE PATHWAY \n");
        }
        if(my_sz->fadefactor>=0.0 && my_sz->fadefactor<=1.0 ){ 
           fprintf(fplog,"|--FADEFACTOR IS %f\n",my_sz->fadefactor);
        }else{
           plumed_error("|--FADEFACTOR VALUE NOT VALID \n");
	}
        switch (my_sz->reset){
           case 1:
             fprintf(fplog,"|--RESETING THE WEIGHTS EVERY BERND'S PATH EVOLUTION\n");
             break; 
           case 0:
             fprintf(fplog,"|--NOT RESETING WEIGHTS DURING BERND PATH'S EVOLUTION\n");
             break; 
        } 
        break; 
      default:
        plumed_error("|--UNKNOWN INDEXING STYLE\n");
        break;
 
   } 



// case of CMAP path

 if(strcmp(my_sz->path_type,"CMAP") == 0){
   colvar.cell_pbc[count]=1; // default is PBC
   iw=seek_word(word,"NOPBC");
   if(iw>=0) {colvar.cell_pbc[count] = 0;}
   iw=seek_word(word,"PBC");
   if(iw>=0) {colvar.cell_pbc[count] = 1;}

   if(colvar.cell_pbc[count]) fprintf(fplog,"|--DISTANCES WITH PBC \n");
   else fprintf(fplog,"|--DISTANCES WITHOUT PBC \n");

   read_sz_map(my_sz,file_maps,file_maps2,file_group,1,fplog);

   colvar.natoms[count]   = my_sz->my_cmap_pack.atoms; 

   snew(colvar.myder[count], colvar.natoms[count]);
   snew(colvar.cvatoms[count], colvar.natoms[count]);

   for(i=0;i<colvar.natoms[count];i++){
     colvar.cvatoms[count][i] = my_sz->my_cmap_pack.list[i];
   }
 
   if(my_sz->umb_on){
     my_sz->umb_map_block=float_3d_array_alloc(my_sz->my_cmap_pack.number+my_sz->my_cmap_pack.gnumber,my_sz->number,my_sz->umbblocksize);
     my_sz->umb_map_avg=float_2d_array_alloc(my_sz->my_cmap_pack.number+my_sz->my_cmap_pack.gnumber,my_sz->number);
     fprintf(fplog,"|---ALLOCATED for UMBRELLA STUFF %lf Kbytes \n",sizeof(real)*((real) ((my_sz->my_cmap_pack.number+
                  my_sz->my_cmap_pack.gnumber)*my_sz->number*my_sz->umbblocksize)/1000));
   }
 }

// RMSD/DRMS case 

 if(strcmp(my_sz->path_type,"MSD") == 0 || strcmp(my_sz->path_type,"DMSD") == 0 ){

   read_sz_rmsd(my_sz,fplog);

   colvar.natoms[count]   = (my_sz->frameset[0])->natoms; 

   snew(colvar.myder[count], colvar.natoms[count]);
   snew(colvar.cvatoms[count], colvar.natoms[count]);

   for(i=0;i<colvar.natoms[count];i++){
     colvar.cvatoms[count][i] = ((my_sz->frameset[0])->atmnum[i])-1; 
   }

   if(my_sz->umb_on){
     my_sz->umb_block=float_4d_array_alloc(3,(*my_sz->frameset[0]).natoms,my_sz->number,my_sz->umbblocksize);
     my_sz->umb_avg=float_3d_array_alloc(3,(*my_sz->frameset[0]).natoms,my_sz->number);
#ifdef PATHREF_FINDIFF 
     fprintf(fplog,"|---ALLOCATING TEST ARRAYS \n");
     my_sz->dpath_dr=float_3d_array_alloc(3,(*my_sz->frameset[0]).natoms,my_sz->number);
     fprintf(fplog,"|---ALLOCATION DONE \n");
#endif
     fprintf(fplog,"|---ALLOCATED for UMBRELLA STUFF %lf Kbytes",sizeof(real)*((real) 3*((*my_sz->frameset[0]).natoms)*(my_sz->number)*(my_sz->umbblocksize)/1000));
   }
 }
// CHIRALITY
 if(strcmp(my_sz->path_type,"CHIRAL") == 0 ){

   read_sz_chiral(my_sz,fplog);
   EXIT();

 //  colvar.natoms[count]   = (my_sz->frameset[0])->natoms; 

//   snew(colvar.myder[count], colvar.natoms[count]);
//   snew(colvar.cvatoms[count], colvar.natoms[count]);

//   for(i=0;i<colvar.natoms[count];i++){
//     colvar.cvatoms[count][i] = ((my_sz->frameset[0])->atmnum[i])-1; 
//   }

 }


// HYBRID STRUCTURE: uses the definition of previous cvs
 if(strcmp(my_sz->path_type,"HYBRID") == 0 ){
    // read which cvs you want to hybridize  
     iw=seek_word(word,"HYBRID");
     if(iw>=0) {
        int c,ii,*t1l,t1s=0, *t2l,t2s,jj;
        ii=0;
        while (1) { 
             //      word[iw+1] is a word ???   
             ii++;
             c=*word[iw+ii];
             if (isalpha(c)!=0 ){
                 //printf("This is a word %s\n",word[iw+ii]); 
                 break; 
             } else {
       //          printf("This is a number %s\n",word[iw+ii]); 
                 if(t1s==0){ t1l=(int *)malloc(sizeof(int));
                             sscanf(word[iw+ii],"%d", &t1l[0]);
                 } else{ 
                        t2s=t1s;t2l=(int *)malloc(t2s*sizeof(int));
                        for(jj=0;jj<t1s;jj++){
                              t2l[jj]=t1l[jj];
                       }
                       free(t1l);
                       t1l=(int *)malloc((t1s+1)*sizeof(int)); 
                       for(jj=0;jj<t1s;jj++){
                              t1l[jj]=t2l[jj];
                       }
                       free(t2l);
                       sscanf(word[iw+ii],"%d", &t1l[t1s]);
                 } 
                 t1s++;
                 //for(jj=0;jj<t1s;jj++){printf("LL %d ",t1l[jj]);}; printf(" NN %d\n",t1s) ; 
             } 
        } 
		 
		//
        // allocate and copy a specific structure
        //
		my_sz->nhybrid=t1s; // number of hybrid functions
        my_sz->lhybrid=(int *)malloc(my_sz->nhybrid*sizeof(int));		// list of the ordinal of hybrid function to be called
        my_sz->lcvhybrid=(int *)malloc(my_sz->nhybrid*sizeof(int));		// list of the ordinal of the cv to be used
        for(ii=0;ii<my_sz->nhybrid;ii++){my_sz->lhybrid[ii]=t1l[ii];}   // list of the type of cv 
		 
        fprintf(fplog,"|--NUMBER OF HYBRID CVS= %d WHICH ARE ",my_sz->nhybrid);
        for(jj=0;jj<my_sz->nhybrid;jj++){
              fprintf(fplog,"CV %d ",my_sz->lhybrid[jj]);my_sz->lhybrid[jj]--;
              if(my_sz->lhybrid[jj]>=count) plumed_error("|--WRONG INPUT: PUT THE HYBRID CV AFTER THE ONES YOU USE TO MAKE THE HYBRID\n");
        };printf("\n"); 
		//
        // set the type of the representation for each cv  
		//
        for(jj=0;jj<my_sz->nhybrid;jj++){
            int kk=my_sz->lhybrid[jj];
            fprintf(fplog,"|- CVTYPE %d\n",colvar.type_s[kk]);my_sz->lcvhybrid[jj]=colvar.type_s[kk];
            // parse the input and check if it is implemented  
		}; 
		 
		//
        // matrix of metrics: some variable count one, some others count more
		// (dummy MSD is one cv per coordinate, just the difference and its derivative are interwined)
		//
		 
     	colvar.natoms[count]=read_sz_hybrid(my_sz,fplog); 
		 
        snew(colvar.myder[count], colvar.natoms[count]);
        snew(colvar.cvatoms[count], colvar.natoms[count]);

        // assuming they are all the same for all the framesets (generally it is) 
        for(i=0;i<colvar.natoms[count];i++){
          colvar.cvatoms[count][i] = my_sz->hbd_frameset[0]->backtable[i]  ;
			fprintf(fplog,"|- read_sz_hybrid: ATOM %d IS INVOLVED \n", my_sz->hbd_frameset[0]->backtable[i]);
        }
		 my_sz->mathybrid_blocks=NULL;
     } else {
		fprintf(fplog,"SOMETHING WENT WRONG \n");
     }  
	 // these are useful only when using an external script to retrieve various information
	 iw=seek_word(word,"INTRAFRAME_DIST");
     if(iw>=0) {
		 if(strcmp(my_sz->path_type,"HYBRID") != 0 ){
			 fprintf(fplog,"|- read_sz_path: INTRAFRAME_DIST is only for hybrid path \n");
		 }
		 // requires a set of matrices that should be equal to the number of frames
		 my_sz->intraframe_dist=1;
		 sscanf(word[iw+1],"%s",my_sz->intraframe_dist_inputfile); // the inputfile
		 sscanf(word[iw+2],"%s",my_sz->intraframe_dist_outputfile); // the outputfile

		 fprintf(fplog,"|- read_sz_path: INTRAFRAME_DIST is on\n");
		 fprintf(fplog,"|- read_sz_path: and the matrix file taken from external source is: %s\n",my_sz->intraframe_dist_inputfile);
		 fprintf(fplog,"|- read_sz_path: and the matrix file printed to : %s\n",my_sz->intraframe_dist_outputfile);
		 fprintf(fplog,"|- read_sz_path: CALLING AND THEN DYING\n");
		 calc_intraframe_dist( my_sz );
		 EXIT();
	 }
	 iw=seek_word(word,"REPARAM");
     if(iw>=0) {
		 if(strcmp(my_sz->path_type,"HYBRID") != 0 ){
			 fprintf(fplog,"|- read_sz_path: INTRAFRAME_DIST is only for hybrid path \n");
		 }
		 // assume only identity matrix as for reparametrization
		 // allocate jacobian
		 
		 // allocate needed structures	
		 init_bernd_evolution(&my_sz,count);

		// find maxiter, smooth and maxerror 
		 int iw3;
		 real smooth, maxerr;
		 int maxiter;
                 char myword[10];  
		 iw3=seek_word(word,"MAXERR");
		 if(iw3>=0) {
			         sscanf(word[iw3+1],"%s", myword);
				 maxerr=atof(myword);
                 } else { maxerr=0.9;}
		 iw3=seek_word(word,"MAXITER");
		 if(iw3>=0) {
			         sscanf(word[iw3+1],"%s", myword);
				 maxiter=atoi(myword);
                 } else { maxiter=30; }
		 iw3=seek_word(word,"SMOOTH");
		 if(iw3>=0) {
			         sscanf(word[iw3+1],"%s", myword);
				 smooth=atof(myword);
                 } else { smooth=1.0; }
		 fprintf(fplog,"|- read_sz_path: REPARAM SMOOTH %f\n",smooth);
		 fprintf(fplog,"|- read_sz_path: REPARAM MAXITER %d\n",maxiter);
		 fprintf(fplog,"|- read_sz_path: REPARAM MAXERR %f\n",maxerr);
	
		 
		 // if you have one matrix per point then do it differently
		 int iw2;
		 iw2=seek_word(word,"MATRICES");
		 
		 if(iw2>=0) {
			 char filenamematrices[200];
			 strcpy(filenamematrices,word[iw2+1]);
			 // call a specific routine that reads in the file, allocate and do the stuff
			 fprintf(fplog,"|- read_sz_path: REPARAM is on\n");
			 fprintf(fplog,"|- read_sz_path: each frame uses its own matrix from file %s\n",filenamematrices);

			 reparam_with_multiple_matrix(my_sz,filenamematrices,maxerr,maxiter,smooth);
			 
			 
			 // this closes the file that was opened by bernd evolution
			 fclose(my_sz->fp_evol);

		 }else{
			 // use a simple diagonal matrix instead 	 
			 int j;
			 
			 fprintf(fplog,"|- read_sz_path: FOUND %d CVS\n",my_sz->hbd_running->hbd_totcvs);
			 for (i=0;i<my_sz->hbd_running->hbd_totcvs;i++){
				 for (j=0;j<my_sz->hbd_running->hbd_totcvs;j++){
					 my_sz->mathybrid[i][j]=0.;
					 if(i==j)my_sz->mathybrid[i][j]=1.;
				 }
			 }
			 
			 reparam_bernd_path(my_sz,maxerr,maxiter,smooth);
			 
			 fclose(my_sz->fp_evol);
			 
		 }

		 
		 EXIT();
		 
	 }
	 iw=seek_word(word,"INTRAFRAME_DIFF");
     if(iw>=0) {
		 my_sz->intraframe_diff=1;
		 if(strcmp(my_sz->path_type,"HYBRID") != 0 ){
			 fprintf(fplog,"|- read_sz_path: INTRAFRAME_DIFF is only for hybrid path \n");EXIT();
		 }
		 sscanf(word[iw+1],"%s",my_sz->intraframe_diff_outputfile); // the outputfile
		 fprintf(fplog,"|- read_sz_path: INTRAFRAME_DIFF is on\n");
		 fprintf(fplog,"|- read_sz_path: and the diff file is printed to : %s\n",my_sz->intraframe_diff_outputfile);
		 fprintf(fplog,"|- read_sz_path: CALLING AND THEN DYING\n");
		 calc_intraframe_diff( my_sz );
		 EXIT();
	 }
	 iw=seek_word(word,"INTERPOLATE");
      	 if(iw>=0) {
		 if(strcmp(my_sz->path_type,"HYBRID") != 0 ){
			 fprintf(fplog,"|- read_sz_path: INTERPOLATE is only for hybrid path \n");
		 }
		 int npoints=atoi(word[iw+1]);
		 // allocate needed structures	
		 init_bernd_evolution(&my_sz,count);
		 
		 // if you have one matrix per point then do it differently
		 int iw2;
		 iw2=seek_word(word,"MATRICES");
		 
		 if(iw2>=0) {
			 char filenamematrices[200];
			 strcpy(filenamematrices,word[iw2+1]);
			 // call a specific routine that reads in the file, allocate and do the stuff
			 fprintf(fplog,"|- read_sz_path: REPARAM is on\n");
			 fprintf(fplog,"|- read_sz_path: each frame uses its own matrix from file %s\n",filenamematrices);

			 interpolate_with_multiple_matrix(my_sz,filenamematrices,npoints);
			 
			 fclose(my_sz->fp_evol);

		 }else{
			 // use a simple diagonal matrix instead 	 
			 int j;
			 
			 fprintf(fplog,"|- read_sz_path: FOUND %d CVS\n",my_sz->hbd_running->hbd_totcvs);
			 for (i=0;i<my_sz->hbd_running->hbd_totcvs;i++){
				 for (j=0;j<my_sz->hbd_running->hbd_totcvs;j++){
					 my_sz->mathybrid[i][j]=0.;
					 if(i==j)my_sz->mathybrid[i][j]=1.;
				 }
			 }
			 
			 interpolate_bernd_path(my_sz,npoints);
			 
			 fclose(my_sz->fp_evol);
			 
		 }
		 EXIT();
		 
	 }
	
	 if(iw>=0) {
		if(strcmp(my_sz->path_type,"HYBRID") != 0 ){
			fprintf(fplog,"|- read_sz_path: INTERPOLATE is only for hybrid path \n");EXIT();
		}
		int first=atoi(word[iw+1]);
		fprintf(fplog,"|- read_sz_path: INTERPOLATE chosen\n");
		fprintf(fplog,"|- read_sz_path: will calculate the interpolation for your frames into %d frames\n",first);
		first--;

		//EXIT();
	 }
	
	 iw=seek_word(word,"TWOFRAMES_DIFF");
	 if(iw>=0) {

		if(strcmp(my_sz->path_type,"HYBRID") != 0 ){
			fprintf(fplog,"|- read_sz_path: INTRAFRAME_DIFF is only for hybrid path \n");EXIT();
		}
		int first=atoi(word[iw+1]);
	    int second=atoi(word[iw+2]);
		int reference=atoi(word[iw+3]);
		fprintf(fplog,"|- read_sz_path: TWOFRAMES_DIFF chosen\n");
		fprintf(fplog,"|- read_sz_path: will calculate the diff beteween frame %d and %d takin frame %d as ref\n",first,second,reference);
		first--;
		second--;
		reference--;
		struct hybrid_frameset *f_first,*f_second,*f_reference,*difference;
		struct hybrid_elem *diff;
		clone_hybrid_frameset(&difference,my_sz->hbd_frameset[0],1,mtd_data.fplog);
		f_first=my_sz->hbd_frameset[first];
		f_second=my_sz->hbd_frameset[second];
		f_reference=my_sz->hbd_frameset[reference];
		calc_diff_twoframes(f_first,f_second,f_reference,difference);
		int j,k;
		for(j=0;j< difference->hbd_nelem; j++){
			 diff=difference->hbd_elem[j];
			 for(k=0;k<diff->ncv;k++){
				 fprintf(mtd_data.fplog,"|- TWOFRAMES DIFFERENCE ELEM %d CV %d IS %f\n",j,k,diff->ref_dist[k]);
			 }
		}

		//EXIT();
	 }
	 iw=seek_word(word,"TWOFRAMES_DIST");
	 if(iw>=0) {
		 
		 if(strcmp(my_sz->path_type,"HYBRID") != 0 ){
			 fprintf(fplog,"|- read_sz_path: INTRAFRAME_DIST is only for hybrid path \n");EXIT();
		 }
		 int first=atoi(word[iw+1]);
		 int second=atoi(word[iw+2]);
		 char matrixfile[200];
		 strcpy(matrixfile,word[iw+3]); //this is the filename for output 
		 fprintf(fplog,"|- read_sz_path: TWOFRAMES_DIST chosen\n");
		 fprintf(fplog,"|- read_sz_path: will calculate the dist beteween frame %d and %d \n",first,second);
		 fprintf(fplog,"|- read_sz_path: input matrix chosen %s\n",matrixfile);
		 first--;
		 second--;

		 calc_twoframe_dist(my_sz,first,second,matrixfile,mtd_data.fplog);

		 
		 //EXIT();
	 }
	 
//     EXIT(); 
 }
         nsz++;

         if(nsz>NMAX_PATH) plumed_error("TOO MANY PATH CVS. INCREASE NMAX_PATH and recompile");

         fprintf(fplog,"\n");
         return colvar.natoms[count];
}


// ------------------------------------------------------------------------------------------------

int PREFIX read_sz_rmsd(struct sz_data *my_sz, FILE *fplog) {

        int l,i,j,k,found;
        char *str,ic[3],str2[100];
        l=0;
/*
  * allocate the pointers
 */
        my_sz->frameset=(struct coordinates_frameset **)malloc((my_sz->number)*sizeof(struct coordinates_frameset *));
        for(i=0;i< my_sz->number;i++){
                my_sz->frameset[i]=(struct coordinates_frameset *)malloc(sizeof(struct coordinates_frameset)) ;
        }

/*
 *  build names
 */
        str=&(my_sz->names[0]);
        for (i=1;i<= my_sz->number ;i++){
                strcpy(str2,my_sz->names);
                if(my_sz->targeted_on==0){
                     if(i<10){
                      ic[0]='0'+i;
                      ic[1]='\0';}
                     else if(i<100) {
                      ic[0]='0'+i/10 ;
                      ic[1]='0'+i%10 ;
                      ic[2]='\0';
                     }
                     else{
                       plumed_error("|--read_sz_input: TOO MANY FRAMES REQUIRED FOR NAME BUILDING!");
                     }
                     strcat(str2,ic);
                     strcat(str2,".pdb");
                }
                fprintf(fplog,"|--%s\n",str2);
                read_sz_coord(str2,my_sz->frameset[i-1],fplog);
        }

/*
 * check over that aligned atoms are the same ( this is requirement for rmsd routine
 * so you don't have to reallocate everything everytime )
 */
       for (i=0;i< my_sz->number-1;i++){
        for (j=i+1;j< my_sz->number;j++){
         if( ( *my_sz->frameset[i] ) .nalign == ( *my_sz->frameset[j] ).nalign ){
            for(k=0;k<(*my_sz->frameset[i]).nalign;k++){
               found=0;
               for(l=0;l<(*my_sz->frameset[j]).nalign;l++){
                  if( (*my_sz->frameset[i]).align_to_frameset[k]==(*my_sz->frameset[j]).align_to_frameset[l] ){found++;}
               }
               if(found==0){fprintf(fplog,"|--ERROR: ATOM %d in frameset %d not found in frameset %d\n",(*my_sz->frameset[i]).align_to_frameset[k],i,j);EXIT();}
               else if(found>1){fprintf(fplog,"|--ERROR: found multiple definition of %d in frameset %d not found in frameset %d\n",(*my_sz->frameset[i]).align_to_frameset[k],i,j);EXIT();}            }
         }
         else{fprintf(fplog,"|--ERROR : ALIGNMENT ATOMS IN THE FRAMESET %d AND %d ARE NOT THE SAME\n",i,j);EXIT();};
        }
       }
// 
// now write the backtable align_to_coord
// 
       for(l=0;l<my_sz->number;l++){// for each frameset in the set
          j=0;
          for(i=0;i<(*my_sz->frameset[l]).natoms;i++){//on all the atoms
             if((*my_sz->frameset[l]).align[i]!=0.){ // only if this atom is used in the alignment
                  (*my_sz->frameset[l]).align_to_coord[j]=i;
                  j++;
             }
             (*my_sz->frameset[l]).frameset_to_coord[i]=i;
           }
       }
       return 0;
};

// ------------------------------------------------------------------------------------------------

int PREFIX read_sz_coord (char *filename, struct coordinates_frameset *p, FILE *fplog){


    char string[400],sm[10],sm1[10];
    char *str,remark[10],end[10],atom[5];
    // new parsing
    char x[12],y[12],z[12],occ[12],beta[12];
    char ind[10], resid[10];
    char name[10],resname[10],chain,hetatm[7] ;	
	
    real tmp0;
    FILE *fp;
    int i,l;

    fp=fopen(filename,"r");
    if (fp == NULL){
       char buf[1024];
       sprintf(buf,"UNSUCCESSFULL OPENING FOR FILE %s",filename);
       plumed_error(buf);
    }

/* PDB
 
         1         2         3         4         5         6         7         8

12345678901234567890123456789012345678901234567890123456789012345678901234567890

ATOM   4150  H   ALA A 431       8.674  16.036  12.858  1.00  0.00

*/
l=0;
p->ndisplace=0;
p->nalign=0;
p->walign=0.;
p->wdisplace=0.;
p->simple=1;
while(1){
  readagain:
  str=fgets(string,100,fp);
  if(str==NULL)break;
  if (feof(fp))break;

  sscanf(str,"%3s",end);if(strstr(end,"END")!=NULL){break;};
  sscanf(str,"%6s",remark);if(strstr(remark,"REMARK")!=NULL){goto readagain;};

  // paste the atom field
  strncpy(atom,string,4);atom[4]='\0';
  strncpy(hetatm,string,6);hetatm[6]='\0';

  if( (strstr(atom,"ATOM")!=NULL) || (strstr(hetatm,"HETATM")!=NULL) ) {
      parse_fixwidth(str,  7, 11, ind);
      parse_fixwidth(str, 13, 16, name);
      parse_fixwidth(str, 18, 20, resname);
      parse_fixwidth(str, 23, 26, resid);
      //parse_fixwidth(str, 22, 22, &chain);
      parse_fixwidth(str, 31, 38, x);
      parse_fixwidth(str, 39, 46, y);
      parse_fixwidth(str, 47, 54, z);
      parse_fixwidth(str, 55, 60, occ);
      parse_fixwidth(str, 61, 66, beta);
      // old buggy parsing
      // sscanf(str,"%s %d %s %s %d %lf %lf %lf %s %s",atom,&(p->atmnum[l]),(p->label[l]),(p->resname[l]),&(p->resid[l]),&x,&y,&z,sm,sm1);
      p->atmnum[l]=atoi(ind);
      p->resid[l]=atoi(resid);
      p->pos[l][0]=atof(x);
      p->pos[l][1]=atof(y);
      p->pos[l][2]=atof(z);
      strcpy(p->label[l],name);
      strcpy(p->resname[l],resname);
  
      #if defined (PLUMED_GROMACS)
      p->pos[l][0] /= 10.;
      p->pos[l][1] /= 10.;
      p->pos[l][2] /= 10.;
      #endif
      
      #if defined (STANDALONE)
      p->pos[l][0] *= mtd_data.ampli;
      p->pos[l][1] *= mtd_data.ampli;
      p->pos[l][2] *= mtd_data.ampli;
      #endif

      #if defined (PLUMED_GENERIC_CXX)
      p->pos[l][0] *= mtd_data.factor_angstrom;
      p->pos[l][1] *= mtd_data.factor_angstrom;
      p->pos[l][2] *= mtd_data.factor_angstrom;
      #endif
  
  //alignment
      tmp0=atof(occ);
      if(tmp0==0.){
        p->align[l]=0.;
      }
      else{
        p->align[l]=tmp0;
        p->align_to_frameset[p->nalign]=l;
        p->nalign++;
        if(p->nalign>=MAXATOMS_RMSD) plumed_error("Number of atoms per frame in path CV exceeded. Increase MAXATOMS_RMSD in metadyn.h and recompile");
        p->walign+=tmp0;
      }
      //displacement
      tmp0=atof(beta);
      if(tmp0==0.){
        p->displace[l]=0;
      }
      else{
        //p->displace[l]=1;
        p->displace[l]=tmp0;
        p->ndisplace++;
        p->wdisplace+=tmp0;
      }

//      fprintf(fplog,"RESID NUM %5d  RESID %5d  RESNAME %4s LABEL %4s PX %12.6f PY %12.6f PZ %12.6f OC %12.6f BE %12.6f\n",p->atmnum[l],p->resid[l],p->resname[l],p->label[l],p->pos[l][0],p->pos[l][1],p->pos[l][2],p->align[l], p->displace[l]);

      if( (p->displace[l]!=p->align[l]) || (p->displace[l]!=1.0) || (p->align[l]!=1.0)  )p->simple=0;
      l++;
      if(l>=MAXATOMS_PATH) plumed_error("Number of atoms per frame in path CV exceeded. Increase MAXATOMS_PATH in metadyn.h and recompile");
   }
}
fclose(fp);
p->natoms=l;
fprintf(fplog,"|---FOUND %d ATOMS FOR DISPLACEMENT\n",p->ndisplace);
fprintf(fplog,"|---TOTAL WEIGHT FOR  DISPLACEMENT %f\n",p->wdisplace);
fprintf(fplog,"|---FOUND %d ATOMS FOR ALIGNMENT\n",p->nalign);
fprintf(fplog,"|---TOTAL WEIGHT   FOR ALIGNMENT %f\n",p->walign);

if(p->nalign==0){
                   fprintf(fplog,"IT SEEMS YOU DONT  WANT TO ALIGN ANY ATOM\n");
                    fprintf(fplog,"Your frameset should look like:\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"                                                         the    the    \n");
                    fprintf(fplog,"                                                         align  mea      \n");
                    fprintf(fplog,"                                                         ment   sure     \n");
                    fprintf(fplog,"                                                          |     |  \n");
                    fprintf(fplog,"                                                          V     V  \n");
                    fprintf(fplog,"ATOM   4150  H   ALA A 431       8.674  16.036  12.858  1.00  0.20\n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                    plumed_error("PluMeD dead with errors: check log file"); 
                }
if(p->ndisplace==0){
                    fprintf(fplog,"IT SEEMS YOU DONT  WANT TO MEASURE THE DISPLACEMENT OF ANY ATOM\n");
                    fprintf(fplog,"Your frameset should look like:\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"                                                         the    the    \n");
                    fprintf(fplog,"                                                         align  mea      \n");
                    fprintf(fplog,"                                                         ment   sure     \n");
                    fprintf(fplog,"                                                          |     |  \n");
                    fprintf(fplog,"                                                          V     V  \n");
                    fprintf(fplog,"ATOM   4150  H   ALA A 431       8.674  16.036  12.858  1.00  0.20\n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                    plumed_error("PluMeD dead with errors: check log file"); 
                   }
// find the back index ( from alignment to frameset )
l=0;
 for(i=0;i< p->natoms;i++){
   if (p->align[i]==0.){p->frameset_to_align[i]=-1 ; } // if negative there's no position in the rmsd frame
   else {p->frameset_to_align[i]=l;l++;}; // if >=0 then provides the index
 }
 if(p->simple){
                    fprintf(fplog,"|-DOING SIMPLE ALIGNMENT: faster\n");
                    fprintf(fplog,"\n");
 }else {
                    fprintf(fplog,"|-DOING DOUBLE ALIGNMENT: slower\n");
                    fprintf(fplog,"\n");
 }
return 0;
};

// ------------------------------------------------------------------------------------------------

void PREFIX read_sz_map(struct sz_data *my_sz, char file_maps[129], char file_maps2[129],
                        char file_group[129], int read_mapfile, FILE *fplog){

 FILE   *fmap;
 int    i,ii,kk,jj,kkk;
 char   stringa[200],end[10],tmpstring[8],tmp2string[8],*cmstring;
 char   buf[1024];
 double r0,cutoff,map,weight;

// open INDEX file
   fmap=fopen(file_maps,"r");
   if(fmap == NULL){
       char buf[1024];
       sprintf(buf,"UNSUCCESSFULL OPENING FOR FILE %s",file_maps);
       plumed_error(buf);
   }

   kk=0;
   ii=0;
   jj=0;
   while(1){
         cmstring=fgets(stringa,200,fmap);
         if(cmstring==NULL){break;}
         sscanf(cmstring,"%3s",end);if(strstr(end,"END")!=NULL){break;};
         sscanf(cmstring,"%7s",tmpstring);if(strstr(tmpstring,"CONTACT")!=NULL){ii=ii+1;}
         sscanf(cmstring,"%5s",tmp2string);if(strstr(tmp2string,"GROUP")!=NULL){jj=jj+1;}

         sscanf(cmstring,"%7s %d %d %d %lf %d %d %lf %lf",tmpstring,&kkk,&(my_sz->my_cmap_pack.index1[kk]),
                &(my_sz->my_cmap_pack.index2[kk]),&r0,&(my_sz->my_cmap_pack.nn[kk]),
                &(my_sz->my_cmap_pack.nd[kk]),&cutoff,&weight);


         my_sz->my_cmap_pack.index1[kk]--;
         my_sz->my_cmap_pack.index2[kk]--;
         my_sz->my_cmap_pack.r0[kk] = (real) r0;
         my_sz->my_cmap_pack.cutoff[kk] = (real) cutoff;
         my_sz->my_cmap_pack.weight[kk] = (real) weight;
         kk=kk+1;
         if(kk>=MAXDIM_CMAP) plumed_error("Number of contacts in path CV exceeded. Increase MAXDIM_CMAP in metadyn.h and recompile");
        }

   my_sz->my_cmap_pack.number=ii;
   my_sz->my_cmap_pack.gnumber=jj;

   fprintf(fplog,"|--%d atomic contacts / %d group contacts \n",my_sz->my_cmap_pack.number,my_sz->my_cmap_pack.gnumber);
   fclose(fmap);

   int iii=0;
   if(read_mapfile==1){
    fmap=fopen(file_maps2,"r");
    if(fmap == NULL){
      sprintf(buf,"UNSUCCESSFULL OPENING FOR FILE %s",file_maps2);
      plumed_error(buf);
    }

    while(1){
    kk=0;
    while(1){
      cmstring=fgets(stringa,200,fmap);
      if(cmstring==NULL){break;}
      sscanf(cmstring,"%3s",end);if(strstr(end,"END")!=NULL){break;};
      sscanf(cmstring,"%d %d %d %lf",&kkk,&ii,&jj,&map);
 
      my_sz->my_cmap_pack.cmap[iii][kk] = (real) map;
 
      kk=kk+1;
    }

    if(cmstring==NULL){break;}
    iii=iii+1;
    }
    if(iii!=my_sz->number){
      sprintf(buf,"NUMBER OF FRAMES FOUND %d IS DIFFERENT FROM EXPECTED %d",iii,my_sz->number);
      plumed_error(buf);
    }
    fclose(fmap);
   }

    if(my_sz->my_cmap_pack.logical_group){
      fmap=fopen(file_group,"r");
      if(fmap == NULL){
       sprintf(buf,"UNSUCCESSFULL OPENING FOR FILE %s",file_group);
       plumed_error(buf); 
      }
      while(1){
       iii=0;
       cmstring=fgets(stringa,200,fmap);
       if(cmstring==NULL){break;}
       char *result = NULL;
       result = strtok( stringa, " \t" );
       while( result != NULL ) {
        if(iii==1) jj = atoi (result);
        if(iii==2) my_sz->my_cmap_pack.group.numatom[jj-1] = atoi (result);
        if(iii>2)  my_sz->my_cmap_pack.group.index[jj-1][iii-3] = atoi (result) - 1; 
        result = strtok( NULL, " \t" );
        iii=iii+1;
       }
      }
      my_sz->my_cmap_pack.group.number=jj;
      fclose(fmap);

      printf("|--Reading GROUP file. Found %d groups \n",  my_sz->my_cmap_pack.group.number);
      for(i=0;i<my_sz->my_cmap_pack.group.number;i++){
       fprintf(fplog,"|--GROUP %d #ofATOMS %d \n",i+1,my_sz->my_cmap_pack.group.numatom[i]);
       for(ii=0;ii<my_sz->my_cmap_pack.group.numatom[i];ii++){
        fprintf(fplog,"|---AT %d \n", my_sz->my_cmap_pack.group.index[i][ii]+1);
       }
      }
     } 
     
// Calculating number of atoms involved and creating the list
// from ATOMIC contacts
            ii=0;
            if(my_sz->my_cmap_pack.number!=0){my_sz->my_cmap_pack.list[ii]=my_sz->my_cmap_pack.index1[0];}

            for(iii=1;iii<my_sz->my_cmap_pack.number;iii++){
             int flag=0;
             for(kk=0;kk<=ii;kk++){
              if(my_sz->my_cmap_pack.index1[iii]==my_sz->my_cmap_pack.list[kk]){flag=1;}
             }
             if(flag==0){
             ii=ii+1;
             my_sz->my_cmap_pack.list[ii]=my_sz->my_cmap_pack.index1[iii];
             }
            }

           for(iii=0;iii<my_sz->my_cmap_pack.number;iii++){
             int flag=0;
             for(kk=0;kk<=ii;kk++){
              if(my_sz->my_cmap_pack.index2[iii]==my_sz->my_cmap_pack.list[kk]){flag=1;}
             }
             if(flag==0){
             ii=ii+1;
             my_sz->my_cmap_pack.list[ii]=my_sz->my_cmap_pack.index2[iii];
             }
            }

            if(my_sz->my_cmap_pack.number!=0){my_sz->my_cmap_pack.atoms=ii+1;}
            else{my_sz->my_cmap_pack.atoms=0;}
// From group

 if(my_sz->my_cmap_pack.logical_group){
           for(i=0;i<my_sz->my_cmap_pack.group.number;i++){
            for(ii=0;ii<my_sz->my_cmap_pack.group.numatom[i];ii++){
             kkk=my_sz->my_cmap_pack.group.index[i][ii];
             int flag=0;
                for(kk=0;kk<my_sz->my_cmap_pack.atoms;kk++){
                 if(my_sz->my_cmap_pack.list[kk]==kkk){
                  my_sz->my_cmap_pack.group.index_to_list[i][ii]=kk;
                  flag=1;
                 }
                }
              if(flag==0){
               my_sz->my_cmap_pack.list[my_sz->my_cmap_pack.atoms]=kkk;
               my_sz->my_cmap_pack.group.index_to_list[i][ii]=my_sz->my_cmap_pack.atoms;
               my_sz->my_cmap_pack.atoms=my_sz->my_cmap_pack.atoms+1;
              }
            }
           }
  }


            fprintf(fplog,"|--Total number of atoms involved %d \n",my_sz->my_cmap_pack.atoms);

// creating connection between my_cmap_pack.list e my_cmap_pack.index_from

           for(iii=0;iii<my_sz->my_cmap_pack.number;iii++){
            for(i=0;i<my_sz->my_cmap_pack.atoms;i++){
             if(my_sz->my_cmap_pack.index1[iii]==my_sz->my_cmap_pack.list[i]){my_sz->my_cmap_pack.index_from1[iii]=i;}
             if(my_sz->my_cmap_pack.index2[iii]==my_sz->my_cmap_pack.list[i]){my_sz->my_cmap_pack.index_from2[iii]=i;}
            }
           }

    return;

}

// ------------------------------------------------------------------------------------------------

void PREFIX cmap_running(int i_c, struct cmap_inpack *inpack, struct cmap_pack *my_cmap_pack){

       int mm,nn,i,j;
       real dist,xp,xq;
       rvec rij;

// atomic contact

       for(i=0;i<my_cmap_pack->number;i++){

        mm=my_cmap_pack->index_from1[i];
        nn=my_cmap_pack->index_from2[i];

// CMAP AND PBC
        if(colvar.cell_pbc[i_c]){
          minimal_image(inpack->r0[mm], inpack->r0[nn], &dist, rij);
        } else {
          dist=sqrt(pow2(inpack->r0[mm][0]-inpack->r0[nn][0])+
                    pow2(inpack->r0[mm][1]-inpack->r0[nn][1])+
                    pow2(inpack->r0[mm][2]-inpack->r0[nn][2]));
         };

/* original implementation

        dist=sqrt(pow2(inpack->r0[mm][0]-inpack->r0[nn][0])+
                  pow2(inpack->r0[mm][1]-inpack->r0[nn][1])+
                  pow2(inpack->r0[mm][2]-inpack->r0[nn][2]));
*/

        if (dist>my_cmap_pack->cutoff[i] || my_cmap_pack->weight[i]==0){
             inpack->cmap[i]=0.;
        }else{
             if (fabs(dist/my_cmap_pack->r0[i]-1.0)<0.00001){
              inpack->cmap[i]=(real) my_cmap_pack->nn[i]/my_cmap_pack->nd[i];
             } else { 
              power(dist/my_cmap_pack->r0[i],my_cmap_pack->nn[i],my_cmap_pack->nd[i],&xp,&xq);
              inpack->cmap[i]=(1.-xp)/(1.-xq)*my_cmap_pack->weight[i];
             }
        }
       }

       if(my_cmap_pack->logical_group){
// group contacts
// evaluating center of mass

        for(i=0;i<my_cmap_pack->group.number;i++){

         my_cmap_pack->group.rcm[i][0]=0.;
         my_cmap_pack->group.rcm[i][1]=0.;
         my_cmap_pack->group.rcm[i][2]=0.;

         for (j=0;j<my_cmap_pack->group.numatom[i];j++){
          mm=my_cmap_pack->group.index_to_list[i][j];
          my_cmap_pack->group.rcm[i][0] += inpack->r0[mm][0];
          my_cmap_pack->group.rcm[i][1] += inpack->r0[mm][1];
          my_cmap_pack->group.rcm[i][2] += inpack->r0[mm][2];
         }

         my_cmap_pack->group.rcm[i][0] /= (real) my_cmap_pack->group.numatom[i];
         my_cmap_pack->group.rcm[i][1] /= (real) my_cmap_pack->group.numatom[i];
         my_cmap_pack->group.rcm[i][2] /= (real) my_cmap_pack->group.numatom[i];

       }


// evaluating contacts

        for(i=0;i<my_cmap_pack->gnumber;i++){

         j=i+my_cmap_pack->number;
         mm=my_cmap_pack->index1[j];
         nn=my_cmap_pack->index2[j];

// CMAP AND PBC
        rvec rij;
        if(colvar.cell_pbc[i_c]){
          minimal_image(my_cmap_pack->group.rcm[mm], my_cmap_pack->group.rcm[nn], &dist, rij);
        } else {
         dist=sqrt(pow2(my_cmap_pack->group.rcm[mm][0]-my_cmap_pack->group.rcm[nn][0])+
                   pow2(my_cmap_pack->group.rcm[mm][1]-my_cmap_pack->group.rcm[nn][1])+
                   pow2(my_cmap_pack->group.rcm[mm][2]-my_cmap_pack->group.rcm[nn][2]));
         };


/* original implementation
         dist=sqrt(pow2(my_cmap_pack->group.rcm[mm][0]-my_cmap_pack->group.rcm[nn][0])+
                   pow2(my_cmap_pack->group.rcm[mm][1]-my_cmap_pack->group.rcm[nn][1])+
                   pow2(my_cmap_pack->group.rcm[mm][2]-my_cmap_pack->group.rcm[nn][2]));
*/ 

         if (dist>my_cmap_pack->cutoff[j] || my_cmap_pack->weight[j]==0){
             inpack->cmap[j]=0.;
         }
         else{
            if(fabs(dist/my_cmap_pack->r0[j]-1.0)<0.00001){
             inpack->cmap[j] = (real) my_cmap_pack->nn[j]/my_cmap_pack->nd[j];
            } else {
              power(dist/my_cmap_pack->r0[j],my_cmap_pack->nn[j],my_cmap_pack->nd[j],&xp,&xq);
              inpack->cmap[j]=(1.-xp)/(1.-xq)*my_cmap_pack->weight[j];
            }
         }

        }
   }
 }

// ------------------------------------------------------------------------------------------------

void PREFIX cmdist_eval(int i_c, int frame,struct cmap_inpack *inpack,struct cmap_outpack *work,
                  struct cmap_pack *my_cmap_pack,int dr1_calc){

       int    jj,k,i,j,ii,jjj,iii;
       int    tot;

       real tmp, dist_r0;
       real tmp4_r0_0,tmp4_r0_1,tmp4_r0_2;
       real tmp1_r0,tmp2_r0,tmp3_r0;
       real tmp1,tmp2,tmp3;
       real pow_P,pow_Q;
       real R01;
       int    P1,Q1;
       rvec rij;

       tmp=0.;
       tot=my_cmap_pack->number+my_cmap_pack->gnumber;

       for(i=0;i<tot;i++){
          tmp=tmp+pow2(my_cmap_pack->cmap[frame][i]-inpack->cmap[i]);
         }
        work->err=tmp;


                                      /* DERIVATIVE CALCULATION:respect to running frame and frameset */

// setting derivatives to zero

        for(k=0;k<my_cmap_pack->atoms;k++){
             work->derr_dr0[0][k]=0.;
             work->derr_dr0[1][k]=0.;
             work->derr_dr0[2][k]=0.;
        }


           for(j=0;j<my_cmap_pack->number;j++){
            if(my_cmap_pack->weight[j]!=0){
              ii=my_cmap_pack->index_from1[j];
              jj=my_cmap_pack->index_from2[j];

// CMAP AND PBC
             if(colvar.cell_pbc[i_c]){
               minimal_image(inpack->r0[ii], inpack->r0[jj], &dist_r0, rij);
             } else {
               rij[0] = inpack->r0[ii][0]-inpack->r0[jj][0];
               rij[1] = inpack->r0[ii][1]-inpack->r0[jj][1];
               rij[2] = inpack->r0[ii][2]-inpack->r0[jj][2];
               dist_r0=sqrt(pow2(inpack->r0[ii][0]-inpack->r0[jj][0])+
                         pow2(inpack->r0[ii][1]-inpack->r0[jj][1])+
                         pow2(inpack->r0[ii][2]-inpack->r0[jj][2]));
             };
/*
              dist_r0=sqrt(pow2(inpack->r0[ii][0]-inpack->r0[jj][0])+pow2(inpack->r0[ii][1]-inpack->r0[jj][1])+pow2(inpack->r0[ii][2]-inpack->r0[jj][2]));
*/
              tmp1_r0=inpack->cmap[j]-my_cmap_pack->cmap[frame][j];
              R01=my_cmap_pack->r0[j];
              P1=my_cmap_pack->nn[j];
              Q1=my_cmap_pack->nd[j];

              if(fabs(dist_r0/R01-1.0)<0.00001){
/* old
               outpack->derr_dr0[0][ii]+=tmp1_r0*(inpack->r0[ii][0]-inpack->r0[jj][0])*P1*(P1-Q1)/Q1;
               outpack->derr_dr0[1][ii]+=tmp1_r0*(inpack->r0[ii][1]-inpack->r0[jj][1])*P1*(P1-Q1)/Q1;
               outpack->derr_dr0[2][ii]+=tmp1_r0*(inpack->r0[ii][2]-inpack->r0[jj][2])*P1*(P1-Q1)/Q1;
               outpack->derr_dr0[0][jj]-=tmp1_r0*(inpack->r0[ii][0]-inpack->r0[jj][0])*P1*(P1-Q1)/Q1;
               outpack->derr_dr0[1][jj]-=tmp1_r0*(inpack->r0[ii][1]-inpack->r0[jj][1])*P1*(P1-Q1)/Q1;
               outpack->derr_dr0[2][jj]-=tmp1_r0*(inpack->r0[ii][2]-inpack->r0[jj][2])*P1*(P1-Q1)/Q1;
*/
// NEWPBC
               work->derr_dr0[0][ii]+=tmp1_r0*rij[0]*P1*(P1-Q1)/Q1;
               work->derr_dr0[1][ii]+=tmp1_r0*rij[1]*P1*(P1-Q1)/Q1;
               work->derr_dr0[2][ii]+=tmp1_r0*rij[2]*P1*(P1-Q1)/Q1;
               work->derr_dr0[0][jj]-=tmp1_r0*rij[0]*P1*(P1-Q1)/Q1;
               work->derr_dr0[1][jj]-=tmp1_r0*rij[1]*P1*(P1-Q1)/Q1;
               work->derr_dr0[2][jj]-=tmp1_r0*rij[2]*P1*(P1-Q1)/Q1;
              }else{
               power(dist_r0/R01,P1,Q1,&pow_P,&pow_Q);

               tmp2_r0=(Q1*pow_Q*(1.-pow_P)-P1*pow_P*(1.-pow_Q))*R01/dist_r0*my_cmap_pack->weight[j];
               tmp3_r0=R01*(1.-pow_Q)*(1.-pow_Q);
/* old
               tmp4_r0_0=(inpack->r0[ii][0]-inpack->r0[jj][0])/dist_r0;
               tmp4_r0_1=(inpack->r0[ii][1]-inpack->r0[jj][1])/dist_r0;
               tmp4_r0_2=(inpack->r0[ii][2]-inpack->r0[jj][2])/dist_r0;
*/
// NEWPBC
               tmp4_r0_0=rij[0]/dist_r0;
               tmp4_r0_1=rij[1]/dist_r0;
               tmp4_r0_2=rij[2]/dist_r0;

               tmp1=2*tmp1_r0*tmp2_r0*tmp4_r0_0/tmp3_r0;
               tmp2=2*tmp1_r0*tmp2_r0*tmp4_r0_1/tmp3_r0;
               tmp3=2*tmp1_r0*tmp2_r0*tmp4_r0_2/tmp3_r0;
               work->derr_dr0[0][ii]+=tmp1;
               work->derr_dr0[1][ii]+=tmp2;
               work->derr_dr0[2][ii]+=tmp3;
               work->derr_dr0[0][jj]-=tmp1;
               work->derr_dr0[1][jj]-=tmp2;
               work->derr_dr0[2][jj]-=tmp3;
              }
           }
          }
// case of group contact

          if(my_cmap_pack->logical_group){
           for(j=0;j<my_cmap_pack->gnumber;j++){

            i=j+my_cmap_pack->number;
            ii=my_cmap_pack->index1[i];
            jj=my_cmap_pack->index2[i];

// CMAP AND PBC
	    if(colvar.cell_pbc[i_c]){
	      minimal_image(my_cmap_pack->group.rcm[ii], my_cmap_pack->group.rcm[jj], &dist_r0, rij);
	    } else {
	      rij[0] = my_cmap_pack->group.rcm[ii][0]-my_cmap_pack->group.rcm[jj][0];
	      rij[1] = my_cmap_pack->group.rcm[ii][1]-my_cmap_pack->group.rcm[jj][1];
	      rij[2] = my_cmap_pack->group.rcm[ii][2]-my_cmap_pack->group.rcm[jj][2];
              dist_r0=sqrt(pow2(my_cmap_pack->group.rcm[ii][0]-my_cmap_pack->group.rcm[jj][0])+
                           pow2(my_cmap_pack->group.rcm[ii][1]-my_cmap_pack->group.rcm[jj][1])+
                           pow2(my_cmap_pack->group.rcm[ii][2]-my_cmap_pack->group.rcm[jj][2]));
	    };

/*
              dist_r0=sqrt(pow2(my_cmap_pack->group.rcm[ii][0]-my_cmap_pack->group.rcm[jj][0])+
                           pow2(my_cmap_pack->group.rcm[ii][1]-my_cmap_pack->group.rcm[jj][1])+
                           pow2(my_cmap_pack->group.rcm[ii][2]-my_cmap_pack->group.rcm[jj][2]));
*/
              tmp1_r0=inpack->cmap[i]-my_cmap_pack->cmap[frame][i];

              R01=my_cmap_pack->r0[i];
              P1=my_cmap_pack->nn[i];
              Q1=my_cmap_pack->nd[i];

              if(fabs(dist_r0/R01-1.0)<0.00001){
               for(jjj=0;jjj<my_cmap_pack->group.numatom[ii];jjj++){
                iii=my_cmap_pack->group.index_to_list[ii][jjj];
/* old
                outpack->derr_dr0[0][iii]+=tmp1_r0*(my_cmap_pack->group.rcm[ii][0]-my_cmap_pack->group.rcm[jj][0])*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[ii]);
                outpack->derr_dr0[1][iii]+=tmp1_r0*(my_cmap_pack->group.rcm[ii][1]-my_cmap_pack->group.rcm[jj][1])*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[ii]);
                outpack->derr_dr0[2][iii]+=tmp1_r0*(my_cmap_pack->group.rcm[ii][2]-my_cmap_pack->group.rcm[jj][2])*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[ii]);
*/
// NEWPBC
                work->derr_dr0[0][iii]+=tmp1_r0*rij[0]*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[ii]);
                work->derr_dr0[1][iii]+=tmp1_r0*rij[1]*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[ii]);
                work->derr_dr0[2][iii]+=tmp1_r0*rij[2]*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[ii]);



               }
               for(jjj=0;jjj<my_cmap_pack->group.numatom[jj];jjj++){
                iii=my_cmap_pack->group.index_to_list[jj][jjj];
/* old
                outpack->derr_dr0[0][iii]-=tmp1_r0*(my_cmap_pack->group.rcm[ii][0]-my_cmap_pack->group.rcm[jj][0])*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[jj]);
                outpack->derr_dr0[1][iii]-=tmp1_r0*(my_cmap_pack->group.rcm[ii][1]-my_cmap_pack->group.rcm[jj][1])*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[jj]);
                outpack->derr_dr0[2][iii]-=tmp1_r0*(my_cmap_pack->group.rcm[ii][2]-my_cmap_pack->group.rcm[jj][2])*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[jj]);
*/
// NEWPBC
                work->derr_dr0[0][iii]-=tmp1_r0*rij[0]*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[jj]);
                work->derr_dr0[1][iii]-=tmp1_r0*rij[1]*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[jj]);
                work->derr_dr0[2][iii]-=tmp1_r0*rij[2]*P1*(P1-Q1)/Q1/((real) my_cmap_pack->group.numatom[jj]);
               }
              }else{
               power(dist_r0/R01,P1,Q1,&pow_P,&pow_Q);

               tmp2_r0=(Q1*pow_Q*(1.-pow_P)-P1*pow_P*(1.-pow_Q))*R01/dist_r0*my_cmap_pack->weight[i];
               tmp3_r0=R01*(1.-pow_Q)*(1.-pow_Q);
/* old
               tmp4_r0_0=(my_cmap_pack->group.rcm[ii][0]-my_cmap_pack->group.rcm[jj][0])/dist_r0;
               tmp4_r0_1=(my_cmap_pack->group.rcm[ii][1]-my_cmap_pack->group.rcm[jj][1])/dist_r0;
               tmp4_r0_2=(my_cmap_pack->group.rcm[ii][2]-my_cmap_pack->group.rcm[jj][2])/dist_r0;
*/
// NEWPBC
               tmp4_r0_0=rij[0]/dist_r0;
               tmp4_r0_1=rij[1]/dist_r0;
               tmp4_r0_2=rij[2]/dist_r0;


               tmp1=2*tmp1_r0*tmp2_r0*tmp4_r0_0/tmp3_r0;
               tmp2=2*tmp1_r0*tmp2_r0*tmp4_r0_1/tmp3_r0;
               tmp3=2*tmp1_r0*tmp2_r0*tmp4_r0_2/tmp3_r0;

               for(jjj=0;jjj<my_cmap_pack->group.numatom[ii];jjj++){
                iii=my_cmap_pack->group.index_to_list[ii][jjj];
                work->derr_dr0[0][iii]+=tmp1/((real) my_cmap_pack->group.numatom[ii]);
                work->derr_dr0[1][iii]+=tmp2/((real) my_cmap_pack->group.numatom[ii]);
                work->derr_dr0[2][iii]+=tmp3/((real) my_cmap_pack->group.numatom[ii]);
               }
               for(jjj=0;jjj<my_cmap_pack->group.numatom[jj];jjj++){
                iii=my_cmap_pack->group.index_to_list[jj][jjj];
                work->derr_dr0[0][iii]-=tmp1/((real) my_cmap_pack->group.numatom[jj]);
                work->derr_dr0[1][iii]-=tmp2/((real) my_cmap_pack->group.numatom[jj]);
                work->derr_dr0[2][iii]-=tmp3/((real) my_cmap_pack->group.numatom[jj]);
               }
              }
           }
          }

           if(dr1_calc){
             for(i=0;i<tot;i++){
              work->derr_dcm[i]=-2.*(inpack->cmap[i]-my_cmap_pack->cmap[frame][i]);
             }
           }

}

// ------------------------------------------------------------------------------------------------

real PREFIX pow2(real x){
return x*x;
}

// ------------------------------------------------------------------------------------------------

void PREFIX power(real x,int p,int q,real *xp,real *xq){
     int i;
     real tot;

     tot=1;
     if(p>=q){
       for(i=1;i<=q;i++){
        tot=tot*x;
        }
        *xq=tot;
       for(i=q+1;i<=p;i++){
        tot=tot*x;
        }
        *xp=tot;}
     else{
       for(i=1;i<=p;i++){
        tot=tot*x;
        }
        *xp=tot;
       for(i=p+1;i<=q;i++){
        tot=tot*x;
        }
        *xq=tot;}
}

// ------------------------------------------------------------------------------------------------


void  PREFIX msd_calculation(struct coordinates_frameset *pframeset,struct cmap_inpack *c_inpack,
                             struct cmap_outpack *c_outpack,real dmsd_dr1[3][MAXATOMS_PATH],int der_frameref_on, int norot, int nocenter){

        int j,l,k,m,n,o,degeneration;

        real tmp0,tmp1,ndisplace,nalign,displace,align,walign,wdisplace;
        real coeff[3][MAXATOMS_PATH];
        real const1,const2; 
        struct rmsd_inpack inpack;
#ifdef RMSD_FULL
	struct rmsd_outpack work;
#else
	struct rmsd_mini_outpack work;
#endif

			/* number of atoms */
        nalign=pframeset->nalign;
	inpack.natoms=pframeset->nalign;
	inpack.totmass=pframeset->walign;
                           /* transfer atoms from the current set */
	for(k=0;k< pframeset->nalign;k++){
                          l=pframeset->align_to_coord[k];
                          inpack.r0[0][k]=c_inpack->r0[l][0];        
                          inpack.r0[1][k]=c_inpack->r0[l][1];        
                          inpack.r0[2][k]=c_inpack->r0[l][2];        
                          //fprintf(mtd_data.fplog,"COORD0 %d  %d %12.6f %12.6f %12.6f\n",k,l,inpack.r0[0][k],inpack.r0[1][k],inpack.r0[2][k]);
        }

                        /* transfer the atoms in the frameset */
        for(k=0;k<pframeset->nalign;k++){
                          l=pframeset->align_to_frameset[k];
                          inpack.r1[0][k]=pframeset->pos[l][0];
                          inpack.r1[1][k]=pframeset->pos[l][1];
                          inpack.r1[2][k]=pframeset->pos[l][2];
                          inpack.mass[k]=pframeset->align[l]; // this contains the weight in mass avg
                          //fprintf(mtd_data.fplog,"COORD1  %12.6f %12.6f %12.6f\n",inpack.r1[0][k],inpack.r1[1][k],inpack.r1[2][k]);
        }
#ifdef RMSD_FULL
        rmsd_pack(inpack,&work,7,1);
#else
       if(norot) { // dont rotate the frameset: fast and easy 
         rmsd_mini_pack_fake(inpack,&work,nocenter,pframeset->simple); 
         if(pframeset->simple==1){
               c_outpack->err=work.err;
               for(k=0;k<pframeset->natoms;k++){
                    for(l=0;l<3;l++){
                        c_outpack->derr_dr0[l][k]=-work.derr_dr0[l][k];   
                        if(der_frameref_on)  dmsd_dr1[l][k]=-work.derr_dr1[l][k];
                        
                    }
               }
               return;
         }
       } else { 
           //  full rmsd through simple derivative  of the eigenvalue   
           if(pframeset->simple==1){
               rmsd_mini_pack(inpack,&work,7,1,1);
               c_outpack->err=work.err*work.err;
               /* DERIVATIVE CALCULATION:respect to running frame */
               
               for(k=0;k<pframeset->natoms;k++){
                    for(l=0;l<3;l++){
                        c_outpack->derr_dr0[l][k]=2.*work.err*work.derr_dr0[l][k];   
                                    if(der_frameref_on){
                                           dmsd_dr1[l][k]=2.*work.err*work.derr_dr1[l][k];   
                                    }
                    }
               }
               // simple RMSD finite difference test system
               //rmsd_findiff_interface(inpack,&outpack);
               return;
            }else{
               //  full rmsd through derivative of rotation matrix  
               degeneration=rmsd_mini_pack(inpack,&work,7,0,1);
               //rmsd_mini_pack_fake(inpack,&outpack,7);
               // simple RMSD finite difference test system
               //rmsd_findiff_interface(inpack,&outpack);
            }
        }	
#endif
        //fprintf(mtd_data.fplog,"ERROR_RMSD_PACK %f\n",outpack.err);	
        //printf("ERROR_RMSD_PACK %f\n",outpack.err);	
	//(*err)=outpack.err;
                            /* check rotation and translation */

        //                for(k=0;k<my_coord0.natoms;k++){
        //                     my_coord0.pos[0][k]-=outpack.cmr0[0];
        //                     my_coord0.pos[1][k]-=outpack.cmr0[1];
        //                     my_coord0.pos[2][k]-=outpack.cmr0[2];
	//		}

        //                ll=init_pdb("current.pdb"); 
        //                plot_pdb(ll,&my_coord0); 

        //                mm=init_pdb("reference.pdb"); 
        //                for(k=0;k<pframeset->natoms;k++){
        //                     pframeset->pos[0][k]-=outpack.cmr1[0];
        //                     pframeset->pos[1][k]-=outpack.cmr1[1];
        //                     pframeset->pos[2][k]-=outpack.cmr1[2];
	//		}
        //                for(k=0;k<pframeset->natoms;k++){

        //                     tmp1=pframeset->pos[0][k]*outpack.d[0][0]+
        //                          pframeset->pos[1][k]*outpack.d[0][1]+
        //                          pframeset->pos[2][k]*outpack.d[0][2];

        //                     tmp2=pframeset->pos[0][k]*outpack.d[1][0]+
        //                          pframeset->pos[1][k]*outpack.d[1][1]+
        //                          pframeset->pos[2][k]*outpack.d[1][2];

        //                     tmp3=pframeset->pos[0][k]*outpack.d[2][0]+
        //                          pframeset->pos[1][k]*outpack.d[2][1]+
        //                          pframeset->pos[2][k]*outpack.d[2][2];

	//	             pframeset->pos[0][k]=tmp1;
        //                     pframeset->pos[1][k]=tmp2;
        //                     pframeset->pos[2][k]=tmp3;
	//		}
        //                plot_pdb(mm,pframeset); 
        //                CkExit();
                            /* REAL MSD CALCULATION */

	tmp0=0.;
	ndisplace=(real) pframeset->ndisplace; 
	walign= pframeset->walign; 
	wdisplace= pframeset->wdisplace; 
	
	for(k=0;k<pframeset->natoms;k++){
		for(l=0;l<3;l++){
			
			displace= pframeset->displace[k]; 
			align= pframeset->align[k]; 
			tmp1=0.;
			
			// contribution from rotated reference frame //
			for(m=0;m<3;m++){
				tmp1-=work.d[l][m]*(pframeset->pos[k][m]-work.cmr1[m]);
			}
			
			// contribution from running centered frame //
			j=pframeset->frameset_to_coord[k]; 
			tmp1+=(c_inpack->r0[j][l]-work.cmr0[l]); 
			
			
			//printf("DISPLACED ATOM %d %f %f \n",k,tmp2,tmp1*tmp1);
			coeff[l][k]=tmp1;// store coefficents for derivative usage// 
			tmp0+=tmp1*tmp1*displace; //squared distance added//
		}
	}  
	tmp0=tmp0/wdisplace;
	
	//printf("ERRR NEW %f \n",tmp0);
	
	c_outpack->err=tmp0;
	
	/* DERIVATIVE CALCULATION:respect to running frame */
	for(k=0;k<pframeset->natoms;k++){
		for(l=0;l<3;l++){
			
			displace= pframeset->displace[k]; 
			align= pframeset->align[k]; 
			
			tmp1 =2.0*coeff[l][k]*displace/wdisplace ;
			
			const1=2.0*align/(walign*wdisplace);
			
			if(const1>0.){
				for(o=0;o<pframeset->natoms;o++){
					tmp1 -=const1*coeff[l][o]*pframeset->displace[o]; 
				} 
			}
			
			j=pframeset->frameset_to_align[k]; //index of the k atom passed to the rmsd routine
			if(j>=0){
				for(m=0;m<pframeset->natoms;m++){
					// displace= pframeset->displace[m]; 
					const1=2.* pframeset->displace[m]/wdisplace ;
					for(n=0;n<3;n++){
						tmp0=0.;
						for(o=0;o<3;o++){
							tmp0+=work.dd_dr0[n][o][l][j]*(pframeset->pos[m][o]-work.cmr1[o]);
						}
						tmp0*=-const1*coeff[n][m];
						tmp1+=tmp0;    
					}
				}
			}
			c_outpack->derr_dr0[l][k]=tmp1;
		}
	}
	
	/* DERIVATIVE CALCULATION:respect to frameset  */
	if(der_frameref_on){
		for(k=0;k<pframeset->natoms;k++){
			
			j=pframeset->frameset_to_align[k]; //index of the k atom passed to the rmsd routine
			
			for(l=0;l<3;l++){
				
				tmp1=0.;
				
				if(j>=0){ // if it is an alignment atom
					for(m=0;m<pframeset->natoms;m++){
						const1=2.* pframeset->displace[m]/wdisplace ;
						for(n=0;n<3;n++){
							tmp0=0.;
							for(o=0;o<3;o++){
								tmp0+=work.dd_dr1[n][o][l][j]*
								(pframeset->pos[m][o]-work.cmr1[o]);
							}
							tmp0*=-const1*coeff[n][m]; 
							tmp1+= tmp0;    
						}
					}
				}
				
				displace= pframeset->displace[k]; 
				align= pframeset->align[k]; 
				
				tmp0=0.;
				for(o=0;o<3;o++){
					tmp0+=coeff[o][k]*work.d[o][l];
				}
				tmp1-=tmp0*2.*displace/wdisplace;
				
				if(j>=0){
					tmp0=0.;
					for(m=0;m<pframeset->natoms;m++){
						for(o=0;o<3;o++){
							tmp0+=coeff[o][m]*work.d[o][l]*pframeset->displace[m];
						}
					}
					tmp1 += tmp0*2.*align/(walign*wdisplace);
				}
				
				dmsd_dr1[l][k]=tmp1;
			}
		}
	}
	return;
};

// ------------------------------------------------------------------------------------------------


int PREFIX rmsd_pack(struct rmsd_inpack inpack,struct rmsd_outpack *work,int iopt,int iopt2)
{
/* declarations */
int i,j,k,l,p,ll,mm,nn,ii,ix,jx;
real rrsq,xx,yy,zz,m[4][4],rr1[4],rr0[4];
//cR double lambda[4],z[4][4],wk[20],s,q[4];
real lambda[4],s,q[4];
real dddq[3][3][4],gamma[3][3][3];
real dm_r1[4][4][3],dm_r0[4][4][3];
real dm_r1_store[4][4][3][MAXATOMS_RMSD];
real dm_r0_store[4][4][3][MAXATOMS_RMSD];
real derr_dr1_tmp[3][MAXATOMS_RMSD];
real derr_dr0_tmp[3][MAXATOMS_RMSD];
real dderr_dr1_dr1_tmp[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
real dderr_dr0_dr0_tmp[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
real dderr_dr1_dr0_tmp[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
real dderr_dr0_dr1_tmp[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
real pi1[3][3],pi0[3][3]; 
real tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
//cR int ier,arg1,arg2,arg3;
real alpha_m1[3][3],alpha_m2[3][3],alpha_m3[3][3],alpha_m4[3][3];
xx=0.;
yy=0.;
zz=0.;
tmp1=0.;
if(iopt==5 || iopt == 7 ){
	for(i=0;i<inpack.natoms;i++){
		xx+=inpack.r0[0][i]*inpack.mass[i];
		yy+=inpack.r0[1][i]*inpack.mass[i];
		zz+=inpack.r0[2][i]*inpack.mass[i];
                tmp1+=inpack.mass[i];
	}
	xx=xx/((real) tmp1);
	yy=yy/((real) tmp1);
	zz=zz/((real) tmp1);
};
work->cmr0[0]=xx;
work->cmr0[1]=yy;
work->cmr0[2]=zz;
for(i=0;i<inpack.natoms;i++){
	work->r0p[0][i]=inpack.r0[0][i]-xx;
	work->r0p[1][i]=inpack.r0[1][i]-yy;
	work->r0p[2][i]=inpack.r0[2][i]-zz;
}
xx=0.;
yy=0.;
zz=0.;
tmp1=0.;
if(iopt==6 || iopt == 7 ){
	for(i=0;i<inpack.natoms;i++){
		xx+=inpack.r1[0][i]*inpack.mass[i];
		yy+=inpack.r1[1][i]*inpack.mass[i];
		zz+=inpack.r1[2][i]*inpack.mass[i];
                tmp1+=inpack.mass[i]; 
	};
	xx=xx/((real) tmp1);
	yy=yy/((real) tmp1);
	zz=zz/((real) tmp1);

};
work->cmr1[0]=xx;
work->cmr1[1]=yy;
work->cmr1[2]=zz;
for(i=0;i<inpack.natoms;i++){
	work->r1p[0][i]=inpack.r1[0][i]-xx;
	work->r1p[1][i]=inpack.r1[1][i]-yy;
	work->r1p[2][i]=inpack.r1[2][i]-zz;
}
// CLEAN M MATRIX
for(i=0;i<4;i++){
	for(j=0;j<4;j++){
          m[i][j]=0.;  
	}
}
// ASSIGN MATRIX ELEMENTS
for(i=0;i<inpack.natoms;i++){
	
        tmp1=sqrt(inpack.mass[i]); 
        rr1[0]=work->r1p[0][i]*tmp1;
        rr1[1]=work->r1p[1][i]*tmp1;
        rr1[2]=work->r1p[2][i]*tmp1;
        rr0[0]=work->r0p[0][i]*tmp1;
        rr0[1]=work->r0p[1][i]*tmp1;
        rr0[2]=work->r0p[2][i]*tmp1;
	
        rrsq=pow(rr0[0],2)+pow(rr0[1],2)+pow(rr0[2],2)+pow(rr1[0],2)+pow(rr1[1],2)+pow(rr1[2],2);
     
        m[0][0] +=  rrsq+2.*(-rr0[0]*rr1[0]-rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[1][1] +=  rrsq+2.*(-rr0[0]*rr1[0]+rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[2][2] +=  rrsq+2.*(+rr0[0]*rr1[0]-rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[3][3] +=  rrsq+2.*(+rr0[0]*rr1[0]+rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[0][1] += 2.*(-rr0[1]*rr1[2]+rr0[2]*rr1[1]);
        m[0][2] += 2.*( rr0[0]*rr1[2]-rr0[2]*rr1[0]);
        m[0][3] += 2.*(-rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][2] -= 2.*( rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][3] -= 2.*( rr0[0]*rr1[2]+rr0[2]*rr1[0]);
        m[2][3] -= 2.*( rr0[1]*rr1[2]+rr0[2]*rr1[1]);

};
m[1][0] = m[0][1];
m[2][0] = m[0][2];
m[2][1] = m[1][2];
m[3][0] = m[0][3];
m[3][1] = m[1][3];
m[3][2] = m[2][3];

// DIAGONALIZE 

ql77_driver(m,lambda);
s=1.0;
if(m[0][0]<0.)s=-1.;//correct for negative values (?)
q[0]=s*m[0][0];
q[1]=s*m[1][0];
q[2]=s*m[2][0];
q[3]=s*m[3][0];
work->err=sqrt(lambda[0]/((real) inpack.natoms));
if(lambda[0]==lambda[1]) plumed_error("DIAGONALIZATION: NON UNIQUE SOLUTION");

if(iopt==0){return 0;}// JUST DIAGONALIZATION REQUIRED 

/*
 * Find the ROTATION matrix
 */
work->d[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3]       ; 
work->d[1][0]=2.0*(q[1]*q[2]-q[0]*q[3]);
work->d[2][0]=2.0*(q[1]*q[3]+q[0]*q[2]);
work->d[0][1]=2.0*(q[1]*q[2]+q[0]*q[3]);
work->d[1][1]=q[0]*q[0]+q[2]*q[2]-q[1]*q[1]-q[3]*q[3];
work->d[2][1]=2.0*(q[2]*q[3]-q[0]*q[1]);
work->d[0][2]=2.0*(q[1]*q[3]-q[0]*q[2]);
work->d[1][2]=2.0*(q[2]*q[3]+q[0]*q[1]);
work->d[2][2]=q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
#ifdef EXTREME_DEBUG
for (i=0;i<3;i++){
printf("D_MATRIX %12.6f %12.6f %12.6f\n",work->d[i][0],work->d[i][1],work->d[i][2]);
}
#endif
/* 
 * first derivative in perturbation theory
 */
dddq[0][0][0]= 2.0*q[0];
dddq[1][0][0]=-2.0*q[3];
dddq[2][0][0]= 2.0*q[2];
dddq[0][1][0]= 2.0*q[3];
dddq[1][1][0]= 2.0*q[0];
dddq[2][1][0]=-2.0*q[1];
dddq[0][2][0]=-2.0*q[2];
dddq[1][2][0]= 2.0*q[1];
dddq[2][2][0]= 2.0*q[0];

dddq[0][0][1]= 2.0*q[1];
dddq[1][0][1]= 2.0*q[2];
dddq[2][0][1]= 2.0*q[3];
dddq[0][1][1]= 2.0*q[2];
dddq[1][1][1]=-2.0*q[1];
dddq[2][1][1]=-2.0*q[0];
dddq[0][2][1]= 2.0*q[3];
dddq[1][2][1]= 2.0*q[0];
dddq[2][2][1]=-2.0*q[1];

dddq[0][0][2]=-2.0*q[2];
dddq[1][0][2]= 2.0*q[1];
dddq[2][0][2]= 2.0*q[0];
dddq[0][1][2]= 2.0*q[1];
dddq[1][1][2]= 2.0*q[2];
dddq[2][1][2]= 2.0*q[3];
dddq[0][2][2]=-2.0*q[0];
dddq[1][2][2]= 2.0*q[3];
dddq[2][2][2]=-2.0*q[2];

dddq[0][0][3]=-2.0*q[3];
dddq[1][0][3]=-2.0*q[0];
dddq[2][0][3]= 2.0*q[1];
dddq[0][1][3]= 2.0*q[0];
dddq[1][1][3]=-2.0*q[3];
dddq[2][1][3]= 2.0*q[2];
dddq[0][2][3]= 2.0*q[1];
dddq[1][2][3]= 2.0*q[2];
dddq[2][2][3]= 2.0*q[3];

#ifdef EXTREME_DEBUG
printf("\n");
for(i=0;i<4;i++){
	for(j=0;j<3;j++){
		printf("MATR %12.6f %12.6f %12.6f\n",dddq[j][0][i],dddq[j][1][i],dddq[j][2][i]);
	}
        printf("\n");
}
#endif
/*
 * Build gamma 3x3x3 matrix
 */
for(i=0;i<3;i++){     //direction 
    for(j=0;j<3;j++){     //direction 
        for(k=0;k<3;k++){     //eigenvector number
            gamma[i][j][k]=0.0;
            for(l=0;l<4;l++){   //components of each eigenvector in pert. series
              if(lambda[0]==lambda[k+1]){
                 plumed_error("FOUND DEGENERACY IN rmsd_pack ROUTINE");
               /*  write(*,*)"FOUND DEGENERACY IN RMSD_ESS ROUTINE "
                 write(*,*)"I'm DYING...."
                 write(*,*)"COPYING STACK HERE "
                 write(*,*)"R0"
                 do ll=1,n
                  write(*,'(f8.3,f8.3,f8.3)')r0(1,ll),r0(2,ll),r0(3,ll)
                 enddo
                 write(*,*)"R"
                 do ll=1,n
                  write(*,'(f8.3,f8.3,f8.3)')r(1,ll),r(2,ll),r(3,ll)
                 enddo
                 stop*/
		 } 
              else{
                gamma[i][j][k]=gamma[i][j][k]+dddq[i][j][l]*m[l][k+1]/(lambda[0]-lambda[k+1]);
	      }
	    }
	}

    }	
}
#ifdef EXTREME_DEBUG
for(i=0;i<3;i++){
	for(j=0;j<3;j++){
		printf("GAMM %12.6f %12.6f %12.6f\n",gamma[j][0][i],gamma[j][1][i],gamma[j][2][i]);
	}
        printf("\n");
}
#endif
/* 
 * Table of Derivative of the quaternion matrix respect to atom position
 */
for(i=0;i<inpack.natoms;i++){

        tmp1=(inpack.mass[i]); 
        rr1[0]=2.*work->r1p[0][i]*tmp1;
        rr1[1]=2.*work->r1p[1][i]*tmp1;
        rr1[2]=2.*work->r1p[2][i]*tmp1;
        rr0[0]=2.*work->r0p[0][i]*tmp1;
        rr0[1]=2.*work->r0p[1][i]*tmp1;
        rr0[2]=2.*work->r0p[2][i]*tmp1;
     

#ifdef EXTREME_DEBUG
        printf("ATOM %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n",rr0[0],rr0[1],rr0[2],rr1[0],rr1[1],rr1[2]);
#endif

        dm_r1 [0][0][0]=(rr1[0]-rr0[0]);
        dm_r1 [0][0][1]=(rr1[1]-rr0[1]);
        dm_r1 [0][0][2]=(rr1[2]-rr0[2]);
                      
        dm_r1 [0][1][0]=0.;
        dm_r1 [0][1][1]= rr0[2];
        dm_r1 [0][1][2]=-rr0[1];
                      
        dm_r1 [0][2][0]=-rr0[2];
        dm_r1 [0][2][1]= 0.;
        dm_r1 [0][2][2]= rr0[0];
                      
        dm_r1 [0][3][0]= rr0[1];
        dm_r1 [0][3][1]=-rr0[0];
        dm_r1 [0][3][2]= 0.;
                      
        dm_r1 [1][1][0]=(rr1[0]-rr0[0]);
        dm_r1 [1][1][1]=(rr1[1]+rr0[1]);
        dm_r1 [1][1][2]=(rr1[2]+rr0[2]);
                      
        dm_r1 [1][2][0]=-rr0[1];
        dm_r1 [1][2][1]=-rr0[0];
        dm_r1 [1][2][2]= 0.;
                      
        dm_r1 [1][3][0]=-rr0[2];
        dm_r1 [1][3][1]= 0.;
        dm_r1 [1][3][2]=-rr0[0];
                      
        dm_r1 [2][2][0]=(rr1[0]+rr0[0]);
        dm_r1 [2][2][1]=(rr1[1]-rr0[1]);
        dm_r1 [2][2][2]=(rr1[2]+rr0[2]);
                      
        dm_r1 [2][3][0]=0.;
        dm_r1 [2][3][1]=-rr0[2];
        dm_r1 [2][3][2]=-rr0[1];
                      
        dm_r1 [3][3][0]=(rr1[0]+rr0[0]);
        dm_r1 [3][3][1]=(rr1[1]+rr0[1]);
        dm_r1 [3][3][2]=(rr1[2]-rr0[2]);
/*
  derivative respec to to the other vector
 */
        dm_r0 [0][0][0]=-(rr1[0]-rr0[0]);
        dm_r0 [0][0][1]=-(rr1[1]-rr0[1]);
        dm_r0 [0][0][2]=-(rr1[2]-rr0[2]);
                      
        dm_r0 [0][1][0]=0.       ;
        dm_r0 [0][1][1]=-rr1[2];
        dm_r0 [0][1][2]=rr1[1];
                      
        dm_r0 [0][2][0]= rr1[2];      
        dm_r0 [0][2][1]= 0.;
        dm_r0 [0][2][2]=-rr1[0];
                      
        dm_r0 [0][3][0]=-rr1[1] ;     
        dm_r0 [0][3][1]= rr1[0];
        dm_r0 [0][3][2]= 0.;
                      
        dm_r0 [1][1][0]=-(rr1[0]-rr0[0]);
        dm_r0 [1][1][1]=(rr1[1]+rr0[1]);
        dm_r0 [1][1][2]=(rr1[2]+rr0[2]);
                      
        dm_r0 [1][2][0]=-rr1[1];
        dm_r0 [1][2][1]=-rr1[0];
        dm_r0 [1][2][2]= 0.;
                      
        dm_r0 [1][3][0]=-rr1[2];
        dm_r0 [1][3][1]= 0.;
        dm_r0 [1][3][2]=-rr1[0];
                      
        dm_r0 [2][2][0]=(rr1[0]+rr0[0]);
        dm_r0 [2][2][1]=-(rr1[1]-rr0[1]);
        dm_r0 [2][2][2]=(rr1[2]+rr0[2]);
                      
        dm_r0 [2][3][0]=0.;
        dm_r0 [2][3][1]=-rr1[2];
        dm_r0 [2][3][2]=-rr1[1];
                      
        dm_r0 [3][3][0]=(rr1[0]+rr0[0]);
        dm_r0 [3][3][1]=(rr1[1]+rr0[1]);
        dm_r0 [3][3][2]=-(rr1[2]-rr0[2]);
/*
 * write the diagonal
 */ 
	for(j=0;j<3;j++){

          dm_r1[1][0][j]=dm_r1[0][1][j];
          dm_r1[2][0][j]=dm_r1[0][2][j];
          dm_r1[3][0][j]=dm_r1[0][3][j];
          dm_r1[2][1][j]=dm_r1[1][2][j];
          dm_r1[3][1][j]=dm_r1[1][3][j];
          dm_r1[3][2][j]=dm_r1[2][3][j];

          dm_r0[1][0][j]=dm_r0[0][1][j];
          dm_r0[2][0][j]=dm_r0[0][2][j];
          dm_r0[3][0][j]=dm_r0[0][3][j];
          dm_r0[2][1][j]=dm_r0[1][2][j];
          dm_r0[3][1][j]=dm_r0[1][3][j];
          dm_r0[3][2][j]=dm_r0[2][3][j];
	  
          for(ll=0;ll<4;ll++){
          	for(mm=0;mm<4;mm++){
          		dm_r0_store[ll][mm][j][i]=dm_r0[ll][mm][j];
          		dm_r1_store[ll][mm][j][i]=dm_r1[ll][mm][j];
		};
	  };
 
	}
#ifdef EXTREME_DEBUG
	for(k=0;k<4;k++){
	for(l=0;l<4;l++){
	 printf("DM_R0 %12.6f %12.6f %12.6f\n",dm_r0[k][l][0],dm_r0[k][l][1],dm_r0[k][l][2]);
	}
        printf("\n"); 
        };
        for(k=0;k<4;k++){
	for(l=0;l<4;l++){
          printf("DM_R1 %12.6f %12.6f %12.6f\n",dm_r1[k][l][0],dm_r1[k][l][1],dm_r1[k][l][2]);
	}
        printf("\n"); 
        };
#endif
/*
 * pi matrix : coefficents in per theory
 */
	for(j=0;j<3;j++){
          pi1[0][j]=0.;
          pi1[1][j]=0.;
          pi1[2][j]=0.;
          pi0[0][j]=0.;
          pi0[1][j]=0.;
          pi0[2][j]=0.;
          work->derr_dr1 [j][i]=0.;
          work->derr_dr0 [j][i]=0.;

          for(k=0;k<4;k++){
            for(l=0;l<4;l++){
              work->derr_dr1[j][i]=work->derr_dr1[j][i]+q[k]*q[l]*dm_r1[l][k][j];
              work->derr_dr0[j][i]=work->derr_dr0[j][i]+q[k]*q[l]*dm_r0[l][k][j];
              for(mm=0;mm<3;mm++){
                pi0[mm][j]+=m[k][mm+1]*dm_r0[l][k][j]*q[l];
                pi1[mm][j]+=m[k][mm+1]*dm_r1[l][k][j]*q[l];  
	      };
	    };
	  };
          work->derr_dr1[j][i]=work->derr_dr1[j][i]/sqrt(4.*inpack.natoms*lambda[0]);
          work->derr_dr0[j][i]=work->derr_dr0[j][i]/sqrt(4.*inpack.natoms*lambda[0]);
	};
	for(j=0;j<3;j++){
		for (k=0;k<3;k++){
			for(l=0;l<3;l++){	    
              		work->dd_dr1[j][k][l][i]=0.;
              		work->dd_dr0[j][k][l][i]=0.;
			for(ii=0;ii<3;ii++){
                  		work->dd_dr1[j][k][l][i]+=gamma[j][k][ii]*pi1[ii][l]; 
                		work->dd_dr0[j][k][l][i]+=gamma[j][k][ii]*pi0[ii][l]; 
				}
			}
		}
	}
}
/*
 * Check arrays
 */
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR0 %12.6f %12.6f %12.6f\n",work->derr_dr0[0][i],work->derr_dr0[1][i],work->derr_dr0[2][i]);
}
printf("\n");
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR1 %12.6f %12.6f %12.6f\n",work->derr_dr1[0][i],work->derr_dr1[1][i],work->derr_dr1[2][i]);
}
for(i=0;i<inpack.natoms;i++){
for(j=0;j<3;j++){
for(k=0;k<3;k++){
printf("DD_DR0 %12.6f %12.6f %12.6f\n",work->dd_dr0[j][k][0][i],work->dd_dr0[j][k][1][i],work->dd_dr0[j][k][2][i]);
}}}
for(i=0;i<inpack.natoms;i++){
for(j=0;j<3;j++){
for(k=0;k<3;k++){
printf("DD_DR1 %12.6f %12.6f %12.6f\n",work->dd_dr1[j][k][0][i],work->dd_dr1[j][k][1][i],work->dd_dr1[j][k][2][i]);
}}}
#endif
/*
 * Second derivative if asked
 *
 *
 */
if(iopt2==2){
/*
 *   dr0 dr
 */

for(i=0;i<3;i++){      //  r0 atom component
	for(k=0;k<inpack.natoms;k++){ // r0 atom index 
		for (j=0;j<3;j++){//  r atom component
        		for(l=0;l<inpack.natoms;l++){// r atom index
            			work->dderr_dr0_dr0[i][k][j][l]=0.;
            			work->dderr_dr0_dr1[i][k][j][l]=0.;
            			work->dderr_dr1_dr0[i][k][j][l]=0.;
            			work->dderr_dr1_dr1[i][k][j][l]=0.;
            			for(p=1;p<4;p++){//eigenvector index 
              				 tmp1=0.;
            				 tmp2=0.;
              				 tmp3=0.;
              				 tmp4=0.;
  		 	       		 for(nn=0;nn<4;nn++){//eigenvector component
  		 	       		 	for(mm=0;mm<4;mm++){//eigenvector component
                    					tmp1+=(m[nn][0]*dm_r0_store[nn][mm][i][k]*m[mm][p]);
                    					tmp2+=(m[nn][0]*dm_r1_store[nn][mm][j][l]*m[mm][p]);
                					tmp3+=(m[nn][0]*dm_r1_store[nn][mm][i][k]*m[mm][p]);
                    					tmp4+=(m[nn][0]*dm_r0_store[nn][mm][j][l]*m[mm][p]);
						}
					 }
              				 work->dderr_dr0_dr1[i][k][j][l]+=2.*tmp1*tmp2/(lambda[0]-lambda[p]);
              				 work->dderr_dr1_dr0[i][k][j][l]+=2.*tmp3*tmp4/(lambda[0]-lambda[p]);
               				 work->dderr_dr1_dr1[i][k][j][l]+=2.*tmp2*tmp3/(lambda[0]-lambda[p]);
               				 work->dderr_dr0_dr0[i][k][j][l]+=2.*tmp1*tmp4/(lambda[0]-lambda[p]);
				};

/*
 *  second order diagonal and semi-diagonal terms
 */
				if(k-l==0){ 
					if(i-j==0){ 
				
						work->dderr_dr1_dr1[i][k][j][l]+=2.;
              					work->dderr_dr0_dr0[i][k][j][l]+=2.;
						if(i==0){
               						tmp5=2.*(-pow(m[0][0],2) - pow(m[1][0],2) 
									+pow(m[2][0],2) +pow(m[3][0],2));
               						tmp6=tmp5;
						}
						if(i==1){
							tmp5=2.0*(-pow(m[0][0],2)+ pow(m[1][0],2) 
									-pow(m[2][0],2) +pow(m[3][0],2));
							tmp6=tmp5;
						}
						
              					if(i==2){
						        tmp5=2.*(-pow(m[0][0],2) + pow(m[1][0],2) 
									+ pow(m[2][0],2) -pow(m[3][0],2));
               						tmp6=tmp5;
						}

					}
					else{

             					if( i==1 && j==0 ){// dy dx 
                					tmp5=4.*(-m[0][0]*m[3][0]-m[1][0]*m[2][0]);
                					tmp6=4.*( m[0][0]*m[3][0]-m[1][0]*m[2][0]);
              					};
             					if( i==2 && j==0 ){// dz dx 
                                                        tmp5=4.*( m[0][0]*m[2][0]-m[1][0]*m[3][0]);
              						tmp6=4.*(-m[0][0]*m[2][0]-m[1][0]*m[3][0]);
						};
             					if( i==0 && j==1 ){// dx dy 
                					tmp5=4.*( m[0][0]*m[3][0]-m[1][0]*m[2][0]);
                					tmp6=4.*(-m[0][0]*m[3][0]-m[1][0]*m[2][0]);
						};
             					if( i==2 && j==1 ){// dz dx 
                					tmp5=4.*(-m[0][0]*m[1][0]-m[2][0]*m[3][0]);
                					tmp6=4.*( m[0][0]*m[1][0]-m[2][0]*m[3][0]);
						};
             					if( i==0 && j==2 ){// dx dz 
                					tmp5=4.*(-m[2][0]*m[0][0]-m[3][0]*m[1][0]);
                					tmp6=4.*( m[2][0]*m[0][0]-m[3][0]*m[1][0]);
              					};
             					if( i==1 && j==2 ){// dy dz 
                					tmp5=4.*( m[1][0]*m[0][0]-m[3][0]*m[2][0]);
                					tmp6=4.*(-m[1][0]*m[0][0]-m[3][0]*m[2][0]);
              					};

					};
            		 		work->dderr_dr1_dr0[i][k][j][l]+=tmp5;
             				work->dderr_dr0_dr1[i][k][j][l]+=tmp6;
				}; 
            			work->dderr_dr0_dr1[i][k][j][l]=
                                work->dderr_dr0_dr1[i][k][j][l]/(2.*sqrt(lambda[0]*((real) inpack.natoms)))
                                 -sqrt(((real) inpack.natoms)/lambda[0])*work->derr_dr0[i][k]*work->derr_dr1[j][l];

    			        work->dderr_dr1_dr0[i][k][j][l]=
          			 work->dderr_dr1_dr0[i][k][j][l]/(2.*sqrt(lambda[0]*((real) inpack.natoms)))
            			-sqrt(((real) inpack.natoms)/lambda[0])*work->derr_dr1[i][k]*work->derr_dr0[j][l];

        		        work->dderr_dr1_dr1[i][k][j][l]=
           			 work->dderr_dr1_dr1[i][k][j][l]/(2.*sqrt(lambda[0]*((real) inpack.natoms)))
            			-sqrt(((real) inpack.natoms)/lambda[0])*work->derr_dr1[i][k]*work->derr_dr1[j][l];

            			work->dderr_dr0_dr0[i][k][j][l]=
          			 work->dderr_dr0_dr0[i][k][j][l]/(2.*sqrt(lambda[0]*((real) inpack.natoms)))
           			 -sqrt(((real) inpack.natoms)/lambda[0])*work->derr_dr0[i][k]*work->derr_dr0[j][l];
			};
		};
	};
};
};
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR0_DR0 %12.6f %12.6f %12.6f\n",work->dderr_dr0_dr0[j][i][0][l],work->dderr_dr0_dr0[j][i][1][l],work->dderr_dr0_dr0[j][i][2][l]);
		}
	}
	
}
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR1_DR1 %12.6f %12.6f %12.6f\n",work->dderr_dr1_dr1[j][i][0][l],work->dderr_dr1_dr1[j][i][1][l],work->dderr_dr1_dr1[j][i][2][l]);
		}
	}
	
}	
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR0_DR1 %12.6f %12.6f %12.6f\n",work->dderr_dr0_dr1[j][i][0][l],work->dderr_dr0_dr1[j][i][1][l],work->dderr_dr0_dr1[j][i][2][l]);
		}
	}
	
}	
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR1_DR0 %12.6f %12.6f %12.6f\n",work->dderr_dr1_dr0[j][i][0][l],work->dderr_dr1_dr0[j][i][1][l],work->dderr_dr1_dr0[j][i][2][l]);
		}
	}
	
}
#endif
/*
 * Now correct for cm - hard part in 2nd derivative
 *
 */
if(iopt==6 || iopt==7){
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
                derr_dr1_tmp[l][i]=work->derr_dr1[l][i];
			for(j=0;j<inpack.natoms;j++){
                        	derr_dr1_tmp[l][i]-=(1./((real) inpack.natoms))*work->derr_dr1[l][j];
			}
		}	
	}
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
			work->derr_dr1[l][i]=derr_dr1_tmp[l][i];
		}
	}	
}
if(iopt==5 || iopt==7){
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
                derr_dr0_tmp[l][i]=work->derr_dr0[l][i];
			for(j=0;j<inpack.natoms;j++){
                        	derr_dr0_tmp[l][i]-=(1./((real) inpack.natoms))*work->derr_dr0[l][j];
			}
		}	
	}
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
			work->derr_dr0[l][i]=derr_dr0_tmp[l][i];
		}
	}	
}
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR0 %12.6f %12.6f %12.6f\n",work->derr_dr0[0][i],work->derr_dr0[1][i],work->derr_dr0[2][i]);
}
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR1_CM %12.6f %12.6f %12.6f\n",work->derr_dr1[0][i],work->derr_dr1[1][i],work->derr_dr1[2][i]);
}
#endif
if(iopt2==2){
	if(iopt==6){//dr1 correction
		for(ix=0;ix<3;ix++){	
			for(jx=0;jx<3;jx++){	
                                alpha_m1[ix][jx]=0.0;
				for(i=0;i<inpack.natoms;i++){	
					for(j=0;j<inpack.natoms;j++){	
                                        alpha_m1[ix][jx]+=work->dderr_dr1_dr1[ix][i][jx][j];
					}
				}
		        	alpha_m1[ix][jx]= alpha_m1[ix][jx]/(((real) inpack.natoms*inpack.natoms));
			}
		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						dderr_dr1_dr1_tmp[ix][i][jx][j]=work->dderr_dr1_dr1[ix][i][jx][j];
						dderr_dr1_dr0_tmp[ix][i][jx][j]=work->dderr_dr0_dr1[ix][i][jx][j];
						dderr_dr0_dr1_tmp[ix][i][jx][j]=work->dderr_dr1_dr0[ix][i][jx][j];
						for(mm=0;mm<inpack.natoms;mm++){
							dderr_dr1_dr0_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*work->dderr_dr1_dr0[ix][i][jx][mm];
							dderr_dr0_dr1_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*work->dderr_dr0_dr1[ix][i][jx][mm];
							dderr_dr1_dr1_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
						         	(work->dderr_dr1_dr1[ix][i][jx][mm]+work->dderr_dr1_dr1[ix][mm][jx][j]);
						}
						dderr_dr1_dr1_tmp[ix][i][jx][j]+=alpha_m1[ix][jx];

					}
				}
			}

		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						work->dderr_dr1_dr1[ix][i][jx][j]=dderr_dr1_dr1_tmp[ix][i][jx][j];
						work->dderr_dr1_dr0[ix][i][jx][j]=dderr_dr1_dr0_tmp[ix][i][jx][j];
						work->dderr_dr0_dr1[ix][i][jx][j]=dderr_dr0_dr1_tmp[ix][i][jx][j];
					}
				}
			}
		}	
	}
	else if(iopt==5){//dr0 correction
		for(ix=0;ix<3;ix++){	
			for(jx=0;jx<3;jx++){	
                                alpha_m1[ix][jx]=0.0;
				for(i=0;i<inpack.natoms;i++){	
					for(j=0;j<inpack.natoms;j++){	
                                        alpha_m1[ix][jx]+=work->dderr_dr0_dr0[ix][i][jx][j];
					}
				}
		        	alpha_m1[ix][jx]= alpha_m1[ix][jx]/(((real) inpack.natoms*inpack.natoms));
			}
		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						dderr_dr0_dr0_tmp[ix][i][jx][j]=work->dderr_dr0_dr0[ix][i][jx][j];
						dderr_dr1_dr0_tmp[ix][i][jx][j]=work->dderr_dr1_dr0[ix][i][jx][j];
						dderr_dr0_dr1_tmp[ix][i][jx][j]=work->dderr_dr0_dr1[ix][i][jx][j];
						for(mm=0;mm<inpack.natoms;mm++){
							dderr_dr1_dr0_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*work->dderr_dr1_dr0[ix][i][jx][mm];
							dderr_dr0_dr1_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*work->dderr_dr0_dr1[ix][i][jx][mm];
							dderr_dr0_dr0_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
						         	(work->dderr_dr0_dr0[ix][i][jx][mm]+work->dderr_dr0_dr0[ix][mm][jx][j]);
						}
						dderr_dr0_dr0_tmp[ix][i][jx][j]+= alpha_m1[ix][jx];
					}
				}
			}
	
		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						work->dderr_dr0_dr0[ix][i][jx][j]=dderr_dr0_dr0_tmp[ix][i][jx][j];
						work->dderr_dr1_dr0[ix][i][jx][j]=dderr_dr1_dr0_tmp[ix][i][jx][j];
						work->dderr_dr0_dr1[ix][i][jx][j]=dderr_dr0_dr1_tmp[ix][i][jx][j];
					}
				}
			}
		}	
	}
	else if(iopt==7){
		for(ix=0;ix<3;ix++){	
			for(jx=0;jx<3;jx++){	
                                alpha_m1[ix][jx]=0.0;
                                alpha_m2[ix][jx]=0.0;
                                alpha_m3[ix][jx]=0.0;
                                alpha_m4[ix][jx]=0.0;
				for(i=0;i<inpack.natoms;i++){	
					for(j=0;j<inpack.natoms;j++){	
                                     		alpha_m1[ix][jx]+=work->dderr_dr0_dr0[ix][i][jx][j];
                           	        	alpha_m2[ix][jx]+=work->dderr_dr1_dr1[ix][i][jx][j];
                               			alpha_m3[ix][jx]+=work->dderr_dr0_dr1[ix][i][jx][j];
                               	       	 	alpha_m4[ix][jx]+=work->dderr_dr1_dr0[ix][i][jx][j];
					}
				}
		        	alpha_m1[ix][jx]= alpha_m1[ix][jx]/(((real) inpack.natoms*inpack.natoms));
		        	alpha_m2[ix][jx]= alpha_m2[ix][jx]/(((real) inpack.natoms*inpack.natoms));
		        	alpha_m3[ix][jx]= alpha_m3[ix][jx]/(((real) inpack.natoms*inpack.natoms));
		        	alpha_m4[ix][jx]= alpha_m4[ix][jx]/(((real) inpack.natoms*inpack.natoms));
			}
		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						dderr_dr0_dr0_tmp[ix][i][jx][j]=work->dderr_dr0_dr0[ix][i][jx][j];
						dderr_dr1_dr1_tmp[ix][i][jx][j]=work->dderr_dr1_dr1[ix][i][jx][j];
						dderr_dr1_dr0_tmp[ix][i][jx][j]=work->dderr_dr1_dr0[ix][i][jx][j];
						dderr_dr0_dr1_tmp[ix][i][jx][j]=work->dderr_dr0_dr1[ix][i][jx][j];
						for(mm=0;mm<inpack.natoms;mm++){
							dderr_dr0_dr0_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
						         	(work->dderr_dr0_dr0[ix][i][jx][mm]+work->dderr_dr0_dr0[ix][mm][jx][j]);
							dderr_dr1_dr1_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
						         	(work->dderr_dr1_dr1[ix][i][jx][mm]+work->dderr_dr1_dr1[ix][mm][jx][j]);
							dderr_dr0_dr1_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
						         	(work->dderr_dr0_dr1[ix][i][jx][mm]+work->dderr_dr0_dr1[ix][mm][jx][j]);
							dderr_dr1_dr0_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
					         		(work->dderr_dr1_dr0[ix][i][jx][mm]+work->dderr_dr1_dr0[ix][mm][jx][j]);
						}
						dderr_dr0_dr0_tmp[ix][i][jx][j]+=alpha_m1[ix][jx];
						dderr_dr1_dr1_tmp[ix][i][jx][j]+=alpha_m2[ix][jx];
						dderr_dr0_dr1_tmp[ix][i][jx][j]+=alpha_m3[ix][jx];
						dderr_dr1_dr0_tmp[ix][i][jx][j]+=alpha_m4[ix][jx];
					}
				}
			}

		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						work->dderr_dr0_dr0[ix][i][jx][j]=dderr_dr0_dr0_tmp[ix][i][jx][j];
						work->dderr_dr1_dr1[ix][i][jx][j]=dderr_dr1_dr1_tmp[ix][i][jx][j];
						work->dderr_dr1_dr0[ix][i][jx][j]=dderr_dr1_dr0_tmp[ix][i][jx][j];
						work->dderr_dr0_dr1[ix][i][jx][j]=dderr_dr0_dr1_tmp[ix][i][jx][j];
					}
				}
			}
		}	

	}
}
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR0_DR0_CM %12.6f %12.6f %12.6f\n",work->dderr_dr0_dr0[j][i][0][l],work->dderr_dr0_dr0[j][i][1][l],work->dderr_dr0_dr0[j][i][2][l]);
		}
	}
	
}
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR1_DR1_CM %12.6f %12.6f %12.6f\n",work->dderr_dr1_dr1[j][i][0][l],work->dderr_dr1_dr1[j][i][1][l],work->dderr_dr1_dr1[j][i][2][l]);
		}
	}
	
}	
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR0_DR1_CM %12.6f %12.6f %12.6f\n",work->dderr_dr0_dr1[j][i][0][l],work->dderr_dr0_dr1[j][i][1][l],work->dderr_dr0_dr1[j][i][2][l]);
		}
	}
	
}	
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR1_DR0_CM %12.6f %12.6f %12.6f\n",work->dderr_dr1_dr0[j][i][0][l],work->dderr_dr1_dr0[j][i][1][l],work->dderr_dr1_dr0[j][i][2][l]);
		}
	}
	
}
#endif

return 0;
}

// ------------------------------------------------------------------------------------------------
//
// iopt=5 reset cm of r0
// iopt=6 reset cm of r1
// iopt=7 reset cm of r0 and r1
// simple=1 correction on the eigenvalue rmsd only
// simple=0 derivative of the rotation matrix 
// simple=2 derivative of the rotation matrix  and eigenvalues
// 

int PREFIX rmsd_mini_pack(struct rmsd_inpack inpack,struct rmsd_mini_outpack *work,int iopt, int simple, int permissive)
{
/* declarations */
int i,j,k,l,ll,mm,ii,iopt2;
real rrsq,xx,yy,zz,m[4][4],rr1[4],rr0[4];
real lambda[4],s,q[4],fact1,fact2;
real dddq[3][3][4],gamma[3][3][3];
real dm_r1[4][4][3],dm_r0[4][4][3];
real dm_r1_store[4][4][3][MAXATOMS_RMSD];
real dm_r0_store[4][4][3][MAXATOMS_RMSD];
real derr_dr1_tmp[3][MAXATOMS_RMSD];
real derr_dr0_tmp[3][MAXATOMS_RMSD];
real pi1[3][3],pi0[3][3],tmp1,dnatoms; 
real dd_dr_temp[3][3][3][MAXATOMS_RMSD];
xx=0.;
yy=0.;
zz=0.;
tmp1=0.;
dnatoms=(inpack.natoms);
if(iopt==5 || iopt == 7 ){
	for(i=0;i<inpack.natoms;i++){
		xx+=inpack.r0[0][i]*inpack.mass[i];
		yy+=inpack.r0[1][i]*inpack.mass[i];
		zz+=inpack.r0[2][i]*inpack.mass[i];
//                tmp1+=inpack.mass[i];
                // no mass weight
//		xx+=inpack.r0[0][i];
//		yy+=inpack.r0[1][i];
//		zz+=inpack.r0[2][i];
//      
	}
//	xx=xx/((real) tmp1);
//	yy=yy/((real) tmp1);
//	zz=zz/((real) tmp1);
        xx=xx/(inpack.totmass);
        yy=yy/(inpack.totmass);
        zz=zz/(inpack.totmass);
//        xx=xx/dnatoms;
//        yy=yy/dnatoms;
//        zz=zz/dnatoms;
};
work->cmr0[0]=xx;
work->cmr0[1]=yy;
work->cmr0[2]=zz;
for(i=0;i<inpack.natoms;i++){
	work->r0p[0][i]=inpack.r0[0][i]-xx;
	work->r0p[1][i]=inpack.r0[1][i]-yy;
	work->r0p[2][i]=inpack.r0[2][i]-zz;
        // additional weighting 
       	work->r0p[0][i]*=inpack.mass[i];
	work->r0p[1][i]*=inpack.mass[i];
	work->r0p[2][i]*=inpack.mass[i];

 
}
xx=0.;
yy=0.;
zz=0.;
tmp1=0.;
if(iopt==6 || iopt == 7 ){
	for(i=0;i<inpack.natoms;i++){
		xx+=inpack.r1[0][i]*inpack.mass[i];
		yy+=inpack.r1[1][i]*inpack.mass[i];
		zz+=inpack.r1[2][i]*inpack.mass[i];
//              tmp1+=inpack.mass[i]; 
//		xx+=inpack.r1[0][i];
//		yy+=inpack.r1[1][i];
//		zz+=inpack.r1[2][i];
//      
	};
//	xx=xx/((real) tmp1);
//	yy=yy/((real) tmp1);
//	zz=zz/((real) tmp1);
        xx=xx/(inpack.totmass);
        yy=yy/(inpack.totmass);
        zz=zz/(inpack.totmass);
//        xx=xx/(dnatoms);
//        yy=yy/(dnatoms);
//        zz=zz/(dnatoms);

};
work->cmr1[0]=xx;
work->cmr1[1]=yy;
work->cmr1[2]=zz;
for(i=0;i<inpack.natoms;i++){
	work->r1p[0][i]=inpack.r1[0][i]-xx;
	work->r1p[1][i]=inpack.r1[1][i]-yy;
	work->r1p[2][i]=inpack.r1[2][i]-zz;
        // additional weighting 
	work->r1p[0][i]*=inpack.mass[i];
	work->r1p[1][i]*=inpack.mass[i];
	work->r1p[2][i]*=inpack.mass[i];

}
// CLEAN M MATRIX
for(i=0;i<4;i++){
	for(j=0;j<4;j++){
          m[i][j]=0.;  
	}
}
// ASSIGN MATRIX ELEMENTS
for(i=0;i<inpack.natoms;i++){
	
        rr1[0]=work->r1p[0][i];
        rr1[1]=work->r1p[1][i];
        rr1[2]=work->r1p[2][i];
        rr0[0]=work->r0p[0][i];
        rr0[1]=work->r0p[1][i];
        rr0[2]=work->r0p[2][i];
	
        rrsq=pow(rr0[0],2)+pow(rr0[1],2)+pow(rr0[2],2)+pow(rr1[0],2)+pow(rr1[1],2)+pow(rr1[2],2);
     
        m[0][0] +=  rrsq+2.*(-rr0[0]*rr1[0]-rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[1][1] +=  rrsq+2.*(-rr0[0]*rr1[0]+rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[2][2] +=  rrsq+2.*(+rr0[0]*rr1[0]-rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[3][3] +=  rrsq+2.*(+rr0[0]*rr1[0]+rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[0][1] += 2.*(-rr0[1]*rr1[2]+rr0[2]*rr1[1]);
        m[0][2] += 2.*( rr0[0]*rr1[2]-rr0[2]*rr1[0]);
        m[0][3] += 2.*(-rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][2] -= 2.*( rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][3] -= 2.*( rr0[0]*rr1[2]+rr0[2]*rr1[0]);
        m[2][3] -= 2.*( rr0[1]*rr1[2]+rr0[2]*rr1[1]);

};
m[1][0] = m[0][1];
m[2][0] = m[0][2];
m[2][1] = m[1][2];
m[3][0] = m[0][3];
m[3][1] = m[1][3];
m[3][2] = m[2][3];

// DIAGONALIZE 
int error;
error=ql77_driver(m,lambda);
s=1.0;
if(m[0][0]<0.)s=-1.;//correct for negative values (?)
q[0]=s*m[0][0];
q[1]=s*m[1][0];
q[2]=s*m[2][0];
q[3]=s*m[3][0];
work->err=sqrt(lambda[0]/(dnatoms));
if (permissive!=1){
	if(lambda[0]==lambda[1]) plumed_error("DIAGONALIZATION: NON UNIQUE SOLUTION");
}else if (permissive==1) {
	if(lambda[0]==lambda[1]) {fprintf(mtd_data.fplog,"|- WARNING: DEGENERACY IN THE DIAGONALIZATION...KEEP ON GOING!\n");} 
}
if(iopt==0){return 0;}// JUST DIAGONALIZATION REQUIRED 

/*
 * Find the ROTATION matrix
 */
work->d[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3]       ; 
work->d[1][0]=2.0*(q[1]*q[2]-q[0]*q[3]);
work->d[2][0]=2.0*(q[1]*q[3]+q[0]*q[2]);
work->d[0][1]=2.0*(q[1]*q[2]+q[0]*q[3]);
work->d[1][1]=q[0]*q[0]+q[2]*q[2]-q[1]*q[1]-q[3]*q[3];
work->d[2][1]=2.0*(q[2]*q[3]-q[0]*q[1]);
work->d[0][2]=2.0*(q[1]*q[3]-q[0]*q[2]);
work->d[1][2]=2.0*(q[2]*q[3]+q[0]*q[1]);
work->d[2][2]=q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
#ifdef EXTREME_DEBUG
for (i=0;i<3;i++){
printf("D_MATRIX %12.6f %12.6f %12.6f\n",work->d[i][0],work->d[i][1],work->d[i][2]);
}
#endif
/* 
 * first derivative in perturbation theory
 */
dddq[0][0][0]= 2.0*q[0];
dddq[1][0][0]=-2.0*q[3];
dddq[2][0][0]= 2.0*q[2];
dddq[0][1][0]= 2.0*q[3];
dddq[1][1][0]= 2.0*q[0];
dddq[2][1][0]=-2.0*q[1];
dddq[0][2][0]=-2.0*q[2];
dddq[1][2][0]= 2.0*q[1];
dddq[2][2][0]= 2.0*q[0];

dddq[0][0][1]= 2.0*q[1];
dddq[1][0][1]= 2.0*q[2];
dddq[2][0][1]= 2.0*q[3];
dddq[0][1][1]= 2.0*q[2];
dddq[1][1][1]=-2.0*q[1];
dddq[2][1][1]=-2.0*q[0];
dddq[0][2][1]= 2.0*q[3];
dddq[1][2][1]= 2.0*q[0];
dddq[2][2][1]=-2.0*q[1];

dddq[0][0][2]=-2.0*q[2];
dddq[1][0][2]= 2.0*q[1];
dddq[2][0][2]= 2.0*q[0];
dddq[0][1][2]= 2.0*q[1];
dddq[1][1][2]= 2.0*q[2];
dddq[2][1][2]= 2.0*q[3];
dddq[0][2][2]=-2.0*q[0];
dddq[1][2][2]= 2.0*q[3];
dddq[2][2][2]=-2.0*q[2];

dddq[0][0][3]=-2.0*q[3];
dddq[1][0][3]=-2.0*q[0];
dddq[2][0][3]= 2.0*q[1];
dddq[0][1][3]= 2.0*q[0];
dddq[1][1][3]=-2.0*q[3];
dddq[2][1][3]= 2.0*q[2];
dddq[0][2][3]= 2.0*q[1];
dddq[1][2][3]= 2.0*q[2];
dddq[2][2][3]= 2.0*q[3];

#ifdef EXTREME_DEBUG
printf("\n");
for(i=0;i<4;i++){
	for(j=0;j<3;j++){
		printf("MATR %12.6f %12.6f %12.6f\n",dddq[j][0][i],dddq[j][1][i],dddq[j][2][i]);
	}
        printf("\n");
}
#endif
int dump_zero;
dump_zero=0;
/*
 * Build gamma 3x3x3 matrix
 */
for(i=0;i<3;i++){     //direction 
    for(j=0;j<3;j++){     //direction 
        for(k=0;k<3;k++){     //eigenvector number
            gamma[i][j][k]=0.0;
            for(l=0;l<4;l++){   //components of each eigenvector in pert. series
              if(lambda[0]==lambda[k+1]){
                if (permissive!=1){		
                 plumed_error("FOUND DEGENERACY IN  rmsd_mini_pack  ROUTINE");
                }else if (permissive==1){
                  // give a signal to put the derivatives to zero
		  dump_zero=1;
                }
              }else{
                gamma[i][j][k]=gamma[i][j][k]+dddq[i][j][l]*m[l][k+1]/(lambda[0]-lambda[k+1]);
	      }
	    }
	}

    }	
}
// if this is set just put derivatives to zero in case of error=0. or degenerate eigenvalues 
if(dump_zero==1){
  fprintf(mtd_data.fplog,"|- WARNING: USING A FAKE VALUE FOR MSD \n"); 
  for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		work->derr_dr1[j][i]=0.;		
		work->derr_dr0[j][i]=0.;		
		for (k=0;k<3;k++){
  			// ders of rotmatrix
			for(l=0;l<3;l++){	    
      				work->dd_dr1[j][k][l][i]=0.; 
      				work->dd_dr0[j][k][l][i]=0.; 
                        } 
                } 
        } 
  } 
  return 1;
}
#ifdef EXTREME_DEBUG
for(i=0;i<3;i++){
	for(j=0;j<3;j++){
		printf("GAMM %12.6f %12.6f %12.6f\n",gamma[j][0][i],gamma[j][1][i],gamma[j][2][i]);
	}
        printf("\n");
}
#endif
/* 
 * Table of Derivative of the quaternion matrix respect to atom position
 */
fact1=sqrt(4.*dnatoms*lambda[0]);
for(i=0;i<inpack.natoms;i++){

        //tmp1=(inpack.mass[i]); 
        tmp1=1.0; 
        rr1[0]=2.*work->r1p[0][i]*tmp1;
        rr1[1]=2.*work->r1p[1][i]*tmp1;
        rr1[2]=2.*work->r1p[2][i]*tmp1;
        rr0[0]=2.*work->r0p[0][i]*tmp1;
        rr0[1]=2.*work->r0p[1][i]*tmp1;
        rr0[2]=2.*work->r0p[2][i]*tmp1;
     

#ifdef EXTREME_DEBUG
        printf("ATOM %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n",rr0[0],rr0[1],rr0[2],rr1[0],rr1[1],rr1[2]);
#endif

        dm_r1 [0][0][0]=(rr1[0]-rr0[0]);
        dm_r1 [0][0][1]=(rr1[1]-rr0[1]);
        dm_r1 [0][0][2]=(rr1[2]-rr0[2]);
                      
        dm_r1 [0][1][0]=0.;
        dm_r1 [0][1][1]= rr0[2];
        dm_r1 [0][1][2]=-rr0[1];
                      
        dm_r1 [0][2][0]=-rr0[2];
        dm_r1 [0][2][1]= 0.;
        dm_r1 [0][2][2]= rr0[0];
                      
        dm_r1 [0][3][0]= rr0[1];
        dm_r1 [0][3][1]=-rr0[0];
        dm_r1 [0][3][2]= 0.;
                      
        dm_r1 [1][1][0]=(rr1[0]-rr0[0]);
        dm_r1 [1][1][1]=(rr1[1]+rr0[1]);
        dm_r1 [1][1][2]=(rr1[2]+rr0[2]);
                      
        dm_r1 [1][2][0]=-rr0[1];
        dm_r1 [1][2][1]=-rr0[0];
        dm_r1 [1][2][2]= 0.;
                      
        dm_r1 [1][3][0]=-rr0[2];
        dm_r1 [1][3][1]= 0.;
        dm_r1 [1][3][2]=-rr0[0];
                      
        dm_r1 [2][2][0]=(rr1[0]+rr0[0]);
        dm_r1 [2][2][1]=(rr1[1]-rr0[1]);
        dm_r1 [2][2][2]=(rr1[2]+rr0[2]);
                      
        dm_r1 [2][3][0]=0.;
        dm_r1 [2][3][1]=-rr0[2];
        dm_r1 [2][3][2]=-rr0[1];
                      
        dm_r1 [3][3][0]=(rr1[0]+rr0[0]);
        dm_r1 [3][3][1]=(rr1[1]+rr0[1]);
        dm_r1 [3][3][2]=(rr1[2]-rr0[2]);
/*
  derivative respec to to the other vector
 */
        dm_r0 [0][0][0]=-(rr1[0]-rr0[0]);
        dm_r0 [0][0][1]=-(rr1[1]-rr0[1]);
        dm_r0 [0][0][2]=-(rr1[2]-rr0[2]);
                      
        dm_r0 [0][1][0]=0.       ;
        dm_r0 [0][1][1]=-rr1[2];
        dm_r0 [0][1][2]=rr1[1];
                      
        dm_r0 [0][2][0]= rr1[2];      
        dm_r0 [0][2][1]= 0.;
        dm_r0 [0][2][2]=-rr1[0];
                      
        dm_r0 [0][3][0]=-rr1[1] ;     
        dm_r0 [0][3][1]= rr1[0];
        dm_r0 [0][3][2]= 0.;
                      
        dm_r0 [1][1][0]=-(rr1[0]-rr0[0]);
        dm_r0 [1][1][1]=(rr1[1]+rr0[1]);
        dm_r0 [1][1][2]=(rr1[2]+rr0[2]);
                      
        dm_r0 [1][2][0]=-rr1[1];
        dm_r0 [1][2][1]=-rr1[0];
        dm_r0 [1][2][2]= 0.;
                      
        dm_r0 [1][3][0]=-rr1[2];
        dm_r0 [1][3][1]= 0.;
        dm_r0 [1][3][2]=-rr1[0];
                      
        dm_r0 [2][2][0]=(rr1[0]+rr0[0]);
        dm_r0 [2][2][1]=-(rr1[1]-rr0[1]);
        dm_r0 [2][2][2]=(rr1[2]+rr0[2]);
                      
        dm_r0 [2][3][0]=0.;
        dm_r0 [2][3][1]=-rr1[2];
        dm_r0 [2][3][2]=-rr1[1];
                      
        dm_r0 [3][3][0]=(rr1[0]+rr0[0]);
        dm_r0 [3][3][1]=(rr1[1]+rr0[1]);
        dm_r0 [3][3][2]=-(rr1[2]-rr0[2]);
/*
 * write the diagonal
 */ 
	for(j=0;j<3;j++){

          dm_r1[1][0][j]=dm_r1[0][1][j];
          dm_r1[2][0][j]=dm_r1[0][2][j];
          dm_r1[3][0][j]=dm_r1[0][3][j];
          dm_r1[2][1][j]=dm_r1[1][2][j];
          dm_r1[3][1][j]=dm_r1[1][3][j];
          dm_r1[3][2][j]=dm_r1[2][3][j];

          dm_r0[1][0][j]=dm_r0[0][1][j];
          dm_r0[2][0][j]=dm_r0[0][2][j];
          dm_r0[3][0][j]=dm_r0[0][3][j];
          dm_r0[2][1][j]=dm_r0[1][2][j];
          dm_r0[3][1][j]=dm_r0[1][3][j];
          dm_r0[3][2][j]=dm_r0[2][3][j];
	  
          for(ll=0;ll<4;ll++){
          	for(mm=0;mm<4;mm++){
          		dm_r0_store[ll][mm][j][i]=dm_r0[ll][mm][j];
          		dm_r1_store[ll][mm][j][i]=dm_r1[ll][mm][j];
		};
	  };
 
	}
#ifdef EXTREME_DEBUG
	for(k=0;k<4;k++){
	for(l=0;l<4;l++){
	 printf("DM_R0 %12.6f %12.6f %12.6f\n",dm_r0[k][l][0],dm_r0[k][l][1],dm_r0[k][l][2]);
	}
        printf("\n"); 
        };
        for(k=0;k<4;k++){
	for(l=0;l<4;l++){
          printf("DM_R1 %12.6f %12.6f %12.6f\n",dm_r1[k][l][0],dm_r1[k][l][1],dm_r1[k][l][2]);
	}
        printf("\n"); 
        };
#endif
/*
 * pi matrix : coefficents in per theory
 */
	for(j=0;j<3;j++){
          pi1[0][j]=0.;
          pi1[1][j]=0.;
          pi1[2][j]=0.;
          pi0[0][j]=0.;
          pi0[1][j]=0.;
          pi0[2][j]=0.;
          work->derr_dr1 [j][i]=0.;
          work->derr_dr0 [j][i]=0.;

          for(k=0;k<4;k++){
            for(l=0;l<4;l++){
              work->derr_dr1[j][i]=work->derr_dr1[j][i]+q[k]*q[l]*dm_r1[l][k][j];
              work->derr_dr0[j][i]=work->derr_dr0[j][i]+q[k]*q[l]*dm_r0[l][k][j];
              for(mm=0;mm<3;mm++){
                pi0[mm][j]+=m[k][mm+1]*dm_r0[l][k][j]*q[l];
                pi1[mm][j]+=m[k][mm+1]*dm_r1[l][k][j]*q[l];  
	      };
	    };
	  };
          work->derr_dr1[j][i]=work->derr_dr1[j][i]/fact1;
          work->derr_dr0[j][i]=work->derr_dr0[j][i]/fact1;

	};
	for(j=0;j<3;j++){
		for (k=0;k<3;k++){
			for(l=0;l<3;l++){	    
              		work->dd_dr1[j][k][l][i]=0.;
              		work->dd_dr0[j][k][l][i]=0.;
			for(ii=0;ii<3;ii++){
                  		work->dd_dr1[j][k][l][i]+=gamma[j][k][ii]*pi1[ii][l]; 
                		work->dd_dr0[j][k][l][i]+=gamma[j][k][ii]*pi0[ii][l]; 
				}
			}
		}
	}
}
/*
 * Check arrays
 */
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR0 %12.6f %12.6f %12.6f\n",work->derr_dr0[0][i],work->derr_dr0[1][i],work->derr_dr0[2][i]);
}
printf("\n");
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR1 %12.6f %12.6f %12.6f\n",work->derr_dr1[0][i],work->derr_dr1[1][i],work->derr_dr1[2][i]);
}
for(i=0;i<inpack.natoms;i++){
for(j=0;j<3;j++){
for(k=0;k<3;k++){
printf("DD_DR0 %12.6f %12.6f %12.6f\n",work->dd_dr0[j][k][0][i],work->dd_dr0[j][k][1][i],work->dd_dr0[j][k][2][i]);
}}}
for(i=0;i<inpack.natoms;i++){
for(j=0;j<3;j++){
for(k=0;k<3;k++){
printf("DD_DR1 %12.6f %12.6f %12.6f\n",work->dd_dr1[j][k][0][i],work->dd_dr1[j][k][1][i],work->dd_dr1[j][k][2][i]);
}}}
#endif
/*
 * Now correct for cm - hard part in 2nd derivative
 *
 */
if(iopt==6 || iopt==7){
        // don't correct for the rot matrix: only for error calc
        if(simple>0){
		for(l=0;l<3;l++){
			for(i=0;i<inpack.natoms;i++){
				fact2=inpack.mass[i]/(inpack.totmass);
                                //--> derr_dr1_tmp[l][i]=outpack->derr_dr1[l][i];
                                derr_dr1_tmp[l][i]=inpack.mass[i]*work->derr_dr1[l][i];
				for(j=0;j<inpack.natoms;j++){
                                	//--> derr_dr1_tmp[l][i]-=(inpack.mass[i]/(inpack.totmass))*outpack->derr_dr1[l][j];
                                	derr_dr1_tmp[l][i]-=inpack.mass[j]*fact2*work->derr_dr1[l][j];
				}
			}	
		}
		for(l=0;l<3;l++){
			for(i=0;i<inpack.natoms;i++){
				work->derr_dr1[l][i]=derr_dr1_tmp[l][i];
			}
		}	
        }
        if(simple==2 || simple==0){
        // correct the rotation matrix
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
			    for(k=0;k<3;k++){
                                for(l=0;l<inpack.natoms;l++){
                                   //--> dd_dr_temp[i][j][k][l]=outpack->dd_dr1[i][j][k][l];
                                   dd_dr_temp[i][j][k][l]=inpack.mass[l]*work->dd_dr1[i][j][k][l];
                                   tmp1=inpack.mass[l]/inpack.totmass; 
                                   for(mm=0;mm <inpack.natoms;mm++){
                                     //--> dd_dr_temp[i][j][k][l]-=outpack->dd_dr1[i][j][k][mm]*tmp1; 
                                     dd_dr_temp[i][j][k][l]-=work->dd_dr1[i][j][k][mm]*tmp1*inpack.mass[mm]; 
                                   }
                                }
                            }	
  	  	       	}	
                }
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
			    for(k=0;k<3;k++){
                                for(l=0;l<inpack.natoms;l++){
                                   work->dd_dr1[i][j][k][l]=dd_dr_temp[i][j][k][l];
                                }
                            }	
  	  	       	}	
                }
        }
	
}
if(iopt==5 || iopt==7){
        if(simple>0){
		for(l=0;l<3;l++){
			for(i=0;i<inpack.natoms;i++){
                                //--> derr_dr0_tmp[l][i]=outpack->derr_dr0[l][i];
				fact2=inpack.mass[i]/(inpack.totmass);
                                derr_dr0_tmp[l][i]=inpack.mass[i]*work->derr_dr0[l][i];
				for(j=0;j<inpack.natoms;j++){
                                	//--> derr_dr0_tmp[l][i]-=(inpack.mass[i]/(inpack.totmass))*outpack->derr_dr0[l][j];
                                	derr_dr0_tmp[l][i]-=inpack.mass[j]*fact2*work->derr_dr0[l][j];
				}
			}	
		}
		for(l=0;l<3;l++){
			for(i=0;i<inpack.natoms;i++){
				work->derr_dr0[l][i]=derr_dr0_tmp[l][i];
			}
		}	
        }
        if(simple==2 || simple==0){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
			    for(k=0;k<3;k++){
                                for(l=0;l<inpack.natoms;l++){
                                   //--> dd_dr_temp[i][j][k][l]=outpack->dd_dr0[i][j][k][l];
                                   dd_dr_temp[i][j][k][l]=inpack.mass[l]*work->dd_dr0[i][j][k][l];
                                   tmp1=inpack.mass[l]/inpack.totmass; 
                                   for(mm=0;mm <inpack.natoms;mm++){
                                     //--> dd_dr_temp[i][j][k][l]-=outpack->dd_dr0[i][j][k][mm]*tmp1; 
                                     dd_dr_temp[i][j][k][l]-=work->dd_dr0[i][j][k][mm]*tmp1*inpack.mass[mm]; 
                                   }
                                }
                            }	
  	  	       	}	
                }
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
			    for(k=0;k<3;k++){
                                for(l=0;l<inpack.natoms;l++){
                                   work->dd_dr0[i][j][k][l]=dd_dr_temp[i][j][k][l];
                                }
                            }	
  	  	       	}	
                }
	}
}
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR0 %12.6f %12.6f %12.6f\n",work->derr_dr0[0][i],work->derr_dr0[1][i],work->derr_dr0[2][i]);
}
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR1_CM %12.6f %12.6f %12.6f\n",work->derr_dr1[0][i],work->derr_dr1[1][i],work->derr_dr1[2][i]);
}
#endif
return 0;
}

// ------------------------------------------------------------------------------------------------

int PREFIX ql77_driver(real m[][4],real* lambda){
int i,j,ok;
double ll[16];
double *lambda2;
lambda2=(double *)malloc(4*sizeof(double));
for(i=0;i<4;i++){
//lambda2[i]=0.;
 for(j=0;j<4;j++){
	 ll[4*i+j]=m[i][j]; }};
#ifdef EXTREME_DEBUG
for(i=0;i<4;i++){
         printf("M_MATRIX %12.6f %12.6f %12.6f %12.6f\n",m[i][0],m[i][1],m[i][2],m[i][3]) ;  
}
#endif
ok=ql77(4,ll,lambda2);
#ifdef EXTREME_DEBUG
printf("EIGV %f %f %f %f\n",lambda2[0],lambda2[1],lambda2[2],lambda2[3]);
#endif
//back to square representation: columns have eigenvectors
        for(j=0;j<4;j++){
            lambda[j]=lambda2[j]; 
	    for(i=0;i<4;i++){
		 m[i][j]=ll[4*j+i]; 
            }
	};
free(lambda2); 
return ok;
};
int PREFIX ql77 (int n,double *x,double *d)
{
  int i,j,k,l,ni;
  double *e,h,g,f,b,s,p,r,c,absp;
  double totwork,work;

  const double eps=7.e-14,tol=1.e-30;

  e=(double *)malloc(n*sizeof(double));

  int danger;
  danger=0;

  totwork = 0;
  for(ni=1; ni<n; ni++)
    totwork += pow((double) n-ni,3);

  work=0;
  for(ni=1; (ni < n); ni++) {
    i=n-ni;
    l=i-1;
    h=0.0;                                                          
    g=x[i+n*(i-1)];
    if (l > 0) {
      for(k=0; (k<l); k++) 
	h=h+x[i+n*k]*x[i+n*k];
      s=h+g*g;
      if (s < tol)
	h=0.0;
      else if (h > 0.0) { 
	l=l+1;
	f=g;
	g=sqrt(s);
	g=-g;
	h=s-f*g;                                                           
	x[i+n*(i-1)]=f-g;
	f=0.0;
	for(j=0; (j<l); j++) {
	  x[j+n*i]=x[i+n*j]/h;
          if ( fabs(h) < tol ) {danger=1; fprintf(mtd_data.fplog,"|- WARNING: PROBLEM WITH DIAGONALIZATION...KEEP ON GOING!\n");} 
	  s=0.0;
	  for(k=0; (k<=j); k++)
	    s=s+x[j+n*k]*x[i+n*k];
	  for(k=j+1; (k<l); k++)
	    s=s+x[k+n*j]*x[i+n*k];
	  e[j]=s/h;
	  f=f+s*x[j+n*i];
	}
	f=f/(h+h);
	for(j=0; (j<l); j++) 
	  e[j]=e[j]-f*x[i+n*j];
	for(j=0; (j<l); j++) {
	  f=x[i+n*j];
	  s=e[j];
	  for(k=0; (k<=j); k++)
	    x[j+n*k]=x[j+n*k]-f*e[k]-x[i+n*k]*s;
	}
      }
    }
    d[i]=h;
    e[i-1]=g;

    work += pow((double) n-ni,3);
  }
  if (danger==1){
    for(i=0; (i<n); i++) {
      d[i] = 0.0;
      for(j=0; (j<n); j++) {
	x[j+n*i] = 0.0 ;
        if(i==j){
	  x[j+n*i] = 1.0 ;
        }  
      }
    }
    free(e);
    return 0;
  }
 

  /*
   *  accumulation of transformation matrix and intermediate d vector
   */
  
  d[0]=x[0];
  x[0]=1.0;

  work=0;
  for(i=1; (i<n); i++) {
    if (d[i] > 0.0) {
      for(j=0; (j<i); j++) {
	s=0.0;
	for(k=0; (k<i); k++) 
	  s=s+x[i+n*k]*x[k+n*j];
	for(k=0; (k<i); k++)
	  x[k+n*j]=x[k+n*j]-s*x[k+n*i];
      }
    }
    d[i]=x[i+n*i];
    x[i+n*i]=1.0;
    for(j=0; (j<i); j++) {
      x[i+n*j]=0.0;
      x[j+n*i]=0.0;
    }
    work += pow((double) i,3);
  }

  /*
   *  ql iterates
   */

  b=0.0;
  f=0.0;
  e[n-1]=0.0;
  totwork += pow((double) n,3);
  work=0;
  for(l=0; (l<n); l++) {
    h=eps*(fabs(d[l])+fabs(e[l]));
    if (h > b) 
      b=h;                                                   
    for(j=l; (j<n); j++) {
      if(fabs(e[j]) <= b) 
	break;
    }
    if (j != l) { 
      do {
	g=d[l];
	p=(d[l+1]-g)*0.5/e[l];
	r=sqrt(p*p+1.0);
	if(p < 0.0)
	  p=p-r;
	else
	  p=p+r;
	d[l]=e[l]/p;
	h=g-d[l];                                                     
	for(i=l+1; (i<n); i++)
	  d[i]=d[i]-h;                                                       
	f=f+h;                                                             
	p=d[j];
	c=1.0;
	s=0.0;
	for(ni=l; (ni<j); ni++) {
	  i=l+j-1-ni;
	  g=c*e[i];
	  h=c*p;
	  if(fabs(p) >= fabs(e[i])) {
	    c=e[i]/p;
	    r=sqrt(c*c+1.0);
	    e[i+1]=s*p*r;
	    s=c/r;
	    c=1.0/r;
	  } else {
	    c=p/e[i];
	    r=sqrt(c*c+1.0);
	    e[i+1]=s*e[i]*r;
	    s=1.0/r;
	    c=c/r;
	  }
	  p=c*d[i]-s*g;
	  d[i+1]=h+s*(c*g+s*d[i]);
	  for(k=0; (k<n); k++) {
	    h=x[k+n*(i+1)];
	    x[k+n*(i+1)]=x[k+n*i]*s+h*c;
	    x[k+n*i]=x[k+n*i]*c-h*s;
	  }
	}
	e[l]=s*p;
	d[l]=c*p;
      } while (fabs(e[l]) > b); 
    }
    d[l]=d[l]+f;

    work += pow((double) n-l,3);
  }

  /*
   *  put eigenvalues and eigenvectors in 
   *  desired ascending order
   */

 
  for(i=0; (i<n-1); i++) {
    k    = i;
    p    = d[i];
    absp = fabs(d[i]);
    for(j=i+1; (j<n); j++) {
      if(fabs(d[j]) < absp) {
	k    = j;
	p    = d[j];
	absp = fabs(d[j]);
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for(j=0; (j<n); j++) {
	p        = x[j+n*i];
	x[j+n*i] = x[j+n*k];
	x[j+n*k] = p;
      }
    }
  }

	  free(e);
 	  return 1;	
	}

//-----------------------------------------------------------------------------------------------------

void PREFIX dmsd_calculation(int i_c,struct coordinates_frameset *pframeset,struct cmap_inpack *c_inpack,
                             struct cmap_outpack *c_outpack,real dmsd_dr1[3][MAXATOMS_PATH]){

    int i,j,ix;
    rvec rij,rij0;
    real mod_rij,mod_rij0, drmsd, fact;

    drmsd = 0.;
    for(i=0;i<colvar.natoms[i_c];i++) for(ix=0;ix<3;ix++) {c_outpack->derr_dr0[ix][i] = dmsd_dr1[ix][i] = 0.;} 

  
    for(i=0;i<colvar.natoms[i_c]-1;i++) {
     for(j=i+1;j<colvar.natoms[i_c];j++) {
      minimal_image(c_inpack->r0[i], c_inpack->r0[j], &mod_rij, rij);
      minimal_image(pframeset->pos[i], pframeset->pos[j], &mod_rij0, rij0);
      drmsd += (mod_rij-mod_rij0)*(mod_rij-mod_rij0);
      for(ix=0;ix<3;ix++) {
        c_outpack->derr_dr0[ix][i] +=  rij[ix]*(mod_rij-mod_rij0)/mod_rij;
        c_outpack->derr_dr0[ix][j] += -rij[ix]*(mod_rij-mod_rij0)/mod_rij;
        dmsd_dr1[ix][i] += -rij0[ix]*(mod_rij-mod_rij0)/mod_rij0;
        dmsd_dr1[ix][j] +=  rij0[ix]*(mod_rij-mod_rij0)/mod_rij0;
      }
     }
    } 

   fact = 2./((real)colvar.natoms[i_c]*((real)colvar.natoms[i_c]-1.));

   c_outpack->err = drmsd * fact; 

   for(i=0;i<colvar.natoms[i_c];i++) {
    c_outpack->derr_dr0[0][i] *= 2.*fact;
    c_outpack->derr_dr0[1][i] *= 2.*fact;
    c_outpack->derr_dr0[2][i] *= 2.*fact; 
    dmsd_dr1[0][i] *= 2.*fact;
    dmsd_dr1[1][i] *= 2.*fact;
    dmsd_dr1[2][i] *= 2.*fact;
   }

}

int PREFIX read_sz_hybrid(struct sz_data *my_sz, FILE *fplog) {
//AVOID

	int size,l,i,j,k,jj,kk,ll,found,ii,ielem;
	char *str,ic[3],str2[100];
	l=0;
	
	/*
	 * Allocate the pointers: fully dynamical
	 */
	
	// allocate a number of pointers to framesets
	my_sz->hbd_frameset=(struct hybrid_frameset **)malloc((my_sz->number)*sizeof(struct hybrid_frameset *));
	
	for(i=0;i< my_sz->number;i++){
		
				// allocate the frame itself
                my_sz->hbd_frameset[i]=(struct hybrid_frameset *)malloc(sizeof(struct hybrid_frameset));
		
				// this array has the size of the real of calls
                my_sz->hbd_frameset[i]->hbd_nelem=my_sz->nhybrid;
                my_sz->hbd_frameset[i]->hbd_elem=(struct hybrid_elem **)malloc((my_sz->hbd_frameset[i]->hbd_nelem)*sizeof(struct hybrid_elem * ));
		
				// reset the number of atoms
                my_sz->hbd_frameset[i]->hbd_totalatoms=0;   
        
			    // nelem=number of real CVs in the hybrid representation
                for (j=0;j<my_sz->hbd_frameset[i]->hbd_nelem;j++){
                        my_sz->hbd_frameset[i]->hbd_elem[j]=(struct hybrid_elem *)malloc(sizeof(struct hybrid_elem));
                        // type of the cv in the list
                        my_sz->hbd_frameset[i]->hbd_elem[j]->cvtype=my_sz->lcvhybrid[j]; 
                        // progressive index in the list of the cvs
                        my_sz->hbd_frameset[i]->hbd_elem[j]->cvindex=my_sz->lhybrid[j]; 
                } 
        }

/*
 *  build names
 */
        str=&(my_sz->names[0]);
	
		/*
		 *	loop over all the frames
		 */
        for (i=1;i<= my_sz->number ;i++){
				strcpy(str2,my_sz->names);
				if(i<10){
					ic[0]='0'+i;
					ic[1]='\0';}
				else if(i<100) {
					ic[0]='0'+i/10 ;
					ic[1]='0'+i%10 ;
					ic[2]='\0';
				}
				else{
					plumed_error("TOO MANY FRAMES REQUIRED FOR NAME BUILDING!");
				}
				strcat(str2,ic);
				strcat(str2,".dat");
				fprintf(fplog,"|-\n");
				fprintf(fplog,"|-********************NOW READING*****************\n");
				fprintf(fplog,"|--%s\n",str2);
				fprintf(fplog,"|-********************NOW READING*****************\n");
				fprintf(fplog,"|-\n");

				FILE *fp;
				// test if the file exists 
				fp=fopen(str2,"r");
				if (fp == NULL) {
					char buf[1024];
					sprintf(buf,"UNSUCCESSFULL OPENING FOR FILE %s",str2);
					plumed_error(buf);
				}
				//
                // loop over the needed cv: open the file with the right reader and allocate accordingly 			
				//
				int start=0;
				int cvindex;
				for(j=0;j<my_sz->nhybrid;j++){ // loop over the calls
					
					// assign the type to the element
					my_sz->hbd_frameset[i-1]->hbd_elem[j]->cvindex=my_sz->lhybrid[j];
					my_sz->hbd_frameset[i-1]->hbd_elem[j]->cvtype=my_sz->lcvhybrid[j];
					
					// put always on for all the variables
					logical.always[my_sz->lhybrid[j]]=1;
					
					// switch off pdb
					my_sz->hbd_frameset[i-1]->hbd_elem[j]->contains_pdb=0;
					
					
					//printf("CVINDEX %d CVTYPE %d\n",my_sz->hbd_frameset[i-1]->hbd_elem[j]->cvindex,my_sz->hbd_frameset[i-1]->hbd_elem[j]->cvtype);

					switch(my_sz->lcvhybrid[j]) { // select the cv
						case 1:  
							my_sz->hbd_frameset[i-1]->hbd_totalatoms+=hbd_read_simple(fp,fplog,"DIST",my_sz->hbd_frameset[i-1]->hbd_elem[j],1); 
							break;
						case 3:  
							my_sz->hbd_frameset[i-1]->hbd_totalatoms+=hbd_read_simple(fp,fplog,"COORD",my_sz->hbd_frameset[i-1]->hbd_elem[j],1); 
							break;
							
						case 4:  
							my_sz->hbd_frameset[i-1]->hbd_totalatoms+=hbd_read_simple(fp,fplog,"ANGLE",my_sz->hbd_frameset[i-1]->hbd_elem[j],1); 
							break;
							
						case 5:  
							my_sz->hbd_frameset[i-1]->hbd_totalatoms+=hbd_read_simple(fp,fplog,"TORSION",my_sz->hbd_frameset[i-1]->hbd_elem[j],1); 
							break;

						case 32:  
							my_sz->hbd_frameset[i-1]->hbd_totalatoms+=hbd_read_simple(fp,fplog,"POSITION",my_sz->hbd_frameset[i-1]->hbd_elem[j],1); 
							break;
	
							
						case 53:  
							my_sz->hbd_frameset[i-1]->hbd_totalatoms+=hbd_read_msd(fp,fplog,my_sz->hbd_frameset[i-1]->hbd_elem[j],1);
							// pass here the info about rototrasl?	
							break;
							
						default: plumed_error("READER NOT IMPLEMENTED"); break;
					} 
				}
                fprintf(fplog,"|- TOTAL NUMBER OF ATOMS IN FRAMESET %d IS %d \n",i,my_sz->hbd_frameset[i-1]->hbd_totalatoms);
			
				my_sz->hbd_frameset[i-1]->hbd_totcvs=0;
				for(ii=0;ii<my_sz->nhybrid ;ii++){
					my_sz->hbd_frameset[i-1]->hbd_totcvs+=my_sz->hbd_frameset[i-1]->hbd_elem[ii]->ncv;
				}
				fprintf(fplog,"|- TOTAL NUMBER OF CVS IN FRAMESET %d IS %d \n",i,my_sz->hbd_frameset[i-1]->hbd_totcvs);
				fclose(fp);
				char string[400];*str;
//
// MATRIX READING: LEAVE THIS ASIDE FOR NOW: take it diagonal
//
//				char string[400],matrix[10],end[10],*str;
//                while(1){
//                     readagain:
//                     str=fgets(string,100,fp);  
//                     if(feof(fp))break;
//                     sscanf(str,"%3s",end);if(strstr(end,"END")!=NULL){break;};
//                     sscanf(str,"%6s",matrix);  
//                     if(strstr(matrix,"MATRIX")!=NULL){
//                            for(jj=0;jj<my_sz->nhybrid ;jj++){ 
//                     		str=fgets(string,100,fp);  
//                                //printf("MATRIX:   %s",string);
//                                char *pp=strtok(string," \t\n");
//                                my_sz->mathybrid[jj][0]=atof(pp);//printf("XXX %d %d %s\n",jj,0,pp);
//                                for(kk=1;kk<my_sz->nhybrid;kk++){
//                                    char *pp=strtok(NULL," \t\n");//printf("XXX %d %d  %s\n",jj,kk,pp);
//                                    my_sz->mathybrid[jj][kk]=atof(pp);
//                                }
//                            }
//     	             }else{
//                             fprintf(mtd_data.fplog,"|- ASSUMING DIAGONAL MATRIX FOR THIS FRAME \n");
//                             break; 
//                     }
//				}
//				                fclose(fp);
//				// print the matrix
//                for(jj=0;jj<my_sz->nhybrid;jj++){
//                   printf("MATRIX "); 
//                   for(kk=0;kk<my_sz->nhybrid;kk++){
//                     	printf(" %12.6f ",my_sz->mathybrid[jj][kk]); 
//				   }  	
//                   printf("\n");
//                }
			
                // is the frameset fixed?
			
                fp=fopen(str2,"r");
                char fixed[100];
                my_sz->hbd_frameset[i-1]->fixed=0; 
                while(1){
                     str=fgets(string,100,fp);  
                     if(feof(fp))break;
                     sscanf(str,"%6s",fixed);  
                     if(strstr(fixed,"FIXED")!=NULL){
                           fprintf(fplog,"|- THIS IS A FIXED FRAME \n"); 
                           my_sz->hbd_frameset[i-1]->fixed=1;break;
                     }
                }
                fclose(fp);
		// now I do not have the jacobian allocated yet
		my_sz->hbd_frameset[i-1]->has_jac_alloc=0;
 
        }
	
	/*
	 * allocate the matrix for the coupling
	 */
	
	fprintf(fplog,"|- ALLOCATING METRIC MATRIX \n");
	i=my_sz->hbd_frameset[0]->hbd_totcvs;
	my_sz->mathybrid=float_2d_array_alloc(i,i);
	my_sz->myinverse=float_2d_array_alloc(i,i);
	
	// final lines should be the MATRIX
	// first set it to default
	
	fprintf(fplog,"|- METRIC MATRIX SET TO DIAGONAL\n");

	for(ii=0;ii<i ;ii++){
		for(jj=0;jj<i ;jj++){
			my_sz->mathybrid[ii][jj]=0.;
		}
		my_sz->mathybrid[ii][ii]=1.;
	}

      
	/*
	 * Now do the backtable: from each element to a global backtable
	 */
	
	fprintf(fplog,"|- WRITING GLOBAL BACKTABLE\n");
	struct  hybrid_elem *ee;
	for (i=0;i< my_sz->number ;i++){
		l=0;
		// a global backtable to keep track of all the index of all atoms
		//snew(my_sz->hbd_frameset[i]->backtable,my_sz->hbd_frameset[i]->hbd_totalatoms);
		my_sz->hbd_frameset[i]->backtable=(int *)malloc(my_sz->hbd_frameset[i]->hbd_totalatoms*sizeof(int));
		//printf("|- CHECK ON TOTAL NUMBER OF ATOMS\n");
		for (j=0;j<my_sz->hbd_frameset[i]->hbd_nelem;j++){
			ee=(my_sz->hbd_frameset[i]->hbd_elem[j]) ; 
			for (k=0;k<ee->nder;k++){
				my_sz->hbd_frameset[i]->backtable[l]= ee->backtable[k];
				//printf("|- FRAMESET %2d BACKTABLE  ID %3d ATOM  %3d \n",i,l,my_sz->hbd_frameset[i]->backtable[l]);
				l++;
			}
		}	
		//printf("|- TOTAL NUMBER OF ATOMS IN FRAMESET %d IS %d \n",i,l);
		//printf("|- ALLOCATING DERIVATIVES FOR %d ELEM \n",my_sz->hbd_frameset[i]->hbd_totalatoms);
		snew(my_sz->hbd_frameset[i]->myder,my_sz->hbd_frameset[i]->hbd_totalatoms); 
	}	
	
	/*
	 * Allocate the running frame: need for comparison with the stored frames
	 */

	fprintf(fplog,"|- ALLOCATING RUNNING STRUCTURE: \n");
	
	clone_hybrid_frameset(&my_sz->hbd_running,my_sz->hbd_frameset[0],0,fplog);
	
	// set the initial backtrack length for jacobian computation to zero

	my_sz->hbd_running->jac_index_length=0;
	
	
	for (j=0;j<my_sz->hbd_running->hbd_nelem;j++){

		// take one frame as reference for allocating the structure: this 
		// impose that all the structs must have the same number of elements
		struct hybrid_elem *elem_run=my_sz->hbd_running->hbd_elem[j];
		struct hybrid_elem *elem_ref=my_sz->hbd_frameset[0]->hbd_elem[j];
		
		int ncv=elem_ref->ncv;
		
		// single elements allocate immediately the tools for the jacobian
		elem_run->has_jac_alloc=1;
		elem_run->natoms_per_cv=(int *)malloc(ncv*sizeof(int));
		for(i=0;i<ncv;i++){
			elem_run->natoms_per_cv[i]=elem_ref->natoms_per_cv[i];	
			//fprintf(fplog,"|- BUILDING RUNNING STRUCTURE: ELEM %d CV %d DEPENDS ON %d atoms\n",j,i,elem_run->natoms_per_cv[i]);
		}
		
		// allocate the list of the atoms
		elem_run->list_per_cv=(int **)malloc(ncv*sizeof(int *));
		for(i=0;i<ncv;i++){
			elem_run->list_per_cv[i]=(int *)malloc(elem_run->natoms_per_cv[i]*sizeof(int));
			if(elem_run->list_per_cv[i]==NULL){fprintf(fplog,"ALLOCATION FAILED"); EXIT();}
			//fprintf(fplog,"|- BUILDING RUNNING STRUCTURE: ELEM %d CV %d DEPENDS ON : ",j,i); 
			for(k=0;k< elem_run->natoms_per_cv[i] ;k++){
				elem_run->list_per_cv[i][k]=elem_ref->list_per_cv[i][k];
				//fprintf(fplog," %d ",elem_run->list_per_cv[i][k]); 
			}
			//fprintf(fplog,"\n");
		}
		
		// the derivative for each cv depends has natoms_per_cv[i]*3 elementzs (rvec is 3 floats)
		elem_run->cvder=(rvec **)malloc(ncv*sizeof(rvec *));// 3n cvs for each structure
		if(elem_run->cvder==NULL){fprintf(fplog,"ALLOCATION FAILED"); EXIT();}

		for(i=0;i<ncv;i++)snew(elem_run->cvder[i],elem_run->natoms_per_cv[i]); // each structure may depend on the other atoms
				
	}
	fprintf(fplog,"|- HBD_RUNNING: NELEM %d NATOM %d\n",my_sz->hbd_running->hbd_nelem, my_sz->hbd_running->hbd_totalatoms);


	
	return my_sz->hbd_frameset[0]->hbd_totalatoms;
	

//AVOID
};

// ------------------------------------------------------------------------------------------------


int PREFIX  hbd_read_simple (  FILE *myfile, FILE *fplog, const char *istring, struct  hybrid_elem *elem, int hasfile){
//AVOID
    FILE *initpos;
    char string[400],remark[10],end[10],*str,type[10];
    initpos=myfile;
	int cvindex,cvtype;
    int done=0;
    int rdone=0;
    double uno;
	int i;
	
	cvindex=elem->cvindex;
	cvtype=elem->cvtype;
	if(hasfile){
		fprintf(fplog,"|----------------------------------------------\n");
		fprintf(fplog,"|- READING %s VARIABLE WITH SIMPLE READER:\n",istring);
		// scan a remark line and check if the type corresponds (zero level of idiotproofing)
		while(1){
		readagain:
			str=fgets(string,100,myfile);
			//printf("LINE %s ",string);
			sscanf(str,"%3s",end);if(strstr(end,"END")!=NULL){break;}; 
			sscanf(str,"%6s",remark);
			if(strstr(remark,"REMARK")!=NULL){
				strcpy(elem->remarkfield,str);
				sscanf(str,"%s %s",remark,type); 
				if(strstr(type,istring)!=NULL){fprintf(fplog,"|- TYPE CONFIRMED: %s \n",istring);rdone=1;}
				goto readagain;
			};
			if(!done){
				snew(elem->ref_dist,1);   
				sscanf(str,"%lf",&uno);done=1; elem->ref_dist[0] = (real) uno; 
				fprintf(fplog,"|- FOUND REFERENCE VALUE: %f\n",elem->ref_dist[0]); 
			}
		} 
		if(rdone==0) plumed_error("MISSING REMARK FIELD IN PATH FRAMES FOR CV CONFIRM");
		if(done==0)  plumed_error("MISSING END FIELD IN PATH FRAMES");
		fprintf(fplog,"|- READ CORRECTLY ENDED\n");
	}else{
		snew(elem->ref_dist,1);
		elem->ref_dist[0] = 0.;
	}
	
    // access colvar for making a correct allocation of the derivative 
    
	elem->nder=colvar.natoms[cvindex]; 
    fprintf(fplog,"|- NATOMS INVOLVED IN CV TOTAL = %d  : ",elem->nder); 
    //snew(elem->myder,elem->nder);   
    snew(elem->backtable,elem->nder);  
	//
	// one element can contain many cv: this is not the case
	// needed for jacobian stuff: only one cv contained
	elem->cvder=(rvec **)malloc(sizeof(rvec *));snew(elem->cvder[0],elem->nder);   
	elem->natoms_per_cv=(int *)malloc(sizeof(int));
	elem->list_per_cv=(int **)malloc(sizeof(int *));
	elem->list_per_cv[0]=(int *)malloc( elem->nder * sizeof(int ));
	// this contains only one cv and one list of atoms
	elem->natoms_per_cv[0]=elem->nder;
	for(i=0;i<elem->natoms_per_cv[0];i++){
		elem->list_per_cv[0][i]=colvar.cvatoms[cvindex][i];
	}
    for(i=0;i<elem->nder;i++){
        // backtable containt the index of the atoms on which you calculate the property
        elem->backtable[i]=colvar.cvatoms[cvindex][i]; 
        fprintf(fplog,"%d ",elem->backtable[i]+1);
    }
    fprintf(fplog,"\n");
	// increment of just one element
	elem->ncv=1;
	elem->has_diff_alloc=0;
	elem->has_jac_alloc=1;
	fprintf(fplog,"|----------------------------------------------\n");
    return elem->nder ;
//AVOID
};

int PREFIX  hbd_read_msd (  FILE *myfile, FILE *fplog, struct hybrid_elem *elem , int hasfile){
//AVOID
	int cvindex,cvtype;
	FILE *initpos;
	fpos_t pos;
    char string[400],remark[10],end[10],*str,type[10];
    initpos=myfile;
    int done=0;
    int rdone=0;
    double uno;
	int i,j;

	
	fprintf(fplog,"|----------------------------------------------\n");
	fprintf(fplog,"|- ENTERED HBD_READ_MSD : \n");
	fprintf(fplog,"|-  \n");
	if(hasfile){
		fprintf(fplog,"|- NOTE : \n");
		fprintf(fplog,"|-  ONLY POSITIONS ARE STORED. THE ALIGNMENT \n");
		fprintf(fplog,"|-  PATTERN IS THAT IN THE MSD \n");
		fprintf(fplog,"|-  VARIABLE DECLARED IN CV %d\n",elem->cvindex+1);
		fprintf(fplog,"|-  \n");
	}


	cvindex=elem->cvindex;
	cvtype=elem->cvtype;
	if(hasfile){
		fgetpos (myfile,&pos); // remember the position
		// scan a remark line and check if the type corresponds (zero level of idiotproofing)
		while(1){
		readagain:
			str=fgets(string,100,myfile);
			if(str==NULL)break;
			//printf("LINE %s ",string);
			sscanf(str,"%3s",end);if(strstr(end,"END")!=NULL){done=1;break;}; 
			sscanf(str,"%6s",remark);
			if(strstr(remark,"REMARK")!=NULL){
				strcpy(elem->remarkfield,str);
				sscanf(str,"%s %s",remark,type); 
				if(strstr(type,"MSD")!=NULL){fprintf(fplog,"|- TYPE CONFIRMED: MSD \n");rdone=1;}
			};
			goto readagain;
		} 
		if(rdone==0) plumed_error("MISSING REMARK FIELD IN PATH FRAMES FOR CV CONFIRM");
		if(done==0)  plumed_error("MISSING END FIELD IN PATH FRAMES");
		fprintf(fplog,"|- READ CORRECTLY ENDED\n");
		//reset to the initial address
		fsetpos (myfile, &pos);
		
		// parse the input through a specific pdb reader
		read_pdb (&(elem->pdb),myfile,fplog);
		// a flag just to tell if you have pdb or not
		elem->contains_pdb=1;
		
		
		// check now if the msd is a dummy variable (I do not want you to mess up with stuff)
		
		if(colvar.rmsd_container[cvindex].dummy==0){
			sprintf(str,"THE %d MSD CV INVOLVED IN HYBRID PATH SHOULD BE DECLARED AS DUMMY \n",cvindex+1);
			plumed_error(str);
			EXIT();
		}else{
			fprintf(fplog,"|- THE MSD CV IS CORRECTLY DECLARED AS DUMMY \n");
		}

		
		// check now if the indexes of the pdb match
		for(i=0;i<elem->pdb->natoms;i++){
			if(elem->pdb->index[i]!=colvar.rmsd_container[cvindex].mypdb->index[i]){
				sprintf(str,"ATOM %d ( index %d ) IN THE FRAME DOES NOT MATCH\n ATOM %d ( index %d ) IN THE MSD CV %d \n",i,elem->pdb->index[i],i,colvar.rmsd_container[cvindex].mypdb->index[i],cvindex+1);
				plumed_error(str);
			}else{
				//	fprintf(fplog,"|- ATOM %d ( index %d ) IN THE FRAME MATCHES \n|- ATOM %d ( index %d ) IN THE MSD CV %d \n",i,elem->pdb->index[i],i,colvar.rmsd_container[cvindex].mypdb->index[i],cvindex+1);
				
			}
		}
	}else{
		// in case of couplingmatrix you  do not need all this information
		elem->contains_pdb=1;
		// copy the pdb from the  frame
		allocate_pdb(&elem->pdb, colvar.natoms[cvindex]);
		copy_pdb(colvar.rmsd_container[cvindex].mypdb, elem->pdb, fplog);
	}
	
	elem->nder=colvar.natoms[cvindex];
	fprintf(fplog,"|- NATOMS INVOLVED IN CV TOTAL = %d  : ",elem->nder); 
    //snew(elem->myder,elem->nder);  
    snew(elem->backtable,elem->nder);   
	

	elem->ncv=3*elem->pdb->natoms;
	
	elem->has_diff_alloc=0;
	elem->has_jac_alloc=1;

	// reference distances
	snew(elem->ref_dist,elem->ncv);  
	
	j=0;
	for (i=0;i<elem->pdb->natoms;i++){
#if defined (PLUMED_GROMACS)	
		elem->ref_dist[j]=elem->pdb->x[i]/10.; /* printf("|- COORD X AT %d IS :%f \n",i,elem->ref_dist[j]);*/ j++;
		elem->ref_dist[j]=elem->pdb->y[i]/10.; /* printf("|- COORD Y AT %d IS :%f \n",i,elem->ref_dist[j]);*/ j++;
		elem->ref_dist[j]=elem->pdb->z[i]/10.; /* printf("|- COORD Z AT %d IS :%f \n",i,elem->ref_dist[j]);*/ j++;
#else
		elem->ref_dist[j]=elem->pdb->x[i]; /* printf("|- COORD X AT %d IS :%f \n",i,elem->ref_dist[j]);*/ j++;
		elem->ref_dist[j]=elem->pdb->y[i]; /* printf("|- COORD Y AT %d IS :%f \n",i,elem->ref_dist[j]);*/ j++;
		elem->ref_dist[j]=elem->pdb->z[i]; /* printf("|- COORD Z AT %d IS :%f \n",i,elem->ref_dist[j]);*/ j++;
#endif
	}
	
	//
	// needed for jacobian stuff
	// one element can contain many cvs
	// and a cv may depend on lots of atoms
	//
	
	// allocate the 3n ders
	elem->cvder=(rvec **)malloc(elem->ncv*sizeof(rvec *));// 3n cvs for each structure
	for(i=0;i<elem->ncv;i++)snew(elem->cvder[i],elem->nder); // each structure may depend on the other atoms

	// allocate the (trivial) number of atoms per cv
	elem->natoms_per_cv=(int *)malloc(elem->ncv*sizeof(int));
	for(i=0;i<elem->ncv;i++)elem->natoms_per_cv[i]=elem->nder; 
	
	// allocate the list of the atoms
	elem->list_per_cv=(int **)malloc(elem->ncv*sizeof(int *));
	for(i=0;i<elem->ncv;i++){
		elem->list_per_cv[i]=(int *)malloc(elem->ncv*sizeof(int));
		for(j=0;j<elem->natoms_per_cv[i];j++)elem->list_per_cv[i][j]=colvar.cvatoms[cvindex][j];
	}

	
    for(i=0;i<elem->nder;i++){
        // backtable containt the index of the atoms on which you calculate the property
        elem->backtable[i]=colvar.cvatoms[cvindex][i]; 
        fprintf(fplog,"%d ",elem->backtable[i]+1);
    }
    fprintf(fplog,"\n");
	
	
	fprintf(fplog,"|- \n"); 
	fprintf(fplog,"|----------------------------------------------\n");
	return elem->nder ;
//AVOID
};

int  PREFIX  hbd_collect_config ( struct hybrid_frameset *running  ) {
//AVOID
	
	int i;
	//fprintf(mtd_data.fplog,"|- HBD_COLLECT_CONFIG: ENTERING\n");

   for (i=0;i< running->hbd_nelem;i++){
           switch (running->hbd_elem[i]->cvtype){  
               case 1: /*fprintf(mtd_data.fplog,"|- HBD_COLLECT_CONFIG: TYPE %d \n", running.hbd_elem[i]->cvtype);*/
				   hbd_copy_simple(running->hbd_elem[i]); break ;
			   case	3: /*fprintf(mtd_data.fplog,"|- HBD_COLLECT_CONFIG: TYPE %d \n", running.hbd_elem[i]->cvtype);*/ 
				   hbd_copy_simple(running->hbd_elem[i]); break ;
			   case 4: /*fprintf(mtd_data.fplog,"|- HBD_COLLECT_CONFIG: TYPE %d \n", running.hbd_elem[i]->cvtype);*/ 
				   hbd_copy_simple(running->hbd_elem[i]); break ;
			   case 5: /*fprintf(mtd_data.fplog,"|- HBD_COLLECT_CONFIG: TYPE %d \n", running.hbd_elem[i]->cvtype);*/
				   hbd_copy_simple(running->hbd_elem[i]); break ;
			   case 32:/*fprintf(mtd_data.fplog,"|- HBD_COLLECT_CONFIG: TYPE %d \n", running.hbd_elem[i]->cvtype);*/ 
				   hbd_copy_simple(running->hbd_elem[i]); break ;
			   case 53:/*fprintf(mtd_data.fplog,"|- HBD_COLLECT_CONFIG: TYPE %d \n", running.hbd_elem[i]->cvtype);*/ 
				   hbd_copy_msd(running->hbd_elem[i]); break ;

                // in case you have a rmsd you just have to retrieve the configurations : dont use chain rules. The metrics does it all  
			   default:fprintf(mtd_data.fplog,"|- HBD_COLLECT_CONFIG: TYPE %d \n", running->hbd_elem[i]->cvtype);  plumed_error("|- NOT IMPLEMENTED: now dying..."); break;
           }
   }
  // fprintf(mtd_data.fplog,"|- HBD_COLLECT_CONFIG: EXITING\n");

   return 1; 
//AVOID
}

/*
 * Just simple copy the value and the derivative
 */
void  PREFIX  hbd_copy_simple (struct hybrid_elem *elem ) {
//AVOID
	
	int i,j,l,k,cvindex;
	cvindex=elem->cvindex;
		
	elem->ref_dist[0]=colvar.ss0[cvindex];

	if(elem->has_jac_alloc){
		for(j=0;j<elem->ncv;j++){
			for(l=0;l<elem->natoms_per_cv[j];l++){
				// cvder has: ncv/atom/coord depencency
				for(k=0;k<3;k++)elem->cvder[j][l][k]=colvar.myder[cvindex][l][k];
				//fprintf(mtd_data.fplog,"DER SIMPLE: %12.6f %12.6f %12.6f \n", elem->cvder[j][l][0],elem->cvder[j][l][1],elem->cvder[j][l][2]);
			}
		}
	}else{
		fprintf(mtd_data.fplog,"NO JACOBIAN DEFINED FOR THIS ELEMENT: CANNOT RETRIEVE THE DERIVATIVE RESPECT ATOMS\n");
		EXIT();
	}
	
//	for(i=0;i<colvar.natoms[cvindex];i++){
//		elem->myder[i][0]=colvar.myder[cvindex][i][0]; 
//		elem->myder[i][1]=colvar.myder[cvindex][i][1]; 
//		elem->myder[i][2]=colvar.myder[cvindex][i][2]; 
//	} 
	 
//AVOID
}
void  PREFIX  hbd_copy_msd ( struct hybrid_elem *elem ) {
//AVOID
	
	int i,j,k,l,cvindex;
	real totalign,wal;
	cvindex=elem->cvindex;
		
	
	// instead of this
	// elem->ref_dist=colvar.ss0[cvindex];
	// copy atoms of the pdb into place
	// by exploiting an unuseful vector in dummy 
	
	int iat;
	for(i=0;i<colvar.natoms[cvindex];i++){ 
		elem->ref_dist[3*i+0]=colvar.myder[cvindex][i][0];
		elem->ref_dist[3*i+1]=colvar.myder[cvindex][i][1];
		elem->ref_dist[3*i+2]=colvar.myder[cvindex][i][2];
	}
	
	// put the element derivative to zero for all the atoms and components
	for(j=0;j<elem->ncv;j++){
		for(l=0;l<elem->natoms_per_cv[j];l++){
			for(k=0;k<3;k++)elem->cvder[j][l][k]=0.;
		}
	}
		
	totalign=colvar.rmsd_container[cvindex].walign;
	
	for(j=0;j<elem->ncv;j++){
		
		// this index runs on the displaced atom
		for(l=0;l<elem->natoms_per_cv[j];l++){
			k=j%3;
			// specific weight of this atom
			wal=(colvar.rmsd_container[cvindex].align[l]);
						
			elem->cvder[j][l][k] = -wal/totalign;
			
			if(l==j/3)elem->cvder[j][l][k]+=1;
			
			
		}

	}	
	return;
	
//AVOID
}


int PREFIX rmsd_mini_pack_fake(struct rmsd_inpack inpack,struct rmsd_mini_outpack *work, int nocenter, int simple)
{
/* declarations */
int i,j,k,l,ll,mm,ii;
real rrsq,xx,yy,zz,m[4][4],rr1[4],rr0[4];
real lambda[4],s,q[4];
real pi1[3][3],pi0[3][3],tmp1,dnatoms; 
xx=0.;
yy=0.;
zz=0.;
tmp1=0.;
//         printf("FAKE \n");
dnatoms=(inpack.natoms);
if(!nocenter) {// if you dont need to center no prob...
	for(i=0;i<inpack.natoms;i++){
		xx+=inpack.r0[0][i]*inpack.mass[i];
		yy+=inpack.r0[1][i]*inpack.mass[i];
		zz+=inpack.r0[2][i]*inpack.mass[i];
	}
        xx=xx/(inpack.totmass);
        yy=yy/(inpack.totmass);
        zz=zz/(inpack.totmass);

}
work->cmr0[0]=xx;
work->cmr0[1]=yy;
work->cmr0[2]=zz;
for(i=0;i<inpack.natoms;i++){
	work->r0p[0][i]=inpack.r0[0][i]-xx;
	work->r0p[1][i]=inpack.r0[1][i]-yy;
	work->r0p[2][i]=inpack.r0[2][i]-zz;
}
xx=0.;
yy=0.;
zz=0.;
tmp1=0.;
if(!nocenter) { // if you dont need to center no prob...
	for(i=0;i<inpack.natoms;i++){
		xx+=inpack.r1[0][i]*inpack.mass[i];
		yy+=inpack.r1[1][i]*inpack.mass[i];
		zz+=inpack.r1[2][i]*inpack.mass[i];
	};
        xx=xx/(inpack.totmass);
        yy=yy/(inpack.totmass);
        zz=zz/(inpack.totmass);
}
work->cmr1[0]=xx;
work->cmr1[1]=yy;
work->cmr1[2]=zz;
for(i=0;i<inpack.natoms;i++){
	work->r1p[0][i]=inpack.r1[0][i]-xx;
	work->r1p[1][i]=inpack.r1[1][i]-yy;
	work->r1p[2][i]=inpack.r1[2][i]-zz;
}
/*
 * Find the ROTATION matrix
 */
work->d[0][0]=1.0 ; 
work->d[1][0]=0.0 ;
work->d[2][0]=0.0 ;
work->d[0][1]=0.0 ;
work->d[1][1]=1.0 ;
work->d[2][1]=0.0 ;
work->d[0][2]=0.0 ;
work->d[1][2]=0.0 ;
work->d[2][2]=1.0 ;

// error
if(simple){
        work->err=0.;     
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
                   work->err+=inpack.mass[i]*(work->r0p[l][i]-work->r1p[l][i])*(work->r0p[l][i]-work->r1p[l][i]);
                }
        }
        work->err/=(inpack.totmass);
//        printf("ERR %f\n",outpack->err);
/*
 * derivative 
 *
 */
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
                        work->derr_dr1[l][i]=-(2.*inpack.mass[i]/(inpack.totmass))*(work->r1p[l][i]-work->r0p[l][i]); 
                        if(!nocenter){
				for(j=0;j<inpack.natoms;j++){
                                       work->derr_dr1[l][i]+=(2.*inpack.mass[j]/(inpack.totmass))*inpack.mass[i]/inpack.totmass*(work->r1p[l][j]-work->r0p[l][j]);
				}
                        }
//                	printf("DER %d DIR %d : %f\n",i,l,outpack->derr_dr1[l][i]);   
		}	
	}
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
                        work->derr_dr0[l][i]=(2.*inpack.mass[i]/(inpack.totmass))*(work->r1p[l][i]-work->r0p[l][i]); 
                        if(!nocenter){
				for(j=0;j<inpack.natoms;j++){
                                       work->derr_dr0[l][i]-=(2.*inpack.mass[j]/(inpack.totmass))*inpack.mass[i]/inpack.totmass*(work->r1p[l][j]-work->r0p[l][j]);
				}
                        }
//                	printf("DER %d DIR %d : %f\n",i,l,outpack->derr_dr0[l][i]);   
		}	
	}
          
}
for(i=0;i<3;i++){
	for(j=0;j<3;j++){
	    for(k=0;k<3;k++){
                for(l=0;l<inpack.natoms;l++){
                   work->dd_dr1[i][j][k][l]=0.;
                   work->dd_dr0[i][j][k][l]=0.;
                }
            }	
       	}	
}
	
	
return 0;
}
int PREFIX rmsd_findiff_interface(struct rmsd_inpack inpack,struct rmsd_mini_outpack *work){
fprintf(mtd_data.fplog,"Entering rmsd finite difference test system\n");
fprintf(mtd_data.fplog,"-------------------------------------------\n");
fprintf(mtd_data.fplog,"TEST1: derivative of the value (derr_dr0/derr_dr1)\n");
// test 1
int i,j,k,l,m;
real step=1.e-7,olderr,delta; 
real derr_dr1[3][MAXATOMS_RMSD];
real derr_dr0[3][MAXATOMS_RMSD];
real dd_dr0[3][3][3][MAXATOMS_RMSD];
real dd_dr1[3][3][3][MAXATOMS_RMSD];
real oldd[3][3];
// get initial value of the error and derivative of it 
rmsd_mini_pack(inpack,work,7,1,0);
fprintf(mtd_data.fplog,"INITIAL ERROR VALUE: %f FOR %d ATOMS\n",work->err,inpack.natoms);
olderr=work->err;
// store the derivative
for(j=0;j<3;j++){
for(i=0;i<inpack.natoms;i++){
derr_dr1[j][i]=work->derr_dr1[j][i];
derr_dr0[j][i]=work->derr_dr0[j][i];
}
}
rmsd_mini_pack(inpack,work,7,0,0);
for(l=0;l<3;l++){
for(m=0;m<3;m++){
oldd[l][m]=work->d[l][m];
for(j=0;j<3;j++){
for(i=0;i<inpack.natoms;i++){
dd_dr1[l][m][j][i]=work->dd_dr1[l][m][j][i];
dd_dr0[l][m][j][i]=work->dd_dr0[l][m][j][i];
}
}
}
}
fprintf(mtd_data.fplog,"TESTING: derr_dr1 \n");
for(j=0;j<3;j++){
   for(i=0;i<inpack.natoms;i++){
       // random displacement
       delta=(drand48()-0.5)*2*step;
       inpack.r1[j][i]+=delta; 
       rmsd_mini_pack(inpack,work,7,2,0);
       inpack.r1[j][i]-=delta; 
       switch(j){
         case 0:
            fprintf(mtd_data.fplog,"TESTING: X  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr1[j][i],(work->err-olderr)/delta,derr_dr1[j][i]-(work->err-olderr)/delta);break;
         case 1:
            fprintf(mtd_data.fplog,"TESTING: Y  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr1[j][i],(work->err-olderr)/delta,derr_dr1[j][i]-(work->err-olderr)/delta);break;
         case 2:
            fprintf(mtd_data.fplog,"TESTING: Z  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr1[j][i],(work->err-olderr)/delta,derr_dr1[j][i]-(work->err-olderr)/delta);break;
      
       }    
   }
}
fprintf(mtd_data.fplog,"TESTING: derr_dr0 \n");
for(j=0;j<3;j++){
   for(i=0;i<inpack.natoms;i++){
       // random displacement
       delta=(drand48()-0.5)*2*step;
       inpack.r0[j][i]+=delta; 
       rmsd_mini_pack(inpack,work,7,2,0);
       inpack.r0[j][i]-=delta; 
       switch(j){
         case 0:
            fprintf(mtd_data.fplog,"TESTING: X  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr0[j][i],(work->err-olderr)/delta,derr_dr0[j][i]-(work->err-olderr)/delta);break;
         case 1:
            fprintf(mtd_data.fplog,"TESTING: Y  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr0[j][i],(work->err-olderr)/delta,derr_dr0[j][i]-(work->err-olderr)/delta);break;
         case 2:
            fprintf(mtd_data.fplog,"TESTING: Z  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr0[j][i],(work->err-olderr)/delta,derr_dr0[j][i]-(work->err-olderr)/delta);break;
       }    
   }
}
fprintf(mtd_data.fplog,"TESTING: dd_dr0 \n");
for(l=0;l<3;l++){
  for(m=0;m<3;m++){
    for(j=0;j<3;j++){
       for(i=0;i<inpack.natoms;i++){
           // random displacement
           delta=(drand48()-0.5)*2*step;
           inpack.r0[j][i]+=delta; 
           rmsd_mini_pack(inpack,work,7,2,0);
           inpack.r0[j][i]-=delta; 
           switch(j){
             case 0:
                fprintf(mtd_data.fplog,"TESTING: DD_DR0 [ %d ][ %d ]:  X %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr0[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr0[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
             case 1:
                fprintf(mtd_data.fplog,"TESTING: DD_DR0 [ %d ][ %d ]:  Y %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr0[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr0[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
             case 2:
                fprintf(mtd_data.fplog,"TESTING: DD_DR0 [ %d ][ %d ]:  Z %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr0[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr0[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;


           }    
       }
    }
  }
}
fprintf(mtd_data.fplog,"TESTING: dd_dr1 \n");
for(l=0;l<3;l++){
  for(m=0;m<3;m++){
    for(j=0;j<3;j++){
       for(i=0;i<inpack.natoms;i++){
           // random displacement
           delta=(drand48()-0.5)*2*step;
           inpack.r1[j][i]+=delta; 
           rmsd_mini_pack(inpack,work,7,2,0);
           inpack.r1[j][i]-=delta; 
           switch(j){
             case 0:
                fprintf(mtd_data.fplog,"TESTING: DD_DR1 [ %d ][ %d ]:  X %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr1[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr1[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
             case 1:
                fprintf(mtd_data.fplog,"TESTING: DD_DR1 [ %d ][ %d ]:  Y %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr1[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr1[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
             case 2:
                fprintf(mtd_data.fplog,"TESTING: DD_DR1 [ %d ][ %d ]:  Z %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr1[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr1[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;


           }    
       }
    }
  }
}
EXIT();
return 0 ;
};
void PREFIX pathref_findiff(int i_c, struct mtd_data_s *mtd_data){
fprintf(mtd_data->fplog,"Entering PATH finite difference test system\n");
fprintf(mtd_data->fplog,"-------------------------------------------\n");
fprintf(mtd_data->fplog,"TEST1: derivative of the value (dpathvar/d_r_ref)\n");
// test 1
int i,j,k,l,m,nneigh,ii;
struct sz_data *pmy_sz;
real step=1.e-9,oldpath,delta,newpath; 
// retrieve and store the derivative
real dpath_dr0[3][MAXATOMS_RMSD][MAXFRAMES_PATH];
pmy_sz=&my_sz_list[ic_to_sz[i_c]];
oldpath=colvar.ss0[i_c];
nneigh=pmy_sz->nneigh;
for(j=0;j<colvar.natoms[i_c];j++){
       for(ii=0;ii<pmy_sz->number;ii++){
            dpath_dr0[0][j][ii]=0.;
            dpath_dr0[1][j][ii]=0.;
            dpath_dr0[2][j][ii]=0.;
       }
       for(ii=0;ii<nneigh;ii++){
            i=pmy_sz->lneigh[ii];
            dpath_dr0[0][j][i]=pmy_sz->dpath_dr[0][j][i];
            dpath_dr0[1][j][i]=pmy_sz->dpath_dr[1][j][i];
            dpath_dr0[2][j][i]=pmy_sz->dpath_dr[2][j][i];
       }
}

// for each frame: 
for(ii=0;ii<nneigh;ii++){
         i=pmy_sz->lneigh[ii];
         // random displacement
         delta=(drand48()-0.5)*2*step;
         // for each atom
         for(k=0;k<pmy_sz->frameset[i]->natoms;k++){         
         for(l=0;l<3;l++){         
// change a bit the reference 
             pmy_sz->frameset[i]->pos[k][l]+=delta; 
             if(colvar.type_s[i_c]==30)spath_restraint(i_c, mtd_data); 
             if(colvar.type_s[i_c]==31)zpath_restraint(i_c, mtd_data); 
             newpath=colvar.ss0[i_c];
// recalculate the variable 
             pmy_sz->frameset[i]->pos[k][l]-=delta; 
             switch(l){ 
                case 0:  fprintf(mtd_data->fplog,"TESTING: NFR %d X %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,k,dpath_dr0[l][k][i],(newpath-oldpath)/delta,dpath_dr0[l][k][i]-(newpath-oldpath)/delta);break;
                case 1:  fprintf(mtd_data->fplog,"TESTING: NFR %d Y %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,k,dpath_dr0[l][k][i],(newpath-oldpath)/delta,dpath_dr0[l][k][i]-(newpath-oldpath)/delta);break;
                case 2:  fprintf(mtd_data->fplog,"TESTING: NFR %d Z %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,k,dpath_dr0[l][k][i],(newpath-oldpath)/delta,dpath_dr0[l][k][i]-(newpath-oldpath)/delta);break;
             }
// print the output 
         } 
         } 
}
fprintf(mtd_data->fplog,"exiting PATH test system\n");
EXIT();
}; 
void  PREFIX sbernd_restraint(int i_c, struct mtd_data_s *mtd_data){
//AVOID
        int iat,i,ii,j;
        real tmp1;
        struct coordinates_frameset *pmy_coord1;
        struct hybrid_frameset *hbd_pmy;
        struct sz_data *pmy_sz;
        struct cmap_inpack inpack;
        struct cmap_outpack work;
        real *save_err; 
        int npts;

        int isign,ipntmin,ipntmin2,ipntmin3,ncv;
        real v1v1,v2v2,v3v3,v1v2;
        real *dv1,*dv2,*dv3,*dv4;
        real *dv1x,*dv2x,*dv3x,*dv4x;

        pmy_sz=&my_sz_list[ic_to_sz[i_c]];
        // suppress the neighlist 
        npts=pmy_sz->number;
        save_err=(real *)malloc(pmy_sz->number*sizeof(real));

        // initialize the vectors for evoulution

        if(mtd_data->istep==0 && logical.debug!=1){
              fprintf(mtd_data->fplog,"|- INITIALIZING ENSING PATHWAY ....\n");
              init_bernd_evolution(&pmy_sz,i_c);  
              fprintf(mtd_data->fplog,"|- INITIALIZING ENSING PATHWAY EVOLUTION DONE!\n");
        }  

        if(strcmp(pmy_sz->path_type,"HYBRID") == 0){ 
            //retrive values of cv and derivatives (one set)
            // assuming they are all the same for all the framesets (generally it is) 
            // in case of distance where the comp weight is on the cv evaluation and not on the 
            // distance acquire the points in CV space 
            // eg: dist between two atoms-> collect the cvval (which is the coordinate) and the  
            // eg: rmsd between two struct-> collect the actual coordinates (which is the coordinate) 
            // note: if D is the distance function from the reference  
            // d D(r_0,r)/dx= dD(r_0,r)/dr dr/dx 
            //   
            // in case of rmsd there is no need for chain rule as the derivative is directly calculated with quaternion
            // d D(r_0,r)/dx via quaternion   (it's like dr/dx=1)
            //
            //
            hbd_collect_config(pmy_sz->hbd_running);
			
			hbd_collect_jacobian(pmy_sz->hbd_running,pmy_sz->mathybrid,pmy_sz->myinverse,mtd_data->fplog,1,1);

			if(pmy_sz->debug_metrics)test_hbd_metrics_new(pmy_sz->hbd_running,pmy_sz->hbd_frameset[0],&work,pmy_sz->mathybrid,mtd_data->fplog);
        }else{ 
            plumed_error("ENSING PATH IS IMPLEMENTED ONLY IN HYBRID MODE ");
        }

        // 
        // calculate all-with-all distance
        // 
        for(i=0;i< npts;i++){
                pmy_sz->lneigh[i]=i;
                if(strcmp(pmy_sz->path_type,"HYBRID") == 0){
                   hbd_pmy=pmy_sz->hbd_frameset[i];
				   hbd_metrics_new(pmy_sz->hbd_running,pmy_sz->hbd_frameset[i],&work,pmy_sz->mathybrid,mtd_data->fplog);
                }

                //fprintf(mtd_data->fplog,"ERR %d %f \n",i,outpack.err);
                //fflush(mtd_data->fplog);
                save_err[i]=work.err;

		}

        // sort the closest 
        //for(i=0;i<npts;i++)printf("BEFORE SORTING %d %f\n",pmy_sz->lneigh[i],save_err[i]);
        realquicksort(save_err,pmy_sz->lneigh,0,npts-1);
        //for(i=0;i<npts;i++)printf("AFTER SORTING %d %f\n",pmy_sz->lneigh[i],save_err[i]);

        // find the indexes and the sign
        ipntmin=pmy_sz->lneigh[0];
        if(pmy_sz->lneigh[1]-pmy_sz->lneigh[0] > 0){ ipntmin2=ipntmin+1 ;}
        else{ipntmin2=ipntmin-1;}  
        isign=ipntmin-ipntmin2;
        ipntmin3=ipntmin+isign;
        //fprintf(mtd_data->fplog,"ipntmin %d ipntmin2 %d inpntmin3 %d\n",ipntmin,ipntmin2,ipntmin3);


        // case hybrid 
        // malloc the vectors
        if(strcmp(pmy_sz->path_type,"HYBRID") == 0){
                 real *v1,*v2,*v3,*p,*s_ipntmin,*s_ipntmin2,*s_ipntmin3;
                 real *dummy1,*dummy2,*dummy3,*dummy4,*dummy5; 
                 // difference of vectors: beware: this influence the whole distance... 
                 ncv=(pmy_sz->hbd_running)->hbd_totcvs;
                 // calculate v1v1 dot prod and derivative on the spot 
                 real *dv1v1_ds1,*dv1v1_dp,*dv1v1_dpref;
				 snew(dv1v1_ds1,ncv);
				 snew(dv1v1_dp,ncv);
				 snew(dummy1,ncv);
				 snew(dummy2,ncv);
				 snew(dummy3,ncv);
				 snew(dummy4,ncv);
     			 snew(dummy5,ncv);


				 // (v_iptnmin-v_run)_v_iptnmin MAT (v_iptnmin-v_run)_v_iptnmin and derivatives:
			     struct hybrid_frameset *f1,*f2,*fref;
			
			     f1=pmy_sz->hbd_frameset[ipntmin];
				 f2=pmy_sz->hbd_running;
				 fref=pmy_sz->hbd_frameset[ipntmin];
			
			     v1v1=hbd_vecmvec_ref(f1,f2,fref,f1,f2,fref,
								      pmy_sz->mathybrid,
									  dv1v1_ds1,dv1v1_dp,dummy1,dummy2,dummy3,dummy4,
									  mtd_data->fplog);

			
                 for(i=0;i<ncv;i++)dv1v1_dp[i]+=dummy3[i]; // correct for the reference alignment
			     for(i=0;i<ncv;i++)dv1v1_ds1[i]+=dummy1[i]+dummy2[i]+dummy4[i];
			
			
                 real *dv3v3_dp,*dv3v3_ds2;
                 // calculate v3v3 dot prod and derivative on the spot
				 snew(dv3v3_dp,ncv);
				 snew(dv3v3_ds2,ncv);
			
				 // dummy vectors are overwritten. No need to realloc.
			
				 f1=pmy_sz->hbd_frameset[ipntmin2];
				 f2=pmy_sz->hbd_running;
				 fref=pmy_sz->hbd_frameset[ipntmin2];
		
				 v3v3=hbd_vecmvec_ref(f1,f2,fref,f1,f2,fref,
								      pmy_sz->mathybrid,
									  dv3v3_ds2,dv3v3_dp,dummy1,dummy2,dummy3,dummy4,
									  mtd_data->fplog);
			
				 for(i=0;i<ncv;i++)dv3v3_dp[i]+=dummy3[i]; // correct for the reference alignment
				 for(i=0;i<ncv;i++)dv3v3_ds2[i]+=dummy1[i]+dummy2[i]+dummy4[i];

                 real *dv1v2_dp,*dv1v2_ds1,*dv1v2_ds3,*dv1v2_ds2;
                 real *dv2v2_ds1,*dv2v2_ds2;
				 snew(dv1v2_ds1,ncv);
		     	 snew(dv1v2_ds2,ncv);
				 snew(dv1v2_ds3,ncv);
			     snew(dv1v2_dp,ncv);


			     snew(dv2v2_ds2,ncv); 
				 snew(dv2v2_ds1,ncv); 
			
                 if((ipntmin3<0) || (ipntmin3>=npts)) {
					 
                    // v2 = s1[ipntmin] - s2 [ipntmin2]
					 
					f1=pmy_sz->hbd_frameset[ipntmin];
					f2=pmy_sz->hbd_frameset[ipntmin2];
					fref=pmy_sz->hbd_frameset[ipntmin]; 
					 

					 
					v2v2=hbd_vecmvec_ref(f1,f2,fref,f1,f2,fref,
									 pmy_sz->mathybrid,
								     dv2v2_ds1,dv2v2_ds2,dummy1,dummy2,dummy3,dummy4,
									 mtd_data->fplog);
					 
					 for(i=0;i<ncv;i++)dv2v2_ds2[i]+=dummy3[i]; // correct for the reference alignment
					 for(i=0;i<ncv;i++)dv2v2_ds1[i]+=dummy1[i]+dummy2[i]+dummy4[i];
				 

					 // this is more tricky: need full expression for projection
					 

					 
					 v1v2=hbd_vecmvec_ref(pmy_sz->hbd_frameset[ipntmin],pmy_sz->hbd_running,pmy_sz->hbd_frameset[ipntmin],
										  pmy_sz->hbd_frameset[ipntmin],pmy_sz->hbd_frameset[ipntmin2],pmy_sz->hbd_frameset[ipntmin],
										  pmy_sz->mathybrid,
										  dv1v2_ds1,dv1v2_dp,dummy1,dummy2,dummy3,dummy4,
										  mtd_data->fplog);
					 
					 // this has derivatives in many respects
					 
					 for(i=0;i<ncv;i++)dv1v2_ds1[i]+=dummy1[i]+dummy2[i]+dummy4[i];  
					 for(i=0;i<ncv;i++)dv1v2_ds2[i]+=dummy3[i];  
					 
					 
                 } else {
                    // v2 = s3[ipntmin3] - s1 [ipntmin]
					 					
					 f1=pmy_sz->hbd_frameset[ipntmin3];
					 f2=pmy_sz->hbd_frameset[ipntmin];
					 fref=pmy_sz->hbd_frameset[ipntmin3]; 
					 
					 v2v2=hbd_vecmvec_ref(f1,f2,fref,f1,f2,fref,
									  pmy_sz->mathybrid,
									  dv2v2_ds1,dv2v2_ds2,dummy1,dummy2,dummy3,dummy4,
									  mtd_data->fplog);
					 
					 for(i=0;i<ncv;i++)dv2v2_ds2[i]+=dummy3[i]; // correct for the reference alignment
					 for(i=0;i<ncv;i++)dv2v2_ds1[i]+=dummy1[i]+dummy2[i]+dummy4[i];
					 
					 v1v2=hbd_vecmvec_ref(pmy_sz->hbd_frameset[ipntmin],pmy_sz->hbd_running,pmy_sz->hbd_frameset[ipntmin],
										  pmy_sz->hbd_frameset[ipntmin3],pmy_sz->hbd_frameset[ipntmin],pmy_sz->hbd_frameset[ipntmin3],
										  pmy_sz->mathybrid,
										  dv1v2_ds1,dv1v2_dp,dummy1,dv1v2_ds3,dummy2,dummy3,
										  mtd_data->fplog);
					 
					 for(i=0;i<ncv;i++)dv1v2_ds1[i]+=dummy1[i]+dummy1[i]+dummy2[i];
					 for(i=0;i<ncv;i++)dv1v2_ds3[i]+=dummy3[i];
					 
					 
                    } 
                 // now I have all the projection needed and derivatives as well 

                 real root,dx,rc; 
                 root=sqrt(v1v2*v1v2-v2v2*(v1v1-v3v3));
                 dx=0.5*((root-v1v2)/v2v2-1.);

                 rc=ipntmin+isign*dx              ; 
                 //fprintf(mtd_data->fplog,"RC %f\n",rc); 
                 colvar.ss0[i_c]=rc; 

                 // some derivatives respect to  each cv 
                 real *dercv=(real *)calloc(ncv,sizeof(real));
                 real prefactor=0.5/(v2v2*root);
			
				 // apply the chain rule
                 //  the total derivative should be the derivative per each cv times the one per atom    

				 int atom_offset,cv_offset,k;
				 struct hybrid_elem *myelem;
					
				 atom_offset=0;
				 cv_offset=0;	
				
			     for(i=0;i<pmy_sz->hbd_running->hbd_totalatoms;i++){
					 colvar.myder[i_c][i][0]=0.;colvar.myder[i_c][i][1]=0.;colvar.myder[i_c][i][2]=0.;
				 }
					
				 for(i=0;i<pmy_sz->hbd_running->hbd_nelem;i++){
						myelem=pmy_sz->hbd_running->hbd_elem[i];
						for(j=0;j<myelem->ncv;j++){

							dercv[cv_offset+j]=prefactor*(2.*v1v2*dv1v2_dp[cv_offset+j]-v2v2*dv1v1_dp[cv_offset+j]+v2v2*dv3v3_dp[cv_offset+j])-dv1v2_dp[cv_offset+j]/v2v2;
							dercv[cv_offset+j]=isign*0.5*dercv[cv_offset+j]; 

							for(k=0;k<myelem->nder;k++){
								//
								//	 fprintf(fplog,"|- hbd_metrics_new: NELEM %d NOW CV %d ATOM %d TOUCHING %d\n",i,j,k,atom_offset+k);
								
								colvar.myder[i_c][atom_offset+k][0]+=myelem->cvder[j][k][0]*dercv[cv_offset+j];	
								colvar.myder[i_c][atom_offset+k][1]+=myelem->cvder[j][k][1]*dercv[cv_offset+j];	
								colvar.myder[i_c][atom_offset+k][2]+=myelem->cvder[j][k][2]*dercv[cv_offset+j];	
							}
						}
						
						atom_offset+=myelem->nder; // sum up the number of cvs
						cv_offset+=myelem->ncv; // sum up the number of cvs
										
				 }
						
                 // evolution: move the frames direction of the gradient with a magnitude 
                 // depending on 
                 if(pmy_sz->ievol>0){
                     real v4v4,v1v4;
                     //fprintf(mtd_data->fplog,"|-CALCULATING FORCE ON THE FRAMES (S)...\n"); 
                     // calculate the force perpendicular to the path  
                     // bernd ensing style distance       
                     // calculate v1v1 dot prod and derivative on the spot 
					 
                     real *dv4v4_ds1,*dv4v4_ds2;
					 snew(dv4v4_ds1,ncv);
					 snew(dv4v4_ds2,ncv);

					 f1=pmy_sz->hbd_frameset[ipntmin];
					 f2=pmy_sz->hbd_frameset[ipntmin2];
					 fref=pmy_sz->hbd_frameset[ipntmin];
					 
					 v4v4=hbd_vecmvec_ref(f1,f2,fref,f1,f2,fref,
										  pmy_sz->mathybrid,
										  dv4v4_ds1,dv4v4_ds2,dummy1,dummy2,dummy3,dummy4,
										  mtd_data->fplog);
					 
					 for(i=0;i<ncv;i++)dv4v4_ds2[i]+=dummy3[i]; // correct for the reference alignment
					 for(i=0;i<ncv;i++)dv4v4_ds1[i]+=dummy1[i]+dummy2[i]+dummy4[i];
                     
					 real *dv1v4_dp,*dv1v4_ds1,*dv1v4_ds2;
					 snew(dv1v4_dp,ncv);
					 snew(dv1v4_ds1,ncv);
					 snew(dv1v4_ds2,ncv);					 
					 
					 v1v4=hbd_vecmvec_ref(pmy_sz->hbd_frameset[ipntmin],pmy_sz->hbd_running,pmy_sz->hbd_frameset[ipntmin],
										  pmy_sz->hbd_frameset[ipntmin2],pmy_sz->hbd_frameset[ipntmin],pmy_sz->hbd_frameset[ipntmin2],
										  pmy_sz->mathybrid,
										  dv1v4_ds1,dv1v4_dp,dummy1,dv1v4_ds2,dummy2,dummy3,
										  mtd_data->fplog);
					 
					 for(i=0;i<ncv;i++)dv1v4_ds1[i]=dummy1[i]+dummy2[i];
					 for(i=0;i<ncv;i++)dv1v4_ds2[i]=dummy3[i];
					 
                     // now I have all the projection needed and derivatives as well 
                     rc=(v1v1+dx*dx*v4v4-2.*dx*v1v4);
                 
                     // some derivatives respect to  each cv 
                     real ddx_dp;
                     //real *dercv;
                     //dercv=(real *)calloc(ncv,sizeof(real));
                     real prefactor=0.5/(v2v2*root);
                     // calculate the weight  
                     //if(dx<-0.5)dx=-0.5; 
                     //if(dx>0)dx=0.; 
                     //fprintf(mtd_data->fplog,"DXS %f RC %f\n",dx,rc); 
                     if(dx<-0.5){
                          tmp1=-0.5;
                     }else if(dx>0.){ 
                          tmp1=0.;
                     }else { tmp1=dx;}
                     real w2=-1.*tmp1;
                     real w1=1.+tmp1;
                     // fade the weight and displacement  away: memory
                     // fade=0: forget at each step
                     // fade=1: always remember everything from past history
                     // should be tuned on the recrossing 
                     pmy_sz->wcount[ipntmin]  *=pmy_sz->fadefactor;
                     pmy_sz->wcount[ipntmin2] *=pmy_sz->fadefactor;
                     for(i=0;i<ncv;i++){
                       pmy_sz->disp[ipntmin][i] *=pmy_sz->fadefactor; 
                       pmy_sz->disp[ipntmin2][i]*=pmy_sz->fadefactor; 
                     }
                     // add the new weight 
                     pmy_sz->wcount[ipntmin] +=w1;
                     pmy_sz->wcount[ipntmin2]+=w2;
                     // stupid workaround for not looping externally on the elements 
                     // displacement[i]=dist*grad(dist)_i/|grad(dist)|
                     //fprintf(mtd_data->fplog,"|-S: DX %f W1 %f W2 %f I %d %d\n",dx,w1,w2,ipntmin,ipntmin2); 
                     real sum=0.;
                     for(i=0;i<ncv;i++){
                           ddx_dp=0.5*prefactor*(2.*v1v2*dv1v2_dp[i]-v2v2*dv1v1_dp[i]+v2v2*dv3v3_dp[i])-0.5*dv1v2_dp[i]/v2v2;
                           dercv[i]=(dv1v1_dp[i]-2.*ddx_dp*v1v4-2.*dx*dv1v4_dp[i]+2.*dx*ddx_dp*v4v4 ); 
                           sum+=dercv[i]*dercv[i];
                     }  
                     sum=sqrt(sum);
                     for(i=0;i<ncv;i++){
                           //fprintf(mtd_data->fplog,"|-S: DERCV %d  %f \n",i,dercv[i]); 
                           pmy_sz->disp[ipntmin][i] +=w1*sqrt(rc)*dercv[i]/sum; 
                           pmy_sz->disp[ipntmin2][i]+=w2*sqrt(rc)*dercv[i]/sum; 
                     }
                     free(dv4v4_ds1);
                     free(dv4v4_ds2);
                     free(dv1v4_ds1);
		     free(dv1v4_ds2);            
                     free(dv1v4_dp);  
					 
		     // experimental an preliminary: just a check
	             if( colvar.it%pmy_sz->ievol==0 && colvar.it!=0 ) do_bernd_evolution(pmy_sz); 
					 
					 
                     // calculate average and make displacement if needed  
                    // if( (pmy_sz->ievol>0) && (mtd_data->istep%pmy_sz->ievol==0) && (mtd_data->istep!=0) )do_bernd_evolution(pmy_sz);  
                     //fprintf(mtd_data->fplog,"|-CALCULATING FORCE ON THE FRAMES (S) DONE!\n"); 
                 }
                 free(dv1v1_ds1); 
                 free(dv1v1_dp); 
                 free(dv3v3_dp); 
                 free(dv3v3_ds2); 
                 free(dv1v2_dp); 
		 free(dv1v2_ds1); 
	         free(dv1v2_ds2); 
		 free(dv1v2_ds3); 
                 free(dv2v2_ds1); 
                 free(dv2v2_ds2); 
                 free(dercv); 
		 free(dummy1);
	         free(dummy2);
		 free(dummy3);   
		 free(dummy4);
	         free(dummy5);
        }
          
        free(save_err);
//        fprintf(mtd_data->fplog,"exiting PATH test system\n");
//        EXIT();
        return;
//AVOID
};
// 
//  (v1-v2)_ref_v3 ' M (v4-v5)_ref_v6  is calculated and derivatives therein 
// 
// allo=0 do not allocate: the vectors should be already there 
// allo=1 allocate the vectors for the derivative
//
real PREFIX hbd_vecmvec_ref(struct hybrid_frameset *v1,  struct hybrid_frameset *v2, 
							struct hybrid_frameset *v3,  struct hybrid_frameset *v4,
							struct hybrid_frameset *v5,  struct hybrid_frameset *v6,							
							real **mat , real *dv1dcv,  real *dv2dcv,  real *dv3dcv,  
							real *dv4dcv, real *dv5dcv,  real *dv6dcv, FILE *fplog){
//AVOID
	real dist;
	int i,j,k,l,m;
	struct hybrid_elem *el1,*el2,*el3,*el4,*el5,*el6,*diff1,*diff2;
	struct hybrid_frameset *d1,*d2;
	
	//fprintf(fplog,"|- ENTERING VECMVEC_REF\n");
		
	// allocate a simple difference vector by cloning one structure

	clone_hybrid_frameset(&d1,v1,1,fplog);
	
	//
	// (v1-v2)_ref_v3
	// and derivatives
	//
	
	// loop over all the elements of the frameset and make the differences:
	// in case of RMSD you make the alignment and store the data
	// (matrix derivatives and so on)
	
	
	for(i=0;i< v1->hbd_nelem; i++){
		
		el1=v1->hbd_elem[i];
		el2=v2->hbd_elem[i];
		el3=v3->hbd_elem[i];
		diff1=d1->hbd_elem[i];

		int incr=0 ;// just an increment position for the derivative vector (being updated)
		switch(el1->cvtype){
				// case simple cv: just do the difference (no pbc or other stuff)
				// note: if angles are sin and cos component then this is good anyway
			case 1:	
				simple_elem_diff(diff1,el1,el2,el3,fplog);
				break;
			case 3:	
				simple_elem_diff(diff1,el1,el2,el3,fplog);
				break;
			case 4:	
				simple_elem_diff(diff1,el1,el2,el3,fplog);
				break;
			case 5:	
				simple_elem_diff(diff1,el1,el2,el3,fplog);
				break;
			case 32:	
				simple_elem_diff(diff1,el1,el2,el3,fplog);
				break;
			case 53:	
				// case msd
				//check_msd_elem_diff(diff1,el1,el2,el3,fplog);
				msd_elem_diff(diff1,el1,el2,el3,fplog,0,-1);
				break;
			default:
				fprintf(fplog,"|- THIS DIFFERENCE IS NOT IMPLEMENTED: CVTYPE %d ",el1->cvtype);
				EXIT();
				break;
		}
		
	}
	
	// now calculate the second term
	
	clone_hybrid_frameset(&d2,v4,1,fplog);
		
	//EXIT();
	
	for(i=0;i< v4->hbd_nelem; i++){
		

		el4=v4->hbd_elem[i];
		el5=v5->hbd_elem[i];
		el6=v6->hbd_elem[i];
		diff2=d2->hbd_elem[i];

		int incr=0 ;// just an increment position for the derivative vector (being updated)
		switch(el4->cvtype){
				// case simple cv: just do the difference (no pbc or other stuff)
				// note: if angles are sin and cos component then this is good anyway
			case 1:	
				simple_elem_diff(diff2,el4,el5,el6,fplog); 
				break;
			case 3:	
				simple_elem_diff(diff2,el4,el5,el6,fplog);
				break;
			case 4:	
				simple_elem_diff(diff2,el4,el5,el6,fplog);
				break;
			case 5:	
				simple_elem_diff(diff2,el4,el5,el6,fplog);
				break;
			case 32:	
				simple_elem_diff(diff2,el4,el5,el6,fplog);
				break;
			case 53:	
				// case msd
				//check_msd_elem_diff(diff2,el4,el5,el6,fplog);
				msd_elem_diff(diff2,el4,el5,el6,fplog,0,-1);
				
				break;
			default:
				fprintf(fplog,"|- THIS DIFFERENCE IS NOT IMPLEMENTED: CVTYPE %d ",el4->cvtype);
				EXIT();
				break;
		}
		
	}
	// now calculate the  projection

	int ii,jj,kk,pos1,pos2,oldpos1;
	struct hybrid_elem *first_elem;
	struct hybrid_elem *second_elem;
	real val;
	pos1=0;
	pos2=0;
	int jstart;
	
	
	dist=0.;
    //fprintf(fplog,"|- hbd_vecmvec: ELEM "); 
	for (i=0;i<v1->hbd_nelem;i++){
		first_elem=v1->hbd_elem[i];
		for (ii=0;ii< first_elem->ncv ;ii++){ // loop over all cvs of the first elem
				pos2=0;
// THE loop here below is for simple metrics and not for general projeciton
//			    pos2=pos1;
//				for (j=i;j<v1->hbd_nelem;j++){
//					second_elem=v1->hbd_elem[j];
//					if(j==i){
//						jstart=ii;
//					}else{
//						jstart=0;
//					}
//					for (jj=jstart;jj< second_elem->ncv ;jj++){ // loop over all the cvs of the second
				for (j=0;j<v1->hbd_nelem;j++){
					second_elem=v1->hbd_elem[j];
					for (jj=0;jj< second_elem->ncv ;jj++){	
						//if(mat[pos1][pos2])fprintf(fplog," %10.6f |",d1->hbd_elem[i]->ref_dist[ii]*mat[pos1][pos2]*d2->hbd_elem[j]->ref_dist[jj]);
	//					fprintf(fplog," %10.6f %f %f  |",d1->hbd_elem[i]->ref_dist[ii],mat[pos1][pos2],d2->hbd_elem[j]->ref_dist[jj]);
						dist+=d1->hbd_elem[i]->ref_dist[ii]*mat[pos1][pos2]*d2->hbd_elem[j]->ref_dist[jj];
                           //fprintf(fplog,"XXX %d %d c1 %f c2 %f D %f \n",pos1,pos2,d1->hbd_elem[i]->ref_dist[ii] ,d2->hbd_elem[j]->ref_dist[jj]  ,d1->hbd_elem[i]->ref_dist[ii]*mat[pos1][pos2]*d2->hbd_elem[j]->ref_dist[jj]); 
						pos2++;
					}
				}
				pos1++;
		}
	}
	//fprintf(fplog,"\n");
	//fprintf(fplog,"|- hbd_vecmvec: DIST IS %f \n",dist);
	//fflush(fplog);
	//EXIT();
	/* now calculate the derivative by using the chain rule
	 *
	 * err= sum_i_j d1_i M_i_j d2_j
	 *
	 * the derivatives respect to each cv are 
	 *
	 * derr/d_cv_k= sum_i_j ( d d1_i/d_cv_k M_i_j d2_j  + d1_i M_i_j d2_j/d_cv_k )
	 *
	 * therefore the output must be a single vector which is long v1->totatoms and has the same ordering of v1->backtable
	 * 
	 * differences in cv are retrieved through the diff1.diff_der1 diff1.diff_der2 diff_der3 (for running reference and reference again)
	 *
	 */
	
	// first element derr/d_cv_k with cv_k that belongs to d1: derr/d_cv_k= sum_i_j  d d1_i/d_cv_k M_i_j d2_j
	// this computes dv1dcv and dv2dcv 
	
	// second element derr/d_cv_k with cv_k that belongs to d2: derr/d_cv_k= sum_i_j d1_i M_i_j d d2_j/d_cv_k
	// this computes dv4dcv and dv5dcv
	
	// third element has both derr/d_cv_k derr/d_cv_k= sum_i_j ( d d1_i/d_cv_k M_i_j d2_j  + d1_i M_i_j d2_j/d_cv_k )
	// this computes dv3dcv and dv6dcv
	
	int startval;
	pos1=0;	
	startval=0;
	// full loop
	for (i=0;i<v1->hbd_nelem;i++){
		first_elem=v1->hbd_elem[i];
		for (ii=0;ii< first_elem->ncv ;ii++){ // loop over all cvs of the first elem
			// reset the value
			dv1dcv[pos1]=0.;
			dv2dcv[pos1]=0.;
			dv3dcv[pos1]=0.;
			dv4dcv[pos1]=0.;
			dv5dcv[pos1]=0.;
			dv6dcv[pos1]=0.;
			// do the inner loop for the summation
			pos2=0;

			for (j=0;j<v1->hbd_nelem;j++){
				second_elem=v1->hbd_elem[j];
				for (jj=0;jj< second_elem->ncv ;jj++){ 
					//fprintf(fplog,"|- hbd_vecmvec:-------------------------------\n");
					for(kk=0;kk<first_elem->ncv;kk++){
						

						dv1dcv[pos1]+= d1->hbd_elem[i]->diff_der1[kk][ii] * mat[startval+kk][pos2] * d2->hbd_elem[j]->ref_dist[jj] ;
						dv2dcv[pos1]+= d1->hbd_elem[i]->diff_der2[kk][ii] * mat[startval+kk][pos2] * d2->hbd_elem[j]->ref_dist[jj] ;
						dv3dcv[pos1]+= d1->hbd_elem[i]->diff_der3[kk][ii] * mat[startval+kk][pos2] * d2->hbd_elem[j]->ref_dist[jj] ;
						
						
						dv4dcv[pos1]+= d2->hbd_elem[i]->diff_der1[kk][ii] * mat[startval+kk][pos2] * d1->hbd_elem[j]->ref_dist[jj] ;
						dv5dcv[pos1]+= d2->hbd_elem[i]->diff_der2[kk][ii] * mat[startval+kk][pos2] * d1->hbd_elem[j]->ref_dist[jj] ;
						dv6dcv[pos1]+= d2->hbd_elem[i]->diff_der3[kk][ii] * mat[startval+kk][pos2] * d1->hbd_elem[j]->ref_dist[jj] ;
						
					}
					pos2++;
				}
			}
			// end of the inner loop for the summation
			pos1++;
		}
		startval+=first_elem->ncv;
	}
	
	destroy_hybrid_frameset(d1,fplog);
	destroy_hybrid_frameset(d2,fplog);

	//fprintf(fplog,"|- EXITING VECMVEC_REF\n");
	return dist;
//AVOID
};

// allocate the things for path evolution
void PREFIX init_bernd_evolution  ( struct sz_data **pmy_sz , int i_c){
//AVOID
              int i, ncv,npts; 
              (*pmy_sz)->iforce=1; 
              //
			  ncv=((*pmy_sz)->hbd_frameset[0])->hbd_totcvs;
              npts=(*pmy_sz)->number;
              (*pmy_sz)->wcount=(real *)calloc(npts,sizeof(real)); 
              (*pmy_sz)->disp=(real **)malloc(npts*sizeof(real *));
              for(i=0;i<npts;i++)(*pmy_sz)->disp[i]=(real *)calloc(ncv,sizeof(real)); 
              if((*pmy_sz)->iforce){
                  char string[100];
                  sprintf(string,"evol_path_%d.dat",i_c+1); // waiting for a better name 
                  strcpy((*pmy_sz)->evol_path,string);
                  (*pmy_sz)->fp_evol=fopen((*pmy_sz)->evol_path,"w"); 
              }
 return;
//AVOID
};
void PREFIX do_bernd_evolution  ( struct sz_data *pmy_sz ){
//AVOID
              //
              int i,j,cv_offset,ii;
	      struct hybrid_elem *myelem;
              fprintf(mtd_data.fplog,"|-EVOLVING FRAMES\n");
              int ncv=(pmy_sz->hbd_frameset[0])->hbd_totcvs;
			  int npts=pmy_sz->number;
              fprintf(mtd_data.fplog,"|-XXXX\n");
              fprintf(mtd_data.fplog,"|-MMMM\n");
              for(i=0;i<npts;i++){
                 fprintf(mtd_data.fplog,"|-XXXX ");
                 if(pmy_sz->wcount[i]>0.){
                   for(j=0;j<ncv;j++){
                     fprintf(mtd_data.fplog," %12.6f ",pmy_sz->disp[i][j]/pmy_sz->wcount[i]);
                   }
                 }else{
                   for(j=0;j<ncv;j++){
                     fprintf(mtd_data.fplog," %6.3f ",pmy_sz->wcount[i]);
                   }
                 }
                 fprintf(mtd_data.fplog,"\n");
                 // before the step 
                 fprintf(mtd_data.fplog,"|-MMMM ");
			     cv_offset=0;	
			     for(ii=0;ii<pmy_sz->hbd_running->hbd_nelem;ii++){
						 myelem=pmy_sz->hbd_frameset[i]->hbd_elem[ii];
						 for(j=0;j<myelem->ncv;j++){
							 fprintf(mtd_data.fplog,"%12.6f ",myelem->ref_dist[j]);
						 }
						 cv_offset+=myelem->ncv; //
			     }
				 fprintf(mtd_data.fplog,"\n");
                 //fprintf(pmy_sz->fp_evol," WW %12.6f ",pmy_sz->wcount[i]);
                 if(!pmy_sz->hbd_frameset[i]->fixed &&  (mtd_data.istep!=0) ){
                      if(pmy_sz->wcount[i]>0.){
						  cv_offset=0;	
						  for(ii=0;ii<pmy_sz->hbd_frameset[i]->hbd_nelem;ii++){
							  myelem=pmy_sz->hbd_frameset[i]->hbd_elem[ii];
							  for(j=0;j<myelem->ncv;j++){
								  // 
								  // this REALLY is the important thing:
								  // here sum up the vector of the displacement.
								  // if you do with GRADIENTS respect to the reference frame
								  // then you should not worry about reference system 
								  // (when rototraslations are involved): that's why I prefer this.
								  // the problem with it is that you cannot explore points which are far
								  // from the initial point
								  // to be implemented: the simple calculation of the average as
								  // bernd implemented it, but including rototraslation
								  //
								  myelem->ref_dist[j]+=(pmy_sz->disp[i][cv_offset+j])/(pmy_sz->wcount[i]);
								  pmy_sz->disp[i][cv_offset+j]=0.; 
							  }
							  cv_offset+=myelem->ncv; //
						  }
						  if(pmy_sz->reset==1){
                             pmy_sz->wcount[i]=0.;
						  }
                      }
                 }                 
                 //fprintf(pmy_sz->fp_evol,"\n"); 
                 ///fprintf(pmy_sz->fp_evol,"%12.6f XXX ", pmy_sz->wcount[i]);
              }   
              // here do the reparametrization!!!! 
              // call bernd reparametrization 
              reparam_bernd_path (pmy_sz ,0.1, 30, 1.);

              fprintf(mtd_data.fplog,"|-EVOLVING FRAMES DONE!\n");
              //
 return;
//AVOID
};
void PREFIX reparam_bernd_path ( struct sz_data *pmy_sz , real maxerr, int maxiter, real smooth){
//AVOID
              int i,j,k,npts,nmiddle;
              fprintf(mtd_data.fplog,"|-REPARAMETRIZATION:  .....\n");
              fprintf(mtd_data.fplog,"|-  SMOOTH %f  \n",smooth);
              fprintf(mtd_data.fplog,"|-  MAXITER %d \n",maxiter);
              fprintf(mtd_data.fplog,"|-  MAXERR %f \n",smooth);
              // calculate the ideal fraction for each of the segments 
              // note: if the boundaries are not fixed: the previous and the 
              //       following points then take the ideal distance from  
              //       the closer fixed segment (defined between two fixed points)  
              //  
              npts=pmy_sz->number;  
              fprintf(mtd_data.fplog,"|-NPTS %d \n",npts);
              // mallocs some useful stuff
              real *len,*sumlen,*dr,*sfrac;  
              real *dummy1,*dummy2,*dummy3,*dummy4; // set of dummies that will be trashed 
              real *dd;
              real *tmp;
              int *closest_fixed,first,second,*start_fixed;
              len    =(real *)calloc(npts,sizeof(real));
              sumlen =(real *)calloc(npts,sizeof(real));
              dr     =(real *)calloc(npts,sizeof(real));
              dd     =(real *)calloc(npts,sizeof(real));
              sfrac  =(real *)calloc(npts,sizeof(real));
              tmp=(real *)malloc(npts*sizeof(real)); 
              closest_fixed=(int *)calloc(npts,sizeof(real)); 
              start_fixed=(int *)calloc(npts,sizeof(real));
	
			  // make the allocation of an array containing the differences
			  struct hybrid_frameset **diff_vec;
	          fprintf(mtd_data.fplog,"|-reparam_bernd_path: ALLOCATING THE DIFFERENCE VECTORS\n");
			  diff_vec=(struct hybrid_frameset **)malloc((npts-1)*sizeof(struct hybrid_frameset*));
	          for(i=0;i<npts-1;i++)clone_hybrid_frameset(&diff_vec[i],pmy_sz->hbd_running,1,mtd_data.fplog);
	          real *d1,*d2,*d3; // derivatives of the vector
			  snew(d1,pmy_sz->hbd_running->hbd_totcvs);
			  snew(d2,pmy_sz->hbd_running->hbd_totcvs);
			  snew(d3,pmy_sz->hbd_running->hbd_totcvs);
	
			  // make some temporary vectors of frameset to be copied between one iteration and the other
			  struct hybrid_frameset **temp_vec;
			  fprintf(mtd_data.fplog,"|-reparam_bernd_path: ALLOCATING THE TEMP VECTORS\n");
			  temp_vec=(struct hybrid_frameset **)malloc(npts*sizeof(struct hybrid_frameset*));
			  for(i=0;i<npts;i++)clone_hybrid_frameset(&temp_vec[i],pmy_sz->hbd_frameset[i],0,mtd_data.fplog);
	
	          struct hybrid_frameset *frame;
	          struct hybrid_elem *elem;
			  real olderror;
			  int now_exit;
	          now_exit=0;
	          olderror=1.e5;
	
              int kk,iter;
              ///iter assignment
              iter=maxiter; 
	
//	dump_frameset_unformatted(pmy_sz);
//    EXIT();
              // find ideal distances: these depend only on fixed points: assign ideality index     
              k=0;
	
			  // find the closest fixed frame. 
			  // closest_fixed contains the index from head and tail from which you have to consider fixed 
			  //  Initialize with nonsense value
              for(k=0;k<npts;k++)closest_fixed[k]=-1;
              for(k=0;k<npts;k++)start_fixed[k]=-1;
              // from head
              for(i=0;i<npts-1;i++){
                     if(pmy_sz->hbd_frameset[i]->fixed)break; // exit asa you find the first fixed
                     for(j=i+1;j<npts;j++){
                        if(pmy_sz->hbd_frameset[j]->fixed){closest_fixed[i]=j;break;}  // got the fixed index!
                     }
              } 
              // from tail : go backward
              for(i=npts-1;i>0;i--){
                     if(pmy_sz->hbd_frameset[i]->fixed)break; // exit asa you find the first fixed
                     for(j=i-1;j>=0;j--){
                        if(pmy_sz->hbd_frameset[j]->fixed){closest_fixed[i]=j; break;} // got the fixed index!
                     }
              } 
              for(i=0;i<npts;i++){
                  fprintf(mtd_data.fplog,"|-CLOSEST FIXED %d\n",closest_fixed[i]);
              }
	
	         // in between two fixed frames: this vector tells which is the closest fixed frame from which you haev to start counting
			 // for ideal segments
             j=0; 
             for(i=0;i<npts;i++){
                 if(closest_fixed[i]<0){ // for all the unassigned points 
                     if(pmy_sz->hbd_frameset[i]->fixed)j=i;  // take the index of the first fixed frame
                     start_fixed[i]=j; // take this
                 }  
             } 
             for(i=0;i<npts;i++){
                 fprintf(mtd_data.fplog,"|-START FIXED %d\n",start_fixed[i]);
             }
             //EXIT(); 
              for(kk=0;kk<iter;kk++){ 
				  fprintf(mtd_data.fplog,"\n|-ITER %d\n",kk);
				  fflush(mtd_data.fplog);
                  for(i=0;i<npts-1;i++){ // loop over all the possible distances
						 
						// if multiple matrices are defined then use them 
						// this let point mathybrid to the right block
						if(pmy_sz->mathybrid_blocks!=NULL)point_to_matrix(i,pmy_sz);
					  
						// the difference here uses the frame i as reference
					    // so diff[0] is the diff betwen frame 1 and 0 centered on 0
							
                        dd[i]=hbd_distanddiff_ref(pmy_sz->hbd_frameset[i],pmy_sz->hbd_frameset[i+1],pmy_sz->hbd_frameset[i],pmy_sz->mathybrid,
												  d1,d2,d3,diff_vec[i],mtd_data.fplog);
						fprintf(mtd_data.fplog,"|-DISTSQUARED %3d TO %3d : %12.6f\n",i,i+1,dd[i]);
						fflush(mtd_data.fplog);

						dd[i]=sqrt(dd[i]);

                  }

                 // assign the ideal segment lenght and direction for each of the points only if enclosed between two fixed points 
                  real prevsum=0.;
                  nmiddle=0; 
                  for(i=0;i<npts;i++){
                      if(closest_fixed[i]<0){ // choose the points which are in between fixed (fixed pts included)
						  if(pmy_sz->hbd_frameset[i]->fixed && sumlen[i-1]>0.){ // got to one fixed point and the sum already cointains some acc value
							  fprintf(mtd_data.fplog,"|-SS I %d NMIDDLE %d  \n",i,nmiddle); 
							  k=1;
							  for(j=i-nmiddle;j<i;j++){
								  dr[j]=sumlen[i-1]/(nmiddle); 
								  
								  sfrac[j]=k*sumlen[i-1]/(nmiddle);k++;

								  fprintf(mtd_data.fplog,"|-SFRAC I %d  : %f\n",j,sfrac[j]); 
							  }
							  sumlen[i]=dd[i]; // at the end of a segment within fixed point reset the value 
							  nmiddle=0;
						  }else{
							  if(i!=start_fixed[i]){
								  sumlen[i]=sumlen[i-1]+dd[i]; // keep on summing up	
							      fprintf(mtd_data.fplog,"|-POINT %d : %f = %f ( %d ) + %f  ( % d )\n",i,sumlen[i],sumlen[i-1],i-1,dd[i],i);
							  }else{
								  sumlen[i]=dd[i]; 
							  }
						  } 
						  nmiddle++; 
                      }
                  }
				  //EXIT();
                  // sumlen in the other case: from the head
                  prevsum=0.;
                  for(i=npts-1;i>=0;i--){
                       if(closest_fixed[i]>=0 && closest_fixed[i]>i){
                                 sumlen[i]=prevsum+dd[i];
                                 prevsum+=dd[i]; 
                                 //fprintf(mtd_data.fplog,"|-HEAD I %d SUMLEN %f\n",i,sumlen[i]); 
                       }  
                  } 
                  // sumlen in the other case: from the tail 
                  prevsum=0.;
                  for(i=0;i<npts;i++){
                       if(closest_fixed[i]>=0 && closest_fixed[i]<i){
                                 sumlen[i]=prevsum+dd[i-1];
                                 prevsum+=dd[i-1]; 
                                 //fprintf(mtd_data.fplog,"|-TAIL I %d SUMLEN %f\n",i,sumlen[i]); 
                       }  
                  } 
				  fflush(mtd_data.fplog);
                  // assign the ideal segment from the closest fixed point
				  // this could be done differently and extending up to 
				  // certain rules for free energy are respected
                  for(i=0;i<npts;i++){
                      if(closest_fixed[i]>=0){
                         if(closest_fixed[i]>i){
                            dr[i]=dr[closest_fixed[i]]; 
                            sfrac[i]=(closest_fixed[i]-i)*dr[closest_fixed[i]];
                         }else{ 
                            dr[i]=dr[closest_fixed[i]-1];  
                            sfrac[i]=(i-closest_fixed[i])*dr[closest_fixed[i]-1];
                         } 
                      }
                  }  


				  fprintf(mtd_data.fplog,"\n|-Summary of ideality indexes: \n|\n");
				   int nerror=0;
				  real error=0.;
                  for(i=0;i<npts;i++){
					  if(closest_fixed[i-1]==i || closest_fixed[i+1]==i )fprintf(mtd_data.fplog,"|\n");
					  fprintf(mtd_data.fplog,"|-SUMLEN-real I %d  : %f SFRAC-ideal %f DR-ideal %f DD-real %f\n",i, sumlen[i],sfrac[i],dr[i],dd[i]);
					 if(dr[i]>1.e-5){error+=pow2((dr[i]-dd[i])/dr[i]); nerror++;}
                  }
				error/=(float) nerror;
				error=sqrt(error);
				fprintf(mtd_data.fplog,"|-AVERAGE ERROR: %f \n",error); 
				//if (error<0.05){
				if (error<maxerr){
				//if (error>olderror){
				          now_exit=1;
				}else{
				          olderror=error;
				}
				  fprintf(mtd_data.fplog,"|-CUMULATIVE DISTANCE: %f\n",error);
				  fflush(mtd_data.fplog);
				  //EXIT();
				  fprintf(mtd_data.fplog,"\n");
                  // calculate the evolution factors 
                  // note: closest_fixed[i]<0 &&  not fixed :   evolve i+1
                  // note: closest_fixed[i]>i               :   evolve i-1
                  // note: closest_fixed[i]<i               :   evolve i+1
                  // and  evolve a-la-stringmethod -> what about periodicity? (at final step?) 
                  for(i=0;i<npts;i++){
					  fflush(mtd_data.fplog);
					  if(closest_fixed[i]<0 ){ //first the points between two fixed points
						  fflush(mtd_data.fplog);
						  int doit=0;
						  // conditions for which the evolution is not to be carried out
						  if(i+1<npts){
							  doit=1;
							  // do not evolve if the next frame is a fixed
							  if(pmy_sz->hbd_frameset[i+1]->fixed)doit=0;
							  // do not evolve if the next frame belongs to tail
							  if(closest_fixed[i+1]>=0)doit=0;
						  }
						  if(doit){
						  fflush(mtd_data.fplog);
							  real mmin=0.;
							  for(j=start_fixed[i];j<npts;j++){
								  if(sfrac[i]<=sumlen[j] &&  sfrac[i]>mmin){
									  break;
								  } 
								  mmin=sumlen[j];
							  }; //now j contains the index: the new frame should be an interpolation
							  // between j and j-1
							  tmp[i]=smooth*(sfrac[i]-mmin)/dd[j]; // 0.9 factor seems to enhance stability of the algo
							  
							  if(pmy_sz->mathybrid_blocks!=NULL)point_to_matrix(j,pmy_sz);
							  
                              fprintf(mtd_data.fplog,"|-MOVING SEGMENT %d ( %f ) between %d and %d ( %f and %f) prefactor %f : Moving frame %d\n",i,sfrac[i],j,j+1,mmin,sumlen[j],tmp[i],i+1);   
						  fflush(mtd_data.fplog);
							  // 
							  // now call the interpolation and put all in a temporary frameset that is replaced at the end of the round 
							  // 
							  // v1=v2+dr(v3-v4)_v2
							  //                      v1                      v2                                             v3                             v4

							  do_step_bernd(temp_vec[i+1],pmy_sz->hbd_frameset[j],pmy_sz->mathybrid,tmp[i],pmy_sz->hbd_frameset[j+1],pmy_sz->hbd_frameset[j]);

							  fflush(mtd_data.fplog);
						  }
					  } else { // points at the extrema
						  if(closest_fixed[i]>i){ // tail
							  tmp[i]=(sfrac[i]-sumlen[i])/dd[i];
							    fprintf(mtd_data.fplog,"EVOLVING %d ",i); 
							  
							  if(pmy_sz->mathybrid_blocks!=NULL)point_to_matrix(i,pmy_sz);

							  do_step_bernd(temp_vec[i],pmy_sz->hbd_frameset[i],pmy_sz->mathybrid,tmp[i],pmy_sz->hbd_frameset[i],pmy_sz->hbd_frameset[i+1]);
						  }else{ //head
							  tmp[i]=(sfrac[i]-sumlen[i])/dd[i-1];
                                 fprintf(mtd_data.fplog,"EVOLVING %d ",i); 
							  
							  if(pmy_sz->mathybrid_blocks!=NULL)point_to_matrix(i,pmy_sz);
							  
							  do_step_bernd(temp_vec[i],pmy_sz->hbd_frameset[i],pmy_sz->mathybrid,tmp[i],pmy_sz->hbd_frameset[i],pmy_sz->hbd_frameset[i-1]);
						  }
					  }
					  // now make a copyback of all the new ref dist into the old frames

					  
                  }
				  //EXIT();							  

				  for(i=0;i<pmy_sz->number;i++){
					  frame=pmy_sz->hbd_frameset[i];
					  for(j=0;j<frame->hbd_nelem;j++){
						  elem=frame->hbd_elem[j];
						  for(k=0;k<elem->ncv;k++){
//							  fprintf(mtd_data.fplog,"|- FRAME %3d ELEM %3d NCV %3d VAl %f\n",i,j,k,elem->ref_dist[k]);
//							  fflush(mtd_data.fplog);
							  elem->ref_dist[k]=temp_vec[i]->hbd_elem[j]->ref_dist[k];
						  }
						  
					  }	
				  }
				  fprintf(mtd_data.fplog,"|- COPIED NEW FRAMES IN PLACE OF THE OLD ONES!\n");
				  fflush(mtd_data.fplog);
				  //EXIT();
				  if(now_exit)break;

              }
              free(len);
              free(sumlen);
              free(dd);
              free(dr);
              free(sfrac);
              free(tmp);
              free(closest_fixed);
              free(start_fixed);
				
	          // destroy the difference frameset
	          for(i=0;i<npts-1;i++)destroy_hybrid_frameset(diff_vec[i], mtd_data.fplog);
	          free(diff_vec);
	
			  // destroy the temp frameset
			  for(i=0;i<npts;i++)destroy_hybrid_frameset(temp_vec[i], mtd_data.fplog);
	          free(temp_vec);
				
			  dump_frameset_formatted(pmy_sz);
			  //dump_frameset_unformatted(pmy_sz);

	
              fprintf(mtd_data.fplog,"|-REPARAMETRIZATION:  DONE!\n");
//              //EXIT();
//AVOID
}
//// v1=v2+dr(v3-v4)_v2
void PREFIX do_step_bernd (  struct hybrid_frameset *v1 , struct hybrid_frameset *v2, real **mat, real dr, struct hybrid_frameset *v3 ,struct hybrid_frameset *v4 ){
//AVOID
	int i,j;
	struct hybrid_elem *myelem;
	struct hybrid_frameset *diff;
	
	// clone the frame
	clone_hybrid_frameset(&diff,v1,1,mtd_data.fplog);
	real *d1,*d2,*d3;
	snew(d1,v1->hbd_totcvs);
	snew(d2,v2->hbd_totcvs);
	snew(d3,v3->hbd_totcvs);
	
	// calculate the difference in the vectors
	hbd_distanddiff_ref(v3,v4,v2,mat,
						d1,d2,d3,diff,mtd_data.fplog);
	
	// here take into consideration periodicity????
	// take into consideration the reference frame for the distance
	for(i=0;i<v1->hbd_nelem;i++){
		myelem=v1->hbd_elem[i];
		for(j=0;j<myelem->ncv;j++){
			v1->hbd_elem[i]->ref_dist[j]=v2->hbd_elem[i]->ref_dist[j]+
			dr*diff->hbd_elem[i]->ref_dist[j];
			//(v3->hbd_elem[i]->ref_dist[j]-v4->hbd_elem[i]->ref_dist[j]);
		}
	}
	
	destroy_hybrid_frameset(diff, mtd_data.fplog);
	free(d1);
	free(d2);
	free(d3);

//AVOID
};
int PREFIX read_msd( char **word, int count, t_plumed_input *input, FILE *fplog ) {
	int iw,i,help;
	FILE *fp;
	struct pdb *mypdb;
	char *filename;
	help=0;
	iw=seek_word(word,"PDB");
	if(iw>=0) {
		fprintf(fplog,"\n%1i-MSD variable: \n",count+1);
		fprintf(fplog,"|- MSD READER: PDB KEYWORD FOUND\n");
		// open the file and check existence
		filename=word[iw+1];
		fp=fopen(filename,"r");
		if (fp == NULL){
			char buf[1024];
			fprintf(fplog,"|- MSD READER: UNSUCCESSFULL OPENING FOR FILE %s\n",filename);
			EXIT();
		}else{
			fprintf(fplog,"|- MSD READER: NOW OPENING  FILE %s \n",filename);
			//
			// now call a simple pdb reader with streamline and retrieve the pdb infos 	
			// 
			read_pdb(&(colvar.rmsd_container[count].mypdb),fp,fplog);
			//
			// call the allocator to the work array
			// take the elements and check if the rmsd work array is sufficient
			//
			readapt_rmsd_work_array(colvar.rmsd_container[count].mypdb->natoms,fplog);

			fprintf(fplog,"|- MSD READER: PDB FILE HAS %d ATOMS \n",colvar.rmsd_container[count].mypdb->natoms);	

			//	
			// call the allocator for the storage: store the pdb into the rmsd_container
			//
			copy_pdb_into_rmsd_container( *(colvar.rmsd_container[count].mypdb), &(colvar.rmsd_container[count]),fplog);
						
			// copy the pdb into the metadynamics work arrays
			
			colvar.natoms[count]   = colvar.rmsd_container[count].mypdb->natoms;
			
			snew(colvar.myder[count], colvar.natoms[count]);
			snew(colvar.cvatoms[count], colvar.natoms[count]);
			
			for(i=0;i<colvar.natoms[count];i++){
				colvar.cvatoms[count][i] = colvar.rmsd_container[count].mypdb->index[i]-1; // from pdb to c indexing
			}
			//count++;
		}

	}else{ fprintf(fplog,"|- MSD READER:  NEEDED PDB KEYWORD FOR MSD IN HYBRID PATH\n"); help=1;}
	// find the sigma
	real sigma;
	iw=seek_word(word,"SIGMA");
	if(iw>=0){ sscanf(word[iw+1],"%lf", &sigma);
		colvar.delta_r[count]  = (real) sigma; }
	
	iw=seek_word(word,"DUMMY");
	if(iw>=0){ 
		fprintf(fplog,"|- MSD READER: ------------------------------------WARNING!!!!-------------------------------------\n");
		fprintf(fplog,"|- MSD READER: DUMMY VARIABLE, WILL NOT MAKE A REAL DYNAMICS BUT KEEP TRACK OF INDEXES (FOR HYBRID)\n");
		fprintf(fplog,"|- MSD READER: WHATEVER DERIVATIVE WILL BE PUT TO ZERO\n");
		fprintf(fplog,"|- MSD READER: ------------------------------------------------------------------------------------\n");

		colvar.rmsd_container[count].dummy=1;
	}else{
		colvar.rmsd_container[count].dummy=0;	
	}

	iw=seek_word(word,"DREF");
	if(iw>=0){ 
                sscanf(word[iw+1],"%i", & colvar.rmsd_container[count].dref_freq);	
		fprintf(fplog,"|- MSD READER: ------------------------------------------------------------------------------------\n");
		fprintf(fplog,"|- MSD READER: DREF: dump the coordinate respect to each component\n");
		fprintf(fplog,"|- MSD READER:       dumping_frequency is %d \n",colvar.rmsd_container[count].dref_freq);
		fprintf(fplog,"|- MSD READER: ------------------------------------------------------------------------------------\n");
	        char dername[200]; 
		char *str,ic[3],str2[100];
                if(count+1<10){
                 ic[0]='0'+count+1;
                 ic[1]='\0';}
                else if(count+1<100) {
                 ic[0]='0'+(count+1)/10 ;
                 ic[1]='0'+(count+1)%10 ;
                 ic[2]='\0';
                }
                else{
                  plumed_error("|--msd reader: TOO MANY FRAMES REQUIRED FOR NAME BUILDING!");
                }
		strcpy(colvar.rmsd_container[count].dreffile,"dref_");
		strcat(colvar.rmsd_container[count].dreffile,ic);
		strcat(colvar.rmsd_container[count].dreffile,".dat");
		colvar.rmsd_container[count].fpdreffile=fopen(colvar.rmsd_container[count].dreffile,"w");
		colvar.rmsd_container[count].dref=1;
	}else{
		colvar.rmsd_container[count].dref=0;	
	}

	iw=seek_word(word,"NOROT");
	if(iw>=0){ 
		colvar.rmsd_container[count].dorot=0;	
	        fprintf(fplog,"|- MSD:\n");
		fprintf(fplog,"|- NO rotation on the reference system\n");
        }else{
		colvar.rmsd_container[count].dorot=1;	
	        fprintf(fplog,"|- MSD:\n");
		fprintf(fplog,"|- rotation enabled on the reference system\n");
        }

	iw=seek_word(word,"NOTRASL");
	if(iw>=0){ 
		colvar.rmsd_container[count].dotrasl=0;	
	        fprintf(fplog,"|- MSD:\n");
		fprintf(fplog,"|- NO traslation on the reference system\n");
        }else{
		colvar.rmsd_container[count].dotrasl=1;	
	        fprintf(fplog,"|- MSD:\n");
		fprintf(fplog,"|- traslation enabled on the reference system\n");
        }

	if(help){
		fprintf(fplog,"|- MSD:\n");
		fprintf(fplog,"|- This is a fake cv useful to point to the atoms which are \n");
		fprintf(fplog,"|- involved into a path CV with rmsd. It does not have any\n");
		fprintf(fplog,"|- meaning\n");
		fprintf(fplog,"|- It requires a PDB file.\n");
		fprintf(fplog,"|- Position are ignored. Only atom indexes, beta and occupancy matter\n");
		fprintf(fplog,"|- \n");
		EXIT();
	}
	if (logical.do_hills){
		if (colvar.delta_r[count]>0){
			fprintf(fplog,"|- MSD READER: SIGMA %f\n",colvar.delta_r[count]);
       		}
	}
	else fprintf(fplog,"\n");
	
	//EXIT();

	return colvar.natoms[count];
};
// the general pdb parser
int PREFIX read_pdb(struct pdb **mypdb, FILE *myfile, FILE *fplog){
	FILE *initpos;
	fpos_t pos;
	int i,j,k;
	char *str,string[200],end[3],remark[6],atom[5],hetatm[7];
	// save the initial address 
	fgetpos (myfile,&pos);
    // now keep on reading
	// goto up to the bottom of the file up to the "END" or "TER" or "EOF"
	// and count the atoms
	i=0;
	while(1){
	readagain:
		str=fgets(string,180,myfile);
        if(str==NULL)break;
		//fprintf(fplog,"LINE %s \n",string);
		strncpy(end,string,3);
		if(strstr(end,"END")!=NULL){break;}; 
		if(strstr(end,"TER")!=NULL){break;}; 
		strncpy(atom,string,4);atom[4]='\0';
		strncpy(hetatm,string,6);hetatm[6]='\0';
		if( (strstr(atom,"ATOM")!=NULL) || (strstr(hetatm,"HETATM")!=NULL) ) {
			//fprintf(fplog,"FOUND ATOM %s \n",string);
			// do first a count of atoms:
			i++;
		}
		goto readagain;
    } 	
	fprintf(fplog,"|- read_pdb: FOUND %d NATOMS\n",i);
	// allocate structures
	allocate_pdb(mypdb,i);
	fprintf(fplog,"|- read_pdb: ALLOCATED NATOMS %d \n",(*mypdb)->natoms);
	//now go back and assign the elements
	//reset to the initial address
    //return 1;
    
	fsetpos (myfile, &pos);
	i=0;
	char x[12],y[12],z[12],occ[12],beta[12];
//	char ind[5], resid[4];
//	char name[4],resname[3],chain ;
	char ind[10], resid[10];
	char name[10],resname[10],chain ;	
	while(1){
	readagain2:
		str=fgets(string,200,myfile);
	//	fprintf(fplog,"NLINE %s \n",string);
		strncpy(end,string,3);
		if(strstr(end,"END")!=NULL){break;}; 
		if(strstr(end,"TER")!=NULL){break;}; 
		//sscanf(str,"%6s",remark);
		strncpy(atom,string,4);atom[4]='\0';
		strncpy(hetatm,string,6);hetatm[6]='\0';
		if( (strstr(atom,"ATOM")!=NULL) || (strstr(hetatm,"HETATM")!=NULL) ) {
			// now do a serious parsing and assign the elements
//			12345678901234567890123456789012345678901234567890123456789012345678901234567890
//			
//			ATOM   4150  H   ALA A 431       8.674  16.036  12.858  1.00  0.00
//
//			COLUMNS        DATA  TYPE    FIELD        DEFINITION
//			-------------------------------------------------------------------------------------
//			1 -  6        Record name   "ATOM  "
//			7 - 11        Integer       serial       Atom  serial number.
//			13 - 16        Atom          name         Atom name.
//			17             Character     altLoc       Alternate location indicator.
//			18 - 20        Residue name  resName      Residue name.
//			22             Character     chainID      Chain identifier.
//			23 - 26        Integer       resSeq       Residue sequence number.
//			27             AChar         iCode        Code for insertion of residues.
//			31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
//			39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
//			47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
//			55 - 60        Real(6.2)     occupancy    Occupancy.
//			61 - 66        Real(6.2)     tempFactor   Temperature  factor.
//			77 - 78        LString(2)    element      Element symbol, right-justified.
//			79 - 80        LString(2)    charge       Charge  on the atom.
//			
			
//			//now implement a fixed width parser since the c parser does not allows it
			parse_fixwidth(str,  7, 11, ind);
			parse_fixwidth(str, 13, 16, name);
			parse_fixwidth(str, 18, 20, resname);
			parse_fixwidth(str, 23, 26, resid);
			//parse_fixwidth(str, 22, 22, &chain);
			parse_fixwidth(str, 31, 38, x);
			parse_fixwidth(str, 39, 46, y);
			parse_fixwidth(str, 47, 54, z);
			parse_fixwidth(str, 55, 60, occ);
			parse_fixwidth(str, 61, 66, beta);

//			fprintf(fplog,"IND %s NAME %s RES %s ID %s  X %s Y %s Z %s O %s B %s\n",ind,name,resname,resid,x,y,z,occ,beta);
//			
			(*mypdb)->index[i]=atoi(ind);
			(*mypdb)->resid[i]=atoi(resid);
			(*mypdb)->x[i]=atof(x);
			(*mypdb)->y[i]=atof(y);
			(*mypdb)->z[i]=atof(z);
			(*mypdb)->beta[i]=atof(beta);
			(*mypdb)->occ[i]=atof(occ);
			strcpy((*mypdb)->name[i],name);
			strcpy((*mypdb)->resname[i],resname);
//			//(*mypdb).chain[i]=&chain;
//		
//			fprintf(fplog,"S_IND %d NAME %s RES %s ID %d  X %f Y %f Z %f O %f B %f\n",(*mypdb)->index[i],(*mypdb)->name[i],(*mypdb)->resname[i],(*mypdb)->resid[i],(*mypdb)->x[i],(*mypdb)->y[i],(*mypdb)->z[i],(*mypdb)->occ[i],(*mypdb)->beta[i]);
//			
			i++;
		}
		goto readagain2;
    } 		
    fprintf(fplog,"|- read_pdb : READ IN %d NATOMS\n",(*mypdb)->natoms);

	//    deallocate_pdb(*mypdb);
//    fprintf(fplog,"|- PDB READER: DEALLOCATED\n");

	// fprintf(fplog,"ALLOCATED NATOMS %d \n",(*mypdb)->natoms);
	// read in
	
	return 1;
};
// this allocate the structure and report its pointer
int PREFIX allocate_pdb(struct pdb **mypdb, int i){
    int j;
	(*mypdb)=(struct pdb *)malloc(sizeof(struct pdb));
	(*mypdb)->natoms=i;
	(*mypdb)->x=(real *)malloc(i*sizeof(real));
	(*mypdb)->y=(real *)malloc(i*sizeof(real));
	(*mypdb)->z=(real *)malloc(i*sizeof(real));
	(*mypdb)->beta=(real *)malloc(i*sizeof(real));
	(*mypdb)->occ=(real *)malloc(i*sizeof(real));
	(*mypdb)->index=(int *)malloc(i*sizeof(int));
	(*mypdb)->resid=(int *)malloc(i*sizeof(int));
	(*mypdb)->name=(char **)malloc(i* sizeof(char *));
	(*mypdb)->resname=(char **)malloc(i* sizeof(char *));
	(*mypdb)->chain=(char **)malloc(i* sizeof(char *));
	for(j=0;j<i;j++){
	  // allocate chars
		(*mypdb)->name[j]=(char *)malloc(4*sizeof(char));
		(*mypdb)->resname[j]=(char *)malloc(4*sizeof(char));
		(*mypdb)->chain[j]=(char *)malloc(4*sizeof(char));
	}
	//fprintf(mtd_data.fplog,"|- allocate_pdb: ALLOCATING %d ATOMS IN PDB\n",(*mypdb)->natoms);
	return 1;
};
int PREFIX deallocate_pdb(struct pdb *mypdb){
    int i;
	// free all the pointed elements
	//fprintf(mtd_data.fplog,"|- deallocate_pdb: DEALLOCATING %d ATOMS IN PDB\n",(*mypdb).natoms);
	free((*mypdb).index);
	free((*mypdb).resid);
	free((*mypdb).x);
	free((*mypdb).y);
	free((*mypdb).z);
	free((*mypdb).occ);
	free((*mypdb).beta);
	for(i=0;i<(*mypdb).natoms;i++){
		free((*mypdb).name[i]);
		free((*mypdb).resname[i]);
		free((*mypdb).chain[i]);
	}
	free((*mypdb).name);
	free((*mypdb).resname);
	free((*mypdb).chain);
	// free the pointer
	free(mypdb);
	return 1;
};
// this is a fixed width reader with indexing starting from  1
void PREFIX parse_fixwidth(char *str,int i, int j, char *str2){
	int ii,jj,kk;
	kk=0;
	for(ii=i-1;ii<=j-1;ii++){
		str2[kk]=str[ii];
		//strncpy(&str2[kk],&str[ii],1);
		kk++;
	}
	str2[kk]='\0';
	//strcpy(&str2[kk],'\0');
};
/*
 * this allocate/reallocate rmsd work array
 * so that rmsd is all dynamcal
 *
 */
void PREFIX readapt_rmsd_work_array(int natoms, FILE *fplog){
	int i,j;
	if(natoms>rmsd_workstruct.maxsize)
	{
		// 
		fprintf(fplog,"|- NEED TO REALLOCATE RMSD WORK ARRAY FROM %d to %d ELEMENTS\n",rmsd_workstruct.natoms,natoms);
		// allocate for the first time
		if(rmsd_workstruct.maxsize!=0){
			// reallocate a suitable array
			// free first
			free(rmsd_workstruct.align);
			free(rmsd_workstruct.lalign);
			free(rmsd_workstruct.displace);
			free(rmsd_workstruct.ldisplace);
			free_2dr_array_alloc(rmsd_workstruct.r0p,0); // the second is a fake for compatibility
			free_2dr_array_alloc(rmsd_workstruct.r1p,0);
			free_2dr_array_alloc(rmsd_workstruct.r0p_rotated,0);
			free_2dr_array_alloc(rmsd_workstruct.r1p_rotated,0);
			free_2dr_array_alloc(rmsd_workstruct.r0,0);
			free_2dr_array_alloc(rmsd_workstruct.r1,0);
			free_2dr_array_alloc(rmsd_workstruct.derr_dr0,0);
			free_2dr_array_alloc(rmsd_workstruct.derr_dr1,0);
			free_4dr_array_alloc(rmsd_workstruct.dd_dr0,3,3,3);
			free_4dr_array_alloc(rmsd_workstruct.dd_dr1,3,3,3);
			free_4dr_array_alloc(rmsd_workstruct.dm_r0_store,4,4,3);
			free_4dr_array_alloc(rmsd_workstruct.dm_r1_store,4,4,3);
			free_2dr_array_alloc(rmsd_workstruct.array_3_n,0);
			free_4dr_array_alloc(rmsd_workstruct.dd_dr_temp,3,3,3);

		}
		rmsd_workstruct.align=(real *)malloc(natoms*sizeof(real));
		rmsd_workstruct.displace=(real *)malloc(natoms*sizeof(real));
		rmsd_workstruct.lalign=(int *)malloc(natoms*sizeof(int));
		rmsd_workstruct.ldisplace=(int *)malloc(natoms*sizeof(int));
		rmsd_workstruct.r0p=float_2d_array_alloc(3,natoms);
		rmsd_workstruct.r1p=float_2d_array_alloc(3,natoms);
		rmsd_workstruct.r0p_rotated=float_2d_array_alloc(3,natoms);
		rmsd_workstruct.r1p_rotated=float_2d_array_alloc(3,natoms);
		rmsd_workstruct.r0=float_2d_array_alloc(3,natoms);
		rmsd_workstruct.r1=float_2d_array_alloc(3,natoms);
		rmsd_workstruct.derr_dr0=float_2d_array_alloc(3,natoms);
		rmsd_workstruct.derr_dr1=float_2d_array_alloc(3,natoms);
		rmsd_workstruct.dd_dr0=float_4d_array_alloc(3,3,3,natoms);
		rmsd_workstruct.dd_dr1=float_4d_array_alloc(3,3,3,natoms);
		
		rmsd_workstruct.dm_r1_store=float_4d_array_alloc(4,4,3,natoms);
		rmsd_workstruct.dm_r0_store=float_4d_array_alloc(4,4,3,natoms);
		rmsd_workstruct.dd_dr_temp=float_4d_array_alloc(3,3,3,natoms);
		rmsd_workstruct.array_3_n=float_2d_array_alloc(3,natoms);
		
		rmsd_workstruct.maxsize=natoms;
	}else{
		fprintf(fplog,"|- NO NEED TO REALLOCATE RMSD WORK ARRAY: REQUESTED SIZE OF  %d IS LOWER THAN ACTUAL %d\n ELEMENTS\n",natoms,rmsd_workstruct.maxsize);
	}
};
void PREFIX readapt_rmsd_work_array_secondder(int natoms, FILE *fplog){
	int i,j;
	// check the first derivative first
	readapt_rmsd_work_array_secondder(natoms,fplog);
	// now do the second der
	if(natoms>rmsd_workstruct.maxsize_secondder)
	{
		// 
		fprintf(fplog,"|- 2nd DER: NEED TO REALLOCATE RMSD WORK ARRAY FROM %d to %d\n ELEMENTS\n",rmsd_workstruct.natoms,natoms);
		// allocate for the first time
		if(rmsd_workstruct.maxsize_secondder!=0){
			// reallocate a suitable array
			// free first
			free_4dr_array_alloc(rmsd_workstruct.dderr_dr1_dr1,3,rmsd_workstruct.maxsize_secondder,3);
			free_4dr_array_alloc(rmsd_workstruct.dderr_dr1_dr0,3,rmsd_workstruct.maxsize_secondder,3);
			free_4dr_array_alloc(rmsd_workstruct.dderr_dr0_dr1,3,rmsd_workstruct.maxsize_secondder,3);
			free_4dr_array_alloc(rmsd_workstruct.dderr_dr0_dr0,3,rmsd_workstruct.maxsize_secondder,3);
		}
		rmsd_workstruct.dderr_dr1_dr1=float_4d_array_alloc(3,natoms,3,natoms);
		rmsd_workstruct.dderr_dr1_dr0=float_4d_array_alloc(3,natoms,3,natoms);
		rmsd_workstruct.dderr_dr0_dr1=float_4d_array_alloc(3,natoms,3,natoms);
		rmsd_workstruct.dderr_dr0_dr0=float_4d_array_alloc(3,natoms,3,natoms);
		rmsd_workstruct.maxsize_secondder=natoms;
	}else{
		fprintf(fplog,"|- 2nd DER: NO NEED TO REALLOCATE RMSD WORK ARRAY: REQUESTED SIZE OF  %d IS LOWER THAN ACTUAL %d\n ELEMENTS\n",natoms,rmsd_workstruct.natoms);
	}
	
	
};
void PREFIX allocate_rmsd_container(struct rmsd_container_t *rmsd_container,int i){
	rmsd_container->align=(real *)malloc(i*sizeof(real));
	rmsd_container->displace=(real *)malloc(i*sizeof(real));
};
void PREFIX copy_pdb_into_rmsd_container(struct pdb mypdb, struct rmsd_container_t *rmsd_container,FILE *fplog){
	int natoms, i,simple,nalign,ndisplace;
	real walign,wdislpace;
	natoms=mypdb.natoms;
	fprintf(fplog,"|- PDB HAS %d ATOMS\n",natoms);

	rmsd_container->align=(real *)malloc(natoms*sizeof(real));
	rmsd_container->displace=(real *)malloc(natoms*sizeof(real));
	rmsd_container->ldisplace=(int *)malloc(natoms*sizeof(int));
	rmsd_container->lalign=(int *)malloc(natoms*sizeof(int));

	rmsd_container->x=(real *)malloc(natoms*sizeof(real));
	rmsd_container->y=(real *)malloc(natoms*sizeof(real));
	rmsd_container->z=(real *)malloc(natoms*sizeof(real));
	rmsd_container->index=(int *)malloc(natoms*sizeof(int));

	
	walign=0.;wdislpace=0.;simple=1;nalign=0;ndisplace=0;
	rmsd_container->natoms=natoms;
	for(i=0;i<natoms;i++){
#if defined (PLUMED_GROMACS)	
		rmsd_container->x[i]=mypdb.x[i]/10.;
		rmsd_container->y[i]=mypdb.y[i]/10.;
		rmsd_container->z[i]=mypdb.z[i]/10.;
#else
		rmsd_container->x[i]=mypdb.x[i];
		rmsd_container->y[i]=mypdb.y[i];
		rmsd_container->z[i]=mypdb.z[i];
#endif
		rmsd_container->index[i]=mypdb.index[i]-1;
		rmsd_container->align[i]=mypdb.occ[i];
		rmsd_container->displace[i]=mypdb.beta[i];
		if(rmsd_container->align[i]!=rmsd_container->displace[i])simple=0;
		// assign various weight function 
		walign+=rmsd_container->align[i];
		wdislpace+=rmsd_container->displace[i];
		if(rmsd_container->align[i]!=0.){
			rmsd_container->lalign[nalign]=i;
			nalign++;
		}
		if(rmsd_container->displace[i]!=0.){
			rmsd_container->ldisplace[ndisplace]=i;
			ndisplace++;
		}
	};
	rmsd_container->ndisplace=ndisplace;
	rmsd_container->nalign=nalign;
	rmsd_container->walign=walign;
	rmsd_container->wdisplace=wdislpace;

	fprintf(fplog,"|- FOUND %d ATOMS FOR ALIGNMENT\n",nalign);
	fprintf(fplog,"|- FOUND %d ATOMS FOR DISPLACEMENT\n",ndisplace);
	if(simple==1){
		fprintf(fplog,"|- ALIGN AND DISPLACE ARE IDENTICAL: WILL USE FASTER SIMPLE ALIGNMENT\n");
	}else{
		fprintf(fplog,"|- ALIGN AND DISPLACE ARE NOT IDENTICAL: WILL USE SLOWER ALIGNMENT\n");	
	}
	
};
void PREFIX msd_restraint(int i_c, struct mtd_data_s *mtd_data){
	int i,j,iat;
	//fprintf(mtd_data->fplog,"|- ENTERED THE MSD ROUTINE: I_C IS %d\n",i_c);
	//
	// get the actual position and copy them into the work array
	// probably faster by copying only the pointers? but then you don't need the allocation
	//
	rmsd_workstruct.natoms=colvar.rmsd_container[i_c].natoms;
	for (i=0;i<colvar.natoms[i_c];i++){
		iat = colvar.cvatoms[i_c][i];
		// the running structure r0
		//fprintf(mtd_data->fplog,"|- ATOM %d ID %d : X %f\n",i,iat,mtd_data->pos[iat][0]);
		rmsd_workstruct.r0[0][i] = mtd_data->pos[iat][0];
		rmsd_workstruct.r0[1][i] = mtd_data->pos[iat][1];
		rmsd_workstruct.r0[2][i] = mtd_data->pos[iat][2];
		//fprintf(mtd_data->fplog,"|- ATOM %d : X %f Y %f Z %f\n",i,rmsd_workstruct.r0[0][i],rmsd_workstruct.r0[1][i],rmsd_workstruct.r0[2][i]);
		// the reference structure r1
		rmsd_workstruct.r1[0][i] = colvar.rmsd_container[i_c].x[i];
		rmsd_workstruct.r1[1][i] = colvar.rmsd_container[i_c].y[i];
		rmsd_workstruct.r1[2][i] = colvar.rmsd_container[i_c].z[i];
		// align and displace
		rmsd_workstruct.align[i]=colvar.rmsd_container[i_c].align[i];
		rmsd_workstruct.displace[i]=colvar.rmsd_container[i_c].displace[i];
		rmsd_workstruct.lalign[i]=colvar.rmsd_container[i_c].lalign[i];
		rmsd_workstruct.ldisplace[i]=colvar.rmsd_container[i_c].ldisplace[i];
	}
	rmsd_workstruct.ndisplace=colvar.rmsd_container[i_c].ndisplace;
	rmsd_workstruct.nalign=colvar.rmsd_container[i_c].nalign;
	rmsd_workstruct.simple=colvar.rmsd_container[i_c].simple;
	rmsd_workstruct.dorot=colvar.rmsd_container[i_c].dorot;
	rmsd_workstruct.dotrasl=colvar.rmsd_container[i_c].dotrasl;
	

	// do the real calculation
	// opts:
	//  -do the derivative respect to the ref struct
	//  -use the rotation matrix
	//  -use the geometric center correction
	//
	// all the output structs are in the rmsd workarray
	//
	if(colvar.rmsd_container[i_c].dummy){
		// exploit this vector that should not be used anyway 
		colvar.ss0[i_c]=0.;
		for(j=0;j<3;j++) {	
			for(i=0;i<colvar.natoms[i_c];i++) {
				iat= colvar.cvatoms[i_c][i];
//				colvar.myder[i_c][i][j] = 0.;
				colvar.myder[i_c][i][j] = mtd_data->pos[iat][j];

			}
		}
	}else{
		
		msd_calculation_dynamic(&rmsd_workstruct,colvar.rmsd_container[i_c].dref);
		
		// copy the result into derivative module
		colvar.ss0[i_c]=rmsd_workstruct.err;
		for(j=0;j<3;j++) {	
			for(i=0;i<colvar.natoms[i_c];i++) {
				colvar.myder[i_c][i][j] = rmsd_workstruct.derr_dr0[j][i];
			}
		}

		if(colvar.rmsd_container[i_c].dref){
			// calculate derivatives now	
			if(colvar.it%colvar.rmsd_container[i_c].dref_freq==0){
			fprintf(colvar.rmsd_container[i_c].fpdreffile,"NEWSET\n");
				for(i=0;i<colvar.natoms[i_c];i++) {
					for(j=0;j<3;j++) {	
						 fprintf(colvar.rmsd_container[i_c].fpdreffile,"%12.6f ",rmsd_workstruct.derr_dr1[j][i]);
					}
					fprintf(colvar.rmsd_container[i_c].fpdreffile,"\n");
				}
			}
		}
	}
// clean the structures (now dummy) 
	clean_rmsd_work_array(&rmsd_workstruct);
	
	//fprintf(mtd_data->fplog,"|- EXITED MSD  ROUTINE\n");
	//EXIT();
};
void  PREFIX msd_calculation_dynamic(struct rmsd_workstruct_s *work,int der_frameref_on){
	//fprintf(mtd_data.fplog,"|- ENTERED THE RMSD_DYNAMIC ROUTINE\n");
	
	//
	// finite difference tests
	//
	//rmsd_dynamic_findiff_interface(work);
	//EXIT();

	int do_rot=work->dorot;
	int do_center=work->dotrasl;
 
	if(do_rot) { 
		// fast version to get simple alignment
		if(work->simple){
			// do the alignment and 
			// -do the weighted centering
			// -do der respct reference frame r1
			// -avoid der respect to rotation matrix 
			// msd_core_dynamic_simple(struct rmsd_workstruct_s *work,int do_center, int do_frameref_der, int do_rotmat_der)
			msd_core_dynamic_simple(work,do_center,der_frameref_on,0);
		}else{
			// slow version: calculate the measure respect to a different set of atoms
			// -do the weighted centering
			// -calculate der respect to rotation matrix
			// -avoid der respect reference frame
			msd_core_dynamic_simple(work,do_center,der_frameref_on,1);
			// then use it into a different structure
			msd_core_dynamic_weighted(work,do_center,der_frameref_on);
		}		
	}else{
		// the simplest case: no rotation needed
		msd_core_dynamic_norot(work,do_center); 
	}
	
	//fprintf(mtd_data.fplog,"|- EXITED THE RMSD_DYNAMIC ROUTINE\n");
};	
void  PREFIX msd_core_dynamic_norot(struct rmsd_workstruct_s *work, int do_center){
	//fprintf(mtd_data.fplog,"|- ENTERED THE RMSD_DYNAMIC_NOROT ROUTINE\n");
	real xx[3], totalign, totdisplace, wdisplace;
	int natoms,i,j,k,l;

	natoms=work->natoms;
	totalign=0.;
	for(i=0;i<work->nalign;i++){
		k=work->lalign[i];
		totalign+=work->align[k];
	}


	wdisplace=0.;
	for(i=0;i<work->ndisplace;i++){
		k=work->ldisplace[i];
		wdisplace+=work->displace[k];
	}
	work->wdisplace=wdisplace;
	
	//
	// recenter r0
	//
	for(j=0;j<3;j++)xx[j]=0.;
	if(do_center) {// if you dont need to center no prob...
		for(j=0;j<3;j++){
			for(i=0;i<work->nalign;i++){
				k=work->lalign[i];
				xx[j]+=work->r0[j][k]*work->align[k];
				//fprintf(mtd_data.fplog,"|-R0 DIR %d AT  %d P %f \n",j,k,work->r0[j][k]);
			}
		}
		for(j=0;j<3;j++)xx[j]=xx[j]/(totalign);
		
	}
	for(j=0;j<3;j++)work->cmr0[j]=xx[j];
	
	// shift the position of all the atoms
	for(i=0;i<natoms;i++){
		for(j=0;j<3;j++)work->r0p[j][i]=work->r0[j][i]-xx[j];
		for(j=0;j<3;j++)work->r0p_rotated[j][i]=work->r0[j][i]-xx[j];
	}
	
	//
	// recenter the r1 structure now
	//	
	for(j=0;j<3;j++)xx[j]=0.;
	if(do_center) {// if you dont need to center no prob...
		for(j=0;j<3;j++){
			for(i=0;i<work->nalign;i++){
				k=work->lalign[i];
				xx[j]+=work->r1[j][k]*work->align[k];
				//fprintf(mtd_data.fplog,"|-R1 DIR %d AT  %d P %f \n",j,k,work->r1[j][k]);
			}
		}
		for(j=0;j<3;j++)xx[j]=xx[j]/(totalign);
		
	}
	for(j=0;j<3;j++)work->cmr1[j]=xx[j];
	
	// shift the position of all the atoms
	for(i=0;i<natoms;i++){
		for(j=0;j<3;j++)work->r1p[j][i]=work->r1[j][i]-xx[j];
		for(j=0;j<3;j++)work->r1p_rotated[j][i]=work->r1[j][i]-xx[j];
	}
	

	/*
	 * Identity in the ROTATION matrix
	 */
	
	work->d[0][0]=1.0 ; 
	work->d[0][1]=0.0 ; 
	work->d[0][2]=0.0 ; 
	work->d[1][0]=0.0 ; 
	work->d[1][1]=1.0 ; 
	work->d[1][2]=0.0 ; 
	work->d[2][0]=0.0 ; 
	work->d[2][1]=0.0 ; 
	work->d[2][2]=1.0 ; 

	work->err=0.;     
	totdisplace=0.;
	for(i=0;i<work->ndisplace;i++){
		k=work->ldisplace[i];
		totdisplace+=work->displace[k];
	}
	
	for(l=0;l<3;l++){
		for(i=0;i<work->ndisplace;i++){
			k=work->ldisplace[i];	
			//fprintf(mtd_data.fplog,"POS1 %d %d : R0P %f R1P %f \n",k,l,work->r0p[l][k],work->r1p[l][k]);
			work->err+=work->displace[k]*(work->r0p[l][k]-work->r1p[l][k])*(work->r0p[l][k]-work->r1p[l][k]);
		}
	}
	work->err/=totdisplace;
	
//	fprintf(mtd_data.fplog,"|- ERR %f\n",work->err);
	/*
	 * derivative 
	 *
	 */
	for(l=0;l<3;l++){
		for(i=0;i<work->natoms;i++){
			work->derr_dr1[l][i]=-(2.*work->displace[i]/(totdisplace))*(work->r0p[l][i]-work->r1p[l][i]); 
			if(do_center){
				for(j=0;j<work->natoms;j++){
					work->derr_dr1[l][i]+=(2.*work->align[i]*work->displace[j]/(totalign*totdisplace))*(work->r0p[l][j]-work->r1p[l][j]);
				}
			}
			//                	printf("DER %d DIR %d : %f\n",i,l,outpack->derr_dr1[l][i]);   
		}	
	}
	
	for(l=0;l<3;l++){
		for(i=0;i<work->natoms;i++){
			work->derr_dr0[l][i]=(2.*work->displace[i]/(totdisplace))*(work->r0p[l][i]-work->r1p[l][i]); 
			if(do_center){
				for(j=0;j<work->natoms;j++){
					work->derr_dr0[l][i]-=(2.*work->align[i]*work->displace[j]/(totalign*totdisplace))*(work->r0p[l][j]-work->r1p[l][j]);
				}
			}
			//                	printf("DER %d DIR %d : %f\n",i,l,outpack->derr_dr0[l][i]);   
		}	
	}
	
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			for(k=0;k<3;k++){
                for(l=0;l<work->natoms;l++){
					work->dd_dr1[i][j][k][l]=0.;
					work->dd_dr0[i][j][k][l]=0.;
                }
            }	
       	}	
	}
	
//	fprintf(mtd_data.fplog,"|- EXITING THE RMSD_DYNAMIC_NOROT ROUTINE\n");
};
//
// core routine that uses the dynamically allocated structure
//
// note that is R1 to be overlapped ON R0 !!!!!!!!!!
//
void  PREFIX msd_core_dynamic_simple(struct rmsd_workstruct_s *work,int do_center, int do_frameref_der, int do_rotmat_der){
	//fprintf(mtd_data.fplog,"|- ENTERED THE RMSD_DYNAMIC_SIMPLE ROUTINE\n");
	real xx[3], totalign, totdisplace, s, tmp1, fact1, fact2;
	//
	// needed workarray for quaternion: the array needed with natom dimensions are
	// allocated in the rmsd_workstruct
	//
	real m[4][4],rr1[4],rr0[4],q[4],lambda[4],dddq[4][4][4],gamma[3][3][3],rrsq;
	real dm_r1[4][4][3],dm_r0[4][4][3];
	real pi1[3][3],pi0[3][3];
	int natoms,i,j,k,l,ii,ll,jj,mm,n,nn,iii;
	
	//
	// recenter r0
	//
	natoms=work->natoms;
	totalign=0.;
	for(i=0;i<work->nalign;i++){
		k=work->lalign[i];
		totalign+=work->align[k];
	}
	for(j=0;j<3;j++)xx[j]=0.;
	if(do_center) {// if you dont need to center no prob...
		for(j=0;j<3;j++){
			for(i=0;i<work->nalign;i++){
				k=work->lalign[i];
				xx[j]+=work->r0[j][k]*work->align[k];
				//fprintf(mtd_data.fplog,"|-R0 DIR %d AT  %d P %f \n",j,k,work->r0[j][k]);
			}
		}
		for(j=0;j<3;j++)xx[j]=xx[j]/(totalign);
	}
	for(j=0;j<3;j++)work->cmr0[j]=xx[j];
	
	// shift the position of all the atoms
	for(i=0;i<natoms;i++){
		for(j=0;j<3;j++)work->r0p[j][i]=work->r0[j][i]-xx[j];
	}
	//fprintf(mtd_data.fplog,"|-NATOMS %d\n",natoms);
//	for(i=0;i<natoms;i++){
//		//for(j=0;j<3;j++)fprintf(mtd_data.fplog,"|-R0P DIR %d AT  %d P %f \n",i,j,work->r0p[j][i]);
//		fprintf(mtd_data.fplog,"ATOM    %3d    C ALA     1      %8.3f%8.3f%8.3f\n",i+1,work->r0p[0][i],work->r0p[1][i],work->r0p[2][i]);
//	}
//	fprintf(mtd_data.fplog,"END\n");
	
	
	//
	// recenter the r1 structure now
	//	
	for(j=0;j<3;j++)xx[j]=0.;
	if(do_center) {// if you dont need to center no prob...
		for(j=0;j<3;j++){
			for(i=0;i<work->nalign;i++){
				k=work->lalign[i];
				xx[j]+=work->r1[j][k]*work->align[k];
				//fprintf(mtd_data.fplog,"|-R1 DIR %d AT  %d P %f \n",j,k,work->r1[j][k]);
			}
		}
		for(j=0;j<3;j++)xx[j]=xx[j]/(totalign);
	}
	for(j=0;j<3;j++)work->cmr1[j]=xx[j];
	
	// shift the position of ALL the atoms
	for(i=0;i<natoms;i++){
		for(j=0;j<3;j++)work->r1p[j][i]=work->r1[j][i]-xx[j];
	}
//	for(i=0;i<natoms;i++){
//		for(j=0;j<3;j++)fprintf(mtd_data.fplog,"|-R1P DIR %d AT  %d P %f  R1  %f\n",i,j,work->r1p[j][i],work->r1[j][i]);
//		fprintf(mtd_data.fplog,"ATOM    %3d    C ALA     1      %8.3f%8.3f%8.3f\n",i+1,work->r0p[0][i],work->r0p[1][i],work->r0p[2][i]);
//
//	}
//	fprintf(mtd_data.fplog,"END\n");

	fflush(mtd_data.fplog);
	
	//
	// start with the real calculation
	//
	
	// CLEAN M MATRIX
	for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			m[i][j]=0.;  
		}
	}
	
	// ASSIGN MATRIX ELEMENTS USING ONLY THE ATOMS INVOLVED IN ALIGNMENT
	for(i=0;i<work->nalign;i++){
		
		k=work->lalign[i];
        tmp1=work->align[k];
		
		// adopt scaled coordinates
		
		rr1[0]=work->r1p[0][k]*tmp1;
        rr1[1]=work->r1p[1][k]*tmp1;
        rr1[2]=work->r1p[2][k]*tmp1;
        rr0[0]=work->r0p[0][k]*tmp1;
        rr0[1]=work->r0p[1][k]*tmp1;
        rr0[2]=work->r0p[2][k]*tmp1;
		
        rrsq=(pow(rr0[0],2)+pow(rr0[1],2)+pow(rr0[2],2)+pow(rr1[0],2)+pow(rr1[1],2)+pow(rr1[2],2));
		
        m[0][0] +=  rrsq+2.*(-rr0[0]*rr1[0]-rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[1][1] +=  rrsq+2.*(-rr0[0]*rr1[0]+rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[2][2] +=  rrsq+2.*(+rr0[0]*rr1[0]-rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[3][3] +=  rrsq+2.*(+rr0[0]*rr1[0]+rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[0][1] += 2.*(-rr0[1]*rr1[2]+rr0[2]*rr1[1]);
        m[0][2] += 2.*( rr0[0]*rr1[2]-rr0[2]*rr1[0]);
        m[0][3] += 2.*(-rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][2] -= 2.*( rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][3] -= 2.*( rr0[0]*rr1[2]+rr0[2]*rr1[0]);
        m[2][3] -= 2.*( rr0[1]*rr1[2]+rr0[2]*rr1[1]);
		
	};
	m[1][0] = m[0][1];
	m[2][0] = m[0][2];
	m[2][1] = m[1][2];
	m[3][0] = m[0][3];
	m[3][1] = m[1][3];
	m[3][2] = m[2][3];
	
	// check if you might expect something strange

	// DIAGONALIZE : minimize the distance respect to the scaled coordinates
	
	ql77_driver(m,lambda);
	s=1.0;
	if(m[0][0]<0.)s=-1.;//correct for negative values (?)
	q[0]=s*m[0][0];
	q[1]=s*m[1][0];
	q[2]=s*m[2][0];
	q[3]=s*m[3][0];
	work->err=lambda[0]/totalign;
	//fprintf(mtd_data.fplog,"|- ERR: %f \n",work->err);
	
	/*
	 * the ROTATION matrix
	 */
	
	work->d[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3]       ; 
	work->d[1][0]=2.0*(q[1]*q[2]-q[0]*q[3]);
	work->d[2][0]=2.0*(q[1]*q[3]+q[0]*q[2]);
	work->d[0][1]=2.0*(q[1]*q[2]+q[0]*q[3]);
	work->d[1][1]=q[0]*q[0]+q[2]*q[2]-q[1]*q[1]-q[3]*q[3];
	work->d[2][1]=2.0*(q[2]*q[3]-q[0]*q[1]);
	work->d[0][2]=2.0*(q[1]*q[3]-q[0]*q[2]);
	work->d[1][2]=2.0*(q[2]*q[3]+q[0]*q[1]);
	work->d[2][2]=q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];

	/* 
	 * first derivative in perturbation theory : derivative of the rotation matrix respect to the 
	 * quternion vectors
	 */

	dddq[0][0][0]= 2.0*q[0];
	dddq[1][0][0]=-2.0*q[3];
	dddq[2][0][0]= 2.0*q[2];
	dddq[0][1][0]= 2.0*q[3];
	dddq[1][1][0]= 2.0*q[0];
	dddq[2][1][0]=-2.0*q[1];
	dddq[0][2][0]=-2.0*q[2];
	dddq[1][2][0]= 2.0*q[1];
	dddq[2][2][0]= 2.0*q[0];
	
	dddq[0][0][1]= 2.0*q[1];
	dddq[1][0][1]= 2.0*q[2];
	dddq[2][0][1]= 2.0*q[3];
	dddq[0][1][1]= 2.0*q[2];
	dddq[1][1][1]=-2.0*q[1];
	dddq[2][1][1]=-2.0*q[0];
	dddq[0][2][1]= 2.0*q[3];
	dddq[1][2][1]= 2.0*q[0];
	dddq[2][2][1]=-2.0*q[1];
	
	dddq[0][0][2]=-2.0*q[2];
	dddq[1][0][2]= 2.0*q[1];
	dddq[2][0][2]= 2.0*q[0];
	dddq[0][1][2]= 2.0*q[1];
	dddq[1][1][2]= 2.0*q[2];
	dddq[2][1][2]= 2.0*q[3];
	dddq[0][2][2]=-2.0*q[0];
	dddq[1][2][2]= 2.0*q[3];
	dddq[2][2][2]=-2.0*q[2];
	
	dddq[0][0][3]=-2.0*q[3];
	dddq[1][0][3]=-2.0*q[0];
	dddq[2][0][3]= 2.0*q[1];
	dddq[0][1][3]= 2.0*q[0];
	dddq[1][1][3]=-2.0*q[3];
	dddq[2][1][3]= 2.0*q[2];
	dddq[0][2][3]= 2.0*q[1];
	dddq[1][2][3]= 2.0*q[2];
	dddq[2][2][3]= 2.0*q[3];
	/*
	 * Build gamma 3x3x3 matrix
	 */
        int dump_zero;
	dump_zero=0;
	for(i=0;i<3;i++){     //direction 
		for(j=0;j<3;j++){     //direction 
			for(k=0;k<3;k++){     //eigenvector number
				gamma[i][j][k]=0.0;
				for(l=0;l<4;l++){   //components of each eigenvector in pert. series
					if(lambda[0]==lambda[k+1] && dump_zero==0){
						fprintf(mtd_data.fplog,"|- FOUND DEGENERACY IN msd_core_dynamic_simple ROUTINE \n");
						fprintf(mtd_data.fplog,"|- I'm DYING....\n");
						fprintf(mtd_data.fplog,"|- LAMBDA_0 %f \n",lambda[0]);
						fprintf(mtd_data.fplog,"|- LAMBDA_K+1 %f \n",lambda[k+1]);
						fprintf(mtd_data.fplog,"|- COPYING STACK HERE \n");
						fprintf(mtd_data.fplog,"|- R0\n");
						for(ll=0;ll<natoms;ll++)fprintf(mtd_data.fplog,"|- %f %f %f \n",work->r0p[0][ll],work->r0p[1][ll],work->r0p[2][ll]);
						fprintf(mtd_data.fplog,"|- R1\n");
						for(ll=0;ll<natoms;ll++)fprintf(mtd_data.fplog,"|- %f %f %f \n",work->r1p[0][ll],work->r1p[1][ll],work->r1p[2][ll]);
						dump_zero=1;
						//plumed_error("FOUND DEGENERACY IN msd_core_dynamic_simple  ROUTINE");
					} 
					else{
						gamma[i][j][k]  +=  dddq[i][j][l]*m[l][k+1]/(lambda[0]-lambda[k+1]);
					}
				}
			}
			
		}	
	}
	// if this is set just put derivatives to zero in case of error=0. or degenerate eigenvalues 
if(dump_zero==1){
  fprintf(mtd_data.fplog,"|- WARNING: USING A FAKE VALUE FOR MSD \n"); 
  for(i=0;i<work->natoms;i++){
	for(j=0;j<3;j++){
                work->derr_dr1[j][i]=0.;
                work->derr_dr0[j][i]=0.;
		for (k=0;k<3;k++){
  			// ders of rotmatrix
			for(l=0;l<3;l++){	    
      				work->dd_dr1[j][k][l][i]=0.; 
      				work->dd_dr0[j][k][l][i]=0.; 
                        } 
                } 
        } 
  } 
  for(i=0;i<natoms;i++){
		work->r1p_rotated[0][i]=work->r1p[0][i];
		work->r1p_rotated[1][i]=work->r1p[0][i];
		work->r1p_rotated[2][i]=work->r1p[0][i];

  }
  for(i=0;i<natoms;i++){
		work->r0p_rotated[0][i]=work->r0p[0][i];
		work->r0p_rotated[1][i]=work->r0p[0][i];
		work->r0p_rotated[2][i]=work->r0p[0][i];
  }	
  for(i=0;i<3;i++){for(j=0;j<3;j++){work->dinv[i][j]=0.;work->d[i][j]=0.;}};
  for(i=0;i<3;i++){work->dinv[i][i]=1.;work->d[i][i]=1.;};

  return ;
}
	//fprintf(mtd_data.fplog,"|- STAGE 1 \n");

	  
	// clean up the rotation matrix for all the atoms
	// probably not needed
	for(j=0;j<3;j++){
		for (k=0;k<3;k++){
			for(l=0;l<3;l++){	  
				for(i=0;i<work->natoms;i++){
					work->dd_dr1[j][k][l][i]=0.;
					work->dd_dr0[j][k][l][i]=0.;
					work->dd_dr_temp[j][k][l][i]=0.;
				}
			}
		}
	}  
	// clean up the derivatives
	for(l=0;l<3;l++){	  
		for(i=0;i<work->natoms;i++){
			work->derr_dr1[l][i]=0.;
			work->derr_dr0[l][i]=0.;
			work->array_3_n[l][i]=0.;
		}
	}
	
		/* 
	 * Table of Derivative of the quaternion matrix respect to atom position: needed only if simple
	 * alignment is required and no correction respect to the rotation matrix is wanted
	 */
		
  for(iii=0;iii<work->nalign;iii++){
	  
	  i=work->lalign[iii];
	  tmp1=work->align[i];

	  // once again: derivative respect to scaled distance
	  
	  rr1[0]=2.*work->r1p[0][i]*tmp1;
	  rr1[1]=2.*work->r1p[1][i]*tmp1;
	  rr1[2]=2.*work->r1p[2][i]*tmp1;
	  rr0[0]=2.*work->r0p[0][i]*tmp1;
	  rr0[1]=2.*work->r0p[1][i]*tmp1;
	  rr0[2]=2.*work->r0p[2][i]*tmp1;
	  
	  
	  dm_r1 [0][0][0]=(rr1[0]-rr0[0]);
	  dm_r1 [0][0][1]=(rr1[1]-rr0[1]);
	  dm_r1 [0][0][2]=(rr1[2]-rr0[2]);
	  
	  dm_r1 [0][1][0]=0.;
	  dm_r1 [0][1][1]= rr0[2];
	  dm_r1 [0][1][2]=-rr0[1];
	  
	  dm_r1 [0][2][0]=-rr0[2];
	  dm_r1 [0][2][1]= 0.;
	  dm_r1 [0][2][2]= rr0[0];
	  
	  dm_r1 [0][3][0]= rr0[1];
	  dm_r1 [0][3][1]=-rr0[0];
	  dm_r1 [0][3][2]= 0.;
	  
	  dm_r1 [1][1][0]=(rr1[0]-rr0[0]);
	  dm_r1 [1][1][1]=(rr1[1]+rr0[1]);
	  dm_r1 [1][1][2]=(rr1[2]+rr0[2]);
	  
	  dm_r1 [1][2][0]=-rr0[1];
	  dm_r1 [1][2][1]=-rr0[0];
	  dm_r1 [1][2][2]= 0.;
	  
	  dm_r1 [1][3][0]=-rr0[2];
	  dm_r1 [1][3][1]= 0.;
	  dm_r1 [1][3][2]=-rr0[0];
	  
	  dm_r1 [2][2][0]=(rr1[0]+rr0[0]);
	  dm_r1 [2][2][1]=(rr1[1]-rr0[1]);
	  dm_r1 [2][2][2]=(rr1[2]+rr0[2]);
	  
	  dm_r1 [2][3][0]=0.;
	  dm_r1 [2][3][1]=-rr0[2];
	  dm_r1 [2][3][2]=-rr0[1];
	  
	  dm_r1 [3][3][0]=(rr1[0]+rr0[0]);
	  dm_r1 [3][3][1]=(rr1[1]+rr0[1]);
	  dm_r1 [3][3][2]=(rr1[2]-rr0[2]);
	  /*
	   derivative respec to to the other vector
	   */
	  dm_r0 [0][0][0]=-(rr1[0]-rr0[0]);
	  dm_r0 [0][0][1]=-(rr1[1]-rr0[1]);
	  dm_r0 [0][0][2]=-(rr1[2]-rr0[2]);
	  
	  dm_r0 [0][1][0]=0.       ;
	  dm_r0 [0][1][1]=-rr1[2];
	  dm_r0 [0][1][2]=rr1[1];
	  
	  dm_r0 [0][2][0]= rr1[2];      
	  dm_r0 [0][2][1]= 0.;
	  dm_r0 [0][2][2]=-rr1[0];
	  
	  dm_r0 [0][3][0]=-rr1[1] ;     
	  dm_r0 [0][3][1]= rr1[0];
	  dm_r0 [0][3][2]= 0.;
	  
	  dm_r0 [1][1][0]=-(rr1[0]-rr0[0]);
	  dm_r0 [1][1][1]=(rr1[1]+rr0[1]);
	  dm_r0 [1][1][2]=(rr1[2]+rr0[2]);
	  
	  dm_r0 [1][2][0]=-rr1[1];
	  dm_r0 [1][2][1]=-rr1[0];
	  dm_r0 [1][2][2]= 0.;
	  
	  dm_r0 [1][3][0]=-rr1[2];
	  dm_r0 [1][3][1]= 0.;
	  dm_r0 [1][3][2]=-rr1[0];
	  
	  dm_r0 [2][2][0]=(rr1[0]+rr0[0]);
	  dm_r0 [2][2][1]=-(rr1[1]-rr0[1]);
	  dm_r0 [2][2][2]=(rr1[2]+rr0[2]);
	  
	  dm_r0 [2][3][0]=0.;
	  dm_r0 [2][3][1]=-rr1[2];
	  dm_r0 [2][3][2]=-rr1[1];
	  
	  dm_r0 [3][3][0]=(rr1[0]+rr0[0]);
	  dm_r0 [3][3][1]=(rr1[1]+rr0[1]);
	  dm_r0 [3][3][2]=-(rr1[2]-rr0[2]);
	  /*
	   * write the diagonal
	   */ 
	  
	  for(j=0;j<3;j++){
		  
		  dm_r1[1][0][j]=dm_r1[0][1][j];
		  dm_r1[2][0][j]=dm_r1[0][2][j];
		  dm_r1[3][0][j]=dm_r1[0][3][j];
		  dm_r1[2][1][j]=dm_r1[1][2][j];
		  dm_r1[3][1][j]=dm_r1[1][3][j];
		  dm_r1[3][2][j]=dm_r1[2][3][j];
		  
		  dm_r0[1][0][j]=dm_r0[0][1][j];
		  dm_r0[2][0][j]=dm_r0[0][2][j];
		  dm_r0[3][0][j]=dm_r0[0][3][j];
		  dm_r0[2][1][j]=dm_r0[1][2][j];
		  dm_r0[3][1][j]=dm_r0[1][3][j];
		  dm_r0[3][2][j]=dm_r0[2][3][j];
		  
		  for(ll=0;ll<4;ll++){
			  for(mm=0;mm<4;mm++){
				  work->dm_r0_store[ll][mm][j][i]=dm_r0[ll][mm][j];
				  work->dm_r1_store[ll][mm][j][i]=dm_r1[ll][mm][j];
			  };
		  };
	  };
	  /*
	   * pi matrix : coefficents in per theory
	   */

	  for(j=0;j<3;j++){
		  pi1[0][j]=0.;
		  pi1[1][j]=0.;
		  pi1[2][j]=0.;
		  pi0[0][j]=0.;
		  pi0[1][j]=0.;
		  pi0[2][j]=0.;
		  work->derr_dr1[j][i]=0.;
		  work->derr_dr0[j][i]=0.;
		  
		  for(k=0;k<4;k++){
			  for(l=0;l<4;l++){
				  work->derr_dr1[j][i]=work->derr_dr1[j][i]+q[k]*q[l]*dm_r1[l][k][j];
				  work->derr_dr0[j][i]=work->derr_dr0[j][i]+q[k]*q[l]*dm_r0[l][k][j];
				  for(mm=0;mm<3;mm++){
					  pi0[mm][j]+=m[k][mm+1]*dm_r0[l][k][j]*q[l];
					  pi1[mm][j]+=m[k][mm+1]*dm_r1[l][k][j]*q[l];  
				  };
			  };
		  };
		  work->derr_dr1[j][i]=work->derr_dr1[j][i]/totalign;
		  work->derr_dr0[j][i]=work->derr_dr0[j][i]/totalign;
		  
	  };
	//fprintf(mtd_data.fplog,"|- STAGE 2c \n");

	  for(j=0;j<3;j++){
		  for (k=0;k<3;k++){
			  for(l=0;l<3;l++){	    
				  work->dd_dr1[j][k][l][i]=0.;
				  work->dd_dr0[j][k][l][i]=0.;
				  for(ii=0;ii<3;ii++){
					  work->dd_dr1[j][k][l][i]+=gamma[j][k][ii]*pi1[ii][l]; 
					  work->dd_dr0[j][k][l][i]+=gamma[j][k][ii]*pi0[ii][l]; 
				  }
			  }
		  }
	  }
  }
	
	
	// end of the calculation of the derivative of the rotation matrix
	
	/*
	 * Now correct for center of mass: only if needed 
	 *
	 */
		//
		// correction for r1 frame
		//
	if(do_frameref_der){
		for(l=0;l<3;l++){
			for(k=0;k<work->nalign;k++){
				i=work->lalign[k];
				work->array_3_n[l][i]=work->align[i]*work->derr_dr1[l][i];
				tmp1=work->align[i]/totalign; 
				if(do_center){		
					for(jj=0;jj<work->nalign;jj++){
						j=work->lalign[jj];
						work->array_3_n[l][i]-=tmp1*work->align[j]*work->derr_dr1[l][j];
					}
				}
				
			}
		}
		for(l=0;l<3;l++){
			for(k=0;k<work->nalign;k++){
				i=work->lalign[k];
				work->derr_dr1[l][i]=work->array_3_n[l][i];
			}
		}	
	}
	//
	// correction for r0 frame
	//
	for(l=0;l<3;l++){
		for(k=0;k<work->nalign;k++){
			i=work->lalign[k];
			work->array_3_n[l][i]=work->align[i]*work->derr_dr0[l][i];
			tmp1=work->align[i]/totalign; 
			if(do_center){		
				for(jj=0;jj<work->nalign;jj++){
					j=work->lalign[jj];
					work->array_3_n[l][i]-=tmp1*work->align[j]*work->derr_dr0[l][j];
				}
			}
		}
	}
	for(l=0;l<3;l++){
		for(k=0;k<work->nalign;k++){
			i=work->lalign[k];
			work->derr_dr0[l][i]=work->array_3_n[l][i];
		}
	}	
	//
	// correction for the rotation matrix: r1 frame
	//
	if(do_frameref_der && do_rotmat_der){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					for(ll=0;ll<work->nalign;ll++){
						l=work->lalign[ll];
						work->dd_dr_temp[i][j][k][l]=work->align[l]*work->dd_dr1[i][j][k][l];
						tmp1=work->align[l]/totalign; 
						if(do_center){		
							for(nn=0;nn<work->nalign;nn++){
								n=work->lalign[nn];
								work->dd_dr_temp[i][j][k][l]-=work->dd_dr1[i][j][k][n]*tmp1*work->align[n]; 
							}
						}
						
					}
				}	
			}	
		}
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					for(ll=0;ll<work->nalign;ll++){
						l=work->lalign[ll];
						work->dd_dr1[i][j][k][l]=work->dd_dr_temp[i][j][k][l];
					}
				}	
			}	
		}
	}
	//
	// correction for the rotation matrix: r0 frame
	//
	if(do_rotmat_der){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					for(ll=0;ll<work->nalign;ll++){
						l=work->lalign[ll];
						work->dd_dr_temp[i][j][k][l]=work->align[l]*work->dd_dr0[i][j][k][l];
						tmp1=work->align[l]/totalign; 
						if(do_center){		
							for(nn=0;nn<work->nalign;nn++){
								n=work->lalign[nn];
								work->dd_dr_temp[i][j][k][l]-=work->dd_dr0[i][j][k][n]*tmp1*work->align[n]; 
							}
						}
					}
				}	
			}	
		}
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					for(ll=0;ll<work->nalign;ll++){
						l=work->lalign[ll];
						work->dd_dr0[i][j][k][l]=work->dd_dr_temp[i][j][k][l];
					}
				}	
			}	
		}
	}
	// write the rotated frames
//	for(i=0;i<natoms;i++){
//		fprintf(mtd_data.fplog,"ATOM %6d  C   ALA     1    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,work->r0p[0][i],work->r0p[1][i],work->r0p[2][i]);
//	}
//	fprintf(mtd_data.fplog,"END\n");
//	
//	for(i=0;i<natoms;i++){
//		fprintf(mtd_data.fplog,"ATOM %6d  C   ALA     2    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,work->r1p[0][i],work->r1p[1][i],work->r1p[2][i]);
//	}
//	fprintf(mtd_data.fplog,"END\n");
//	
	for(i=0;i<natoms;i++){
		work->r1p_rotated[0][i]=work->d[0][0]*work->r1p[0][i]+
								work->d[0][1]*work->r1p[1][i]+
								work->d[0][2]*work->r1p[2][i];
		work->r1p_rotated[1][i]=work->d[1][0]*work->r1p[0][i]+
								work->d[1][1]*work->r1p[1][i]+
								work->d[1][2]*work->r1p[2][i];
		work->r1p_rotated[2][i]=work->d[2][0]*work->r1p[0][i]+
								work->d[2][1]*work->r1p[1][i]+
								work->d[2][2]*work->r1p[2][i];
//		fprintf(mtd_data.fplog,"ATOM %6d  C   ALA     2    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,work->r1p_rotated[0][i],work->r1p_rotated[1][i],work->r1p_rotated[2][i]);

	}
//	fprintf(mtd_data.fplog,"END\n");
//	EXIT();
	// invert the matrix
	real det;
	det= work->d[0][0]*(work->d[2][2]*work->d[1][1]-work->d[2][1]*work->d[1][2])
		-work->d[1][0]*(work->d[2][2]*work->d[0][1]-work->d[2][1]*work->d[0][2])
		+work->d[2][0]*(work->d[1][2]*work->d[0][1]-work->d[1][1]*work->d[0][2]);
	
						
	work->dinv[0][0]= (work->d[2][2]*work->d[1][1]-work->d[2][1]*work->d[1][2])/det;//a22a11-a21a12
	work->dinv[0][1]=-(work->d[2][2]*work->d[0][1]-work->d[2][1]*work->d[0][2])/det;//-(a22a01-a21a02)
	work->dinv[0][2]= (work->d[1][2]*work->d[0][1]-work->d[1][1]*work->d[0][2])/det;//a12a01-a11a02
	
	work->dinv[1][0]=-(work->d[2][2]*work->d[1][0]-work->d[2][0]*work->d[1][2])/det;  //-(a22a10-a20a12)
	work->dinv[1][1]= (work->d[2][2]*work->d[0][0]-work->d[2][0]*work->d[0][2])/det;	// 	a22a00-a20a02
	work->dinv[1][2]=-(work->d[1][2]*work->d[0][0]-work->d[1][0]*work->d[0][2])/det;	// 	-(a12a00-a10a02)
			
	work->dinv[2][0]= (work->d[2][1]*work->d[1][0]-work->d[2][0]*work->d[1][1])/det;  //  a21a10-a20a11
	work->dinv[2][1]=-(work->d[2][1]*work->d[0][0]-work->d[2][0]*work->d[0][1])/det;	//	-(a21a00-a20a01)		
	work->dinv[2][2]= (work->d[1][1]*work->d[0][0]-work->d[1][0]*work->d[0][1])/det;  //a11a00-a10a01
	
//	for(i=0;i<3;i++){
//		fprintf(mtd_data.fplog,"MAT ");
//		for(j=0;j<3;j++){
//			real v=0.;
//			for(ii=0;ii<3;ii++){
//				//for(jj=0;jj<3;jj++){
//				v+=work->d[i][ii]*work->dinv[ii][j];
//				//}
//			}
//			fprintf(mtd_data.fplog," %12.6f ",v);
//
//		}
//		fprintf(mtd_data.fplog,"\n");
//	}
	for(i=0;i<natoms;i++){
		work->r0p_rotated[0][i]=work->dinv[0][0]*work->r0p[0][i]+
								work->dinv[0][1]*work->r0p[1][i]+
								work->dinv[0][2]*work->r0p[2][i];
		work->r0p_rotated[1][i]=work->dinv[1][0]*work->r0p[0][i]+
								work->dinv[1][1]*work->r0p[1][i]+
								work->dinv[1][2]*work->r0p[2][i];
		work->r0p_rotated[2][i]=work->dinv[2][0]*work->r0p[0][i]+
								work->dinv[2][1]*work->r0p[1][i]+
								work->dinv[2][2]*work->r0p[2][i];
		// ATOM      1  C   ALA     2      -1.132   0.547  -0.390  0.50  0.50
		//fprintf(mtd_data.fplog,"ATOM %6d  C   ALA     1    %8.3f%8.3f%8.3f  1.00  1.00\n",i+1,work->r0p_rotated[0][i],work->r0p_rotated[1][i],work->r0p_rotated[2][i]);
				
	}
	//fprintf(mtd_data.fplog,"END\n");
	
	//fprintf(mtd_data.fplog,"|- EXITED THE RMSD_DYNAMIC_SIMPLE ROUTINE\n");
};
void PREFIX msd_core_dynamic_weighted(struct rmsd_workstruct_s *work,int do_center, int do_frameref_der){
	//fprintf(mtd_data.fplog,"|- ENTERING THE RMSD_DYNAMIC_WEIGHTED ROUTINE\n");
	
	real tmp0,tmp1,walign,wdisplace,ndisplace,const1;
	int  i,j,k,l,m,n,o,oo,mm,kk;
	
	// found various weight
	ndisplace=(real) work->ndisplace; 
	walign=0.;
	for(i=0;i<work->nalign;i++){
		k=work->lalign[i];
		walign+=work->align[k];
	}
	work->walign=walign;

	wdisplace=0.;
	for(i=0;i<work->ndisplace;i++){
		k=work->ldisplace[i];
		wdisplace+=work->displace[k];
	}
	work->wdisplace=wdisplace;
	
	tmp0=0.;
	for(kk=0;kk<work->ndisplace;kk++){
		k=work->ldisplace[kk];
		for(l=0;l<3;l++){
			
			tmp1=0.;
			// contribution from rotated reference frame //
			for(m=0;m<3;m++){
				tmp1-=work->d[l][m]*work->r1p[m][k];
			}
			
			// contribution from running centered frame //
			tmp1+= work->r0p[l][k];  
			//fprintf(mtd_data.fplog,"WEIGHTED %3d COMP %1d VAL %f\n",k,l,tmp1*sqrt(work->displace[k]/wdisplace));
			
			work->array_3_n[l][k]=tmp1; // store coefficents for derivative usage// 
			tmp0+=tmp1*tmp1*work->displace[k]; //squared distance added//
		}
	}  
	tmp0=tmp0/wdisplace;
	
    //fprintf(mtd_data.fplog,"|- ERRR NEW %f \n",tmp0);
	
	work->err=tmp0;

	/* DERIVATIVE CALCULATION:respect to running frame */
	for(k=0;k<work->natoms;k++){
		for(l=0;l<3;l++){
			
			tmp1 =2.0*work->array_3_n[l][k]*work->displace[k]/wdisplace ;
			
			const1=2.0*work->align[k]/(walign*wdisplace);
			
			if(const1!=0.){
				for(oo=0;oo<work->ndisplace;oo++){
					o=work->ldisplace[oo];
					tmp1 -=const1*work->array_3_n[l][o]*work->displace[o]; 
				} 
			}
			
				for(mm=0;mm<work->ndisplace;mm++){
					m=work->ldisplace[mm];
					const1=2.* work->displace[m]/wdisplace ;
					for(n=0;n<3;n++){
						tmp0=0.;
						for(o=0;o<3;o++){
							tmp0+=work->dd_dr0[n][o][l][k]*work->r1p[o][m];
						}
						tmp0*=-const1*work->array_3_n[n][m];
						
						tmp1+=tmp0;
					}
				}
			
			work->derr_dr0[l][k]=tmp1;
			
		}
	}

	if(do_frameref_der){
		for(k=0;k<work->natoms;k++){
			for(l=0;l<3;l++){
				/////////////////////////////////////				
				tmp1=0.;
				for(mm=0;mm<work->ndisplace;mm++){
					m=work->ldisplace[mm];
					const1=2.* work->displace[m]/wdisplace ;
					for(n=0;n<3;n++){
						tmp0=0.;
						for(o=0;o<3;o++){
							tmp0+=work->dd_dr1[n][o][l][k]*work->r1p[o][m];
						}
						tmp0*=-const1*work->array_3_n[n][m];
						tmp1+= tmp0;  
					}
					
				}
				
				tmp0=0.;
				for(o=0;o<3;o++){
							tmp0+=work->array_3_n[o][k]*work->d[o][l];
				}
				tmp1+=-tmp0*2.*work->displace[k]/wdisplace;
				
				
				tmp0=0.;
				for(mm=0;mm<work->ndisplace;mm++){
					m=work->ldisplace[mm];
					for(o=0;o<3;o++){
						tmp0+=work->array_3_n[o][m]*work->d[o][l]*work->displace[m];
					}
				}
				tmp1 += tmp0*2.*work->align[k]/(walign*wdisplace);
				
				work->derr_dr1[l][k]=tmp1;
				
				/////////////////////////////////////
			}
		}
	}

	//fprintf(mtd_data.fplog,"|- EXITED THE RMSD_DYNAMIC_WEIGHTED ROUTINE\n");
}
void PREFIX rmsd_dynamic_findiff_interface(struct rmsd_workstruct_s *work){
	fprintf(mtd_data.fplog,"|- ENTERED THE RMSD_DYNAMIC_FINDIFF ROUTINE\n");
	fprintf(mtd_data.fplog,"Entering rmsd finite difference test system\n");
	fprintf(mtd_data.fplog,"-------------------------------------------\n");
	fprintf(mtd_data.fplog,"TEST1: derivative of the value (derr_dr0/derr_dr1)\n");
	// test 1
	int i,j,k,l,m, weighted;
	real step=1.e-7,olderr,delta; 
	real derr_dr1[3][MAXATOMS_RMSD];
	real derr_dr0[3][MAXATOMS_RMSD];
	real dd_dr0[3][3][3][MAXATOMS_RMSD];
	real dd_dr1[3][3][3][MAXATOMS_RMSD];
	real oldd[3][3];
	weighted=1;
	// get initial value of the error and derivative of it 
	msd_core_dynamic_simple(work,1,1,1);
	if(weighted)msd_core_dynamic_weighted(work,1,1);
	//msd_core_dynamic_norot(work,1); 

	fprintf(mtd_data.fplog,"INITIAL ERROR VALUE: %f FOR %d ATOMS\n",work->err,work->natoms);
	olderr=work->err;
	// store the derivative
	for(j=0;j<3;j++){
		for(i=0;i<work->natoms;i++){
			derr_dr1[j][i]=work->derr_dr1[j][i];
			derr_dr0[j][i]=work->derr_dr0[j][i];
		}
	}
	for(l=0;l<3;l++){
		for(m=0;m<3;m++){
			oldd[l][m]=work->d[l][m];
			for(j=0;j<3;j++){
				for(i=0;i<work->natoms;i++){
					dd_dr1[l][m][j][i]=work->dd_dr1[l][m][j][i];
					dd_dr0[l][m][j][i]=work->dd_dr0[l][m][j][i];
				}
			}
		}
	}
	fprintf(mtd_data.fplog,"TESTING: derr_dr1 \n");
	for(j=0;j<3;j++){
		for(i=0;i<work->natoms;i++){
			// random displacement
			delta=(drand48()-0.5)*2*step;
			work->r1[j][i]+=delta; 
			msd_core_dynamic_simple(work,1,1,1);
			if(weighted)msd_core_dynamic_weighted(work,1,1);
			//msd_core_dynamic_norot(work,1); 

			work->r1[j][i]-=delta; 
			switch(j){
				case 0:
					fprintf(mtd_data.fplog,"TESTING: X  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr1[j][i],(work->err-olderr)/delta,derr_dr1[j][i]-(work->err-olderr)/delta);break;
				case 1:
					fprintf(mtd_data.fplog,"TESTING: Y  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr1[j][i],(work->err-olderr)/delta,derr_dr1[j][i]-(work->err-olderr)/delta);break;
				case 2:
					fprintf(mtd_data.fplog,"TESTING: Z  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr1[j][i],(work->err-olderr)/delta,derr_dr1[j][i]-(work->err-olderr)/delta);break;
					
			}    
		}
	}
	//EXIT();
	fprintf(mtd_data.fplog,"TESTING: derr_dr0 \n");
	for(j=0;j<3;j++){
		for(i=0;i<work->natoms;i++){
			// random displacement
			delta=(drand48()-0.5)*2*step;
			work->r0[j][i]+=delta; 
			msd_core_dynamic_simple(work,1,1,1);
			if(weighted)msd_core_dynamic_weighted(work,1,1);
			work->r0[j][i]-=delta; 
			switch(j){
				case 0:
					fprintf(mtd_data.fplog,"TESTING: X  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr0[j][i],(work->err-olderr)/delta,derr_dr0[j][i]-(work->err-olderr)/delta);break;
				case 1:
					fprintf(mtd_data.fplog,"TESTING: Y  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr0[j][i],(work->err-olderr)/delta,derr_dr0[j][i]-(work->err-olderr)/delta);break;
				case 2:
					fprintf(mtd_data.fplog,"TESTING: Z  %4d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",i,derr_dr0[j][i],(work->err-olderr)/delta,derr_dr0[j][i]-(work->err-olderr)/delta);break;
					
			}    
		}
	}
	fprintf(mtd_data.fplog,"TESTING: dd_dr0 \n");
	for(l=0;l<3;l++){
		for(m=0;m<3;m++){
			for(j=0;j<3;j++){
				for(i=0;i<work->natoms;i++){
					// random displacement
					delta=(drand48()-0.5)*2*step;
					work->r0[j][i]+=delta; 
					msd_core_dynamic_simple(work,1,1,1);
					if(weighted)msd_core_dynamic_weighted(work,1,1);
					work->r0[j][i]-=delta; 
					switch(j){
						case 0:
							fprintf(mtd_data.fplog,"TESTING: DD_DR0 [ %d ][ %d ]:  X %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr0[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr0[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
						case 1:
							fprintf(mtd_data.fplog,"TESTING: DD_DR0 [ %d ][ %d ]:  Y %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr0[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr0[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
						case 2:
							fprintf(mtd_data.fplog,"TESTING: DD_DR0 [ %d ][ %d ]:  Z %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr0[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr0[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
							
							
					}    
				}
			}
		}
	}
	fprintf(mtd_data.fplog,"TESTING: dd_dr1 \n");
	for(l=0;l<3;l++){
		for(m=0;m<3;m++){
			for(j=0;j<3;j++){
				for(i=0;i<work->natoms;i++){
					// random displacement
					delta=(drand48()-0.5)*2*step;
					work->r1[j][i]+=delta; 
					msd_core_dynamic_simple(work,1,1,1);
					if(weighted)msd_core_dynamic_weighted(work,1,1);
					work->r1[j][i]-=delta; 
					switch(j){
						case 0:
							fprintf(mtd_data.fplog,"TESTING: DD_DR1 [ %d ][ %d ]:  X %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr1[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr1[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
						case 1:
							fprintf(mtd_data.fplog,"TESTING: DD_DR1 [ %d ][ %d ]:  Y %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr1[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr1[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
						case 2:
							fprintf(mtd_data.fplog,"TESTING: DD_DR1 [ %d ][ %d ]:  Z %d ANAL %18.9f NUMER %18.9f DELTA %18.9f\n",l,m,i,dd_dr1[l][m][j][i],(work->d[l][m]- oldd[l][m])/delta,dd_dr1[l][m][j][i]-(work->d[l][m]- oldd[l][m])/delta);break;
							
					}    
				}
			}
		}
	}
	EXIT();
	
	fprintf(mtd_data.fplog,"|- EXITED THE RMSD_DYNAMIC_FINDIFF ROUTINE\n");
};
// after you use the struct make sure that you delete the structs that are not overwritten
void PREFIX clean_rmsd_work_array(struct rmsd_workstruct_s *work){
	int i,j,k,l;
	// always overwritten
	// r0p, r1p,cmr0,cmr1,align,displace,lalign,ldisplace,dd_dr1,dd_dr0,derr_dr1,derr_dr0
	//for (i=0;i<work->maxsize;i++){
	
	//}

};
// need already allocated pdb
void PREFIX copy_pdb (struct pdb *source, struct pdb *sink, FILE *fplog){
	int i,j,k;
	//deallocate_pdb(sink);
	//allocate_pdb(&sink,source->natoms);
	sink->natoms=source->natoms;
	for(i=0;i<source->natoms;i++){
		sink->x[i]=source->x[i];
		sink->y[i]=source->y[i];
		sink->z[i]=source->z[i];
		sink->beta[i]=source->beta[i];
		sink->occ[i]=source->occ[i];
		sink->index[i]=source->index[i];
		sink->resid[i]=source->resid[i];
		strcpy(sink->name[i],source->name[i]);
		strcpy(sink->resname[i],source->resname[i]);
	}

};
int  PREFIX hbd_collect_jacobian ( struct hybrid_frameset *running, real ** mathybrid ,real ** myinverse, FILE *fplog , int absval, int diag){
//AVOID
	int i,j,k,l,cvindex;
	int ii,iii,jj,jjj,n1,n2,nn,pos1,pos2;
	int random_matrix=0;
	//random_matrix=1;
	int unitary_matrix;
	struct hybrid_elem *first_elem;
	struct hybrid_elem *second_elem;		
	int jstart;
	
	// collect config has to be called before thhis function to store the derivatives
	// this only calculates the matrix and make the overlaps
	
	// for each cv calculate its vector of contributions (need a specific caller...)
	
	if (running->has_jac_alloc==0){
		fprintf(fplog,"|- hbd_collect_jacobian: ALLOCATING THE JACOBIAN\n");
		// now proceed with allocation
		running->has_jac_alloc=1;
		
		// now prepare the structure that tells which elements should be coupled
		// 
		// screen first how many elements you have to couple
		// 

		

		nn=0;
		for (i=0;i<running->hbd_nelem;i++){
			first_elem=running->hbd_elem[i];
			for (ii=0;ii< first_elem->ncv ;ii++){ // loop over all cvs of the first elem
				for(iii=0;iii<first_elem->natoms_per_cv[ii];iii++){
					n1=first_elem->list_per_cv[ii][iii];
					for (j=i;j<running->hbd_nelem;j++){
						second_elem=running->hbd_elem[j];
						if(j==i){
							jstart=ii;
						}else{
							jstart=0;
						}						
						for (jj=jstart;jj< second_elem->ncv ;jj++){ // loop over all the cvs of the second
							for(jjj=0;jjj<second_elem->natoms_per_cv[jj];jjj++){
								n2=second_elem->list_per_cv[jj][jjj];
								if(n1==n2){
									nn++;
								}
							}
						}
					}
				}
			}
		}
		
		// allocate the structure
		running->jac_index=(struct jac_index_s *)malloc(sizeof(struct jac_index_s)*nn);
		running->jac_index_length=nn;
		
		fprintf(fplog,"|- hbd_collect_jacobian: A TOTAL OF %d COMPONENTS OVERLAP\n",running->jac_index_length);
		
		// write the elements with the same loop
		
		nn=0;
		pos1=0;
		pos2=0;
		for (i=0;i<running->hbd_nelem;i++){
			first_elem=running->hbd_elem[i];
			for (ii=0;ii< first_elem->ncv ;ii++){ // loop over all cvs of the first elem
				pos2=pos1;
				for (j=i;j<running->hbd_nelem;j++){
					second_elem=running->hbd_elem[j];
					if(j==i){
						jstart=ii;
					}else{
						jstart=0;
					}
					for (jj=jstart;jj< second_elem->ncv ;jj++){ // loop over all the cvs of the second
						for(iii=0;iii<first_elem->natoms_per_cv[ii];iii++){
							n1=first_elem->list_per_cv[ii][iii];
							for(jjj=0;jjj<second_elem->natoms_per_cv[jj];jjj++){
								n2=second_elem->list_per_cv[jj][jjj];
								if(n1==n2){
									running->jac_index[nn].elem1=i;
									running->jac_index[nn].elem2=j;								
									running->jac_index[nn].cv1=ii;
									running->jac_index[nn].cv2=jj;
									running->jac_index[nn].natom1=iii;
									running->jac_index[nn].natom2=jjj;
									running->jac_index[nn].pos1=pos1;  // these last number represent the global indexing over the jacobian
									running->jac_index[nn].pos2=pos2;  // matrix
									//fprintf(fplog,"ELEMENT %d CV %d OVERLAPS WITH ELEMENT %d CV %d : ATOM %d  POS1 %d POS2 %d\n",i,ii,j,jj,n1,pos1,pos2);
									nn++;
								}
							}
							
						}
						pos2++;
					}
				}
				pos1++;
			}
		}
		
	}
	
	//EXIT();
	// calculate n*(n-1)/2 coupling element with cv
	// each term is calculated out of various indexes
	// in mathybrid
	
	// cleanup
	
	for (i=0;i<running->hbd_totcvs;i++){
		for (j=0;j<running->hbd_totcvs;j++){
			mathybrid[i][j]=0.;
		}
	}
	//fprintf(fplog,"MATHYBRID DIM IS %d ",running.hbd_totcvs);
	// looping over the overlaps
	//fprintf(fplog,"JAC LENGTH IS %d \n",running->jac_index_length);

	real v;
	for(nn=0;nn<running->jac_index_length;nn++){
		pos1=running->jac_index[nn].pos1;
		pos2=running->jac_index[nn].pos2;
		i=running->jac_index[nn].elem1;
		j=running->jac_index[nn].elem2;
		ii=running->jac_index[nn].cv1;
		jj=running->jac_index[nn].cv2;
		iii=running->jac_index[nn].natom1;
		jjj=running->jac_index[nn].natom2;
		v=0.;
		v+=running->hbd_elem[i]->cvder[ii][iii][0]*running->hbd_elem[j]->cvder[jj][jjj][0];
		v+=running->hbd_elem[i]->cvder[ii][iii][1]*running->hbd_elem[j]->cvder[jj][jjj][1];
		v+=running->hbd_elem[i]->cvder[ii][iii][2]*running->hbd_elem[j]->cvder[jj][jjj][2];
		//fprintf(fplog,"JAC TERMS ARE %d %d \n",pos1,pos2);
		if(pos1<=pos2){
			mathybrid[pos1][pos2]+=v	;
		}else{
			mathybrid[pos2][pos1]+=v	;		
		}
	}
	// unitary matrix is only for debug purposes
	if(logical.debug){
		pos1=0;
		int seed=0;
		for (i=0;i<running->hbd_nelem;i++){
			first_elem=running->hbd_elem[i];
			for (ii=0;ii< first_elem->ncv ;ii++){ // loop over all cvs of the first elem
				pos2=pos1;
				for (j=i;j<running->hbd_nelem;j++){
					second_elem=running->hbd_elem[j];
					if(j==i){
						jstart=ii;
					}else{
						jstart=0;
					}
					for (jj=jstart;jj< second_elem->ncv ;jj++){ // loop over all the cvs of the second
						mathybrid[pos1][pos2]=1.;
						pos2++;
					}
				}
				pos1++;
			}
		}
	}
	// this is to check the whole matrix with random assignment
	if(random_matrix){
		pos1=0;
		int seed=0;
		for (i=0;i<running->hbd_nelem;i++){
			first_elem=running->hbd_elem[i];
			for (ii=0;ii< first_elem->ncv ;ii++){ // loop over all cvs of the first elem
				pos2=pos1;
				for (j=i;j<running->hbd_nelem;j++){
					second_elem=running->hbd_elem[j];
					if(j==i){
						jstart=ii;
					}else{
						jstart=0;
					}
					for (jj=jstart;jj< second_elem->ncv ;jj++){ // loop over all the cvs of the second
						mathybrid[pos1][pos2]=rando(&seed)*100.0;
						pos2++;
					}
				}
				pos1++;
			}
		}
	}
	for (i=0;i<running->hbd_totcvs;i++){
		for (j=i;j<running->hbd_totcvs;j++){
			mathybrid[j][i]=mathybrid[i][j];
			//fprintf(fplog,"XXX %12.6f \n",mathybrid[i][j]);
		}
	}
	// just give back the dcv/dx matrix
		// now do the inverse (I know it is stupid! )
		for (i=0;i<running->hbd_totcvs;i++){
			for (j=0;j<running->hbd_totcvs;j++){
				if (absval==1){
					if(mathybrid[i][j]!=0.)mathybrid[i][j]=1./sqrt(pow(mathybrid[i][j],2));
				}else{
					if(mathybrid[i][j]!=0.)mathybrid[i][j]=1./mathybrid[i][j];
					
				}
			}
		}
	// keep only the diagonal
        if(diag==1){
          	for (i=0;i<running->hbd_totcvs;i++){
          		for (j=0;j<running->hbd_totcvs;j++){
          			if(i!=j)mathybrid[i][j]=0;
          		}
          	}
        }
	
	// calculate the jacobian terms
	if(random_matrix){
		fprintf (fplog,"|- RANDOM MATRIX!!!! \n");
		for (i=0;i<running->hbd_totcvs;i++){
			fprintf (fplog,"|-ELEM %4d ",i);
			for (j=0;j<running->hbd_totcvs;j++){
				fprintf(fplog,"%6.2f ",mathybrid[i][j]);
			}
			fprintf (fplog,"\n");
		}
	}
	
	// make its inverse for the weighting (not needed in the current scheme)	
	//matinverse(mathybrid,myinverse,running.hbd_totcvs,fplog);
	
	//fprintf(fplog,"|- hbd_collect_jacobian: CALCULATION DONE! \n");
	//EXIT();
	return 0;
//AVOID
};
// finite difference test for the general metrics 
void PREFIX check_hbd_vecmvec_ref(struct hybrid_frameset *r1,struct hybrid_frameset *r2,struct hybrid_frameset *r3,struct hybrid_frameset *r4,struct hybrid_frameset *r5,struct hybrid_frameset *r6,real **matrix, FILE  *fplog){
//AVOID
	int i,j,k;
	real oldist,newdist;
	real *dv1dcv,*dv2dcv,*dv3dcv,*dv4dcv,*dv5dcv,*dv6dcv;
	fprintf(fplog,"|- check_hbd_vecmvec_ref: ENTERING\n");
	
	// allocate the temporary derivative vectors
	fprintf(fplog,"|- check_hbd_vecmvec_ref: ALLOCATING TEMP VECTORS: %d DIM\n",r1->hbd_totcvs);
	
	snew(dv1dcv,r1->hbd_totcvs);
	snew(dv2dcv,r1->hbd_totcvs);
	snew(dv3dcv,r1->hbd_totcvs);
	snew(dv4dcv,r1->hbd_totcvs);
	snew(dv5dcv,r1->hbd_totcvs);
	snew(dv6dcv,r1->hbd_totcvs);
	
	oldist=hbd_vecmvec_ref(r1,r2,r3,r4,r5,r6,matrix,dv1dcv,dv2dcv,dv3dcv,dv4dcv,dv5dcv,dv6dcv,fplog);
	
	real oldval,fact,delta;
	int pos,seed;
	seed=0;
	fact=1.e-5;
	// start with r1
	
	pos=0;	
	fprintf(fplog,"|- check_hbd_vecmvec_ref: IN   (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ------------------------------\n");
	fprintf(fplog,"|- now checking                d[ (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ]/d_v1 --------------------\n",i);
	for(i=0;i<r1->hbd_nelem;i++){
		fprintf(fplog,"|- check_hbd_vecmvec_ref: ELEM %d CHECK STARTED------------------------------\n",i);
		struct hybrid_elem *myelem;
		myelem=r1->hbd_elem[i];
		for(j=0;j<myelem->ncv;j++){ 
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: CV %d CHECK: \n",j);
			oldval=myelem->ref_dist[j];
			delta=fact*(2.*rando(&seed)-1.);
			myelem->ref_dist[j]+=delta;
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: OLD %f NEW %f\n",oldval,myelem->ref_dist[j]); 
			newdist=hbd_vecmvec_ref(r1,r2,r3,r4,r5,r6,matrix,dv1dcv,dv2dcv,dv3dcv,dv4dcv,dv5dcv,dv6dcv,fplog);
			fprintf(fplog,"|- check_hbd_vecmvec_ref: FINDIFF  NUM %12.6f ANAL %12.6f DIFF %12.6f\n",(newdist-oldist)/delta,dv1dcv[pos],((newdist-oldist)/delta)-dv1dcv[pos]); 
			myelem->ref_dist[j]=oldval;
			pos++;
		}
	}
	
	pos=0;
	fprintf(fplog,"|- check_hbd_vecmvec_ref: IN   (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ------------------------------\n");
	fprintf(fplog,"|- now checking                d[ (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ]/d_v2 --------------------\n",i);
	for(i=0;i<r2->hbd_nelem;i++){
		fprintf(fplog,"|- check_hbd_vecmvec_ref: ELEM %d CHECK STARTED------------------------------\n",i);
		struct hybrid_elem *myelem;
		myelem=r2->hbd_elem[i];
		for(j=0;j<myelem->ncv;j++){ 
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: CV %d CHECK: \n",j);
			oldval=myelem->ref_dist[j];
			delta=fact*(2.*rando(&seed)-1.);
			myelem->ref_dist[j]+=delta;
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: OLD %f NEW %f\n",oldval,myelem->ref_dist[j]); 
			newdist=hbd_vecmvec_ref(r1,r2,r3,r4,r5,r6,matrix,dv1dcv,dv2dcv,dv3dcv,dv4dcv,dv5dcv,dv6dcv,fplog);
			fprintf(fplog,"|- check_hbd_vecmvec_ref: FINDIFF  NUM %12.6f ANAL %12.6f DIFF %12.6f\n",(newdist-oldist)/delta,dv2dcv[pos],((newdist-oldist)/delta)-dv2dcv[pos]); 
			myelem->ref_dist[j]=oldval;
			pos++;
		}
	}
	
	pos=0;
	fprintf(fplog,"|- check_hbd_vecmvec_ref: IN   (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ------------------------------\n");
	fprintf(fplog,"|- now checking                d[ (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ]/d_v3 --------------------\n",i);
	for(i=0;i<r3->hbd_nelem;i++){
		fprintf(fplog,"|- check_hbd_vecmvec_ref: ELEM %d CHECK STARTED------------------------------\n",i);
		struct hybrid_elem *myelem;
		myelem=r3->hbd_elem[i];
		for(j=0;j<myelem->ncv;j++){ 
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: CV %d CHECK: \n",j);
			oldval=myelem->ref_dist[j];
			delta=fact*(2.*rando(&seed)-1.);
			myelem->ref_dist[j]+=delta;
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: OLD %f NEW %f\n",oldval,myelem->ref_dist[j]); 
			newdist=hbd_vecmvec_ref(r1,r2,r3,r4,r5,r6,matrix,dv1dcv,dv2dcv,dv3dcv,dv4dcv,dv5dcv,dv6dcv,fplog);
			fprintf(fplog,"|- check_hbd_vecmvec_ref: FINDIFF  NUM %12.6f ANAL %12.6f DIFF %12.6f\n",(newdist-oldist)/delta,dv3dcv[pos],((newdist-oldist)/delta)-dv3dcv[pos]); 
			myelem->ref_dist[j]=oldval;
			pos++;
		}
	}
	
	pos=0;
	fprintf(fplog,"|- check_hbd_vecmvec_ref: IN   (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ------------------------------\n");
	fprintf(fplog,"|- now checking                d[ (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ]/d_v4 --------------------\n",i);
	for(i=0;i<r4->hbd_nelem;i++){
		fprintf(fplog,"|- check_hbd_vecmvec_ref: ELEM %d CHECK STARTED------------------------------\n",i);
		struct hybrid_elem *myelem;
		myelem=r4->hbd_elem[i];
		for(j=0;j<myelem->ncv;j++){ 
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: CV %d CHECK: \n",j);
			oldval=myelem->ref_dist[j];
			delta=fact*(2.*rando(&seed)-1.);
			myelem->ref_dist[j]+=delta;
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: OLD %f NEW %f\n",oldval,myelem->ref_dist[j]); 
			newdist=hbd_vecmvec_ref(r1,r2,r3,r4,r5,r6,matrix,dv1dcv,dv2dcv,dv3dcv,dv4dcv,dv5dcv,dv6dcv,fplog);
			fprintf(fplog,"|- check_hbd_vecmvec_ref: FINDIFF  NUM %12.6f ANAL %12.6f DIFF %12.6f\n",(newdist-oldist)/delta,dv4dcv[pos],((newdist-oldist)/delta)-dv4dcv[pos]); 
			myelem->ref_dist[j]=oldval;
			pos++;
		}
	}
	
	pos=0;
	fprintf(fplog,"|- check_hbd_vecmvec_ref: IN   (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ------------------------------\n");
	fprintf(fplog,"|- now checking                d[ (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ]/d_v5 --------------------\n",i);	
	for(i=0;i<r5->hbd_nelem;i++){
		fprintf(fplog,"|- check_hbd_vecmvec_ref: ELEM %d CHECK STARTED------------------------------\n",i);
		struct hybrid_elem *myelem;
		myelem=r5->hbd_elem[i];
		for(j=0;j<myelem->ncv;j++){ 
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: CV %d CHECK: \n",j);
			oldval=myelem->ref_dist[j];
			delta=fact*(2.*rando(&seed)-1.);
			myelem->ref_dist[j]+=delta;
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: OLD %f NEW %f\n",oldval,myelem->ref_dist[j]); 
			newdist=hbd_vecmvec_ref(r1,r2,r3,r4,r5,r6,matrix,dv1dcv,dv2dcv,dv3dcv,dv4dcv,dv5dcv,dv6dcv,fplog);
			fprintf(fplog,"|- check_hbd_vecmvec_ref: FINDIFF  NUM %12.6f ANAL %12.6f DIFF %12.6f\n",(newdist-oldist)/delta,dv5dcv[pos],((newdist-oldist)/delta)-dv5dcv[pos]); 
			myelem->ref_dist[j]=oldval;
			pos++;
		}
	}
	
	pos=0;
	fprintf(fplog,"|- check_hbd_vecmvec_ref: IN   (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ------------------------------\n");
	fprintf(fplog,"|- now checking                d[ (v1-v2)_ref_v3 M (v4-v5)_ref_v6 ]/d_v6 --------------------\n",i);	
	for(i=0;i<r6->hbd_nelem;i++){
		fprintf(fplog,"|- check_hbd_vecmvec_ref: ELEM %d CHECK STARTED------------------------------\n",i);
		struct hybrid_elem *myelem;
		myelem=r6->hbd_elem[i];
		for(j=0;j<myelem->ncv;j++){ 
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: CV %d CHECK: \n",j);
			oldval=myelem->ref_dist[j];
			delta=fact*(2.*rando(&seed)-1.);
			myelem->ref_dist[j]+=delta;
			//fprintf(fplog,"|- check_hbd_vecmvec_ref: OLD %f NEW %f\n",oldval,myelem->ref_dist[j]); 
			newdist=hbd_vecmvec_ref(r1,r2,r3,r4,r5,r6,matrix,dv1dcv,dv2dcv,dv3dcv,dv4dcv,dv5dcv,dv6dcv,fplog);
			fprintf(fplog,"|- check_hbd_vecmvec_ref: FINDIFF  NUM %12.6f ANAL %12.6f DIFF %12.6f\n",(newdist-oldist)/delta,dv6dcv[pos],((newdist-oldist)/delta)-dv6dcv[pos]); 
			myelem->ref_dist[j]=oldval;
			pos++;
		}
	}
	
	
	
	
	free(dv1dcv);
	free(dv2dcv);
	free(dv3dcv);
	free(dv4dcv);
	free(dv5dcv);
	free(dv6dcv);
	
	
	fprintf(fplog,"|- check_hbd_vecmvec_ref: EXITING\n");

	EXIT();
	
//AVOID
};
int PREFIX hbd_metrics_new ( struct hybrid_frameset *running , struct hybrid_frameset *reference , struct cmap_outpack *outpack, real **matrix, FILE  *fplog){
//AVOID
	int i,j,k;
	real dist;
	real *dv1dcv,*dv2dcv,*dv3dcv,*dv4dcv,*dv5dcv,*dv6dcv;
//	fprintf(fplog,"|- hbd_metrics_new: ENTERING\n");
	
	// allocate the temporary derivative vectors
	//fprintf(fplog,"|- hbd_metrics_new: ALLOCATING TEMP VECTORS: %d DIM\n",running->hbd_totcvs);

	snew(dv1dcv,running->hbd_totcvs);
	snew(dv2dcv,running->hbd_totcvs);
	snew(dv3dcv,running->hbd_totcvs);
	snew(dv4dcv,running->hbd_totcvs);
	snew(dv5dcv,running->hbd_totcvs);
	snew(dv6dcv,running->hbd_totcvs);
	
	//
	// make the difference between two vectors 
	// by keeping alignment respect to a third one:
	// this would be eventually important 
	// when doing a dot product or difference 
	// in bernd formulation or neb/sm
	//	
	
	// here call the vectorial oper   (v1-v2)_v3_align * matrix * (v4-v5)_v6_align
	
	outpack->err=hbd_vecmvec_ref(running,reference,reference,running,reference,reference, matrix,dv1dcv,dv2dcv,dv3dcv,dv4dcv,dv5dcv,dv6dcv,fplog);
	
	// dump error derivatives respect to all the atoms involved (in outpack)
	
	// in this specific implementation only running has derivatives respect to the atoms
	// namely needs to sum up dv1dcv and dv4dcv
	
	int pos;
	pos=0;	
	for(i=0;i<running->hbd_nelem;i++){		
		// sum up over all the elements
		struct hybrid_elem *myelem;
		myelem=running->hbd_elem[i];
		for(j=0;j<myelem->ncv;j++){ // sum up over all the cvs contained in each element
			dv1dcv[pos]+=dv4dcv[pos]; // put all the needed in dv1dcv
			dv2dcv[pos]+=dv3dcv[pos]+dv5dcv[pos]+dv6dcv[pos];
			pos++;
		}
	}
	// clean up
	for(j=0;j<running->hbd_totalatoms;j++){
		outpack->derr_dr0[0][j]=0.;
		outpack->derr_dr0[1][j]=0.;
		outpack->derr_dr0[2][j]=0.;
	}

	
	// now apply the chain rule
	
	int atom_offset,cv_offset;
	struct hybrid_elem *myelem;

	atom_offset=0;
	cv_offset=0;	
	pos=0;
	
	for(i=0;i<running->hbd_nelem;i++){
		
		//fprintf(fplog,"|- hbd_metrics_new: NOW ELEM %d :  ATOFFSET %d CVOFFSET %d\n", i,atom_offset,cv_offset);
		myelem=running->hbd_elem[i];
	
		for(j=0;j<myelem->ncv;j++){
			for(k=0;k<myelem->natoms_per_cv[j];k++){
				outpack->derr_dr0[0][atom_offset+k]=0.;
				outpack->derr_dr0[1][atom_offset+k]=0.;
				outpack->derr_dr0[2][atom_offset+k]=0.;
			}
		}
						
		// every element has its own atom list:
			
		
		for(j=0;j<myelem->ncv;j++){
			for(k=0;k<myelem->natoms_per_cv[j];k++){
				//
				//	 fprintf(fplog,"|- hbd_metrics_new: NELEM %d NOW CV %d ATOM %d TOUCHING %d\n",i,j,k,atom_offset+k);
				
				outpack->derr_dr0[0][atom_offset+k]+=myelem->cvder[j][k][0]*dv1dcv[cv_offset+j];	
				outpack->derr_dr0[1][atom_offset+k]+=myelem->cvder[j][k][1]*dv1dcv[cv_offset+j];	
				outpack->derr_dr0[2][atom_offset+k]+=myelem->cvder[j][k][2]*dv1dcv[cv_offset+j];	
			}
		}
			
		atom_offset+=myelem->nder; // sum up the number of cvs
		cv_offset+=myelem->ncv; // sum up the number of cvs

	}
			
	// define outpack.derr_dr1 (rvec totalatoms)
	// the derivatives respect to all the other cvs 
	// are a bit different and do not have ders respect to 
	// position of the atoms (they are only parameters)
	
	
	// deallocate the  temporary derivative vectors
	
	free(dv1dcv);
	free(dv2dcv);
	free(dv3dcv);
	free(dv4dcv);
	free(dv5dcv);
	free(dv6dcv);
	
//	fprintf(fplog,"|- hbd_metrics_new: EXITING\n");
//	EXIT();
	
//AVOID
}; 
//
// a simple utility that performs matrix inversion
//
int PREFIX matinverse (real **a, real **id, int n, FILE *fplog) 
{
	int   i, j,k;                    
	int   s;                       // index for elimination
	int   pindex;                  // Pivotindex
	int   error = 0;              // errorflag
	real f;                      // multiply fact
	const real Epsilon = 0.0000001;   // limit
	real Maximum;                // max
	int pivot = 1;
	
	real** asave;
	asave=float_2d_array_alloc(n,n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
		{
			asave[i][j] =a[i][j]; 
		}
	}
	
	// build id matrix
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
		{
			id[i][j] = 0.0;
			if (i == j)
				id[i][j] = 1.0;
		}
	}
	s = 0;
	do {
		Maximum = fabs(a[s][s]);
		if (pivot) // find the highest in the coulumn
		{
			pindex = s ; 
			for (i = s+1; i < n; i++)
				if (fabs(a[i][s]) > Maximum) {
					Maximum = fabs(a[i][s]) ;
					pindex = i;
				}
		}
		error = (Maximum < Epsilon);
		
		if (error) break;     
		
		if (pivot)
		{
			if (pindex != s)  
			{ double h; // swap two rows with the pivot
				for (j = s ; j < n; j++) {
					h = a[s][j];
					a[s][j] = a[pindex][j];
					a[pindex][j]= h;
				}
				for (j = 0 ; j < n; j++) {
					h = id[s][j];
					id[s][j] = id[pindex][j];
					id[pindex][j]= h;						
				}
			}
		}
		// normalize to one the chosen row
		f = a[s][s];
		for (j = s; j < n; j++){
			a[s][j]  = a[s][j] / f;
		}
		for (j = 0; j < n; j++){
			id[s][j] = id[s][j] / f;
		}
		
		//subtract the pivot row* a factor 
		for (i = 0; i < n; i++ ) {
			if (i != s) 
			{
				f = -a[i][s];                 
				for (j = s; j < n ; j++) {   
					a[i][j] += f*a[s][j];
				}
				for (j = 0; j < n ; j++) {  
					id[i][j] += f*id[s][j];
				}
			}
		}
		s++;
	} while ( s < n ) ;
	
	if (error) 
	{
		fprintf(fplog,"Inverse: Matrix is singular\n");
		return 0; 
	}
	
//	fprintf(fplog,"TEST \n");
//	for(i=0;i<n;i++){
//		for(j=0;j<n;j++){
//			f=0.;
//			for(k=0;k<n;k++)f+=asave[i][k]*id[k][j];
//			if(f<0.000000001){f=0;}
//			fprintf(fplog,"%6.2f ",f);
//
//		}
//		fprintf(fplog,"\n");$
//	}
//	fprintf(fplog,"INVMAT \n");
//
//	for(i=0;i<n;i++){
//		for(j=0;j<n;j++){
//			fprintf(fplog,"%6.2f ",id[i][j]);}
//		fprintf(fplog,"\n");
//	}		
	return 1;  
};

/*
 * a cloning routine for hybrid frameset and all the stuff in there
 */
void PREFIX clone_hybrid_frameset(struct hybrid_frameset **sink,struct hybrid_frameset *source,  int need_alloc_diff, FILE *fplog){
//AVOID
	

	int i;
	struct hybrid_elem *elem;
	//fprintf(fplog,"|- clone_hybrid_frameset: ENTERING\n");
	
	// first allocate the block of memory
	(*sink)=(struct hybrid_frameset *)malloc(sizeof(struct hybrid_frameset));
	
	//copy some basics integer stuff
	(*sink)->hbd_nelem=source->hbd_nelem;
	(*sink)->hbd_totalatoms=source->hbd_totalatoms;
	(*sink)->hbd_totcvs=source->hbd_totcvs;
	(*sink)->fixed=source->fixed;

	
	// allocate an array of pointers for all the elements
	
	(*sink)->hbd_elem=(struct hybrid_elem **)malloc((*sink)->hbd_nelem*sizeof(struct hybrid_elem * ));
	
	// loop over all the elements 
	
	for(i=0;i<(*sink)->hbd_nelem;i++){
		elem=source->hbd_elem[i];
	//	fprintf(fplog,"|- clone_hybrid_frameset: ELEM %d HAS PDB %d \n",i,source->hbd_elem[i]->contains_pdb);
		// now clone each element (allocation of the pointer is done inside the routine)
		clone_hybrid_elem(&((*sink)->hbd_elem[i]),elem,need_alloc_diff, fplog);

	}

	// allocate the derivative
	
	snew((*sink)->myder,(*sink)->hbd_totalatoms); 
	
	// allocate the backtable
	
	snew((*sink)->backtable,(*sink)->hbd_totalatoms);

//	fprintf(fplog,"|- clone_hybrid_frameset: TOTATOM %d\n",(*sink)->hbd_totalatoms);
	for(i=0;i<(*sink)->hbd_totalatoms;i++){
		//fprintf(fplog,"|- clone_hybrid_frameset: %3d   %3d \n",i,source->backtable[i]);
		(*sink)->backtable[i]=source->backtable[i];	
	}
	
	// avoid the structures for jacobian
	
	(*sink)->has_jac_alloc=source->has_jac_alloc;
	if((*sink)->has_jac_alloc==1){
		(*sink)->jac_index_length=source->jac_index_length;
		(*sink)->jac_index=(struct jac_index_s *)malloc(sizeof(struct jac_index_s)*(*sink)->jac_index_length);
		//alloc the jacobian stuff and copy it
		for(i=0;i<(*sink)->jac_index_length;i++){
			(*sink)->jac_index[i].elem1=source->jac_index[i].elem1;
			(*sink)->jac_index[i].elem2=source->jac_index[i].elem2;
			(*sink)->jac_index[i].cv1=source->jac_index[i].cv1;
			(*sink)->jac_index[i].cv2=source->jac_index[i].cv2;
			(*sink)->jac_index[i].natom1=source->jac_index[i].natom1;
			(*sink)->jac_index[i].natom2=source->jac_index[i].natom2;
			(*sink)->jac_index[i].pos1=source->jac_index[i].pos1;
			(*sink)->jac_index[i].pos2=source->jac_index[i].pos2;
		}
		//loop over the elements to copy each portion of the jacobian structures
		
	}
	
	
	//fprintf(fplog,"|- clone_hybrid_frameset: EXITING\n");
//AVOID
};
void PREFIX destroy_hybrid_frameset(struct hybrid_frameset *sink, FILE *fplog){
//AVOID
	int i;
	free(sink->myder);
	free(sink->backtable);
	for(i=0;i<sink->hbd_nelem;i++){
		destroy_hybrid_elem(sink->hbd_elem[i], fplog);
	}
	if(sink->has_jac_alloc){
		free(sink->jac_index);
	};
	free(sink->hbd_elem);
	free(sink);
//AVOID
};
void PREFIX clone_hybrid_elem(struct hybrid_elem **sink,struct hybrid_elem *source, int need_alloc_diff, FILE *fplog){
//AVOID
	int i,j,k ;
	//fprintf(fplog,"|- clone_hybrid_elem: ENTERING\n");
	//EXIT();

	(*sink)=(struct hybrid_elem *)malloc(sizeof(struct hybrid_elem));
	(*sink)->cvtype=source->cvtype;
	(*sink)->cvindex=source->cvindex;
	(*sink)->ncv=source->ncv;
	(*sink)->has_jac_alloc=source->has_jac_alloc;
        strcpy((*sink)->remarkfield,source->remarkfield);
	snew((*sink)->ref_dist,(*sink)->ncv);
	for(i=0;i<(*sink)->ncv;i++){
		(*sink)->ref_dist[i]=source->ref_dist[i];
	}
	(*sink)->nder=source->nder;
	// pdb 
	(*sink)->contains_pdb=source->contains_pdb;
	if((*sink)->contains_pdb==1){
		allocate_pdb(&((*sink)->pdb),(*sink)->nder);   
		copy_pdb(source->pdb,(*sink)->pdb, fplog);
	}
	// now clone all the vectorial stuff
	// backtable
	(*sink)->backtable=(int *)malloc((*sink)->nder*sizeof(int));
	for(i=0;i<(*sink)->nder;i++){
		(*sink)->backtable[i]=source->backtable[i];
	}
	// myder
	snew((*sink)->myder,colvar.natoms[(*sink)->cvindex]);

	if(need_alloc_diff==1){
		//fprintf(fplog,"|- clone_hybrid_elem: ALLOCATING A %d X %d DIFFERENCE VECTOR\n",(*sink)->ncv,(*sink)->ncv);
		(*sink)->diff_der1=float_2d_array_alloc((*sink)->ncv,(*sink)->ncv);
		(*sink)->diff_der2=float_2d_array_alloc((*sink)->ncv,(*sink)->ncv);
		(*sink)->diff_der3=float_2d_array_alloc((*sink)->ncv,(*sink)->ncv);
		(*sink)->has_diff_alloc=1;
	}else{
		(*sink)->has_diff_alloc=0;
	}
	(*sink)->has_jac_alloc=source->has_jac_alloc;
	if((*sink)->has_jac_alloc==1){
		(*sink)->natoms_per_cv=(int *)malloc((*sink)->ncv*sizeof(int));
		(*sink)->list_per_cv=(int **)malloc((*sink)->ncv*sizeof(int *));
		(*sink)->cvder=(rvec **)malloc((*sink)->ncv*sizeof(rvec *));
		for(i=0;i<(*sink)->ncv;i++){
			(*sink)->natoms_per_cv[i]=source->natoms_per_cv[i];
			(*sink)->cvder[i]=(rvec *)malloc((*sink)->natoms_per_cv[i]*sizeof(rvec));	
			(*sink)->list_per_cv[i]=(int *)malloc((*sink)->natoms_per_cv[i]*sizeof(int));
			for(j=0;j<(*sink)->natoms_per_cv[i];j++){
				(*sink)->list_per_cv[i][j]=source->list_per_cv[i][j];
			}
		}
	
	}
	
	//fprintf(fplog,"|- clone_hybrid_elem: EXITING\n");
	
//AVOID
};
void PREFIX destroy_hybrid_elem(struct hybrid_elem *sink, FILE *fplog){
//AVOID
	int i;
	free(sink->ref_dist);
	free(sink->backtable);
	free(sink->myder);
	if(sink->contains_pdb==1)deallocate_pdb(sink->pdb);
	if(sink->has_jac_alloc==1){
		for(i=0;i<sink->ncv;i++)
		{
			free(sink->list_per_cv[i]);
			free(sink->cvder[i]);
		}
		free(sink->cvder);
		free(sink->list_per_cv);
		free(sink->natoms_per_cv);
	}
	if(sink->has_diff_alloc==1){
		free_2dr_array_alloc(sink->diff_der1,sink->ncv);
		free_2dr_array_alloc(sink->diff_der2,sink->ncv);
		free_2dr_array_alloc(sink->diff_der3,sink->ncv);
	}
	free(sink);
//AVOID
};
/*
 * This simply makes the difference between a single cv element and returns the value and the derivative in the
 * difference vector allocated 
 */
void PREFIX simple_elem_diff(struct hybrid_elem *diff,struct hybrid_elem  *e1,struct hybrid_elem  *e2,struct hybrid_elem  *e3, FILE *fplog){
//AVOID
	int i,j,k;
	// consistency check: is this a truly simple element?
	if(e1->ncv!=1){
		fprintf(fplog,"|- simple_elem_diff: ERROR IN DIFFERENCE: THIS IS SUPPOSED TO BE A SIMPLE ELEMENT. NOW ABORTING\n");
	}
	// now do the calculation
	//fprintf(fplog,"|- simple_elem_diff: E1VAL %f E2VAL %f \n",e1->ref_dist[0],e2->ref_dist[0]);
	// just put this in the difference
	diff->ref_dist[0]=e1->ref_dist[0]-e2->ref_dist[0]; // e3 is not needed
	// the derivatives of the difference respect to the cv (not respect to
	diff->diff_der1[0][0]= 1.;
	diff->diff_der2[0][0]=-1.;
	diff->diff_der3[0][0]= 0.;
//AVOID
};
void PREFIX msd_elem_diff(struct hybrid_elem *diff,struct hybrid_elem  *e1,struct hybrid_elem  *e2,struct hybrid_elem  *e3, FILE *fplog, int dorandom,int print_rot_mat){
//AVOID
	// make the difference by using the alignment tools
	int i,j,i_c;
	i_c=e1->cvindex;
	
	// copy in the structure the general infos for alignment
	for (i=0;i<colvar.natoms[i_c];i++){
		// align and displace
		rmsd_workstruct.align[i]=colvar.rmsd_container[i_c].align[i];	
		rmsd_workstruct.displace[i]=colvar.rmsd_container[i_c].displace[i];
		rmsd_workstruct.lalign[i]=colvar.rmsd_container[i_c].lalign[i];
		rmsd_workstruct.ldisplace[i]=colvar.rmsd_container[i_c].ldisplace[i];
	}
	rmsd_workstruct.ndisplace=colvar.rmsd_container[i_c].ndisplace;
	rmsd_workstruct.wdisplace=colvar.rmsd_container[i_c].wdisplace;
	rmsd_workstruct.nalign=colvar.rmsd_container[i_c].nalign;
	rmsd_workstruct.simple=colvar.rmsd_container[i_c].simple;
	rmsd_workstruct.natoms=colvar.natoms[i_c];
	rmsd_workstruct.dotrasl=colvar.rmsd_container[i_c].dotrasl;
	rmsd_workstruct.dorot=colvar.rmsd_container[i_c].dorot;

	
	/*
	 *   e1 aligned to e3
	 *
	 */
	
	/* dirty workaround for identical structures: make a tiny random change */
	
	real delta;

	//int dorandom;
	//dorandom=0;
	//if(e1==e3 || e2==e3){dorandom=1;/*fprintf(fplog,"DOING RANDOM PERTURBATION \n"); */
	//}else{dorandom=0;}
	
	// copy the specific structure in general workstructure (should be already properly allocated at this point)
	j=0;
	for (i=0;i<e1->ncv;i=i+3){
		
		if(dorandom){
			delta=(drand48()-0.5)*2*1.e-4;
		}else{delta=0.;}
		
		rmsd_workstruct.r1[0][j] = e1->ref_dist[i+0];
		rmsd_workstruct.r1[1][j] = e1->ref_dist[i+1];
		rmsd_workstruct.r1[2][j] = e1->ref_dist[i+2];
		
		rmsd_workstruct.r0[0][j] = e3->ref_dist[i+0]+delta;
		rmsd_workstruct.r0[1][j] = e3->ref_dist[i+1]+delta;
		rmsd_workstruct.r0[2][j] = e3->ref_dist[i+2]+delta;
		
		// align and displace
		rmsd_workstruct.align[j]=colvar.rmsd_container[i_c].align[j];
		rmsd_workstruct.displace[j]=colvar.rmsd_container[i_c].displace[j];
		rmsd_workstruct.lalign[j]=colvar.rmsd_container[i_c].lalign[j];
		rmsd_workstruct.ldisplace[j]=colvar.rmsd_container[i_c].ldisplace[j];
		
		j++;
	}
	
	// perform the alignment
	
	msd_calculation_dynamic(&rmsd_workstruct,1);	

//	fprintf(fplog,"ERRR %f \n",rmsd_workstruct.err);
	
	real ****dd_de1,****dd_de2,****dd_de3;
	real **e1_p,**e2_p,**e3_p;
	int comp_i,comp_j,at_i,at_j;

	dd_de1=(rmsd_workstruct.dd_dr1);
	dd_de3=(rmsd_workstruct.dd_dr0);
	e1_p=rmsd_workstruct.r1p;
	e3_p=rmsd_workstruct.r0p;

	
	// now calculate the error
	// the component for the diff vector
	for(i=0;i<e1->ncv;i++){
		at_i=i/3;
		comp_i=i%3;
		diff->ref_dist[i]=rmsd_workstruct.r1p_rotated[comp_i][at_i];
		
	//	fprintf(fplog,"|- DIFF1 at %d comp %d : %f WW %f \n",at_i,comp_i, diff->ref_dist[i], rmsd_workstruct.wdisplace); 

		// the components along the vector
		for(j=0;j<e1->ncv;j++){	
			at_j=j/3;
			comp_j=j%3;
			
			//fprintf(fplog,"COMPONENT AT %d COMP %d \n",at_j,comp_j); 
			// first part of the alignment
						
			diff->diff_der1[i][j]=dd_de1[comp_i][0][comp_j][at_j]*e1_p[0][at_i]+
					      dd_de1[comp_i][1][comp_j][at_j]*e1_p[1][at_i]+
					      dd_de1[comp_i][2][comp_j][at_j]*e1_p[2][at_i] ;
			
			if(at_i==at_j)diff->diff_der1[i][j]+=rmsd_workstruct.d[comp_i][comp_j];
			
			diff->diff_der1[i][j]-=rmsd_workstruct.d[comp_i][comp_j]*rmsd_workstruct.align[at_j]/rmsd_workstruct.walign;
			
			diff->diff_der1[i][j]*=sqrt(rmsd_workstruct.displace[at_i]/rmsd_workstruct.wdisplace);
			
			diff->diff_der3[i][j]=dd_de3[comp_i][0][comp_j][at_j]*e1_p[0][at_i]+
					      dd_de3[comp_i][1][comp_j][at_j]*e1_p[1][at_i]+
					      dd_de3[comp_i][2][comp_j][at_j]*e1_p[2][at_i] ;
			
			
		}

	}


	
//	for (i=0;i<e1->ncv;i=i+3){
//		fprintf(fplog,"COMPONENT X E1: %f   E2: %f  E3: %f\n",e1->ref_dist[i+0],e2->ref_dist[i+0],e3->ref_dist[i+0]);
//		fprintf(fplog,"COMPONENT Y E1: %f   E2: %f  E3: %f\n",e1->ref_dist[i+1],e2->ref_dist[i+1],e3->ref_dist[i+1]);
//		fprintf(fplog,"COMPONENT Z E1: %f   E2: %f  E3: %f\n",e1->ref_dist[i+2],e2->ref_dist[i+2],e3->ref_dist[i+2]);
//    }

	
	/*
	 *   e2 aligned to e3
	 *
	 */
	
	j=0;
	for (i=0;i<e1->ncv;i=i+3){
		
		rmsd_workstruct.r1[0][j] = e2->ref_dist[i+0];
		rmsd_workstruct.r1[1][j] = e2->ref_dist[i+1];
		rmsd_workstruct.r1[2][j] = e2->ref_dist[i+2];
	
		j++;
	}
	
	msd_calculation_dynamic(&rmsd_workstruct,1);
	
	//fprintf(fplog,"ERRR %f \n",rmsd_workstruct.err);
	// if thisnumber is larger or equal to zero then open a file and print the 
	// rotation matrices with a label that depends on print_rot_mat 
	if(print_rot_mat>=0){
		FILE *fp;
	        char buf[1024];
		sprintf(buf,"rot_mat_%d.dat",print_rot_mat);
		fp=fopen(buf,"w");
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){	
				fprintf(fp,"%20.10f ", rmsd_workstruct.d[i][j]);
			}	
		}
		fclose(fp);
	}	
	

	dd_de2=(rmsd_workstruct.dd_dr1);
	dd_de3=(rmsd_workstruct.dd_dr0);
	e2_p=rmsd_workstruct.r1p;
	e3_p=rmsd_workstruct.r0p;
	
	for(i=0;i<e1->ncv;i++){
		at_i=i/3;
		comp_i=i%3;
		
		diff->ref_dist[i]-=rmsd_workstruct.r1p_rotated[comp_i][at_i];
	//	 fprintf(fplog,"|- DIFF2 at %d comp %d : %f %f %f \n",at_i,comp_i, rmsd_workstruct.r1p_rotated[comp_i][at_i],rmsd_workstruct.displace[at_i],rmsd_workstruct.wdisplace); 
		
		diff->ref_dist[i]*=sqrt(rmsd_workstruct.displace[at_i]/rmsd_workstruct.wdisplace);
	//	 fprintf(fplog,"|- DIFF at %d comp %d : %f \n",at_i,comp_i, diff->ref_dist[i]); 
		
		// the components along the vector
		for(j=0;j<e1->ncv;j++){	
			at_j=j/3;
			comp_j=j%3;
			diff->diff_der2[i][j]=-dd_de1[comp_i][0][comp_j][at_j]*e2_p[0][at_i]
					      -dd_de1[comp_i][1][comp_j][at_j]*e2_p[1][at_i]
					      -dd_de1[comp_i][2][comp_j][at_j]*e2_p[2][at_i] ;
			
			if(at_i==at_j)diff->diff_der2[i][j]-=rmsd_workstruct.d[comp_i][comp_j];
			
			diff->diff_der2[i][j]+=rmsd_workstruct.d[comp_i][comp_j]*rmsd_workstruct.align[at_j]/rmsd_workstruct.walign;

			diff->diff_der2[i][j]*=sqrt(rmsd_workstruct.displace[at_i]/rmsd_workstruct.wdisplace);
			
			diff->diff_der3[i][j] -= dd_de3[comp_i][0][comp_j][at_j]*e2_p[0][at_i]+
					         dd_de3[comp_i][1][comp_j][at_j]*e2_p[1][at_i]+
						 dd_de3[comp_i][2][comp_j][at_j]*e2_p[2][at_i] ;
			diff->diff_der3[i][j]*=sqrt(rmsd_workstruct.displace[at_i]/rmsd_workstruct.wdisplace);
			
		}
		
	}
	
//	for(i=0;i<e1->ncv;i++){
//		fprintf(fplog,"|- msd_elem_diff: ELEM %4d  VAL  %12.6f\n",i,diff->ref_dist[i]);
//	}
	
//AVOID
};
void PREFIX check_msd_elem_diff(struct hybrid_elem *diff,struct hybrid_elem *v1, struct hybrid_elem *v2, struct hybrid_elem *v3, FILE *fplog){
//AVOID
	int i,j,k;
	real oldist,newdist;
	struct hybrid_elem *oldiff;
	fprintf(fplog,"|- check_msd_elem_diff: ENTERING\n");
	
	// generate a new difference vector
	
	clone_hybrid_elem(&oldiff,diff,1,fplog);
	
	// collect the initial stuff and store in oldiff
	msd_elem_diff(oldiff,v1,v2,v3,fplog,0,-1);
	
	real oldval,fact,delta;
	int pos,seed;
	seed=0;
	fact=1.e-5;
	
	// start with v1:  calculate d (delta) / d ref_dist
	
	// choose the i-th component of difference vector
	for(i=0;i<v1->ncv;i++){
		fprintf(fplog,"|- check_msd_elem_diff: V1:   CV %d CHECK STARTED------------------------------\n",i);
		for(j=0;j<v1->ncv;j++){// choose the derivative respect a cv
			oldval=v1->ref_dist[j];
			delta=fact*(2.*rando(&seed)-1.);
			v1->ref_dist[j]+=delta;
			msd_elem_diff(diff,v1,v2,v3,fplog,0,-1);
			fprintf(fplog,"|- check_msd_elem_diff: FINDIFF  NUM %12.6f ANAL %12.6f DIFF %12.6f\n",
					(diff->ref_dist[i]-oldiff->ref_dist[i])/delta,diff->diff_der1[i][j],(diff->ref_dist[i]-oldiff->ref_dist[i])/delta-diff->diff_der1[i][j]);
			v1->ref_dist[j]=oldval;
		}
	}
	
	for(i=0;i<v2->ncv;i++){
		fprintf(fplog,"|- check_msd_elem_diff: V2:   CV %d CHECK STARTED------------------------------\n",i);
		for(j=0;j<v2->ncv;j++){// choose the derivative respect a cv
			oldval=v2->ref_dist[j];
			delta=fact*(2.*rando(&seed)-1.);
			v2->ref_dist[j]+=delta;
			msd_elem_diff(diff,v1,v2,v3,fplog,0,-1);
			fprintf(fplog,"|- check_msd_elem_diff: FINDIFF  NUM %12.6f ANAL %12.6f DIFF %12.6f\n",
					(diff->ref_dist[i]-oldiff->ref_dist[i])/delta,diff->diff_der2[i][j],(diff->ref_dist[i]-oldiff->ref_dist[i])/delta-diff->diff_der2[i][j]);
			v2->ref_dist[j]=oldval;
		}
	}
	
	for(i=0;i<v3->ncv;i++){
		fprintf(fplog,"|- check_msd_elem_diff: V3:   CV %d CHECK STARTED------------------------------\n",i);
		for(j=0;j<v3->ncv;j++){// choose the derivative respect a cv
			oldval=v3->ref_dist[j];
			delta=fact*(2.*rando(&seed)-1.);
			v3->ref_dist[j]+=delta;
			msd_elem_diff(diff,v1,v2,v3,fplog,0,-1);
			fprintf(fplog,"|- check_msd_elem_diff: FINDIFF  NUM %12.6f ANAL %12.6f DIFF %12.6f\n",
					(diff->ref_dist[i]-oldiff->ref_dist[i])/delta,diff->diff_der3[i][j],(diff->ref_dist[i]-oldiff->ref_dist[i])/delta-diff->diff_der3[i][j]);
			v3->ref_dist[j]=oldval;
		}
	}
	
	fprintf(fplog,"|- check_msd_elem_diff: EXITING\n");
	
	EXIT();
	
//AVOID
};
// needs a special debugging routine since the depency on overlapping atoms should be taken into account
// as independent
void PREFIX test_hbd_metrics_new(struct hybrid_frameset *r1 , struct hybrid_frameset *r2 , struct cmap_outpack *outpack, real **matrix, FILE  *fplog){
//AVOID
	fprintf(fplog,"|- test_hbd_metrics_new: ENTERING\n");
	real fact,delta,oldist,newdist;
	int i,j,k,atindex,cvindex,index;
	int seed;
	struct hybrid_elem *elem;
	
	// update the config
	for(i=0;i<r1->hbd_nelem;i++){
		elem=r1->hbd_elem[i];
		cvindex=elem->cvindex;
		switch(elem->cvtype){
			case 1:
				dist_restraint(cvindex,&mtd_data);
				break;
			case 3:
				coord_restraint(cvindex,&mtd_data);
				break;
			case 4:
				angle_restraint(cvindex,&mtd_data);
				break;
			case 5:
				torsion_restraint(cvindex,&mtd_data);
				break;
			case 32:
				position_restraint(cvindex,&mtd_data);
				break;
			case 53:
				msd_restraint(cvindex,&mtd_data);
				break;
		}
		// collect the derivatives and the config
	}
	hbd_collect_config(r1);
	// save old err value
	hbd_metrics_new(r1,r2,outpack,matrix,fplog);
	oldist=outpack->err;

	seed=0;
	fact=1.e-5;
	index=0;
	
	// loop over each atom of the element
	for(i=0;i<r1->hbd_nelem;i++){
		elem=r1->hbd_elem[i];
		cvindex=elem->cvindex;
		fprintf(fplog,"|-FINDIFF ELEM ---------------------------------------------------------------------------------------\n");
		//if(i==7)EXIT();
		for(j=0;j<elem->nder;j++){
			for(k=0;k<3;k++){
				atindex=elem->backtable[j];
				// change a bit
				delta=fact*(2.*rando(&seed)-1.);
				//delta=0.;
				// change position
				mtd_data.pos[atindex][k]+=delta;
				// update only the cv and ders that belongs to that element
				switch(elem->cvtype){
					case 1:
						dist_restraint(cvindex,&mtd_data);
						break;
					case 3:
						coord_restraint(cvindex,&mtd_data);
						break;
					case 4:
						angle_restraint(cvindex,&mtd_data);
						break;
					case 5:
						torsion_restraint(cvindex,&mtd_data);
						break;
					case 32:
						position_restraint(cvindex,&mtd_data);
						break;
					case 53:
						msd_restraint(cvindex,&mtd_data);
						break;
				}
				// collect the derivatives and the config
				hbd_collect_config(r1);
				// call metrics
				hbd_metrics_new(r1,r2,outpack,matrix,fplog);
				newdist=outpack->err;
			
				//fprintf(fplog,"|-FINDIFF OLDIST %f NEWDIST %f DELTA %e\n",oldist,newdist,delta);
				// do findiff
				fprintf(fplog,"|-FINDIFF ELEM %3d AT %5d COMP %1d ANAL %12.6f NUM %12.6f DELTA %12.6f DELTAPERC %12.6f\n",
						i,j,k,outpack->derr_dr0[k][index],(newdist-oldist)/delta,(outpack->derr_dr0[k][index]-(newdist-oldist)/delta),
						100*	(outpack->derr_dr0[k][index]-(newdist-oldist)/delta)/outpack->derr_dr0[k][index]);
				
				// paste back the atom position
				mtd_data.pos[atindex][k]-=delta;
				// recall the right function (so to reset the value of the cv to previous value)
				switch(elem->cvtype){
					case 1:
						dist_restraint(cvindex,&mtd_data);
						break;
					case 3:
						coord_restraint(cvindex,&mtd_data);
						break;
					case 4:
						angle_restraint(cvindex,&mtd_data);
						break;
					case 5:
						torsion_restraint(cvindex,&mtd_data);
						break;
					case 32:
						position_restraint(cvindex,&mtd_data);
						break;
					case 53:
						msd_restraint(cvindex,&mtd_data);
						break;
				}	
				hbd_collect_config(r1);

			}
			index++;
		}
	}
	fprintf(fplog,"|- test_hbd_metrics_new: EXITING\n");

	EXIT();
//AVOID
};

// 
//  (v1-v2)_ref_v3 ' M (v1-v2)_ref_v3   is calculated 
//  you should already provide an allocated distance frameset (use "clone" clone_hybrid_frameset(&d1,v1,1,fplog) with specific flag )
//
real PREFIX hbd_distanddiff_ref(struct hybrid_frameset *v1,  struct hybrid_frameset *v2, struct hybrid_frameset *v3, 
							real **mat , real *dv1dcv,  real *dv2dcv,  real *dv3dcv,struct hybrid_frameset *d1,  
							FILE *fplog){
//AVOID
	real dist;
	int i,j,k,l,m;
	struct hybrid_elem *el1,*el2,*el3,*diff1;
	
	//fprintf(fplog,"|- hbd_distanddiff_ref : ENTERING\n");
	
	//
	// (v1-v2)_ref_v3
	// and derivatives
	//
	
	// loop over all the elements of the frameset and make the differences:
	// in case of RMSD you make the alignment and store the data
	// (matrix derivatives and so on)
	
	
	for(i=0;i< v1->hbd_nelem; i++){
		
		el1=v1->hbd_elem[i];
		el2=v2->hbd_elem[i];
		el3=v3->hbd_elem[i];
		diff1=d1->hbd_elem[i];
		
		int incr=0 ;// just an increment position for the derivative vector (being updated)
		switch(el1->cvtype){
				// case simple cv: just do the difference (no pbc or other stuff)
				// note: if angles are sin and cos component then this is good anyway
			case 1:	
				simple_elem_diff(diff1,el1,el2,el3,fplog);
				break;
			case 3:	
				simple_elem_diff(diff1,el1,el2,el3,fplog);
				break;
			case 4:	
				simple_elem_diff(diff1,el1,el2,el3,fplog);
				break;
			case 5:	
				simple_elem_diff(diff1,el1,el2,el3,fplog);
				break;
			case 32:	
				simple_elem_diff(diff1,el1,el2,el3,fplog);
				break;
			case 53:	
				// case msd
				//check_msd_elem_diff(diff1,el1,el2,el3,fplog);
				msd_elem_diff(diff1,el1,el2,el3,fplog,1,i);
				break;
			default:
				fprintf(fplog,"|- THIS DIFFERENCE IS NOT IMPLEMENTED: CVTYPE %d ",el1->cvtype);
				EXIT();
				break;
		}
		
	}
	
	// now calculate the  projection
	int ii,jj,kk,pos1,pos2,oldpos1;
	struct hybrid_elem *first_elem;
	struct hybrid_elem *second_elem;
	real val;
	pos1=0;
	pos2=0;
	int jstart;

	dist=0.;
	int startval;
	startval=0;
	// full loop: can be optimized ofcourse...
	for (i=0;i<v1->hbd_nelem;i++){
		first_elem=v1->hbd_elem[i];
		for (ii=0;ii< first_elem->ncv ;ii++){ // loop over all cvs of the first elem
			// reset the value
			dv1dcv[pos1]=0.;
			dv2dcv[pos1]=0.;
			dv3dcv[pos1]=0.;
			// do the inner loop for the summation
			pos2=0;
			
			for (j=0;j<v1->hbd_nelem;j++){
				second_elem=v1->hbd_elem[j];
				for (jj=0;jj< second_elem->ncv ;jj++){ 
					
					dist+=d1->hbd_elem[i]->ref_dist[ii]*mat[pos1][pos2]*d1->hbd_elem[j]->ref_dist[jj];
					
					//fprintf(fplog,"|-COMPONENT %d %d %d %d %f\n",i,ii,j,jj,d1->hbd_elem[i]->ref_dist[ii]*mat[pos1][pos2]*d1->hbd_elem[j]->ref_dist[jj]);
					//fprintf(fplog,"|- hbd_vecmvec:-------------------------------\n");
					for(kk=0;kk<first_elem->ncv;kk++){
						
						dv1dcv[pos1]+= 2*d1->hbd_elem[i]->diff_der1[kk][ii] * mat[startval+kk][pos2] * d1->hbd_elem[j]->ref_dist[jj] ;
						dv2dcv[pos1]+= 2*d1->hbd_elem[i]->diff_der2[kk][ii] * mat[startval+kk][pos2] * d1->hbd_elem[j]->ref_dist[jj] ;
						dv3dcv[pos1]+= 2*d1->hbd_elem[i]->diff_der3[kk][ii] * mat[startval+kk][pos2] * d1->hbd_elem[j]->ref_dist[jj] ;
						
					}
					pos2++;
				}
			}
			// end of the inner loop for the summation
			pos1++;
		}
		startval+=first_elem->ncv;
	}
	
	
	//
	//fprintf(fplog,"|- hbd_distanddiff_ref : EXITING\n");

	return dist;
//AVOID
};
void PREFIX dump_frameset_formatted  ( struct sz_data *pmy_sz ){
//AVOID
	fprintf(mtd_data.fplog,"|- dump_frameset: ENTERING\n");
	int i,j;
	struct hybrid_elem *elem;
	struct hybrid_frameset *frame;

	// introduce a progressive number for the evolution
	
	// the writable file is already open and in    pmy_sz->fp_evol
	// dump first a separator line
	
	fprintf(pmy_sz->fp_evol,"NEWSET\n");
	
	//check that it does not exist
	
	for(i=0;i<pmy_sz->number;i++){
		frame=pmy_sz->hbd_frameset[i];
		fprintf(pmy_sz->fp_evol,"NEWFRAME\n");
		for(j=0;j<frame->hbd_nelem;j++){
			elem=frame->hbd_elem[j];
			//fprintf(mtd_data.fplog,"|- dump_frameset: %s\n",elem->remarkfield);
			switch (elem->cvtype){  
				case 1: 					
					hbd_dump_simple(elem,pmy_sz->fp_evol); break ;
				case 3:
					hbd_dump_simple(elem,pmy_sz->fp_evol); break ;
				case 4: 
					hbd_dump_simple(elem,pmy_sz->fp_evol); break ;
				case 5:
					hbd_dump_simple(elem,pmy_sz->fp_evol); break ;
				case 32:
					hbd_dump_simple(elem,pmy_sz->fp_evol); break ;
				case 53:
					hbd_dump_msd(elem,pmy_sz->fp_evol); break ;
				default:fprintf(mtd_data.fplog,"|- HBD_COLLECT_CONFIG: TYPE %d \n", elem->cvtype);  plumed_error("|- NOT IMPLEMENTED: now dying..."); break;
			}
			
		}
		// is fixed? dump the label then!
		if(frame->fixed){
			//fprintf(mtd_data.fplog,"|- dump_frameset: FIXED\n");
			fprintf(pmy_sz->fp_evol,"FIXED\n");
		}
	}
	fflush(pmy_sz->fp_evol);
	
	fprintf(mtd_data.fplog,"|- dump_frameset: EXITING\n");
//AVOID
};
void PREFIX hbd_dump_simple( struct hybrid_elem *elem , FILE *fp ){
//AVOID
	// simply dump the remark line, the value and END
	fprintf(fp,"%s",elem->remarkfield);
	fprintf(fp,"%12.6f\n",elem->ref_dist[0]);
	fprintf(fp,"END\n");
//AVOID
}
void PREFIX hbd_dump_msd( struct hybrid_elem *elem , FILE *fp ){
//AVOID
	int i;
	struct pdb *mypdb;
	real occ, beta;
	mypdb=elem->pdb;
	int j=0;
	for(i=0;i<mypdb->natoms;i++){
		mypdb->x[i]=elem->ref_dist[j] ;j++;
		mypdb->y[i]=elem->ref_dist[j] ;j++;
		mypdb->z[i]=elem->ref_dist[j] ;j++;
	}
	// dump the whole pdb 
	fprintf(fp,"%s",elem->remarkfield);
	for(i=0;i<mypdb->natoms;i++){
		occ=colvar.rmsd_container[elem->cvindex].align[i];
		beta=colvar.rmsd_container[elem->cvindex].displace[i];
		fprintf(fp,"ATOM  %5d %4s %3s %5d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",mypdb->index[i],mypdb->name[i],mypdb->resname[i],
			   mypdb->resid[i],mypdb->x[i],mypdb->y[i],mypdb->z[i],occ,beta);
		
	}
	fprintf(fp,"END\n");
//AVOID
};
void PREFIX dump_frameset_unformatted  ( struct sz_data *pmy_sz ){
//AVOID
	
	fprintf(mtd_data.fplog,"|- dump_frameset: ENTERING\n");
	int i,j,jj;
	struct hybrid_elem *elem;
	struct hybrid_frameset *frame;
	
	// introduce a progressive number for the evolution
	
	// the writable file is already open and in    pmy_sz->fp_evol
	// dump first a separator line
	
			
	fprintf(pmy_sz->fp_evol,"%10d \n",mtd_data.istep);
	for(i=0;i<pmy_sz->number;i++){
		frame=pmy_sz->hbd_frameset[i];
		for(j=0;j<frame->hbd_nelem;j++){
			elem=frame->hbd_elem[j];
			fprintf(pmy_sz->fp_evol," X ");
			for(jj=0;jj<elem->ncv;jj++){
					fprintf(pmy_sz->fp_evol,"%12.6f ",elem->ref_dist[jj] );
			}
		}
		// is fixed? dump the label then!
		if(frame->fixed){
			//fprintf(mtd_data.fplog,"|- dump_frameset: FIXED\n");
			fprintf(pmy_sz->fp_evol," F ");
		}
		fprintf(pmy_sz->fp_evol,"\n"); 
		
	}
	fflush(pmy_sz->fp_evol);
	
	fprintf(mtd_data.fplog,"|- dump_frameset: EXITING\n");

	
//AVOID
}
void PREFIX read_couplingmatrix  ( char **word, t_plumed_input *input,FILE *fplog ){
//AVOID
	int i,j,iw, help;
	struct sz_data *my_sz;
	// set up a minimal and fake path
	my_sz=&(couplingmatrix.mysz);
	my_sz->number = 1;
	couplingmatrix.is_on=1;
	help=0;
	// this reads the dumping of the couplingmatrix
	fprintf(fplog,"|-COUPLINGMATRIX BEING READED:\n");
	// read the cvs involved
	iw=seek_word(word,"CV");
	if(iw>=0) {
        int c,ii,*t1l,t1s=0, *t2l,t2s,jj;
        ii=0;
        while (1) { 
			//      word[iw+1] is a word ???   
			ii++;
			c=*word[iw+ii];
			if (isalpha(c)!=0 ){
				//printf("This is a word %s\n",word[iw+ii]); 
				break; 
			} else {
				//          printf("This is a number %s\n",word[iw+ii]); 
				if(t1s==0){ t1l=(int *)malloc(sizeof(int));
					sscanf(word[iw+ii],"%d", &t1l[0]);
				} else{ 
					t2s=t1s;t2l=(int *)malloc(t2s*sizeof(int));
					for(jj=0;jj<t1s;jj++){
						t2l[jj]=t1l[jj];
					}
					free(t1l);
					t1l=(int *)malloc((t1s+1)*sizeof(int)); 
					for(jj=0;jj<t1s;jj++){
						t1l[jj]=t2l[jj];
					}
					free(t2l);
					sscanf(word[iw+ii],"%d", &t1l[t1s]);
				} 
				t1s++;
				//for(jj=0;jj<t1s;jj++){printf("LL %d ",t1l[jj]);}; printf(" NN %d\n",t1s) ; 
			} 
        } 
		my_sz->nhybrid=t1s; // number of hybrid functions
        my_sz->lhybrid=(int *)malloc(my_sz->nhybrid*sizeof(int));		// list of the ordinal of hybrid function to be called
        my_sz->lcvhybrid=(int *)malloc(my_sz->nhybrid*sizeof(int));		// list of the ordinal of the cv to be used
        for(ii=0;ii<my_sz->nhybrid;ii++){my_sz->lhybrid[ii]=t1l[ii];}   // list of the type of cv 
		
        fprintf(fplog,"|--NUMBER OF HYBRID CVS= %d WHICH ARE ",my_sz->nhybrid);
		for(jj=0;jj<my_sz->nhybrid;jj++){
			fprintf(fplog,"CV %d ",my_sz->lhybrid[jj]);my_sz->lhybrid[jj]--;
        };printf("\n"); 
		//
        // set the type of the representation for each cv  
		//
        for(jj=0;jj<my_sz->nhybrid;jj++){
            int kk=my_sz->lhybrid[jj];
            fprintf(fplog,"|- CVTYPE %d\n",colvar.type_s[kk]);my_sz->lcvhybrid[jj]=colvar.type_s[kk];
            // parse the input and check if it is implemented  
		}; 
	}else{
		fprintf(fplog,"|- MISSING DEFINITION OF THE CVS INVOLVED\n");
		help=1;
		EXIT();
	}	
	
	// read the frequency of dumping
	iw=seek_word(word,"DUMPFREQ");
	if(iw>=0) {
		sscanf(word[iw+1],"%d", &couplingmatrix.dumpfreq);
		fprintf(fplog,"|- DUMPFREQ IS:%d \n",couplingmatrix.dumpfreq);
	}else{
		couplingmatrix.dumpfreq=1;
		fprintf(fplog,"|- DUMPFREQ MISSING: ASSUMING %d\n",couplingmatrix.dumpfreq);
	}
	
	iw=seek_word(word,"FILENAME");
	if(iw>=0) {
		//sscanf(word[iw+1],"%s", &couplingmatrix.filename);
		strcpy(couplingmatrix.filename,word[iw+1]);
		fprintf(fplog,"|- FILENAME IS: %s \n",couplingmatrix.filename);
	}else{
		strcpy(couplingmatrix.filename,"matrix.dat");
		fprintf(fplog,"|- FILENAME MISSING: ASSUMING %s\n",couplingmatrix.filename);
		
	}
	// read some other options (frequency/displacement)
	
	// create a fake hybrid path to allocate the elements correctly
	// this is a bit reduntant with read_sz_hybrid 
	
	my_sz->hbd_frameset=(struct hybrid_frameset **)malloc(sizeof(struct hybrid_frameset *)); //only one frameset
	struct hybrid_frameset *myframe;
	// allocate the frame itself
	my_sz->hbd_frameset[0]=(struct hybrid_frameset *)malloc(sizeof(struct hybrid_frameset));
	myframe=my_sz->hbd_frameset[0];
	
	// this array has the size of the real of calls
	myframe->hbd_nelem=my_sz->nhybrid;
	myframe->hbd_elem=(struct hybrid_elem **)malloc((myframe->hbd_nelem)*sizeof(struct hybrid_elem * ));
	
	// reset the number of atoms
	myframe->hbd_totalatoms=0;   
	
	// nelem=number of real CVs in the hybrid representation
	for (j=0;j<myframe->hbd_nelem;j++){
		myframe->hbd_elem[j]=(struct hybrid_elem *)malloc(sizeof(struct hybrid_elem));
		// type of the cv in the list
		myframe->hbd_elem[j]->cvtype=my_sz->lcvhybrid[j]; 
		// progressive index in the list of the cvs
		myframe->hbd_elem[j]->cvindex=my_sz->lhybrid[j];
		
		// put always on for all the variables
		logical.always[my_sz->lhybrid[j]]=1;
		
		// switch off pdb
		myframe->hbd_elem[j]->contains_pdb=0;
		
		FILE *fp; //fake file pointer
		fp=NULL;
		switch(my_sz->lcvhybrid[j]) { // select the cv
			case 1:  
				myframe->hbd_totalatoms+=hbd_read_simple(fp,fplog,"DIST",myframe->hbd_elem[j],0); 
				break;
			case 3:  
				myframe->hbd_totalatoms+=hbd_read_simple(fp,fplog,"COORD",myframe->hbd_elem[j],0); 
				break;
				
			case 4:  
				myframe->hbd_totalatoms+=hbd_read_simple(fp,fplog,"ANGLE",myframe->hbd_elem[j],0); 
				break;
				
			case 5:  
				myframe->hbd_totalatoms+=hbd_read_simple(fp,fplog,"TORSION",myframe->hbd_elem[j],0); 
				break;

			case 32:  
				myframe->hbd_totalatoms+=hbd_read_simple(fp,fplog,"POSITION",myframe->hbd_elem[j],0); 
				break;
				
			case 53:  
				myframe->hbd_totalatoms+=hbd_read_msd(fp,fplog,myframe->hbd_elem[j],0);
				break;
				
			default: plumed_error("READER NOT IMPLEMENTED"); break;
		}
	} 
	
	// now allocate the structures according to the number of atoms
	
	myframe->hbd_totcvs=0;
	for(i=0;i<my_sz->nhybrid ;i++){
		myframe->hbd_totcvs+=myframe->hbd_elem[i]->ncv;
	}
	
	fprintf(fplog,"|- TOTAL NUMBER OF CVS IN THE FRAMESET IS %d \n",myframe->hbd_totcvs);
	
	fprintf(fplog,"|- ALLOCATING METRIC MATRIX \n");
	i=myframe->hbd_totcvs;
	my_sz->mathybrid=float_2d_array_alloc(i,i);
	my_sz->myinverse=float_2d_array_alloc(i,i);
	
	
	fprintf(fplog,"|- WRITING GLOBAL BACKTABLE\n");
	struct  hybrid_elem *ee;
	int l,k;
	for (i=0;i< my_sz->number ;i++){
		l=0;
		myframe->backtable=(int *)malloc(myframe->hbd_totalatoms*sizeof(int));
		for (j=0;j<myframe->hbd_nelem;j++){
			ee=(myframe->hbd_elem[j]) ; 
			for (k=0;k<ee->nder;k++){
				myframe->backtable[l]= ee->backtable[k];
				l++;
			}
		}	
		snew(myframe->myder,myframe->hbd_totalatoms); 
	}	
	
	
	myframe->jac_index_length=0;
	myframe->has_jac_alloc=0;

	
	for (j=0;j<myframe->hbd_nelem;j++){
		
		// take one frame as reference for allocating the structure: this 
		// impose that all the structs must have the same number of elements
		struct hybrid_elem *elem=myframe->hbd_elem[j];
				
		// single elements allocate immediately the tools for the jacobian
		elem->has_jac_alloc=1;

		elem->cvder=(rvec **)malloc(elem->ncv*sizeof(rvec *));// 3n cvs for each structure
		if(elem->cvder==NULL){fprintf(fplog,"ALLOCATION FAILED"); EXIT();}
		
		for(i=0;i<elem->ncv;i++)snew(elem->cvder[i],elem->natoms_per_cv[i]); // each structure may depend on the other atoms
		
	}
	
	// print some help
	if(help){
		fprintf(fplog,"|- COUPLINGMATRIX");
		fprintf(fplog,"|- allows to calculate the projection matrix for some gradients of interest");
		fprintf(fplog,"|- SYNTAX: ");
		fprintf(fplog,"|- COUPLINGMATRIX CV 1 5 7 4 DUMPFREQ 10 FILENAME mymatrix.dat ");
		fprintf(fplog,"|- ");
		fprintf(fplog,"|- ");
		fprintf(fplog,"|- ");
		fprintf(fplog,"|- ");
		EXIT();
	}
	// open the file
	couplingmatrix.fp=fopen(couplingmatrix.filename,"w");
	// done!
	//EXIT();
//AVOID
}
/*
 * this performs the calculation of the couplingmatrix
 */
void PREFIX calc_couplingmatrix  ( int istep ){
//AVOID
	
	if(istep%couplingmatrix.dumpfreq==0){
		
		fprintf(mtd_data.fplog,"|- CALC_COUPLINGMATRIX: entering\n");		
		struct hybrid_frameset *hbd_pmy;
		struct sz_data *pmy_sz;
		struct hybrid_frameset *frame;
		pmy_sz=&couplingmatrix.mysz;
		
		
		frame=pmy_sz->hbd_frameset[0];
		
		hbd_collect_config(frame);  // this collect cv values and needed derivatives
		
		hbd_collect_jacobian(frame,pmy_sz->mathybrid,pmy_sz->myinverse,mtd_data.fplog,0,0);
		
		// dump the value
		fprintf(couplingmatrix.fp,"NEWSET\n");
		int i,j;
		
		for (i=0;i<frame->hbd_totcvs;i++){
			for (j=0;j<frame->hbd_totcvs;j++){
				fprintf(couplingmatrix.fp,"%12.6e ",pmy_sz->mathybrid[i][j]);
			}
			fprintf(couplingmatrix.fp,"\n");
		}
		fprintf(mtd_data.fplog,"|- CALC_COUPLINGMATRIX: exiting\n");

	}

	
    //EXIT();
//AVOID
}
// this routine calculates the maragliano vanden-eijnden projector
void PREFIX calc_projector_test  ( struct sz_data *pmy_sz ){
//AVOID
	fprintf(mtd_data.fplog,"|-calc_projector_test: CALCULATING THE PROJECTOR\n");
	int i,j,k;
	struct hybrid_frameset *frame1,*frame2,*difference;
	struct hybrid_elem *el1,*el2,*el3,*diff;
	struct cmap_outpack outpack;
	real dist;
	real *distance;
	
	distance=(real *)malloc(pmy_sz->number*sizeof(real));
	
	clone_hybrid_frameset(&difference,pmy_sz->hbd_frameset[0],1,mtd_data.fplog);

	// loop over all the frames and use one single matrix as template
	for(i=0;i<pmy_sz->number-1;i++){
		fprintf(mtd_data.fplog,"|-calc_projector_test: calculating projector between frame %d and %d\n",i,i+1);
		frame1=pmy_sz->hbd_frameset[i];
		frame2=pmy_sz->hbd_frameset[i+1];
		hbd_metrics_new(frame1,frame2,&outpack,pmy_sz->mathybrid,mtd_data.fplog);
		distance[i]=outpack.err;
		fprintf(mtd_data.fplog,"|-calc_projector_test: ERR %f\n",distance[i]);
		
		for(j=0;j< frame1->hbd_nelem; j++){
			
			el1=frame1->hbd_elem[j];
			el2=frame2->hbd_elem[j];
			el3=frame1->hbd_elem[j]; // <-- reference system! to be changed accordingly
			diff=difference->hbd_elem[j];
			int incr=0 ;// just an increment position for the derivative vector (being updated)
			switch(el1->cvtype){
					// case simple cv: just do the difference (no pbc or other stuff)
					// note: if angles are sin and cos component then this is good anyway
				case 1:	
					simple_elem_diff(diff,el1,el2,el3,mtd_data.fplog);
					break;
				case 3:	
					simple_elem_diff(diff,el1,el2,el3,mtd_data.fplog);
					break;
				case 4:	
					simple_elem_diff(diff,el1,el2,el3,mtd_data.fplog);
					break;
				case 5:	
					simple_elem_diff(diff,el1,el2,el3,mtd_data.fplog);
					break;
				case 32:	
					simple_elem_diff(diff,el1,el2,el3,mtd_data.fplog);
					break;
				case 53:	
					// case msd
					msd_elem_diff(diff,el1,el2,el3,mtd_data.fplog,0,-1);
					break;
				default:
					fprintf(mtd_data.fplog,"|- THIS DIFFERENCE IS NOT IMPLEMENTED: CVTYPE %d ",el1->cvtype);
					EXIT();
					break;
			}
			for(k=0;k<el1->ncv;k++){
				fprintf(mtd_data.fplog,"|- DIFFERENCE ELEM %d CV %d IS %f\n",j,k,diff->ref_dist[k]);
			}
		}
		
	}
	
	fprintf(mtd_data.fplog,"|-calc_projector_test: EXITING\n");
	EXIT();
//AVOID
}
void PREFIX calc_intraframe_dist( struct sz_data *pmy_sz ){
//AVOID
	fprintf(mtd_data.fplog,"|-calc_intraframe_dist: ENTERING\n");
	// allocating the matrix
	real ***matrices;
	int i,j,k,n,totcvs;
	n=pmy_sz->number;
	totcvs=pmy_sz->hbd_frameset[0]->hbd_totcvs;
	matrices=(real ***)malloc(n*sizeof( real**));
	for(i=0;i<n;i++){
		matrices[i]=(real **)malloc(totcvs*sizeof( real *));
		for(j=0;j<totcvs;j++){
				matrices[i][j]=(real *)malloc(totcvs*sizeof( real));
		}
	}
	// reading in
	fprintf(mtd_data.fplog,"|-calc_intraframe_dist: READING MATRICES\n");
	FILE *fp;
	fp=fopen(pmy_sz->intraframe_dist_inputfile,"r");
	char *stringa,*cmstring;
	int char_per_number;
	char_per_number=20;
	stringa=(char *)malloc(totcvs*char_per_number*sizeof(char)); // a high limit for stringa
	i=-1;
	while(1){
		cmstring=fgets(&(stringa[0]),totcvs*char_per_number,fp);
		if(cmstring==NULL){break;}
		// is it a newset?
		char newset[7];
		sscanf(cmstring,"%6s",newset);
		if(strstr(newset,"NEWSET")!=NULL){
			 //fprintf(mtd_data.fplog,"|-calc_intraframe_dist: %d %s\n",i,newset);
			 i++;j=0;k=0;
			 if(i>n){
				 fprintf(mtd_data.fplog,"|-calc_intraframe_dist: EXCEEDING NUMBER OF FRAMES\n");EXIT();
			 }
			
		}else{
			char *result = NULL;
			result=strtok( stringa, " " );
			while(result!=NULL){
				//fprintf(mtd_data.fplog,"|-calc_intraframe_dist: %d %d %d %s\n",i,j,k,result);
				matrices[i][j][k]=atof(result);
				result=strtok (NULL, " ");
				if(k==totcvs){
					fprintf(mtd_data.fplog,"|-calc_intraframe_dist: EXCEEDING NUMBER OF CVS\n");EXIT();
				}
				k++;

			}
			j++;k=0;
		}
		
	}
	if(i!=n-1){
		fprintf(mtd_data.fplog,"|-calc_intraframe_dist: number of frames is not matching with input \n");EXIT();
	}
	// close
	fclose(fp);
	fprintf(mtd_data.fplog,"|-calc_intraframe_dist: READING MATRICES DONE!\n");
	
	// calc the distance
	
	struct hybrid_frameset *frame1,*frame2,*difference;

	real *distance;
	
	struct cmap_outpack outpack;
	
	distance=(real *)malloc(pmy_sz->number*sizeof(real));
	
	clone_hybrid_frameset(&difference,pmy_sz->hbd_frameset[0],1,mtd_data.fplog);
	
	// loop over all the frames and use one single matrix as template
	for(i=0;i<pmy_sz->number-1;i++){
		fprintf(mtd_data.fplog,"|-calc_projector_test: calculating projector between frame %d and %d\n",i,i+1);
		frame1=pmy_sz->hbd_frameset[i];
		frame2=pmy_sz->hbd_frameset[i+1];
		hbd_metrics_new(frame1,frame2,&outpack,matrices[i],mtd_data.fplog);
		distance[i]=outpack.err;
		fprintf(mtd_data.fplog,"|-calc_projector_test: ERR %f\n",distance[i]);
	}
	fp=fopen(pmy_sz->intraframe_dist_outputfile,"w");
	for(i=0;i<n-1;i++){
		fprintf(fp,"%12.6f\n",distance[i]);
	}
	fclose(fp);
	
	// printout

	fprintf(mtd_data.fplog,"|-calc_intraframe_dist: EXITING\n");	
//AVOID
}
void PREFIX calc_intraframe_diff( struct sz_data *pmy_sz ){
//AVOID
	fprintf(mtd_data.fplog,"|-calc_intraframe_diff: ENTERING\n");
	int i,j,k;
	struct hybrid_frameset *frame1,*frame2,*difference;
	struct hybrid_elem *el1,*el2,*el3,*diff;
		
	clone_hybrid_frameset(&difference,pmy_sz->hbd_frameset[0],1,mtd_data.fplog);
	
	FILE *fp;
	fp=fopen(pmy_sz->intraframe_diff_outputfile,"w");
	// loop over all the frames and use one single matrix as template
	for(i=0;i<pmy_sz->number-1;i++){
		
		
		fprintf(mtd_data.fplog,"|-calc_projector_test: calculating projector between frame %d and %d\n",i,i+1);
		frame1=pmy_sz->hbd_frameset[i];
		frame2=pmy_sz->hbd_frameset[i+1];
		
		calc_diff_twoframes(frame1,frame2,frame1,difference);

		fprintf(fp,"NEWSET\n");
		for(j=0;j< frame1->hbd_nelem; j++){
			diff=difference->hbd_elem[j];
			for(k=0;k<diff->ncv;k++){
				fprintf(mtd_data.fplog,"|- DIFFERENCE ELEM %d CV %d IS %f\n",j,k,diff->ref_dist[k]);
				fprintf(fp,"%12.6f\n",diff->ref_dist[k]);
			}
		}
		
	}
	fclose(fp);
	
	fprintf(mtd_data.fplog,"|-calc_intraframe_dist: EXITING\n");		
//AVOID
}

void PREFIX calc_diff_twoframes(struct hybrid_frameset *v1,  struct hybrid_frameset *v2, struct hybrid_frameset *v3,struct hybrid_frameset *d1 ){
//AVOID
	int i,j,k,l,m;
	struct hybrid_elem *el1,*el2,*el3,*diff1;
	
	for(i=0;i< v1->hbd_nelem; i++){
		el1=v1->hbd_elem[i];
		el2=v2->hbd_elem[i];
		el3=v3->hbd_elem[i];
		diff1=d1->hbd_elem[i];
		
		int incr=0 ;// just an increment position for the derivative vector (being updated)
		switch(el1->cvtype){
				// case simple cv: just do the difference (no pbc or other stuff)
				// note: if angles are sin and cos component then this is good anyway
			case 1:	
				simple_elem_diff(diff1,el1,el2,el3,mtd_data.fplog);
				break;
			case 3:	
				simple_elem_diff(diff1,el1,el2,el3,mtd_data.fplog);
				break;
			case 4:	
				simple_elem_diff(diff1,el1,el2,el3,mtd_data.fplog);
				break;
			case 5:	
				simple_elem_diff(diff1,el1,el2,el3,mtd_data.fplog);
				break;
			case 32:	
				simple_elem_diff(diff1,el1,el2,el3,mtd_data.fplog);
				break;
			case 53:	
				// case msd
                                //check_msd_elem_diff(diff1,el1,el2,el3,mtd_data.fplog);
				// last flag means that you want the rot_mat_$i.dat file
				msd_elem_diff(diff1,el1,el2,el3,mtd_data.fplog,0,i);
				break;
			default:
				fprintf(mtd_data.fplog,"|- THIS DIFFERENCE IS NOT IMPLEMENTED: CVTYPE %d ",el1->cvtype);
				EXIT();
				break;
		}
		
	}
	
//AVOID
}
void PREFIX calc_twoframe_dist( struct sz_data *pmy_sz , int first, int second,  char *matrixfile, FILE *fplog){
//AVOID
	int i,j,k;
	int totcvs=pmy_sz->hbd_frameset[0]->hbd_totcvs;
	// alloc the matrix
	real **matrix;
	matrix=(real **)malloc(totcvs*sizeof( real*));
	for(j=0;j<totcvs;j++){
		matrix[j]=(real *)malloc(totcvs*sizeof( real));
	}
	// read it
	FILE *fp;
	fp=fopen(matrixfile,"r");
	char *stringa,*cmstring;
	int char_per_number;
        int len;	
	char_per_number=20;
	j=0;k=0;
	stringa=(char *)malloc(totcvs*char_per_number*sizeof(char)); // a high limit for stringa
	while(1){
		cmstring=fgets(&(stringa[0]),totcvs*char_per_number,fp);
		if(cmstring==NULL){break;}
                len = strlen(cmstring); 
                if( cmstring[len-1] == '\n' )cmstring[len-1] = 0;
		// is it a newset?
		char newset[7];
		sscanf(cmstring,"%6s",newset);
		if(strstr(newset,"NEWSET")!=NULL){
			j=0;k=0;			
		}else{
			char *result = NULL;
			result=strtok( stringa, " " );
			while(result!=NULL){
				//fprintf(mtd_data.fplog,"|-calc_twoframe_dist: %d %d XX%sXX\n",j,k,result);
				matrix[j][k]=atof(result);
				if(k==totcvs){
					fprintf(fplog,"|-calc_intraframe_dist: EXCEEDING NUMBER OF CVS. I was expecting a max of %d\n",totcvs);EXIT();
				}
				result=strtok (NULL, " ");
				k++;
			}
			j++;k=0;
		}
	}
	// close
	fclose(fp);
	//fprintf(fplog,"|-calc_twoframe_dist: READING MATRIX DONE!\n");
	
	struct hybrid_frameset *frame1,*frame2;
	
	real distance;
	
	struct cmap_outpack outpack;
		
	frame1=pmy_sz->hbd_frameset[first];
	frame2=pmy_sz->hbd_frameset[second];
	hbd_metrics_new(frame1,frame2,&outpack,matrix,fplog);
	fprintf(mtd_data.fplog,"|-calc_twoframe_dist: TWOFRAMES_DIST:(squared) %f\n",outpack.err);
	EXIT();
	
//AVOID
}
void PREFIX reparam_with_multiple_matrix( struct sz_data *pmy_sz , char *filename,   real maxerr, int maxiter, real smooth){
//AVOID
	// malloc
	int i,j,k,n,totcvs;
	n=pmy_sz->number;
	totcvs=pmy_sz->hbd_frameset[0]->hbd_totcvs;
	pmy_sz->mathybrid_blocks=(real ***)malloc(n*sizeof( real**));
	for(i=0;i<n;i++){
		pmy_sz->mathybrid_blocks[i]=(real **)malloc(totcvs*sizeof( real *));
		for(j=0;j<totcvs;j++){
			pmy_sz->mathybrid_blocks[i][j]=(real *)malloc(totcvs*sizeof( real));
		}
	}
	// reading in
	fprintf(mtd_data.fplog,"|-reparam_with_multiple_matrix: READING MATRICES from %s\n",filename);
	FILE *fp;
	fp=fopen(filename,"r");
	char *stringa,*cmstring;
	int char_per_number;
	char_per_number=20;
	stringa=(char *)malloc(totcvs*char_per_number*sizeof(char)); // a high limit for stringa
	i=-1;
	while(1){
		cmstring=fgets(&(stringa[0]),totcvs*char_per_number,fp);
		if(cmstring==NULL){break;}
		
		// is it a newset?
		char newset[7];
		sscanf(cmstring,"%6s",newset);
		if(strstr(newset,"NEWSET")!=NULL){
			//fprintf(mtd_data.fplog,"|-calc_intraframe_dist: %d %s\n",i,newset);
			i++;j=0;k=0;
			if(i>n){
				fprintf(mtd_data.fplog,"|-reparam_with_multiple_matrix: EXCEEDING NUMBER OF FRAMES\n");EXIT();
			}
			
		}else{
			char *result = NULL;
			result=strtok( stringa, " " );
			while(result!=NULL){
				//fprintf(mtd_data.fplog,"|-reparam_with_multiple_matrix: %d %d %d %s\n",i,j,k,result);
				pmy_sz->mathybrid_blocks[i][j][k]=atof(result);
				result=strtok (NULL, " ");
				if(k==totcvs){
					fprintf(mtd_data.fplog,"|-reparam_with_multiple_matrix: EXCEEDING NUMBER OF CVS\n");EXIT();
				}
				k++;
				
			}
			j++;k=0;
		}
		
	}
	if(i!=n-1){
		fprintf(mtd_data.fplog,"|-reparam_with_multiple_matrix: number of frames is not matching with input \n");EXIT();
	}
	// close
	fclose(fp);
	fprintf(mtd_data.fplog,"|-reparam_with_multiple_matrix: READING MATRICES DONE!\n");
	// now do the reparam with these matrix
	
	reparam_bernd_path(pmy_sz,maxerr,maxiter, smooth);
	
	// call reparametrization
//AVOID
};
void PREFIX point_to_matrix(int i,struct sz_data *my_sz){
//AVOID
	my_sz->mathybrid=my_sz->mathybrid_blocks[i];
//AVOID
}
void PREFIX interpolate_bernd_path( struct sz_data *pmy_sz , int target_npts ) {
//AVOID
              int i,j,k,npts,nmiddle;
              fprintf(mtd_data.fplog,"|-INTERPOLATION:  .....\n");
              //  
              npts=pmy_sz->number;  
              fprintf(mtd_data.fplog,"|-NPTS %d \n",npts);
              // mallocs some useful stuff
              real *len,*sumlen,*dr,*sfrac;  
              real *dummy1,*dummy2,*dummy3,*dummy4; // set of dummies that will be trashed 
              real *dd;
              real *tmp;
              int *closest_fixed,first,second,*start_fixed;
              len    =(real *)calloc(npts,sizeof(real));
              sumlen =(real *)calloc(npts,sizeof(real));
              dr     =(real *)calloc(npts,sizeof(real));
              dd     =(real *)calloc(npts,sizeof(real));
              sfrac  =(real *)calloc(npts,sizeof(real));
              tmp=(real *)malloc(npts*sizeof(real)); 
	
	      // make the allocation of an array containing the differences
	      struct hybrid_frameset **diff_vec;
	      fprintf(mtd_data.fplog,"|-interpolate_bernd_path: ALLOCATING THE DIFFERENCE VECTORS\n");
	      diff_vec=(struct hybrid_frameset **)malloc((npts-1)*sizeof(struct hybrid_frameset*));
	      for(i=0;i<npts-1;i++)clone_hybrid_frameset(&diff_vec[i],pmy_sz->hbd_running,1,mtd_data.fplog);
	      real *d1,*d2,*d3; // derivatives of the vector
	      snew(d1,pmy_sz->hbd_running->hbd_totcvs);
	      snew(d2,pmy_sz->hbd_running->hbd_totcvs);
	      snew(d3,pmy_sz->hbd_running->hbd_totcvs);
	
	      // make some temporary vectors of frameset to be copied between one iteration and the other
	      struct hybrid_frameset **temp_vec;
	      fprintf(mtd_data.fplog,"|-interpolate_bernd_path: ALLOCATING THE TEMP VECTORS\n");

	      temp_vec=(struct hybrid_frameset **)malloc(target_npts*sizeof(struct hybrid_frameset*));
	      for(i=0;i<target_npts;i++)clone_hybrid_frameset(&temp_vec[i],pmy_sz->hbd_running,0,mtd_data.fplog);
	
	      struct hybrid_frameset *frame;
	      struct hybrid_elem *elem;
              int kk;
	
              for(i=0;i<npts-1;i++){ // loop over all the possible distances
						 
		       // if multiple matrices are defined then use them 
		       // this let point mathybrid to the right block
		       if(pmy_sz->mathybrid_blocks!=NULL)point_to_matrix(i,pmy_sz);
					  
		       // the difference here uses the frame i as reference
		       // so diff[0] is the diff betwen frame 1 and 0 centered on 0
							
                        dd[i]=hbd_distanddiff_ref(pmy_sz->hbd_frameset[i],pmy_sz->hbd_frameset[i+1],pmy_sz->hbd_frameset[i],pmy_sz->mathybrid,
												  d1,d2,d3,diff_vec[i],mtd_data.fplog);
			fprintf(mtd_data.fplog,"|-DISTSQUARED %3d TO %3d : %12.6f\n",i,i+1,dd[i]);
			fflush(mtd_data.fplog);

			dd[i]=sqrt(dd[i]);

              }

	      sumlen[0]=dd[0]; // at the end of a segment within fixed point reset the value 
              for(i=1;i<npts-1;i++){
			  sumlen[i]=sumlen[i-1]+dd[i]; // keep on summing up	
              }

              // now make the real calculation		
 	      double mystep;
	      mystep= sumlen[npts-2]/(double)(target_npts-1);
	      fprintf(mtd_data.fplog,"|- STEP SHOULD BE TAKEN OF %f \n",mystep);
              for(i=0;i<target_npts;i++){		
		double this_step;			
		this_step=i*mystep;
		int myindex;
		myindex=0;
                for(j=0;j<npts-1;j++){		
			if(this_step<=sumlen[j]) {myindex=j; break;}	
		}
		fprintf(mtd_data.fplog,"|-FOR STEP %4d LENGTH %12.6f TAKE INTERVAL FROM %3d AND %3d (MAX DIST %12.6f ) \n",i,this_step,myindex,myindex+1,sumlen[myindex]);   
		// now do the interpolation in the real sense
		double delta;
		if (myindex>0){delta=(this_step-sumlen[myindex-1])/dd[myindex];}
		else{ delta=(this_step-sumlen[myindex-1])/dd[myindex]; }
                if(pmy_sz->mathybrid_blocks!=NULL)point_to_matrix(myindex,pmy_sz);

		// print some info so that copying of starting points is simple
		if ( pow(this_step-sumlen[myindex],2)>pow(this_step-sumlen[myindex+1],2)   ){
			fprintf(mtd_data.fplog,"|- COPYINFO: closer_config_to_frame %d is_from_frame %d \n",i+1,myindex+2); 
		}else{
			fprintf(mtd_data.fplog,"|- COPYINFO: closer_config_to_frame %d is_from_frame %d \n",i+1,myindex+1); 
		}

		//fprintf(mtd_data.fplog,"|-DELTA %f \n",delta);
		do_step_bernd(temp_vec[i],pmy_sz->hbd_frameset[myindex],pmy_sz->mathybrid,delta,pmy_sz->hbd_frameset[myindex+1],pmy_sz->hbd_frameset[myindex]);
		// adjust the "fixed" label	
		temp_vec[i]->fixed=0;
		if ( delta<1.e-4 && pmy_sz->hbd_frameset[myindex]->fixed ) {
			temp_vec[i]->fixed=1;
		}
		if ( delta>(1-1.e-4) && pmy_sz->hbd_frameset[myindex+1]->fixed ) {
			temp_vec[i]->fixed=1;
		}
	
		fflush(mtd_data.fplog);	
              } 
	      // now dump results and exit 
	      pmy_sz->number=target_npts;	
              pmy_sz->hbd_frameset=temp_vec;		
	      dump_frameset_formatted(pmy_sz);
	

//AVOID
}
void PREFIX interpolate_with_multiple_matrix( struct sz_data *pmy_sz , char *filename, int npts) {
//AVOID
	// malloc
	int i,j,k,n,totcvs;
	n=pmy_sz->number;
	totcvs=pmy_sz->hbd_frameset[0]->hbd_totcvs;
	pmy_sz->mathybrid_blocks=(real ***)malloc(n*sizeof( real**));
	for(i=0;i<n;i++){
		pmy_sz->mathybrid_blocks[i]=(real **)malloc(totcvs*sizeof( real *));
		for(j=0;j<totcvs;j++){
			pmy_sz->mathybrid_blocks[i][j]=(real *)malloc(totcvs*sizeof( real));
		}
	}
	// reading in
	fprintf(mtd_data.fplog,"|-interpolate_with_multiple_matrix: READING MATRICES from %s\n",filename);
	FILE *fp;
	fp=fopen(filename,"r");
	char *stringa,*cmstring;
	int char_per_number;
	char_per_number=20;
	stringa=(char *)malloc(totcvs*char_per_number*sizeof(char)); // a high limit for stringa
	i=-1;
	while(1){
		cmstring=fgets(&(stringa[0]),totcvs*char_per_number,fp);
		if(cmstring==NULL){break;}
		
		// is it a newset?
		char newset[7];
		sscanf(cmstring,"%6s",newset);
		if(strstr(newset,"NEWSET")!=NULL){
			//fprintf(mtd_data.fplog,"|-calc_intraframe_dist: %d %s\n",i,newset);
			i++;j=0;k=0;
			if(i>n){
				fprintf(mtd_data.fplog,"|-interpolate_with_multiple_matrix: EXCEEDING NUMBER OF FRAMES\n");EXIT();
			}
			
		}else{
			char *result = NULL;
			result=strtok( stringa, " " );
			while(result!=NULL){
				//fprintf(mtd_data.fplog,"|-reparam_with_multiple_matrix: %d %d %d %s\n",i,j,k,result);
				pmy_sz->mathybrid_blocks[i][j][k]=atof(result);
				result=strtok (NULL, " ");
				if(k==totcvs){
					fprintf(mtd_data.fplog,"|-interpolate_with_multiple_matrix: EXCEEDING NUMBER OF CVS\n");EXIT();
				}
				k++;
				
			}
			j++;k=0;
		}
		
	}
	if(i!=n-1){
		fprintf(mtd_data.fplog,"|-interpolate_with_multiple_matrix: number of frames is not matching with input \n");EXIT();
	}
	// close
	fclose(fp);
	fprintf(mtd_data.fplog,"|-interpolate_with_multiple_matrix: READING MATRICES DONE!\n");
	// now do the reparam with these matrix
	
	interpolate_bernd_path(pmy_sz,npts);
	
//AVOID
}

int PREFIX read_sz_chiral(struct sz_data *my_sz, FILE *fplog){
	fprintf(fplog,"|- CHIRALITY READER: ENTERING \n");  
        int l,i,j,k,found;
        char *str,ic[3],str2[100],stringa[300];
        l=0;
/*
  * allocate the pointers
 */
        my_sz->chiral_frameset=(struct chirality_frameset **)malloc((my_sz->number)*sizeof(struct chirality_frameset *));
        for(i=0;i< my_sz->number;i++){
                my_sz->chiral_frameset[i]=(struct chirality_frameset *)malloc(sizeof(struct chirality_frameset)) ;
        }

/*
 *  build names
 */
        str=&(my_sz->names[0]);
        for (i=1;i<= my_sz->number ;i++){
                strcpy(str2,my_sz->names);
                if(my_sz->targeted_on==0){
                     if(i<10){
                      ic[0]='0'+i;
                      ic[1]='\0';}
                     else if(i<100) {
                      ic[0]='0'+i/10 ;
                      ic[1]='0'+i%10 ;
                      ic[2]='\0';
                     }
                     else{
                       plumed_error("|--read_sz_input: TOO MANY FRAMES REQUIRED FOR NAME BUILDING!");
                     }
                     strcat(str2,ic);
                     strcat(str2,".G");
                }
                fprintf(fplog,"|--%s\n",str2);
		// this is not particularly complex: try to read it from here	
		FILE *fp;
		fp=fopen(str2,"r");
                if (fp == NULL){
                         char buf[1024];
                         sprintf(buf,"UNSUCCESSFULL OPENING FOR FILE %s",str2);
                         plumed_error(buf);
                }
		// check how many fields are there: last one is the value, all the previous are for combination 
		while(1){
			  readagain:
  			  str=fgets(stringa,200,fp);
                          if(str==NULL)break;
                          if (feof(fp))break;	
		          char *result = NULL;
                          int iii,jjj; iii=0; 
		          real val;		
                          result = strtok( stringa, " \t" );	 		
			  while( result != NULL ) {
				if (iii<15){jjj=atoi(result);fprintf(fplog,"|-INDEX %d \n",jjj);}
				else{val=atof(result);fprintf(fplog,"|-VAL %f \n",val);}
				result = strtok( NULL, " \t");
			 	iii++;	
			  }
                }

 	}

	fprintf(fplog,"|- CHIRALITY READER: EXITING \n");  
}

