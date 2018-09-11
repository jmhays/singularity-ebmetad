#include <vector>
#include <valarray>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

extern "C" {
#include "metadyn.h"
;}
#include "recon_types.h"
#include "recon_utils.h"
#include "recon_basins.h"

extern "C" {

class structureObj{
private:    
public:
    rp::rVector x,y,z;
    void storeCoords( int *natoms, double *pos );
};

extern void readhistoryframe_(int* trjfile,int* natoms, char* ele, double* box, double* pos, int* ios);
extern void readpdbframe_(int* trjfile, int* natoms, char* ele, char* resname, int* resno, double* box, double* pos, double* lenunits, int* safe_read );

void search_traj_( int* ftype, int* natoms, int* imcon, double* lenunits, int* trjfile, char* metainp ) {

   std::cout<<"Check on name "<<metainp<<std::endl;

   // Make some fake mass and charge arrays
   double mass[(*natoms)], charge[(*natoms)];
   rp::index nat; nat=(*natoms);
   for(rp::index i=0;i<nat;++i){ mass[i]=1.0; charge[i]=0.0;}
   
   // Make a fake imcon  
   int pbc=0; double box[3];
   if( (*imcon)>0 && (*ftype)==1 ){ pbc=1; }
   
   // Call plumed to read in input file (and get CVs)  
   int ll; double h,s,w; int ncolvars; double per_dat[nconst_max]; 
   double tstep=1.0; // This is a fake timestep (maybe fix this in the future although is not that important)
   init_metadyn_( natoms, &tstep, mass, charge, &pbc, box, metainp, &ncolvars, per_dat, &w, &h, &s, ll );

   rp::index nbasins; rp::index ncv; ncv=ncolvars; rp::rVector period(ncv);
   for(rp::index i=0;i<ncv;++i) period[i]=per_dat[i];
   rp::real sizeParam, height, width; height=h; width=w; sizeParam=s;

   std::cout<<"Found "<<ncolvars<<" collective coordinates and parameters "<<height<<" and "<<sizeParam<<" in input"<<std::endl;

   // First read in basins
   std::vector<rp::basinObj> basinList; rp::index periodic=0;
   if( period[0]!=0. ) periodic=1; 
   readBasins( ncv, periodic, sizeParam, width, height, basinList );
   nbasins=basinList.size();

   // Open cvfile
   std::string filename; filename="TRAJ_COLVARS";
   FILE* cvfile; cvfile=fopen(filename.c_str(),"w+");

   double pos[3*(*natoms)]; int safe_read=0; 
   std::vector<std::string> resString((*natoms)), elString((*natoms)); 
   char ele[5*(*natoms)]; char resname[4*(*natoms)]; int resno[(*natoms)];
   // Declare stuff for calculating colvars
   rp::rVector colvars(ncv); double cvout[ ncv ]; rp::index firsttime=1; 
   std::vector<structureObj> structureList(nbasins); std::vector<double> mineng( nbasins ); double energy;
   int nn=0;
   while (safe_read==0) {

       nn++;
       std::cout<<"Frame "<<nn<<std::endl;

       if( (*ftype)==1 ){ readhistoryframe_( trjfile, natoms, ele, box, pos, &safe_read ); }
       else if( (*ftype)==2 ){ readpdbframe_( trjfile, natoms, ele, resname, resno, box, pos, lenunits, &safe_read ); }

       if(safe_read!=0) break;

       // Call plumed to calculate the collective coordinates
       cv_calculation_(box, pos, &ncolvars, cvout);
       fprintf(cvfile,"%d ",nn);
 
       // Transfer CVs to rVector
       for(rp::index i=0;i<ncv;++i){fprintf(cvfile,"%f ",cvout[i]); colvars[i]=cvout[i];}
       fprintf(cvfile,"\n");  

       if (firsttime==1) {
          // First sort out the residue names and atom names
          std::istringstream ele_in (ele); std::istringstream res_in (resname);
          for(rp::index i=0;i<(*natoms);++i){
             ele_in.seekg(i*5); ele_in>>elString[i];
             if( (*ftype)==1 ){resno[i]=0; resString[i]="DCB";}
             else if( (*ftype)==2 ){res_in.seekg(i*4); res_in>>resString[i]; }
          }
          // Now store all the basins
          for(rp::index i=0;i<nbasins;++i){ 
              mineng[i]=basinList[i].dist( period, colvars );
              structureList[i].storeCoords( natoms, pos ); 
          }
          firsttime=0;
       }
       else{
          energy=basinList[0].dist( period, colvars );
          for(rp::index i=0;i<nbasins;++i){
            energy=basinList[i].dist( period, colvars );
            if( energy<mineng[i] ){ structureList[i].storeCoords( natoms, pos ); mineng[i]=energy; }
          }
       }
   }
   fclose(cvfile);

   // Now write out the structures
   filename="basins.pdb";
   FILE* outfile; outfile=fopen(filename.c_str(),"w+");
   for(rp::index i=0;i<nbasins;++i){
      std::cout<<"Final energy "<<i<<" "<<mineng[i]<<std::endl;
      for(rp::index j=0;j<(*natoms);++j){
         fprintf(outfile,"ATOM  %5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                  int(j),elString[j].c_str(),resString[j].c_str(),resno[j],
                  (*lenunits)*structureList[i].x[j],(*lenunits)*structureList[i].y[j],(*lenunits)*structureList[i].z[j],1.0,1.0);
      }
      fprintf(outfile,"END\n");
   }
   fclose(outfile);

   return; 
}

void structureObj::storeCoords( int *natoms, double *pos ){
  if ( x.size()!=(*natoms) ){ x.resize((*natoms)); y.resize((*natoms));z.resize((*natoms)); } 
  for(rp::index i=0;i<(*natoms);++i){ x[i]=pos[i]; y[i]=pos[i+(*natoms)]; z[i]=pos[i+2*(*natoms)]; } 
}

}
