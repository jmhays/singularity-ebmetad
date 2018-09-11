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
#include "recon_bespoke.h"

// Routine to read the plumed file
rp::index read_plumed_file( char* metainp, rp::rVector& period, rp::real& sizeParam );
// New CV read in routine
void read_colvars( const rp::index ncv, char* cvfilename, std::vector<rp::rVector>& cvs );

// This extern C required so that we can call from fortran 
extern "C" {
void isomap_( int* num_neighbours, int* expo, int* outdim, char* metainp, char* cvfile ) {

   std::cout<<"Check on names "<<metainp<<" "<<cvfile<<std::endl;
   // Read the plumed file to get the number of CVS
   rp::rVector period; rp::index ncv; rp::real sizeParam;
   ncv=read_plumed_file( metainp, period, sizeParam );
   std::cout<<"Found "<<ncv<<" collective coordinates in input"<<std::endl;

   // Transfer input pointers to normal variables
   rp::index nneigh, dimout, eps; nneigh=(*num_neighbours); dimout=(*outdim); eps=(*expo);

   // Setup a bespoke CV object - this actually does the isomap
   rp::bespokeCVObj bespokeCV( ncv, dimout ); bespokeCV.embedPoints( eps, nneigh, period ); 

   // Setup the alpha helix vector
   rp::rVector alphaHelix(ncv);
   for(rp::index i=0;i<ncv;++i){
      if(i%2==0){alphaHelix[i]=-R_PI/3.;} else{alphaHelix[i]=-R_PI/4.;}
   }

   // Open a file to output the embedded basins to
   std::ofstream nodefile; nodefile.open( "basin.pops" );
   std::ofstream cvdefs; cvdefs.open( "bespoke.defs" ); 

   // Get the embedded vectors and the basin centers from the bespoke CVs 
   std::vector<rp::rVector> outVecs, centers; std::vector<rp::basinObj> basinList;
   bespokeCV.getEmbeddedVectors( outVecs ); bespokeCV.getCenters( centers ); bespokeCV.getBasins( basinList );

   // Now output this data
   rp::index nbasins; nbasins=outVecs.size();

   cvdefs<<"NBASINS "<<nbasins<<"\n \n"; cvdefs<<"BASINS>"<<std::endl;

   for(rp::index i=0;i<nbasins;++i){
      nodefile<<outVecs[i]<<" "<<calcDistance( period, centers[i], alphaHelix )<<std::endl;
      cvdefs<<basinList[i]<<std::endl;
   }
   nodefile.close();

   cvdefs<<"BASINS<"<<"\n \n"; cvdefs<<"POINTS>"<<std::endl;
   for(rp::index i=0;i<nbasins;++i){ cvdefs<<outVecs[i]<<std::endl; }
   cvdefs<<"POINTS<"<<std::endl; cvdefs.close();

   // Now read in the colvars
   std::vector<rp::rVector> all_colvars; read_colvars( ncv, cvfile, all_colvars );
   rp::index ncolvars; ncolvars=all_colvars.size();

   // Open a file for outputting the embedded points
   std::ofstream pfile; pfile.open( "embedded_cvs.dat" );

   // And embed the colvars in the low dimensionality space
   rp::rVector outVec( dimout ); rp::rMatrix derivatives( dimout, ncv );
   for(rp::index i=0;i<ncolvars;++i){
      bespokeCV.embedPoint( all_colvars[i], outVec, derivatives ); pfile<<outVec<<std::endl;
   }
   pfile.close();

   return;
}
}

void read_colvars( const rp::index ncv, char* cvfilename, std::vector<rp::rVector>& cvs ){
   // Open file containing collective coordinates
   std::ifstream cvfile; cvfile.open(cvfilename);
   if (!cvfile.good()) ERROR("Something bad happened while opening collective coordinates file for restart");

   std::string line; getline(cvfile, line);     // Read in the first line of crap

   rp::real time; rp::rVector colvars(ncv); rp::index nnodes=0;
   while ( getline(cvfile, line) ) {
      std::istringstream cvin (line);
      cvin>>time;      // Get the line and read the timestep
      for(rp::index i=0;i<ncv;++i){ cvin>>colvars[i]; }
      nnodes++; cvs.push_back( colvars );
   }
   std::cout<<"Read in "<<nnodes<<" points from COLVAR file"<<std::endl;
   return;
}

rp::index read_plumed_file( char* metainp, rp::rVector& period, rp::real& sizeParam ){

  // Make some fake mass and charge arrays
   int natoms=10000;
   double mass[natoms], charge[natoms];
   for(rp::index i=0;i<natoms;++i){ mass[i]=1.0; charge[i]=0.0;}

   // Make a fake imcon  
   int pbc=1; double box[3];

   // Call plumed to read in input file (and get CVs)  
   int ll; double h,s,w; int ncolvars; double per_dat[nconst_max];
   double tstep=1.0; // This is a fake timestep (maybe fix this in the future although is not that important)
   init_metadyn_( &natoms, &tstep, mass, charge, &pbc, box, metainp, &ncolvars, per_dat, &w, &h, &s, ll );
   // Set value of size parameter from what was read in
   sizeParam=s;

   std::cout<<"I have read in the plumed file and found "<<ncolvars<<" colvars"<<std::endl;

   rp::index ncv; ncv=ncolvars; period.resize( ncv );
   for(rp::index i=0;i<ncv;++i){ period[i]=per_dat[i]; }

   return ncv;

}
