#ifdef RECONMETAD 
#ifdef CVS
#include "recon_bespoke.h"
//  #include "recon_cgminimize.h"

namespace rp {

long search1(const real& x, const rVector& xlist, long jold=0);

// This is f(x) = x
inline void bfunction::dist( const real& x, real& f, real& df ) const
{ f=x; df=1; }

// This is f(x) = 1.0 - ( 1 / ( 1 + (x/sigma)^2 ) )
inline void bfunction::sigmoid(const real &x, real& f, real& df) const
{ double sx=x*rsigma; sx=1.0/(1.0+sx*sx); f=1.0-sx; df=x*(sx*sx)*dsigma_sq; }

// This is f(x) = 1.0 - ( 1 / ( 1 + x/sigma) )
inline void bfunction::compress(const real& x, real& f, real& df) const
{ double sx=x*rsigma; sx=1.0/(1.0+sx); f=1.0-sx; df=(sx*sx)*rsigma; }

// This is f(x) = 1.0 - ( 1 + ( 2^(a/b) - 1 )(x/sigma)^a )^-(b/a)
inline void bfunction::general(const real& x, real& f, real& df) const
{ double sx=x*rsigma; sx=c*pow( sx, a ); f=pow( 1.0 + sx, d ); df=b*sx/x*f/(1.0+sx); f=1.0-f; }

void bfunction::set_mode( const real& aa, const real& bb, const real& sig ){
   if( sig<0 ){ 
      pfdf = &bfunction::dist; 
   } else if ( aa==1 && bb==1 ){
      rsigma = 1.0 / sig;
      pfdf = &bfunction::compress;  
   } else if ( aa==2 && bb==2 ){
      rsigma = 1.0 / sig; dsigma_sq = 2.0*rsigma*rsigma;
      pfdf = &bfunction::sigmoid;
   } else {
      a = aa; b = bb; rsigma = 1.0 / sig; 
      c = pow(2, a/b) - 1; d = - b / a;
      pfdf = &bfunction::general;
   }
}

void bespokeCVObj::setup( const std::string efilename, const rp::rVector& p, const real& smoo, const std::vector<rp::index>& nsp, const std::vector<rp::index>& ngr , FILE* fplog ){
#ifdef DEBUG
   if( p.size()!=ncv ) ERROR("SIZE MISMATCH FOR PERIODS");
#endif

  // std::pout<<"|--BESPOKE CV DETAILS, ";

  kt=smoo;  // Copy smoothing parameter
  // Copy number of points in grid and spline
  for(index i=0;i<dimout;++i){ ngrid[i]=ngr[i]; nspline[i]=nsp[i]; }

  // Setup everything for interpolation
  interpolator.setupGrid( nspline[0], nspline[1] ); gx.resize( nspline[0] ); gy.resize( nspline[1] ); 
  sgridU.resize( nspline[0], nspline[1] ); sgridDU.resize( 3, nspline[0], nspline[1], dimout ); 

  // Setup cv derivative interpolators
  for(index i=0;i<ncv;++i){
     dinterpolator[i].setupGrid( nspline[0], nspline[1] );
     dsgridU[i].resize( nspline[0], nspline[1] ); 
     dsgridDU[i].resize( 3, nspline[0], nspline[1], dimout );
  }

  // This reads in a file containing all the information on the embedding
  std::ifstream ifile; ifile.open( efilename.c_str() );
  if (!ifile.good()) ERROR("Something bad happened while opening embedding data file for restart"); 

  // Start output to file
  fprintf(fplog,"\n|--BESPOKE CV DETAILS \n");

  // Read in the input file 
  std::string line; rVector lowP(dimout), highP(ncv); real w, total_weight=0.0;
  index n, nbas=0, np=0, nl=0, nw=0, flag=0; 
  bool read_weights=false, p_read=false, bas_read=false, l_read=false, w_read=false; 
  while ( getline( ifile, line ) ){
     std::istringstream input (line);

     // Read the number of basins in the basin file
     if(line.find("NLANDMARKS")!=std::string::npos && flag==0){
        std::string tmpstr;
        while( input.good() ){ input>>tmpstr; if(tmpstr.find("NLANDMARKS")!=std::string::npos){ flag=1; break; } }
        if(flag==1){ if( ( input>>nland ).fail() ) ERROR("FAILED TO READ NLANDMARKS SAFELY"); }
        else{ ERROR("FAILED TO READ NUMBER OF BASINS FROM INPUT FILE"); }
        // std::pout<<"NUMBER OF LANDMARKS "<<nland<<", VALUE OF KT IN BOLTZMAN AVERAGES "<<kt<<std::endl;; 
        fprintf(fplog,"|- NUMBER OF LANDMARKS %d, VALUE OF KT IN BOLTZMANN AVERAGES %f \n", nland, kt);
        // Resize everything that is used in definition of bespoke CVs 
        bstress.resize( nland, ncv, dimout ); flag=2;
        // Copy period data to bespoke stress object
        bstress.setPeriods( p );
     }
     else if(line.find("NLANDMARKS")!=std::string::npos){ 
        ERROR("FOUND NLANDMARKS MORE THAN ONCE IN INPUT"); 
     }
     // These find the instructions for the switching function in the low and high D spaces
     else if(line.find("LOW_D_FUNCTION")!=std::string::npos) {
        std::string tmpstr, fname; real sigma=-1.0, a=0, b=0;
        while( input.good()) { 
           input>>tmpstr; 
           if(tmpstr.find("TYPE")!=std::string::npos){ if( (input>>fname).fail() ) ERROR("FAILED TO READ LOW_D FUNCTION NAME CORRECTLY"); }
           else if(tmpstr.find("SIGMA")!=std::string::npos){ if( (input>>sigma).fail() ) ERROR("FAILED TO READ LOW_D FUNCTION SIGMA CORRECTLY"); }
           else if(tmpstr.find("POWERS")!=std::string::npos){ if( (input>>a>>b).fail() ) ERROR("FAILED TO READ LOW_D FUNCTION POWERS CORRECTLY"); }
        }
        bstress.set_function( "LOW_D", fname, a, b, sigma, fplog );          
     } 
     else if(line.find("HIGH_D_FUNCTION")!=std::string::npos) {
        std::string tmpstr, fname; real sigma=-1.0, a=0, b=0;
        while( input.good()) { 
           input>>tmpstr;
           if(tmpstr.find("TYPE")!=std::string::npos){ if( (input>>fname).fail() ) ERROR("FAILED TO READ HIGH_D FUNCTION NAME CORRECTLY"); }
           else if(tmpstr.find("SIGMA")!=std::string::npos){ if( (input>>sigma).fail() ) ERROR("FAILED TO READ HIGH_D FUNCTION SIGMA CORRECTLY"); }
           else if(tmpstr.find("POWERS")!=std::string::npos){ if( (input>>a>>b).fail() ) ERROR("FAILED TO READ HIGH_D FUNCTION POWERS CORRECTLY"); }
        }
        bstress.set_function( "HIGH_D", fname, a, b, sigma, fplog );
     }
     // This finds the ends of the blocks of basins and embedded points
     else if(line.find("HIGH_D<")!=std::string::npos && flag==2 ){
        if( !bas_read ) ERROR("FOUND END OF HIGH_D BLOCK BEFORE START");
        bas_read=false; if(nbas!=nland) ERROR("MISMATCH FOR NUMBER OF LANDMARKS (HIGH D) IN INPUT");
     }
     else if(line.find("LOW_D<")!=std::string::npos && flag==2){
        if( !p_read ) ERROR("FOUND END OF LOW_D BLOCK BEFORE START");
        p_read=false; if(np!=nland) ERROR("MISMATCH FOR NUMBER OF LANDMARKS (LOW D) IN INPUT");
     }
     else if(line.find("WEIGHTS<")!=std::string::npos){
        if( !w_read ) ERROR("FOUND END OF WEIGHTS BLOCK BEFORE START");
        w_read=false; if(nw!=nland) ERROR("WRONG NUMBER OF WEIGHTS FOR BESPOKE CVS");
     }
     else if(line.find("LIMITS<")!=std::string::npos){
        if( !l_read ) ERROR("FOUND END OF LIMITS BLOCK BEFORE START");
        l_read=false; if(nl!=dimout) ERROR("WRONG NUMBER OF LIMITS FOR BESPOKE CVS - MUST BE "<<dimout);  
     }

     // These read the blocks of embedded points and basins
     else if( bas_read && flag==2 ){ 
       if( (input>>highP).fail() ){ ERROR("FAILED TO READ IN HIGH_D LANDMARK SAFELY"); } 
       bstress.setHighPoint(nbas, highP); nbas++;
     }
     else if( p_read && flag==2){ 
        if( (input>>lowP).fail() ){ ERROR("FAILED TO READ IN LOW_D LANDMARK SAFELY"); } 
        bstress.setLowPoint(np, lowP); np++;
     }
     else if( w_read && flag==2){
        if( (input>>w).fail() ){ ERROR("FAILED TO READ IN WEIGHT OF LANDMARK SAFELY"); }
        bstress.setWeight(nw, w); total_weight+=w; nw++;
     }
     // This reads the limts for the integration
     else if( l_read ){ 
        if( (input>>limlow[nl]>>limhigh[nl]).fail() ){ ERROR("FAILED TO READ IN LIMITS FOR INTEGRATION SAFELY"); } nl++; 
     }

     // These find the beginings of the blocks of embedded points and basins
     else if( line.find("HIGH_D>")!=std::string::npos && flag==2 ){ bas_read=true; }
     else if( line.find("HIGH_D>")!=std::string::npos ){ ERROR("FOUND HIGH_D BLOCK IN BESPOKE DEFS BEFORE FINDING NUMBER OF LANDMARKS"); }

     else if( line.find("LOW_D>")!=std::string::npos && flag==2 ){ p_read=true; }
     else if( line.find("LOW_D>")!=std::string::npos ){ ERROR("FOUND LOW_D BLOCK IN BESPOKE DEFS BEFORE FINDING NUMBER OF LANDMARKS"); }

     else if( line.find("WEIGHTS>")!=std::string::npos && flag==2 ){ w_read=true; read_weights=true; }
     else if( line.find("WEIGHTS>")!=std::string::npos ){ ERROR("FOUND WEIGHTS BLOCK IN BESPOKE DEFS BEFORE FINDING NUMBER OF LANDMARKS"); }

     else if( line.find("LIMITS>")!=std::string::npos ){ l_read=true; }

     else if( w_read && p_read ) {ERROR("INPUT FOR BESPOKE CVS IS NONSENSICAL"); }
     else if( w_read && bas_read ) {ERROR("INPUT FOR BESPOKE CVS IS NONSENSICAL"); }
     else if( w_read && l_read ) {ERROR("INPUT FOR BESPOKE CVS IS NONSENSICAL"); }
     else if( bas_read && p_read ){ ERROR("INPUT FOR BESPOKE CVS IS NONSENSICAL"); }
     else if( bas_read && l_read ){ ERROR("INPUT FOR BESPOKE CVS IS NONSENSICAL"); }
     else if( p_read && l_read ){ ERROR("INPUT FOR BESPOKE CVS IS NONSENSICAL"); }
  } 
  ifile.close(); //std::pout<<"\n";

  if( flag==0 ) ERROR("NEVER FOUND NBASINS IN BESPOKE CV INPUT");
  if( nbas!=nland ) ERROR("FOUND NO HIGH DIMENSIONALITY INFORMATION IN BESPOKE CV INPUT");
  if( np!=nland ) ERROR("FOUND NO LOW DIMENSIONALITY INFORMATION IN BESPOKE CV INPUT");

  if( !read_weights ){
     WARNING("NO WEIGHT DATA IN BESPOKE FILE SO SETTING ALL WEIGHTS EQUAL");
     w = 1.0 / (double) nland; for(index i=0;i<nland;++i){ bstress.setWeight(i, w ); }
  } else { bstress.normalizeWeights( total_weight ); }

  // for(index i=0;i<dimout;++i){ std::pout<<"|--LIMITS FOR BESPOKE CV "<<i+1<<" "<<limlow[i]<<" "<<limhigh[i]<<std::endl; }
  for(index i=0;i<dimout;++i){ fprintf(fplog,"|--LIMITS FOR BESPOKE CV %d : %f %f \n", i+1, limlow[i], limhigh[i] ); } 
  fprintf(fplog,"\n");
 
  //std::pout<<"\n";
}

void bespokeCVObj::embedPoint( const rVector& colvars, rVector& embededP, rMatrix& derivatives ) {
#ifdef DEBUG 
  if(colvars.size()!=ncv) ERROR("EXPECTED DIFFERENT NUMBER OF COLVARS");
  if(embededP.size()!=dimout) ERROR("MISMATCH FOR EXPECTED NUMBER OF BESPOKE CVS");
  if(derivatives.rows()!=dimout || derivatives.cols()!=ncv) ERROR("DERIVATIVE MATRIX IS WRONG SIZE");
#endif
  rVector pos( dimout ), dx( dimout ), dX( ncv ), lb( dimout ), ub( dimout ); 
  std::valarray<index> box( dimout ); rMatrix dXx( ncv, dimout ); 

  // First calculate the distace from all high dimensionality points
  error=bstress.calcDistances( colvars );

  // compute the value of the stress on a grid 
  real minchi=-1;
  for(index i=0;i<nspline[0];++i){
     pos[0] = gx[i] = limlow[0] + i*( ( limhigh[0] - limlow[0] ) / ( nspline[0] - 1.0 ) ); 
     for(index j=0;j<nspline[1];++j){
         pos[1] = gy[j] = limlow[1] + j*( ( limhigh[1] - limlow[1] ) / ( nspline[1] - 1.0 ) );
         bstress.computeStress( pos, sgridU(i,j), dx, dX, dXx  );    
         if(minchi<0 || minchi>sgridU(i,j)) minchi=sgridU(i,j);
         for(index k=0;k<dimout;++k){ sgridDU(i, j, k) = dx[k]; }
         for(index h=0;h<ncv;++h){ 
            dsgridU[h](i,j)=dX[h];
            for(index k=0;k<dimout;++k){ dsgridDU[h](i, j, k) = dXx(h,k); }
         }
     }
  }

  // Setup our interpolator
  interpolator.set_table2d( gx, gy, sgridU, sgridDU );
  // Setup derivative interpolators
  for(index h=0;h<ncv;++h){ dinterpolator[h].set_table2d( gx, gy, dsgridU[h], dsgridDU[h] ); }

  // This is for testing - outputs the interpolated stress function
  //std::ofstream interpol("interpolated");

  // This now takes averages over the grid of interpolated points
  real f, echi, dechi, tr, tx, ty, tt; tr=tx=ty=tt=0; box[0]=0;
  rVector df(ncv), dtr(ncv), dtx(ncv), dty(ncv), dtt(ncv);
  df=0; dtr=0; dtx=0; dty=0; dtt=0; 

#if defined(MPI)
  int rank, size, ia, ib;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  ia=rank*(ngrid[0]/size); ib=(rank+1==size?ngrid[0]:(rank+1)*(ngrid[0]/size));
 //  std::cerr<<"Running on "<<rank<<" of "<<size<<": from "<<ia<<" to "<<ib<<" tt being " <<tt<<"\n";
  box[0]=search1(limlow[0] + ia*( ( limhigh[0] - limlow[0] ) / ( ngrid[0] - 1.0 ) ),gx,0);
  for(index i=ia;i<ib;++i){
#else
  box[0]=0;
  for(index i=0;i<ngrid[0];++i){
#endif
     //MPI_Barrier(MPI_COMM_WORLD); if (i<ia || i>=ib) continue; 
     pos[0] = limlow[0] + i*( ( limhigh[0] - limlow[0] ) / ( ngrid[0] - 1.0 ) );
     if( pos[0]>gx[box[0]+1] ){ box[0]=box[0]+1; }

     lb[0]=gx[box[0]]; ub[0]=gx[box[0]+1]; box[1]=0;
     for(index j=0;j<ngrid[1];++j){
         pos[1] = limlow[1] + j*( ( limhigh[1] - limlow[1] ) / ( ngrid[1] - 1.0 ) );
         if( pos[1]>gy[box[1]+1] ){ box[1]=box[1]+1; }
         lb[1]=gy[box[1]]; 
         ub[1]=gy[box[1]+1]; 

         interpolator.get_fdf( box, lb, ub , pos, f, dx ); 
         for(index h=0;h<ncv;++h){ dinterpolator[h].get_fdf( box, lb, ub , pos, df[h], dx ); }       

         //interpol<<pos[0]<<" "<<pos[1]<<" "<<f<<std::endl;

         // Accumulate all our averages
         echi=exp(-(f-minchi)/kt); tt+=echi; 
         tx+=pos[0]*echi; ty+=pos[1]*echi; 
         tr+=sqrt(pos[0]*pos[0]+pos[1]*pos[1])*echi;

         // Accumulate averages for derivatives
         for(index h=0;h<ncv;++h){
            dechi= ( -1.0 / kt ) * echi * df[h];
            dtt[h]+=dechi;
            dtx[h]+=pos[0]*dechi; dty[h]+=pos[1]*dechi;
            dtr[h]+=sqrt(pos[0]*pos[0]+pos[1]*pos[1])*dechi;
         } 
         //std::cerr<<rank<<":"<<i<<","<<j<<"("<<pos[0]<<","<<pos[1]<< ") " <<f<<","<<echi<<","<<tt<<"\n";
     }
     //interpol<<std::endl;
  }

// Gather the information on derivatives and sums together from all the nodes
#if defined(MPI)
  rVector buffer(ncv);
  MPI_Datatype mpireal=(sizeof(real)==4?MPI_FLOAT:MPI_DOUBLE);
  //std::cerr<<"Partial on "<<rank<<" of "<<size<<": tt "<<tt<<" dtr "<<dtr[ncv-1]<<"\n";
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce((void*) &tt, (void*) &buffer[0], 1, mpireal, MPI_SUM, MPI_COMM_WORLD); tt=buffer[0]; MPI_Barrier(MPI_COMM_WORLD); 
  MPI_Allreduce((void*) &tr, (void*) &buffer[0], 1, mpireal, MPI_SUM, MPI_COMM_WORLD); tr=buffer[0]; MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce((void*) &tx, (void*) &buffer[0], 1, mpireal, MPI_SUM, MPI_COMM_WORLD); tx=buffer[0]; MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce((void*) &ty, (void*) &buffer[0], 1, mpireal, MPI_SUM, MPI_COMM_WORLD); ty=buffer[0]; MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce((void*) &dtt[0], (void*) &buffer[0], ncv, mpireal, MPI_SUM, MPI_COMM_WORLD); dtt=buffer; MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce((void*) &dtr[0], (void*) &buffer[0], ncv, mpireal, MPI_SUM, MPI_COMM_WORLD); dtr=buffer; MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce((void*) &dtx[0], (void*) &buffer[0], ncv, mpireal, MPI_SUM, MPI_COMM_WORLD); dtx=buffer; MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce((void*) &dty[0], (void*) &buffer[0], ncv, mpireal, MPI_SUM, MPI_COMM_WORLD); dty=buffer; MPI_Barrier(MPI_COMM_WORLD);

  //std::cerr<<"Total on "<<rank<<" of "<<size<<": tt "<<tt<<" dtr "<<dtr[ncv-1]<<"\n";

#endif
  //interpol.close();

  // First step in completion of calculation of derivatives
  for(index h=0;h<ncv;++h) { 
     dtr[h] = ( dtr[h]*tt - tr*dtt[h] ); 
     dtx[h] = ( dtx[h]*tt - tx*dtt[h] ); 
     dty[h] = ( dty[h]*tt - ty*dtt[h] ); 
  }
  // Complete the calculation of averages
  tr*=1.0/tt; tx*=1.0/tt; ty*=1.0/tt;
  // Complete the calculation of the derivatives
  dtr*=1.0/(tt*tt); dtx*=1.0/(tt*tt); dty*=1.0/(tt*tt);

  
  // Calculate the cvs
  embededP[0] = tr*cos( atan2(ty,tx) ); embededP[1] = tr*sin( atan2(ty,tx) );
  //std::cout<<"BESPOKE DATA : r "<<tr<<" theta "<<atan2(ty,tx)<<" (x,y) "<<embededP<<std::endl;
  // And the derivatives
  for(index h=0;h<ncv;++h){
      derivatives(0,h)=cos( atan2(ty,tx) ) * dtr[h] +
        tr* ( ty)/(tx*tx+ty*ty) * sin( atan2(ty,tx) ) * dtx[h] +
        tr* (-tx)/(tx*tx+ty*ty) * sin( atan2(ty,tx) ) * dty[h] ;
      derivatives(1,h)=sin( atan2(ty,tx) ) * dtr[h] +
        tr* (-ty)/(tx*tx+ty*ty) * cos( atan2(ty,tx) ) * dtx[h] +
        tr* ( tx)/(tx*tx+ty*ty) * cos( atan2(ty,tx) ) * dty[h] ;
  }
  //std::cout<<"BESPOKE CV DATA: <R> "<<tr<<" <theta> "<<atan2(ty,tx)<<" p "<<embededP[0]<<" "<<embededP[1]<<std::endl;
  return;
}

void bespokeStress::set_function( const std::string& highlow, const std::string& fname, real& a, real& b, real& sigma, FILE* fplog ) {
  // These are sanity checks
  if(fname.find("distance")!=std::string::npos) {
      sigma=-1.0; a=b=0; 
      fprintf(fplog,"|--USING DISTANCE FOR STRESS IN %s SPACE \n", highlow.c_str() );
      // std::pout<<"USING DISTANCE FOR STRESS IN "<<highlow<<std::endl;
  } else if(fname.find("compress")!=std::string::npos) {
      if(sigma<=0.0) ERROR("INVALID CHOICE FOR SIGMA WITH COMPRESS FUNCTION - "<<highlow);
      a=b=1; 
      fprintf(fplog,"|--USING COMPRESS FUNCTION WITH SIGMA EQUALS %f FOR STRESS IN %s SPACE \n", sigma, highlow.c_str() );
      // std::pout<<"USING COMPRESS FUNCTION WITH SIGMA EQUALS "<<sigma<<" FOR STRESS IN "<<highlow<<std::endl;
  } else if(fname.find("sigmoid")!=std::string::npos) {
      if(sigma<=0.0) ERROR("INVALID CHOICE FOR SIGMA WITH SIGMOID FUNCTION - "<<highlow);
      a=b=2; 
      fprintf(fplog,"|--USING SIGMOID FUNCTION WITH SIGMA EQUALS %f FOR STRESS IN %s SPACE \n", sigma, highlow.c_str() );
      // std::pout<<"USING SIGMOID FUNCTION WITH SIGMA EQUALS "<<sigma<<" FOR STRESS IN "<<highlow<<std::endl;
  } else if(fname.find("general")!=std::string::npos) {
      if(sigma<=0.0) ERROR("INVALID CHOICE FOR SIGMA WITH GENERAL FUNCTION - "<<highlow);
      if(a<=0 || b<=0 ) ERROR("INVALID CHOICE FOR POWERS WITH GENERAL FUNCTION - "<<highlow);
      fprintf(fplog,"|--USING GENERAL SIGMIOD FUNCTION WITH SIGMA EQUALS %f AND POWERS EQUAL TO %f AND %f FOR STRESS IN %s SPACE \n", sigma, a, b, highlow.c_str() ); 
      // std::pout<<"USING GENERAL SIGMIOD FUNCTION WITH SIGMA EQUALS "<<sigma<<" AND POWERS EQUAL TO "<<a<<" AND "<<b<<" FOR STRESS IN "<<highlow<<std::endl;
  } else {
     ERROR("INVALID FUNCTION TYPE FOR "<<highlow);
  }

  // Now set the function
  if(highlow.find("LOW_D")!=std::string::npos) { f_lowd.set_mode( a, b, sigma ); } 
  else if(highlow.find("HIGH_D")!=std::string::npos) { f_highd.set_mode( a, b, sigma ); } 
  else { ERROR("BIZZARE ERROR IN INPUT - FUNCITON IS NOT FOR HIGH OR LOW D SPACE"); }
}

void bespokeStress::computeStress( const rVector& pos, real& stress, rVector& dxstress, rVector& dXstress, rMatrix& dXxstress ) const {
#ifdef DEBUG
  if( pos.size()!=dimlow ) ERROR("SIZE MISMATCH FOR POS");
  if( dxstress.size()!=dimlow ) ERROR("SIZE MISMATCH FOR dxstress");
  if( dXstress.size()!=dimhigh) ERROR("SIZE MISMATCH FOR dXstress");
  if( dXxstress.rows()!=dimhigh || dXxstress.cols()!=dimlow ) ERROR("SIZE MISMATCH FOR dXxstress");
#endif

  // Note to self small x is for low dimensionality big X is for high dimensionality

  rVector pdstress( dimlow ), pdXstress( dimhigh ); 
  real w, dd, ddj, FonD, fond, distf, dist, fdLD, dfdLD, t1; 
  index i,j,k; 
  real *spdxx=&dXxstress(0,0), *pdxx, ds0, ds1;
  stress=0.0; dxstress=0; dXstress=0; dXxstress=0;
  
#if defined(MPI)
  int rank, size, ia, ib;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  ia=rank*(npoints/size); ib=(rank+1==size?npoints:(rank+1)*(npoints/size));
  for (i=ia; i<ib; ++i) {
#else
  for(i=0;i<npoints;++i) {
#endif
     // Calculate the distance from the point
     dist=0.0; 
     for(j=0;j<dimlow;++j) { pdstress[j]=pos[j] - lowP(i,j); dist+=pdstress[j]*pdstress[j]; }
     if (dist<=0.0) continue;  
     // Compute the value of the switching function in the low dimensionality space
     dist=sqrt(dist); f_lowd.fdf( dist, fdLD, dfdLD ); 
     // The derivative of the switching function in the low d space
     fond=dfdLD / dist; 
     // The derivative of the swtiching function in the high d space
     FonD=sHD[i].df;
     // The weight of this point in the stress
     w=weights[i];
     // The difference between the values of the switching function in the low and high d space ( this is used to accumulate stress )
     distf=( sHD[i].f - fdLD );  
     // The cross derivative (a second derivative?)
     dd = -2.0 * w * FonD * fond;
     // Accumulate the total stress
     stress+=w*distf*distf;  
     // The derivative stress with respect to the high D distance
     distf+=distf; t1=w*distf*FonD;

     if(FonD!=0.0){
        diffHD.getRow( i, pdXstress );
 
        // These are the cross terms
        if(fond==0.0){ 
           for(j=0;j<dimhigh;++j) dXstress[j]+=pdXstress[j]*t1; 
        } else{
          pdxx=spdxx; 
          ds0=pdstress[0]; ds1=pdstress[1];
          for(j=0;j<dimhigh;++j){
             dXstress[j]+=pdXstress[j]*t1;
             ddj=pdXstress[j]*dd; 
             // User pointers for fast processing (N.B. This relies on dimlow being 2)
             *pdxx += ds0*ddj; pdxx++;
             *pdxx += ds1*ddj; pdxx++;
          }
        }
     }
     if (fond==0.0) continue;
     t1=-w*distf*fond; 
     for(k=0;k<dimlow;++k) dxstress[k]+=pdstress[k]*t1; 
  }
#if defined(MPI)
  std::valarray<real> buffer(dimhigh*dimlow);
  MPI_Datatype mpireal=(sizeof(real)==4?MPI_FLOAT:MPI_DOUBLE);
  MPI_Allreduce((void*) &stress, (void*) &buffer[0], 1, mpireal, MPI_SUM, MPI_COMM_WORLD);
  stress=buffer[0];  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce((void*) &dxstress[0], (void*) &buffer[0], dimlow, mpireal, MPI_SUM, MPI_COMM_WORLD);
  for(i=0;i<dimlow;++i) dxstress[i]=buffer[i];  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce((void*) &dXstress[0], (void*) &buffer[0], dimhigh, mpireal, MPI_SUM, MPI_COMM_WORLD);
  for(j=0;j<dimhigh;++j) dXstress[j]=buffer[j];  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce((void*) &dXxstress[0], (void*) &buffer[0], dimhigh*dimlow, mpireal, MPI_SUM, MPI_COMM_WORLD);
  k=0; for(j=0;j<dimhigh;++j) for(i=0;i<dimlow;++i) dXxstress(j,i)=buffer[k++];  MPI_Barrier(MPI_COMM_WORLD);
#endif

//  Weights now automatically divide by the total weight 
//  stress /= real(npoints);          // Complete the function value
//  dxstress /= real(npoints);        // Complete the derivative (low d)
//  dXstress /= real(npoints);        // Complete the derivative (high d)
//  dXxstress /= real(npoints);       // Complete the cross derivative
}

real bespokeStress::calcDistances( const rVector& colvars ) {
#ifdef DEBUG
  if( colvars.size()!=dimhigh ) ERROR("SIZE MISMATCH FOR COLVARS");
  if( diffHD.rows()!=npoints || diffHD.cols()!=dimhigh ) ERROR("SIZE MISMATCH FOR diffHD "<<diffHD.size()<<" "<<npoints);
#endif

  real mdist=-1.0; rVector dHD(npoints); 
  for(index i=0;i<npoints;++i){ 
    // Calculate the distance between this center and the point
    dHD[i]=0.0;
    for(index j=0;j<dimhigh;++j){ 
      diffHD(i,j)=pbc( ( colvars[j] - highP(i,j) ) , periods[j] );  dHD[i]+=diffHD(i,j)*diffHD(i,j);
    }
    dHD[i]=sqrt(dHD[i]); 
    if ( dHD[i]<mdist || mdist<0 )  mdist=dHD[i];
    // Compute the value of the high dimensionality switching function
    f_highd.fdf( dHD[i], sHD[i].f, sHD[i].df ); sHD[i].df/=dHD[i];  
  }

  return mdist;
}

void InterpolateBicubic::set_table2d( const rVector& gx, const rVector& gy, const rMatrix& fx, const rTensor& dfx ){
#ifdef DEBUG
  if( gx.size()!=ngrid[0] ) ERROR("SIZE MISMATCH FOR X GRID");
  if( gy.size()!=ngrid[1] ) ERROR("SIZE MISMATCH FOR Y GRID");
  if( fx.rows()!=ngrid[0] || fx.cols()!=ngrid[1] ) ERROR("SIZE MISMATCH FOR FUNCTION VALUES");
  if( dfx.rank()!=3 ) ERROR("DERIVATIVES SHOULD BE A THIRD RANK TENSOR");
  if( dfx.dim(0)!=ngrid[0] || dfx.dim(1)!=ngrid[1] || dfx.dim(2)!=2 ) ERROR("SIZE MISMATCH FOR DERIVATIVES");
#endif

  // This computes the (numerical) second derivatives
  dcross*=0.0; //on the boundaries, assume no curvature
  for (index i=1; i<ngrid[0]-1; ++i){ 
      for (index j=1; j<ngrid[1]-1; ++j){
          dcross(i,j) = ( fx( i+1, j+1 ) + fx( i-1, j-1 ) - fx( i+1, j-1) - fx( i-1, j+1 ) ) /
                        ( ( gx[i+1] - gx[i-1] ) * ( gy[j+1] - gy[j-1] ) );
      }
  }

  rMatrix tc(4,4); real d1, d2;
  rVector y(4), dy1(4), dy2(4), d2y12(4);

#if defined(MPI)
  int rank, size, ik=0, tsz=(ngrid[0]-1)*(ngrid[1]-1)*16; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  double *pclist=&clist(0,0,0,0);
  for (index i=0; i<tsz; ++i) pclist[i]=0.0;
#endif

  
  for (index i=0;i<ngrid[0]-1;++i){
      d1 = gx[i+1] - gx[i];                                     
      for (index j=0; j<ngrid[1]-1;++j){
         d2 = gy[j+1] - gy[j];
#if defined(MPI)
         if ((ik++)%size!=rank) continue; 
#endif         
         y[0] = fx( i, j ); y[1] = fx( i+1 ,j ); y[2] = fx( i+1, j+1 ); y[3] = fx( i, j+1 );
         dy1[0] = dfx( i, j, 0 ); dy1[1] = dfx( i+1, j, 0 ); dy1[2] = dfx( i+1, j+1, 0 ); dy1[3] = dfx( i, j+1, 0 );
         dy2[0] = dfx( i, j, 1 ); dy2[1] = dfx( i+1, j, 1 ); dy2[2] = dfx( i+1, j+1, 1 ); dy2[3] = dfx( i, j+1, 1 );
         d2y12[0] = dcross( i, j ); d2y12[1] = dcross( i+1, j ); d2y12[2] = dcross( i+1, j+1 ); d2y12[3] = dcross( i, j+1 );

         IBicCoeff( y, dy1, dy2, d2y12, d1, d2, tc);    

         for(index k=0; k<4; ++k){ for(index n=0; n<4; ++n){ clist(i,j,k,n)=tc(k,n); } }     
      }
  }
#if defined(MPI)
  rVector buffer(tsz);
  MPI_Datatype mpireal=(sizeof(real)==4?MPI_FLOAT:MPI_DOUBLE);
  MPI_Allreduce((void*) pclist, (void*) &buffer[0], tsz, mpireal, MPI_SUM, MPI_COMM_WORLD); 
  for (index i=0; i<tsz; ++i) pclist[i]=buffer[i];
#endif
}

void InterpolateBicubic::IBicCoeff( const rVector& y, const rVector& dy1, const rVector& dy2, const rVector& d2y12, const real& d1, const real& d2, rMatrix& c ){
    real xx, d1d2=d1*d2;
    rVector cl(16), x(16); 

    for(index i=0;i<4;i++){ x[i] = y[i]; x[i+4] = dy1[i]*d1; x[i+8] = dy2[i]*d2; x[i+12] = d2y12[i]*d1d2; }
    for(index i=0;i<16;i++){ xx=0.0; for(index k=0;k<16;k++){ xx += wt[i][k]*x[k]; } cl[i]=xx; }
    index l=0; for(index i=0;i<4;i++){ for(index j=0;j<4;j++){ c(i,j)=cl[l++]; } }
}

long search1(const real& x, const rVector& xlist, long jold) {
    long jl=jold,ju,jm,inc=1,n=xlist.size();
    if (n<2) ERROR("Invalid size of interpolation arrays");
    if (jl<0 || jl>=n) { jl=0; ju=n-1; }
    else {
        if(x>=xlist[jl]){
            while(true){
                ju=jl+inc; 
                if(ju>=n-1){ju=n-1; break;}
                else if(x<xlist[ju]) break;
                jl=ju; inc+=inc;
            }
        }
        else{
            ju=jl; 
            while(true){
                jl=jl-inc;
                if(jl<=0){jl=0; break;}
                else if(x>=xlist[jl]) break;
                ju=jl; inc+=inc;
            }
        }
    }
    while(ju-jl>1){
        jm=(ju+jl)/2;
        if (x>xlist[jm]) jl=jm; else ju=jm;
    }
    jold=jl;
    return jl;
}

void InterpolateBicubic::get_fdf( const std::valarray<index>& box, const rVector& lb, const rVector& ub, const rVector& pos, real& f, rVector& df ){

    if( ub[0] == lb[0] || ub[1] == lb[1] ) ERROR("Wrong boundaries detected in interpolation");

    real d1 = ub[0] - lb[0], d2 = ub[1] - lb[1];
    real t = (pos[0] - lb[0]) / d1, u = (pos[1] - lb[1]) / d2;

    //faster access by pointer arithmetic (dirty dirty dirty)
    real *cbase=&clist(box[0],box[1],3,3), *c3, *c2, *c1, *c0; 

    real ansy=0.0;
    for (int i=3; i>=0; i--) {    // Note to self - this has to be an int as index cannot be less than zero
        c3=cbase; c2=c3-1; c1=c2-1; c0=c1-1; cbase=c0-1;  
        ansy= t*ansy + ( ( (*c3)*u + (*c2) )*u + (*c1) )*u + (*c0);
    }
    f=ansy;
} 

void bespokeStress::output() const {

  std::cout<<"THE HIGH DIMENSIONALITY POINTS"<<std::endl;
  for(index i=0;i<npoints;++i){ 
     std::cout<<"POINT "<<i<<" ";
     for(index j=0;j<dimhigh;++j){ std::cout<<highP(i,j)<<" "; }
     std::cout<<std::endl;
  }
  std::cout<<"THE LOW DIMENSIONALITY POINTS"<<std::endl;
  for(index i=0;i<npoints;++i){
     std::cout<<"POINT "<<i<<" ";
     for(index j=0;j<dimlow;++j){ std::cout<<lowP(i,j)<<" "; }
     std::cout<<std::endl;
  }
}

}  // END OF NAMESPACE
#endif
#endif
