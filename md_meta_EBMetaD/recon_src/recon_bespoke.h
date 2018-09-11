#ifndef __rp_bespoke
#define __rp_bespoke

#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <numeric>
#include <valarray>
#include <iostream>
#include <string>
#include <fstream>

#include "recon_types.h"
#include "recon_utils.h"
// #include "recon_basins.h"

namespace rp {

class rTensor{
private:
  index order, sz;
  std::valarray<index> dims; 
  std::valarray<real> data;
public:
  rTensor(const index ord=0, ...) {
     order=ord; dims.resize(order);
     if(order==0) return;

     va_list ap; va_start(ap, ord); sz=1;
     for(index i=0;i<order;++i){ dims[i]=va_arg(ap, index); sz*=dims[i]; }
     va_end(ap); data.resize(sz); 
  }
  rTensor(const rTensor& t) : order(t.order), sz(t.sz), dims(t.dims), data(t.data) {}
  rTensor& operator=(const rTensor& t){
     if (&t==this) return *this;
#ifdef DEBUG
     if( t.order!=order ) ERROR("Tensor order mismatch");
     for(index i=0;i<order;++i){ 
         if( t.dims[i]!=dims[i] ) ERROR("Tensor size mismatch");
     }
#endif
     data=t.data;
     return *this;
  }
  index rank() const { return order; } index dim(const index i) const { return dims[i]; }
  void resize( const index ord, ...){ 
     order=ord; dims.resize(order);
     
     va_list ap; va_start(ap, ord); sz=1;
     for(index i=0;i<order;++i){ dims[i]=va_arg(ap, index); sz*=dims[i]; }
     va_end(ap); data.resize(sz);
  }
  // variadic access
  inline real& operator () (const index first, ...){
     index pos=first;
     if (order<=1) return data[pos];

     va_list ap; va_start(ap, first);
     for (index i=1; i<order; ++i) pos=pos*dims[i]+va_arg(ap, index);
     va_end(ap);
     return data[pos];
  }
  inline const real& operator()(const index first, ...) const {
     index pos=first;
     if (order<=1) return data[pos];

     va_list ap; va_start(ap, first);
     for (index i=1; i<order; ++i) pos=pos*dims[i]+va_arg(ap, index);
     va_end(ap);
     return data[pos];
  }
};

class InterpolateBicubic{
private:
  index ndim;
  std::valarray<index> ngrid;
  std::valarray< std::valarray<int> > wt; 
  rMatrix dcross;
  rTensor clist;  
public:
  InterpolateBicubic(const index& nd=2) : ndim(nd), ngrid(ndim), wt(0), dcross(0,0), clist(ndim*ndim,0,0,0,0) {} 
  InterpolateBicubic(const InterpolateBicubic& b) : ndim(b.ndim), ngrid(b.ngrid), wt(b.wt), dcross(b.dcross), clist(b.clist) {}
  void resize(const index& nd){
     ndim=nd; ngrid.resize(ndim); dcross.resize(0,0); clist.resize(ndim*ndim,0,0,0,0);
  }
  void setupGrid(const index first, ...){
     static int wt_d[16*16]=
     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
     2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
     0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
     -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
     9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
     -6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
     2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
     -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
     4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1};

     // This is to set up the coefficient matrix
     index l=0; wt.resize(16); 
     for (index i=0;i<16;i++){ 
          wt[i].resize(16);
          for (index j=0;j<16;j++){ wt[i][j]=wt_d[l++]; } 
     }

     ngrid[0]=first;
     va_list ap; va_start(ap, first); 
     for(index i=1;i<ndim;++i){ ngrid[i]=va_arg(ap, index); }
     va_end(ap);
     
     if(ndim==2){
        dcross.resize( ngrid[0], ngrid[1] ); clist.resize( 4, ngrid[0]-1, ngrid[1]-1, 4, 4 );
     } else {
        ERROR("INTERPOLATOR ONLY WORKS WITH 2 BESPOKE COLLECTIVE COORDINATES");
     }
  }
  void set_table2d( const rVector& gx, const rVector& gy, const rMatrix& fx, const rTensor& dfx );
  void IBicCoeff( const rVector& y, const rVector& dy1, const rVector& dy2, const rVector& d2y12, const real& d1, const real& d2, rMatrix& c );
  void get_fdf( const std::valarray<index>& box, const rVector& lb, const rVector& ub, const rVector& pos, real& f, rVector& df );
};    

typedef struct { real f, df; } bespokeStorage; 

class bfunction {
  typedef void (bfunction::*fpoint)( const real&, real&, real& ) const;
private:
  real rsigma, dsigma_sq; // dsigma_sq = 2*sigma*sigma, rsigma = 1/sigma
  real a, b, c, d; // c = ( 2^(a/b) - 1 ), d = - b / a 
  fpoint pfdf;

  void dist( const real& x, real& f, real& df ) const ; 
  void sigmoid( const real& x, real& f, real& df ) const ;
  void compress( const real& x, real& f, real& df ) const ;
  void general( const real& x, real& f, real& df ) const ;  
public:
  bfunction( index aa=0, index bb=0, real sig=-1.0){ set_mode( aa, bb, sig); }
  bfunction( const bfunction& nf);
  bfunction& operator=(const bfunction& nf);

  void set_mode( const real& aa , const real& bb, const real& sig);
  inline void fdf ( const real& x, real& f, real& df ) const { (this->*pfdf)(x, f, df); }
};

// This class defines a set of mappings from a high to a low dimensional space
class bespokeStress{
private:
  index npoints;
  index dimhigh, dimlow;
  rVector periods, weights;
  rMatrix highP, lowP, diffHD;
  std::valarray<bespokeStorage> sHD;    
  bfunction f_lowd, f_highd;
public:
  bespokeStress( index n=0, index d=0 ) : npoints(0), dimhigh(n), dimlow(d), periods(n), weights(0), highP(0,0), lowP(0,0), diffHD(0,0), sHD(0) {}
  bespokeStress( const bespokeStress& m ) : npoints(m.npoints), dimhigh(m.dimhigh), dimlow(m.dimlow), periods(m.periods), weights(m.weights), highP(m.highP), lowP(m.lowP), diffHD(m.diffHD), sHD(m.sHD) {} 
  void resize (const index& np ,const index& n, const index& d){ 
     npoints=np; dimhigh=n; dimlow=d; periods.resize(dimhigh);  
     weights.resize(npoints); highP.resize(npoints,dimhigh); lowP.resize(npoints,dimlow); 
     diffHD.resize(npoints, dimhigh);
     sHD.resize(npoints);
  }
  void set_function( const std::string& highlow, const std::string& fname, real& a, real& b, real& sigma, FILE *fplog );
  void setPeriods( const rVector& p ){ periods=p; }
  void setHighPoint( const index& i, const rVector& p ){ highP.setRow(i,p); }
  void setLowPoint( const index& i, const rVector& p ){ lowP.setRow(i,p); }
  void setWeight( const index& i, const real& w ){ weights[i]=w; }
  void normalizeWeights( const real& total ){ weights/=total; }
  real calcDistances( const rVector& colvars );
  void computeStress( const rVector& pos, real& stress, rVector& dxstress, rVector& dXstress, rMatrix& dXxstress ) const;

  // This is for debugging readin
  void output () const;
};

// Class that will be used to calculate the CVs
class bespokeCVObj{
private:
  index ncv, dimout;
  real kt;
  std::vector<index> ngrid, nspline;
  rVector limlow, limhigh;
  index nland;
  bespokeStress bstress;
  rVector gx, gy;
  // The interpolator for the function value
  InterpolateBicubic interpolator;
  rMatrix sgridU;
  rTensor sgridDU;
  // The interpolators for the derivatives 
  std::valarray<InterpolateBicubic> dinterpolator;
  std::valarray<rMatrix> dsgridU;
  std::valarray<rTensor> dsgridDU;
  real error;
public:
  bespokeCVObj( index ncolvars=0, index d=0 ) : ncv(ncolvars), dimout(d), kt(0), ngrid(dimout), nspline(dimout), limlow(dimout), limhigh(dimout), nland(0), bstress(0,0), gx(0), gy(0), interpolator(), sgridU(0,0), sgridDU(dimout+1,0,0,0), dinterpolator(ncv), dsgridU(ncv), dsgridDU(ncv), error(0) {} 
  void setup( const std::string efilename, const rp::rVector& p, const real& smoo, const std::vector<rp::index>& nsp, const std::vector<rp::index>& ngr, FILE *fplog ); 
  void embedPoint( const rVector& colvars, rVector& embededP, rMatrix& derivatives ) ;
//  real reconstructError( const rVector& colvars ) const { return error; };
  real reconstructError() const { return error; };
};

};    // END OF NAMESPACE

#endif
