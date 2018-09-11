#include "metadyn.h"

void PREFIX intracontact_restraint(int i_c, struct mtd_data_s *mtd_data) 
{
  int npe,rank; 
  int firstAtom, secondAtom, ix, i, j, l, indice;
  int nn = colvar.nn[i_c], mm = colvar.mm[i_c];
  real *tmpx,*tmpy,*tmpz; 
  rvec rij;
  real mod_rij, func, dfunc, rdist, r6dist, r12dist, num, iden, cont, threshold;

  func = 0.;
  indice = 0;
  threshold=pow(0.00001,1./(nn-mm));
  snew(tmpx,colvar.natoms[i_c]); 
  snew(tmpy,colvar.natoms[i_c]); 
  snew(tmpz,colvar.natoms[i_c]); 
  for(i=0;i<colvar.natoms[i_c];i++) for(ix=0;ix<3;ix++) colvar.myder[i_c][i][ix] = 0.;
  if(logical.parallel_hills){
    npe=plumed_comm_size(mtd_data);
    rank=plumed_comm_rank(mtd_data);
    indice=rank*(colvar.natoms[i_c]-colvar.type[i_c]);
    for(i=0;i<rank;i++) indice-=i;
  }else{
    npe=1;
    rank=0;
  };

  for(i=rank;i<colvar.natoms[i_c];i+=npe){
    firstAtom = colvar.cvatoms[i_c][i];
    for(j=i+colvar.type[i_c];j<colvar.natoms[i_c];j++) {
      secondAtom = colvar.cvatoms[i_c][j];
      minimal_image(mtd_data->pos[firstAtom], mtd_data->pos[secondAtom], &mod_rij, rij);
      rdist = mod_rij/colvar.r_0[i_c];
      if(rdist>threshold) {
        cont = 0.;
        dfunc = 0.;
      } else if(rdist>0.999999 && rdist<1.000001) {
        cont = nn/mm;
        dfunc = 0.5*nn*(nn-mm)/mm;
      } else if(rdist<=0.){
        cont = 1.;
        dfunc = 0.;
      } else {
        r6dist = pow(rdist, nn-1.);
        r12dist = pow(rdist, mm-1.);
        num = 1.-(r6dist*rdist);
        iden = 1./(1.-(r12dist*rdist));
        cont = num*iden;
        dfunc = ((-nn*r6dist*iden)+(cont*(iden*mm)*r12dist))/(mod_rij*colvar.r_0[i_c]);
      }
      func += (colvar.map0[i_c][indice]-cont)*(colvar.map0[i_c][indice]-cont);
      dfunc *= (colvar.map0[i_c][indice]-cont);
      for(ix=0;ix<3;ix++) {
        colvar.myder[i_c][i][ix] += -rij[ix]*dfunc;
        colvar.myder[i_c][j][ix] += rij[ix]*dfunc;
      } 
      indice++;
    }
    indice += ((colvar.natoms[i_c]-i-colvar.type[i_c])*(npe-1));
    for(l=0;l<npe;l++) indice-=l;
  }
  if(logical.parallel_hills){ 
    for(i=0;i<colvar.natoms[i_c];i++){
      tmpx[i]=colvar.myder[i_c][i][0];
      tmpy[i]=colvar.myder[i_c][i][1];
      tmpz[i]=colvar.myder[i_c][i][2];
    }
    plumed_sum(mtd_data,1,&func); 
    plumed_sum(mtd_data,colvar.natoms[i_c],tmpx);
    plumed_sum(mtd_data,colvar.natoms[i_c],tmpy);
    plumed_sum(mtd_data,colvar.natoms[i_c],tmpz);
    for (i=0;i<colvar.natoms[i_c];i++){
      colvar.myder[i_c][i][0]=tmpx[i];
      colvar.myder[i_c][i][1]=tmpy[i];
      colvar.myder[i_c][i][2]=tmpz[i];
    }
  }
  func = sqrt(func);
  for(i=0;i<colvar.natoms[i_c];i++) for(ix=0;ix<3;ix++) colvar.myder[i_c][i][ix] /= func;

  colvar.ss0[i_c] = func;

  free(tmpx);
  free(tmpy);
  free(tmpz);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int PREFIX read_intracontact(char **word, int count, t_plumed_input *file, int *iline, FILE *fplog)
{
  int i, resnr, atnr, j, indice = 0, iw;
  double r_0, delta, db1, db2, db3;
  char string[400], c2, str;
  rvec rij, *refpos;
  real mod_rij, rdist, r6dist, r12dist, num, iden;
  real threshold;

  colvar.type[count]=3;
  iw = seek_word(word,"NATOMS");
  if(iw>=0) sscanf(word[iw+1],"%i", &colvar.natoms[count]);
  iw=seek_word(word,"SIGMA");
  if(iw>=0) sscanf(word[iw+1],"%lf", &delta);
  iw=seek_word(word,"NN");
  if(iw>=0) sscanf(word[iw+1],"%i", &colvar.nn[count]);
  iw=seek_word(word,"MM");
  if(iw>=0) sscanf(word[iw+1],"%i", &colvar.mm[count]);
  iw=seek_word(word,"R_0");
  if(iw>=0) sscanf(word[iw+1],"%lf", &r_0);
  iw=seek_word(word,"JUMP");
  if(iw>=0) sscanf(word[iw+1],"%i", &colvar.type[count]);

  colvar.delta_r[count]  = (real) delta;
  colvar.r_0[count]  = (real) r_0;

  fprintf(fplog, "\n%1i-CONTACT MAP: ATOMS %i SCALE %f TYPE %i\n", count+1,colvar.natoms[count], colvar.delta_r[count], colvar.type[count]);
  fprintf(fplog, "|--FUNCTIONAL FORM: n = %i m = %i r_0 = %f \n", colvar.nn[count], colvar.mm[count], r_0);
  fflush(fplog);
  snew(colvar.cvatoms[count], colvar.natoms[count]);
  snew(colvar.myder[count], colvar.natoms[count]);
  snew(colvar.map0[count], 0.5*colvar.natoms[count]*colvar.natoms[count]);
  snew(refpos, colvar.natoms[count]);

  for(i=0;i<colvar.natoms[count];i++){
    (*iline)++;
    sscanf(file->words[*iline][2],"%d", &atnr);
    sscanf(file->words[*iline][3],"%lf", &db1);
    sscanf(file->words[*iline][4],"%lf", &db2);
    sscanf(file->words[*iline][5],"%lf", &db3);
    atnr--;
    colvar.cvatoms[count][i] = atnr;
    fprintf(fplog, "|--%7i%15lf%15lf%15lf\n", atnr, db1, db2, db3);
    refpos[i][0] = (real) db1;
    refpos[i][1] = (real) db2;
    refpos[i][2] = (real) db3;
  }
  threshold=pow(0.00001,1./(colvar.nn[count]-colvar.mm[count]));
  indice=0;
  for(i=0;i<colvar.natoms[count];i++){
    for(j=i+colvar.type[count];j<colvar.natoms[count];j++) {
      rij[0] = refpos[i][0]-refpos[j][0];
      rij[1] = refpos[i][1]-refpos[j][1];
      rij[2] = refpos[i][2]-refpos[j][2];
      mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
      rdist = mod_rij/colvar.r_0[count];
      if(rdist<=0.){
        colvar.map0[count][indice] = 1.;
      } if(rdist>0.999999 && rdist<1.000001) {
        colvar.map0[count][indice] = colvar.nn[count]/colvar.mm[count];
      } else if(rdist>threshold) {
        colvar.map0[count][indice] = 0.;
      } else {
        r6dist = pow(rdist, colvar.nn[count]);
        r12dist = pow(rdist, colvar.mm[count]);
        num = 1. - r6dist;
        iden = 1./(1.- r12dist);
        colvar.map0[count][indice] = num*iden;
      }
      indice++;
    }
  }
  sfree(refpos);

  return colvar.natoms[count];
}
