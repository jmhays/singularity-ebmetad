#include "solute.h"

/* GROMACS */
/*
    {
        double io = compute_io(ir,top_global->natoms,groups,mdebin->ebin->nener,1);
        if ((io > 2000) && MASTER(cr))
            fprintf(stderr,
                    "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                    io);
    }
// I put it here because I need to modify the topology before the local setup
    if(repl_ex_nst>0&&rex.t0>0.) {
      init_re_solute_temp(cr, fr, top_global, ir->opts.ref_t[0], fplog);
      rex.rest=1;
    }
// Carlo

...

        bExchanged = FALSE;
        if ((repl_ex_nst > 0) && (step > 0) && !bLastStep &&
            do_per_step(step,repl_ex_nst)) 
        {
            real VpmePP=0.;
            if(rex.rest&&MASTER(cr)) 
              ener_rest(cr, enerd, replica_exchange_get_repl(repl_ex), 
                        replica_exchange_get_nrepl(repl_ex), VpmePP);

            bExchanged = replica_exchange(fplog,cr,repl_ex,
                                          state_global,enerd->term,
                                          state,step,t);
*/


void init_re_solute_temp(const t_commrec *mcr, t_forcerec *fr, const gmx_mtop_t *top, real rteio, FILE *fplog)
{
  int aj, j, k, atnr, natoms, cond;
  int molsol=-1, soltype=-1;
  int types[100], ntype=0;

  rex.q = rex.c6 = rex.c12 = sqrt(rteio/rex.t0);
  if(MASTER(mcr)) {
    if(rex.q<=0.) { fprintf(stderr, "SOMETHING WRONG WITH THE RESCALING FACTOR OF REST (%lf)\n", rex.q); exit(0); }
    fprintf(fplog, "\n|-REST: there are %i molecule types\n", top->nmoltype);
  }

  /* all solvent molecule types are identified */
  for(j=0;j<top->nmoltype;j++) { if(!strcmp(*top->moltype[j].name,"SOL")) {molsol=j; break;} }
  if(molsol==-1&&MASTER(mcr)) {fprintf(stderr, "Not found any solvent molecules to apply solute tempering!!!\n"); exit(0);}
  rex.rest=1;
  types[0] = top->moltype[molsol].atoms.atom[0].type;
  ntype=1;
  for(aj=1;aj<top->moltype[molsol].atoms.nr;aj++) {
    /* all the charges of the molecules molsol are rescaled */
    top->moltype[molsol].atoms.atom[aj].q *= rex.q;
    /* all the atom types are saved */
    soltype = top->moltype[molsol].atoms.atom[aj].type;
    cond=1;
    for(k=0; k<ntype; k++) if(soltype==types[k]) { cond=0; }
    if(cond) { types[ntype]=soltype; ntype++; }
  }
  if(MASTER(mcr)) fprintf(fplog, "|-REST: IN SOL %i atom types found\n", ntype);
  /* Here you can add any kind of molecules */
  /*
  molsol=-1;
  for(j=0;j<top->nmoltype;j++) { if(!strcmp(*top->moltype[j].name,"DMPC")) {molsol=j; break;} }
  if(molsol!=-1) {
    types[ntype] = top->moltype[molsol].atoms.atom[0].type;
    ntype++;
    for(aj=1;aj<top->moltype[molsol].atoms.nr;aj++) {
      // all the charges of the molecules molsol are rescaled 
      top->moltype[molsol].atoms.atom[aj].q *= rex.q;
      // all the atom types are saved 
      soltype = top->moltype[molsol].atoms.atom[aj].type;
      cond=1;
      for(k=0; k<ntype; k++) if(soltype==types[k]) { cond=0; }
      if(cond) { types[ntype]=soltype; ntype++; }
    }
  if(MASTER(mcr)) fprintf(fplog, "|-REST: AFTER DMPC %i atom types found\n", ntype);
  }
  */

  /* rescaling VdW */
  atnr = fr->ntype;

  for(j=0; (j<atnr); j++){
    for(aj=0; (aj<ntype); aj++){
      if(types[aj]!=j){ 
        C6(fr->nbfp, atnr, types[aj], j)  *= rex.c6;
        C12(fr->nbfp, atnr, types[aj], j) *= rex.c12;
      } else {
        C6(fr->nbfp, atnr, types[aj], j)  *= rex.c6*rex.c6;
        C12(fr->nbfp, atnr, types[aj], j) *= rex.c12*rex.c12;
      }
    }
  }

  if(MASTER(mcr)) { fprintf(fplog, "|-REST: SOLUTE TEMPERING INITIALIZED! WITH FACTOR %lf\n\n", rex.q); fflush(fplog); }
}

void rest_print_table(t_forcerec *fr, FILE *fplog)
{
  int j, aj;
  int atnr = fr->ntype;

  fprintf(fplog, "|-REST FINAL TABLE\n");
  for(j=0; (j<atnr); j++){
    for(aj=0; (aj<atnr); aj++){
      fprintf(fplog, "|-REST: i %i j %i C6 %lf C12 %lf\n", aj, j,  C6(fr->nbfp, atnr, aj, j), C12(fr->nbfp, atnr, aj, j));
    }
  }
  fflush(fplog);

}

void ener_rest(const t_commrec *mcr, gmx_enerdata_t *grps, int repl, int nrepl, real VpmePP)
{
    int i;
    static int count=0;

    if(!count) {
      if(grps->grpp.nener<2) {
        fprintf(stderr, "SOLUTE TEMPERING ERROR: in the mdp file you must define two ener groups for solute and solvent!!!\n"); 
        exit(0);
      }
      snew(rex.Epot, nrepl);
      snew(rex.FRpot, nrepl);
    }

    for(i=0;i<nrepl;i++) { rex.Epot[i] = 0.; rex.FRpot[i] = 0.; }

    rex.Epot[repl] = grps->grpp.ener[egCOULSR][0] + grps->grpp.ener[egLJSR][0] +
                     grps->grpp.ener[egLJ14][0] + grps->grpp.ener[egCOUL14][0] +
                     grps->grpp.ener[egBHAMSR][0] + grps->grpp.ener[egBHAMLR][0] +
                     VpmePP;

    rex.FRpot[repl] = grps->grpp.ener[egCOULSR][1] + grps->grpp.ener[egLJSR][1] +
                      grps->grpp.ener[egLJ14][1] + grps->grpp.ener[egCOUL14][1] +
                      grps->grpp.ener[egBHAMSR][1] + grps->grpp.ener[egBHAMLR][1];

    gmx_sum_sim(nrepl, rex.FRpot, mcr->ms);
    gmx_sum_sim(nrepl, rex.Epot, mcr->ms);

    //fprintf("fplog, REPL %i EPOT %lf FRPOT %lf\n", repl, rex.Epot[repl], rex.FRpot[repl]); 

    count++;
}

real delta_res(real betaA, real betaB, int a, int b)
{
  real delta_prot, delta_pw, delta;
  real beta_ref = 1./(rex.t0*BOLTZ);

  delta_prot = (betaA-betaB)*(rex.Epot[b]-rex.Epot[a]);
  delta_pw = (sqrt(beta_ref*betaA)-sqrt(beta_ref*betaB))*(rex.FRpot[b]-rex.FRpot[a]);
  delta = delta_prot + delta_pw;

  return delta;
}
