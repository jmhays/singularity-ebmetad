#include <stdio.h>
#include <time.h>
#include <string.h>
#include "config.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "network.h"
#include "repl_ex.h"
#include "nrnb.h"

void init_re_solute_temp(const t_commrec *mcr, t_forcerec *fr, const gmx_mtop_t *top, real rteio, FILE *fplog);
void rest_print_table(t_forcerec *fr, FILE *fplog);
real delta_res(real betaA, real betaB, int a, int b);
void ener_rest(const t_commrec *mcr, gmx_enerdata_t *grps, int repl, int nrepl, real VpmePP);

struct rex_s
{
  int    rest;
  real   t0;
  real   c6;
  real   c12;
  real   q;
  real   stpme;
  real   *Epot;
  real   *FRpot;
}rex;
