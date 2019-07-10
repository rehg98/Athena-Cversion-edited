#include "copyright.h"
/*============================================================================*/
/*! \file blast.c
 *  \brief Problem generator for spherical blast wave problem.
 *
 * PURPOSE: Problem generator for spherical blast wave problem.
 *
 * REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
 *   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.     */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"


static Real total_rad_mom(const GridS *pG, const int i, const int j, const int k);
static Real internal_E(const GridS *pG, const int i, const int j, const int k);
static Real bubble_vol(const GridS *pG, const int i, const int j, const int k);

Real cs_ambient;

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  Prim1DS W;
  Cons1DS U1d;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real pressure,drat,prat,rad,pa,da,x1,x2,x3, vol=0.0, dvol=pGrid->dx1*pGrid->dx2*pGrid->dx3; 
  Real b0=0.0,Bx=0.0,rin;
  double theta;

  rin = par_getd("problem","radius");
  pa  = par_getd("problem","pamb");
  da  = par_getd_def("problem","damb",1.0);
  drat = par_getd_def("problem","drat",1.0);
  prat = par_getd("problem","prat");
#ifdef MHD
  b0 = par_getd("problem","b0");
  theta = (PI/180.0)*par_getd("problem","angle");
#endif

	cs_ambient = sqrt(Gamma * pa/da);

/* setup uniform ambient medium with spherical over-pressured region */

  W.d = da;
  W.Vx = 0.0;
  W.Vy = 0.0;
  W.Vz = 0.0;
#ifdef MHD
  Bx = b0*cos(theta);
  W.By = b0*sin(theta);
  W.Bz = 0.0;
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	rad = sqrt(x1*x1 + x2*x2 + x3*x3);
#ifndef ISOTHERMAL
        W.P = pa;
	if (rad < rin) {
      W.P = prat*pa;
      vol += dvol;	
    }
#endif
        W.d = da;
	if (rad < rin) W.d = drat*da;

        U1d = Prim1D_to_Cons1D(&(W),&Bx);

	pGrid->U[k][j][i].d  = U1d.d;
	pGrid->U[k][j][i].M1 = U1d.Mx;
	pGrid->U[k][j][i].M2 = U1d.My;
	pGrid->U[k][j][i].M3 = U1d.Mz;

	

#ifndef ISOTHERMAL
	pGrid->U[k][j][i].E  = U1d.E;
#endif
#ifdef MHD
	pGrid->B1i[k][j][i] = Bx;
	pGrid->B2i[k][j][i] = U1d.By;
	pGrid->B3i[k][j][i] = U1d.Bz;
	pGrid->U[k][j][i].B1c = Bx;
	pGrid->U[k][j][i].B2c = U1d.By;
	pGrid->U[k][j][i].B3c = U1d.Bz;
	if (i == ie && ie > is) pGrid->B1i[k][j][i+1] = Bx;
	if (j == je && je > js) pGrid->B2i[k][j+1][i] = U1d.By;
	if (k == ke && ke > ks) pGrid->B3i[k+1][j][i] = U1d.Bz;
#endif /* MHD */
      }
    }
  }
#ifdef RESISTIVITY 
  eta_Ohm = 0.0;
  Q_AD    = par_getd("problem","Q_AD");
  Q_Hall  = 0.0;
  d_ind   = 0.0;
#endif

dump_history_enroll(total_rad_mom, "total_rad_mom");
dump_history_enroll(internal_E, "internal_E");
dump_history_enroll(bubble_vol, "bubble_vol");
}


/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
	dump_history_enroll(total_rad_mom, "total_rad_mom");
	dump_history_enroll(internal_E, "internal_E");
	dump_history_enroll(bubble_vol, "bubble_vol");
}

ConsFun_t get_usr_expr(const char *expr)
{
	if(strcmp(expr,"total_rad_mom")==0) return total_rad_mom;
	if (strcmp(expr, "internal_E")==0) return internal_E;
	return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef RESISTIVITY
void get_eta_user(GridS *pG, int i, int j, int k,
                             Real *eta_O, Real *eta_H, Real *eta_A)
{

  *eta_O = 0.0;
  *eta_H = 0.0;
  *eta_A = 0.0;

  return;
}
#endif


void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}

static Real total_rad_mom(const GridS *pG, const int i, const int j, const int k)
{
	Real x1, x2, x3;
	Real M1 = pG->U[k][j][i].M1;
	Real M2 = pG->U[k][j][i].M2;
	Real M3 = pG->U[k][j][i].M3;

	cc_pos(pG,i,j,k,&x1,&x2,&x3);
	Real rad = sqrt(x1*x1 + x2*x2 + x3*x3);

	return 1.0 / rad * (M1 * x1 + M2 * x2 + M3 * x3);
}


static Real internal_E(const GridS *pG, const int i, const int j, const int k)
{
    ConsS U = pG->U[k][j][i];
    PrimS W = Cons_to_Prim(&U);
	
	Real cs = sqrt(Gamma * W.P/W.d); 
	
	Real rmom = total_rad_mom(pG, i, j, k);
	Real vr = rmom / W.d;

	if (fabs(vr) > 0.01 || cs > 1.1*cs_ambient) return W.P / Gamma_1;

	return 0;
}


static Real bubble_vol(const GridS *pG, const int i, const int j, const int k)
{
	ConsS U = pG->U[k][j][i];
    PrimS W = Cons_to_Prim(&U);
	
	Real cs = sqrt(Gamma * W.P/W.d); 

	Real rmom = total_rad_mom(pG, i, j, k);
	Real vr = rmom / W.d;

	if (fabs(vr) > 0.01 || cs > 1.1*cs_ambient) return 1;

	return 0;
}


