/*********************************************************************/
/*                                                                   */
/*   fi_lib  --- A fast interval library (Version 1.2)               */
/*                                                                   */
/*  Authors:                                                         */
/*  --------                                                         */
/*  Werner Hofschuster, Walter Kraemer                               */
/*  Wissenschaftliches Rechnen/Softwaretechnologie                   */
/*  Universitaet Wuppertal, Germany                                  */
/*                                                                   */ 
/*  Copyright:                                                       */
/*  ----------                                                       */
/*  Copyright (C) 1997-2000 Institut fuer Wissenschaftliches Rechnen */
/*                          und Mathematische Modellbildung (IWRMM)  */
/*                                           and                     */
/*                          Institut fuer Angewandte Mathematik      */
/*                          Universitaet Karlsruhe, Germany          */
/*            (C) 2000-2005 Wiss. Rechnen/Softwaretechnologie        */
/*                          Universitaet Wuppertal, Germany          */
/*                                                                   */
/*  This library is free software; you can redistribute it and/or    */
/*  modify it under the terms of the GNU Library General Public      */
/*  License as published by the Free Software Foundation; either     */
/*  version 2 of the License, or (at your option) any later version. */
/*                                                                   */
/*  This library is distributed in the hope that it will be useful,  */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of   */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.             */
/*  See the GNU Library General Public License for more details.     */
/*                                                                   */
/*  You should have received a copy of the GNU Library General Public*/
/*  License along with this library; if not, write to the Free       */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston,            */
/*  MA  02111-1307  USA                                              */
/*                                                                   */
/*********************************************************************/

#include <stdio.h>
#include <math.h>
#include <float.h>

/*********************************************************************/
/*                    CONFIGURE FI_LIB.H here !!!                    */
/*********************************************************************/

#define INTEL 1      /* Version for PC ("little endian") must have 1,*/
                     /* otherwise ("big endian") 0                   */
#define LINT_ARGS 1  /* new style prototyping compilers must have 1, */
                     /* old one 0                                    */

/*********************************************************************/
/*                     end of configuring section                    */
/*********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/*********************************************************************/
/* Interval data type                                                */
/*********************************************************************/

#ifndef CXSC_INCLUDE
typedef struct interval { double INF, SUP;} interval ;
#endif

/*********************************************************************/
/* PROTOTYPE manages the difficulties in prototyping of K & R        */
/* and new style ANSI compilers                                      */
/*********************************************************************/

#ifndef _PROTOTYPE
#ifdef LINT_ARGS
#define _PROTOTYPE(function,params) function params
#else 
#define _PROTOTYPE(function,params) function()
#endif
#endif

/*********************************************************************/
/* preprocessing (global changes)                                    */
/*********************************************************************/

#define local
#define r_succ(params) q_succ(params)
#define r_pred(params) q_pred(params)

/*********************************************************************/
/* replacement for function ldexp from math.h (libm.a)               */
/*********************************************************************/

#define POWER2(x,k) {if (x!=0) {a_diee *test=(a_diee *) &x; test->ieee.expo+=k;}}               /* only for normalised numbers ! */

/*********************************************************************/
/* replacement for function frexp from math.h (libm.a)               */
/*********************************************************************/

#define FREXPO(x,k) {if (x!=0) {a_diee *test=(a_diee *) &x; k=test->ieee.expo;} else k=0;}                          /* only for normalised numbers ! */
 
/*********************************************************************/
/* first 24 bits of a double value ( cut / round )                   */
/*********************************************************************/

#define CUT24(x) (float)(x)

/*********************************************************************/
/* assignment double -> long int                                     */
/*********************************************************************/

#define CUTINT(x) (int)(x)

/*********************************************************************/
/* data type union to access sign, mantisse, exponent                */
/*********************************************************************/

typedef union 
	{
		double f;
		struct 
		{
#if INTEL
			unsigned int mant1 :32;
			unsigned int mant0 :20;
			unsigned int expo  :11;
			unsigned int sign  : 1;
#else
			unsigned int sign  : 1;
			unsigned int expo  :11;
			unsigned int mant0 :20;
			unsigned int mant1 :32;
#endif
		} ieee;
	} a_diee;

/*********************************************************************/
/* NAN (not a number) test                                           */
/*********************************************************************/

#define NANTEST(x) ((x!=x)) 
/* #define NANTEST(x) ((((a_diee) x).ieee.expo == 2047)&&((((a_diee) x).ieee.mant0!=0)||(((a_diee) x).ieee.mant1!=0))) */

/*********************************************************************/
/* error handling                                                    */
/*********************************************************************/

#define INV_ARG 1            /* for error handling: invalid argument */
#define OVER_FLOW 2          /* for error handling: overflow         */
#define DIV_ZERO 3           /* for error handling: division by zero */

/*********************************************************************/
/* constants for the elementary functions                            */
/*********************************************************************/

extern double q_ln2h;
extern double q_l10i;
extern double q_l2i;
extern double q_l10;
extern double q_l2;
extern double q_p2h;
extern double q_p2mh;
extern double q_mine;
extern double q_minr;
extern double q_pi;
extern double q_piha;
extern double q_pih[7];
extern double q_pi2i;
extern double q_pi2p[3];
extern double q_sqra;
extern double q_ctht;

extern double q_ext1;
extern double q_ext2;
extern double q_ex2a;
extern double q_ex2b;
extern double q_ex2c;
extern double q_ext3;
extern double q_ext4;
extern double q_ext5;
extern double q_extm;
extern double q_extn;
extern double q_exil;
extern double q_exl1;
extern double q_exl2;
extern double q_e10i;
extern double q_e1l1;
extern double q_e1l2;
extern double q_exa[5];
extern double q_exb[9];
extern double q_exc[7];
extern double q_exd[7];
extern double q_exld[32];
extern double q_extl[32];

extern double q_exem;
extern double q_exep;
extern double q_exmm;
extern double q_exmp;

extern double q_snhm;
extern double q_snhp;
extern double q_cshm;
extern double q_cshp;
extern double q_cthm;
extern double q_cthp;
extern double q_tnhm;
extern double q_tnhp;

extern double q_lgt1;
extern double q_lgt2;
extern double q_lgt3;
extern double q_lgt4;
extern double q_lgt5;
extern double q_lgt6;
extern double q_lgb[2];
extern double q_lgc[4];
extern double q_lgld[129];
extern double q_lgtl[129];

extern double q_logm;
extern double q_logp;
extern double q_lgpm;
extern double q_lgpp;
extern double q_sqtm;
extern double q_sqtp;

extern double q_atna[7];
extern double q_atnb[8];
extern double q_atnc[7];
extern double q_atnd[6];
extern double q_atnt;

extern double q_asnm;
extern double q_asnp;
extern double q_acsm;
extern double q_acsp;
extern double q_actm;
extern double q_actp;
extern double q_atnm;
extern double q_atnp;

extern double q_csnm;
extern double q_csnp;
extern double q_ccsm;
extern double q_ccsp;
extern double q_cctm;
extern double q_cctp;
extern double q_ctnm;
extern double q_ctnp;

extern double q_sinc[6];
extern double q_sins[6];
extern double q_sint[5];

extern double q_sinm;
extern double q_sinp;
extern double q_cosm;
extern double q_cosp;
extern double q_cotm;
extern double q_cotp;
extern double q_tanm;
extern double q_tanp;

extern double q_lg2m;
extern double q_lg2p;
extern double q_l10m;
extern double q_l10p;
extern double q_e2em;
extern double q_e2ep;
extern double q_e10m;
extern double q_e10p;

extern double q_at3i;

extern double q_erfm;
extern double q_erfp;
extern double q_efcm;
extern double q_efcp;
extern double q_expz[28];
extern double q_epA2[5];
extern double q_eqA2[5];
extern double q_epB1[7];
extern double q_eqB1[7];
extern double q_epB2[6];
extern double q_eqB2[7];
extern double q_epB3[5];
extern double q_eqB3[5];
extern double q_erft[7];

/*********************************************************************/
/* prototypes for interval arithmetic (basic operations)             */
/*********************************************************************/

_PROTOTYPE(double q_min,(double x,double y));
_PROTOTYPE(double q_max,(double x,double y));
_PROTOTYPE(double q_abs,(double x));    
_PROTOTYPE(interval j_abs,(interval x));    
_PROTOTYPE(interval add_ii,(interval x, interval y));
_PROTOTYPE(interval add_di,(double x, interval y));
_PROTOTYPE(interval add_id,(interval x, double y));
_PROTOTYPE(interval sub_ii,(interval x, interval y));
_PROTOTYPE(interval sub_id,(interval x, double y));
_PROTOTYPE(interval sub_di,(double x, interval y));
_PROTOTYPE(interval eq_ii,(interval y));
_PROTOTYPE(interval eq_id,(double y));
_PROTOTYPE(interval mul_ii,(interval x, interval y));
_PROTOTYPE(interval mul_id,(interval x, double y));
_PROTOTYPE(interval mul_di,(double x, interval y));
_PROTOTYPE(interval div_ii,(interval x, interval y));
_PROTOTYPE(interval div_di,(double x, interval y));
_PROTOTYPE(interval div_id,(interval x, double y));
_PROTOTYPE(int in_di,(double x, interval y));
_PROTOTYPE(int in_ii,(interval x, interval y));
_PROTOTYPE(int ieq_ii,(interval x, interval y));
_PROTOTYPE(int is_ii,(interval x, interval y));
_PROTOTYPE(int ig_ii,(interval x, interval y));
_PROTOTYPE(int ise_ii,(interval x, interval y));
_PROTOTYPE(int ige_ii,(interval x, interval y));
_PROTOTYPE(int dis_ii,(interval x, interval y));
_PROTOTYPE(double q_mid,(interval x));
_PROTOTYPE(interval hull,(interval x, interval y));
_PROTOTYPE(interval intsec,(interval x, interval y));
_PROTOTYPE(double q_diam,(interval x));

/*********************************************************************/
/* functions for internal use only                                   */
/*********************************************************************/

_PROTOTYPE(double q_abortr1,(int n, double *x, int fctn));
_PROTOTYPE(interval q_abortr2,(int n, double *x1, double *x2, int fctn));
_PROTOTYPE(double q_abortnan,(int n, double *x, int fctn));
_PROTOTYPE(double q_abortdivd,(int n, double *x));
_PROTOTYPE(interval q_abortdivi,(int n, double *x1, double *x2));
_PROTOTYPE(double q_ep1,(double x));
_PROTOTYPE(double q_epm1,(double x));
_PROTOTYPE(double q_cth1,(double x));
_PROTOTYPE(double q_log1,(double x));
_PROTOTYPE(double q_l1p1,(double x));
_PROTOTYPE(double q_atn1,(double x));
_PROTOTYPE(double q_sin1,(double x, long int k));
_PROTOTYPE(double q_cos1,(double x, long int k));
_PROTOTYPE(double q_rtrg,(double x, long int k));
_PROTOTYPE(double q_r2tr,(double r, long int k));

/*********************************************************************/
/* functions for IO                                                  */
/*********************************************************************/

_PROTOTYPE(double scandown,());
_PROTOTYPE(double scanup,());
_PROTOTYPE(interval scanInterval,());
_PROTOTYPE(void printdown,(double x));
_PROTOTYPE(void printup,(double x));
_PROTOTYPE(void printInterval,(interval x));

/*********************************************************************/
/*               prototypes for library functions                    */
/*********************************************************************/

_PROTOTYPE(double q_sqr,(double x));
_PROTOTYPE(double q_sqrt,(double x));
_PROTOTYPE(double q_exp,(double x));
_PROTOTYPE(double q_expm,(double x));
_PROTOTYPE(double q_sinh,(double x));
_PROTOTYPE(double q_cosh,(double x));
_PROTOTYPE(double q_coth,(double x));
_PROTOTYPE(double q_tanh,(double x));
_PROTOTYPE(double q_log,(double x));
_PROTOTYPE(double q_lg1p,(double x));
_PROTOTYPE(double q_asnh,(double x));
_PROTOTYPE(double q_acsh,(double x));
_PROTOTYPE(double q_acth,(double x));
_PROTOTYPE(double q_atnh,(double x));
_PROTOTYPE(double q_asin,(double x));
_PROTOTYPE(double q_acos,(double x));
_PROTOTYPE(double q_acot,(double x));
_PROTOTYPE(double q_atan,(double x));
_PROTOTYPE(double q_sin,(double x));
_PROTOTYPE(double q_cos,(double x));
_PROTOTYPE(double q_cot,(double x));
_PROTOTYPE(double q_tan,(double x));
_PROTOTYPE(double q_exp2,(double x));
_PROTOTYPE(double q_ex10,(double x));
_PROTOTYPE(double q_log2,(double x));
_PROTOTYPE(double q_lg10,(double x));
_PROTOTYPE(double q_erf,(double x));
_PROTOTYPE(double q_erfc,(double x));

_PROTOTYPE(double q_pred,(double x));
_PROTOTYPE(double q_succ,(double x));
_PROTOTYPE(double q_comp,(int s, double m, int e));
_PROTOTYPE(double q_cmps,(double m, int e));
_PROTOTYPE(int q_sign,(double x));
_PROTOTYPE(double q_mant,(double x));
_PROTOTYPE(double q_mnts,(double x));
_PROTOTYPE(int q_expo,(double x));

_PROTOTYPE(interval j_exp,(interval x));
_PROTOTYPE(interval j_expm,(interval x));
_PROTOTYPE(interval j_sinh,(interval x));
_PROTOTYPE(interval j_cosh,(interval x));
_PROTOTYPE(interval j_coth,(interval x));
_PROTOTYPE(interval j_tanh,(interval x));
_PROTOTYPE(interval j_log,(interval x));
_PROTOTYPE(interval j_lg1p,(interval x));
_PROTOTYPE(interval j_sqrt,(interval x));
_PROTOTYPE(interval j_sqr,(interval x));
_PROTOTYPE(interval j_asnh,(interval x));
_PROTOTYPE(interval j_acsh,(interval x));
_PROTOTYPE(interval j_acth,(interval x));
_PROTOTYPE(interval j_atnh,(interval x));
_PROTOTYPE(interval j_asin,(interval x));
_PROTOTYPE(interval j_acos,(interval x));
_PROTOTYPE(interval j_acot,(interval x));
_PROTOTYPE(interval j_atan,(interval x));
_PROTOTYPE(interval j_sin,(interval x));
_PROTOTYPE(interval j_cos,(interval x));
_PROTOTYPE(interval j_cot,(interval x));
_PROTOTYPE(interval j_tan,(interval x));
_PROTOTYPE(interval j_exp2,(interval x));
_PROTOTYPE(interval j_ex10,(interval x));
_PROTOTYPE(interval j_log2,(interval x));
_PROTOTYPE(interval j_lg10,(interval x));
_PROTOTYPE(interval j_erf,(interval x));
_PROTOTYPE(interval j_erfc,(interval x));

#ifdef __cplusplus
}
#endif
