/*********************************************************************/
/*                                                                   */
/*   Example: hornerc.c  (Interval Horner's scheme in ANSI-C)        */
/*          (For copyright and info's see file "fi_lib.h")           */
/*                                                                   */
/*********************************************************************/

#include<stdio.h>
#include<string.h>
#include"filib/fi_lib.h"     /* use library fi_lib */

/* --- main program ------------------------------------------------------ */

int main()
{
  interval coeff[16];
  interval p, x, res;
  int i;

  /* --- Computation of the coefficients -------------------------------- */
  coeff[0] = eq_id(1.0);
  coeff[1] = eq_id(1.0);
  p = eq_id(1.0);
  for (i=2; i<=15; i++){
    p = mul_id(p,(double)i);
    coeff[i] = div_di(1.0,p);
  }
    
  printf("Interval Horner's scheme in ANSI-C with fi_lib\n");
  printf("==============================================\n\n");

  printf("Computation of the polynom (sum_i=0^15 1/i! x^i), \n");
  printf(" that means the first terms of the taylor series \n");
  printf(" of the exponential function.\n\n");
  printf("Enclosures for the polynom coefficients:\n");

  for (i=0; i<=15; i++){
    printf("   coeff[%2d] = ",i);
    printInterval(coeff[i]);
    printf("\n");
  }
  printf("\n");

  printf("Now, you can choose an interval argument (e.g. 'x = 1 1', \n");
  printf("'x = 1.01 1.02', 'x= -1 -1' or 'x= -2.0 -1.99' ...): \n\n");
  printf("x = ");
  x = scanInterval();

  /* --- interval Horner's scheme ---------------------------------------- */
  res = coeff[15];
  for (i=14; i>=0; i--) {
    res = mul_ii(res,x);
    res = add_ii(res,coeff[i]);
  }

  printf("Result: "); 
  printInterval(res);
  printf("\n"); 

  return 0;
}

