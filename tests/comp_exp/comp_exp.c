/*********************************************************************/
/*                                                                   */
/*   Example: comp_exp.c  (Computation of the exponential function)  */
/*          (For copyright and info's see file "fi_lib.h")           */
/*                                                                   */
/*********************************************************************/

#include<stdio.h>
#include<string.h>
#include"filib/fi_lib.h"     /* use library fi_lib */

/* --- main program ------------------------------------------------------ */

int main()
{
  interval x;
    
  printf("\n");
  printf("Computation of the exponential function in ANSI-C with fi_lib\n");
  printf("=============================================================\n\n");

  printf("Insert an interval argument (e.g. 'x = 1 1' or 'x = 1.01 1.02') \n");
  printf("x = ");
  x = scanInterval();

  printf("Argument x = "); 
  printInterval(x);
  printf("\n");
 
  printf("    exp(x) = "); 
  printInterval( j_exp(x) );
  printf("\n\n"); 

  return 0;
}

