/*********************************************************************/
/*       fi_lib  --- A fast interval library (Version 1.2)           */
/*        (For copyright and info`s see file "fi_lib.h")             */
/*********************************************************************/

#include "fi_lib.h" 



double scandown()
{
  double x;
  scanf("%lf",&x);
  if (fabs(x)>(1e17-1)*1e27) return q_pred(q_pred(x));
  else return q_pred(x);
}

double scanup()
{
  double x;
  scanf("%lf",&x);
  if (fabs(x)>(1e17-1)*1e27) return q_succ(q_succ(x));
  else return q_succ(x);
}

interval scanInterval()
{ 
  interval dummy;
  dummy.INF = scandown();
  dummy.SUP = scanup();
  return dummy;
}
