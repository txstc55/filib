/*********************************************************************/
/*       fi_lib  --- A fast interval library (Version 1.2)           */
/*        (For copyright and info`s see file "fi_lib.h")             */
/*********************************************************************/

#include "fi_lib.h" 



#ifdef LINT_ARGS
local interval j_cos(interval x)
#else
local interval j_cos(x)

interval x;
#endif
{
  interval res;
  double erg1,erg2,y1,y2;
  long int k1,k2,q1,q2;

  if (x.INF==x.SUP)              /* point interval */
    { 
      if ((x.INF<-q_sint[2])||(x.SUP>q_sint[2]))
        { res.INF=-1.0; res.SUP=1.0; }
      else 
        { 
          res.INF=q_cos(x.INF);
          if (res.INF<0)
            {
              res.SUP=res.INF*q_cosm;
              res.INF*=q_cosp;
            }
          else
            {
              res.SUP=res.INF*q_cosp;
              res.INF*=q_cosm;
            }
        }
    }
  else
    {
      if ((x.SUP-x.INF)>=(2.0*q_pi)) 
        { res.INF=-1.0; res.SUP=1.0; }
      else if ((x.INF<-q_sint[2])||(x.SUP>q_sint[2]))
        { res.INF=-1.0; res.SUP=1.0; }
      else 
        {
          /* argument reduction infimum */
          erg1=x.INF*q_pi2i; 
          if (erg1>0) {k1=CUTINT(erg1+0.5);
                       q1=((CUTINT(erg1))+1)%4; }
          else        {k1=CUTINT(erg1-0.5);
                       q1=(CUTINT(erg1))%4; }
          y1=q_rtrg(x.INF,k1);
            
          /* argument reduction supremum */
          erg2=x.SUP*q_pi2i; 
          if (erg2>0) {k2=CUTINT(erg2+0.5);
                       q2=((CUTINT(erg2))+1)%4;}
          else        {k2=CUTINT(erg2-0.5);
                       q2=(CUTINT(erg2))%4;}
          y2=q_rtrg(x.SUP,k2);
          if (q1<0) q1+=4;
          if (q2<0) q2+=4;
         
          if (q1==q2) 
          {
            if ((x.SUP-x.INF)>=q_pi) 
              {
                res.INF=-1.0; res.SUP=1.0;
              }
            else if ((q1==1) || (q1==2))
              {
                res.INF=q_cos1(y2,k2);
                if (res.INF<0) res.INF*=q_sinp; else res.INF*=q_sinm;
                res.SUP=q_cos1(y1,k1);
                if (res.SUP<0) res.SUP*=q_sinm; else res.SUP*=q_sinp; 
              }
            else
              {
                res.INF=q_cos1(y1,k1);
                if (res.INF<0) res.INF*=q_sinp; else res.INF*=q_sinm;
                res.SUP=q_cos1(y2,k2);
                if (res.SUP<0) res.SUP*=q_sinm; else res.SUP*=q_sinp; 
              }
           }
           else  /* now we have q1<>q2 */
           { 
             if (q1==0)
             {
               if (q2==1)
               {
                 erg1=q_cos1(y1,k1);
                 erg2=q_cos1(y2,k2);
                 if (erg1<erg2) res.INF=erg1*q_sinm; else res.INF=erg2*q_sinm;
                 res.SUP=1.0;
               }
               else if (q2==2)
               {
                 res.INF=q_cos1(y2,k2)*q_sinp;
                 res.SUP=1.0;
               }
               else  /* q2==3 */
               {
                 res.INF=-1.0;
                 res.SUP=1.0;
               }
             }
             else if (q1==1)
             {
               if (q2==0)
               {
                 res.INF=-1.0;
                 erg1=q_cos1(y1,k1);
                 erg2=q_cos1(y2,k2);
                 if (erg1>erg2) res.SUP=erg1*q_sinp; else res.SUP=erg2*q_sinp;
               }
               else if (q2==2)
               {
                 res.INF=q_cos1(y2,k2)*q_sinp;
                 res.SUP=q_cos1(y1,k1)*q_sinp;
               }
               else  /* q2==3 */
               {
                 res.INF=-1.0;
                 res.SUP=q_cos1(y1,k1)*q_sinp;
               }
             }
             else if (q1==2)
             {
               if (q2==0)
               {
                 res.INF=-1.0;
                 res.SUP=q_cos1(y2,k2)*q_sinp;
               }
               else if (q2==1)
               {
                 res.INF=-1.0;
                 res.SUP=1.0;
               }
               else  /* q2==3 */
               {
                 res.INF=-1.0;
                 erg1=q_cos1(y1,k1);
                 erg2=q_cos1(y2,k2);
                 if (erg1>erg2) res.SUP=erg1*q_sinm; else res.SUP=erg2*q_sinm;
               }
             }
             else  /* now we have q1==3 */            
             {
               if (q2==0)
               {
                 res.INF=q_cos1(y1,k1)*q_sinp;
                 res.SUP=q_cos1(y2,k2)*q_sinp;
               }
               else if (q2==1)
               {
                 res.INF=q_cos1(y1,k1)*q_sinp;
                 res.SUP=1.0;
               }
               else  /* q2==2 */
               {
                 erg1=q_cos1(y1,k1);
                 erg2=q_cos1(y2,k2);
                 if (erg1<erg2) res.INF=erg1*q_sinp; else res.INF=erg2*q_sinp;
                 res.SUP=1.0;
               }
             }
           }
        }
    }

  if (res.INF<-1.0) res.INF=-1.0;
  if (res.SUP>1.0) res.SUP=1.0;
  return(res);
}
