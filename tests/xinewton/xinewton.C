/*********************************************************************/
/*                                                                   */
/*   Example: xinewton.c  (Compute all zeros of a function)          */
/*       (For copyright and info's see file "fi_lib.h")              */
/*                                                                   */
/*********************************************************************/

#include "xinterval.hpp" 
#include <iostream> 
#include <stdio.h>         // for function sprintf

using namespace std;

// --------------------------------------------------------------------
// ---   Initialisation, data type list	                            ---
// --------------------------------------------------------------------

struct list                  // Data type for resulting list
{
  interval intval;
  int info;
  list* next;
  
  list() {intval=_interval(0.0); next=NULL;}
};

static int MaxZeroNo=10;
const  int NoError = 0,    // Error constants 
       IllMaxZeroNo = 1, 
       NotAllZeros = 2;
const  int MaxCount = 10000;

//----------------------------------------------------------------------------
// Definition of the function f(x) and the derivate df(x)
//     Test function f = cosh(x) + 10 * x^2 * sin(x)^2 - 34                    
//----------------------------------------------------------------------------

interval f(interval x) {
  return ( cosh(x)+10*sqr(x)*sqr(sin(x))-34 );
}

interval df(interval x) {
  return( sinh(x)+20*x*sqr(sin(x))+20*sqr(x)*sin(x)*cos(x) );
}

//----------------------------------------------------------------------------
// main function xinewton
//----------------------------------------------------------------------------

void xinewton(interval fy, interval y, double eps, int yUnique, 
                               list* &zerolist, int& zerono)
{
   if (zerono > MaxZeroNo) return;
   interval* V;
   xinterval z;
   double c;
   interval ci;
   int i;
   int absorbed;
   list* pointer;
   pointer = zerolist;

   if (!in(0.0,fy)) return;

   c = mid(y);
   ci=_interval(c);
   z = c - f(ci)%df(y);
   V = y & z;

   if (V[1] == y) {                      // bisection
      V[1] = _interval(y.INF,c);
      V[2] = _interval(c,y.SUP);
   }

   if ( (V[1] != EmptyIntval() ) && (V[2] == EmptyIntval()))
      yUnique = yUnique || in(V[1],y);
   else
      yUnique = 0;

   for (i=1;i<=2;i++) {
      if (V[i] == EmptyIntval())  continue; 
      if (drel(V[i])<=eps) {
         fy=f(V[i]);
         if (in(0.0,fy)) {
            zerono ++;
            if (zerono > MaxZeroNo) return;

            absorbed=0;
            pointer=zerolist;
            //- Try to absorb V[i] by an already computed element of list -
            if (zerolist != NULL) {
               while ((pointer != NULL) && (!absorbed)) {
                  absorbed = !( (sup(V[i])<inf(pointer->intval)) ||
                                (sup(pointer->intval)<inf(V[i]) ));
                  if (absorbed) { 
                     pointer->intval = ((pointer->intval) | V[i]);
                     pointer->info = 0;
                     zerono --;
                  }
                  pointer = pointer -> next;
               }
            }

            if (!absorbed) {   // Store x in the resulting list
               pointer=zerolist;
               if (zerolist != NULL) {
                  while (pointer->next != NULL) pointer = pointer->next;
                  list* insert = new list;
                  pointer->next = insert;
                  insert->next = NULL;
                  insert->intval = V[i];
                  insert->info = yUnique; 
               } else {
                  zerolist = new list;
                  zerolist->next = NULL;
                  zerolist->intval = V[i];
                  zerolist->info = yUnique;         
               }
            }

         }
       } else {	
         xinewton(f(V[i]),V[i],eps,yUnique,zerolist,zerono); 
       }
   }
}

// --------------------------------------------------------------------
// ---   function AllZerosErrMsg (error message)                    ---
// --------------------------------------------------------------------

char* AllZerosErrMsg (int Err)
{ 
   static char Msg[80] = "";

   switch (Err) {
      case NoError:      break;
      case IllMaxZeroNo: sprintf(Msg,"Error: Parameter for maximum number of zeros must lie in 1,...,%1d!", MaxCount); 
                         break;
      case NotAllZeros:  sprintf(Msg,"Warning: Not all zeros found due to the user limit of %1d zero(s)!",MaxZeroNo);
                         break;
      default:           sprintf(Msg,"Error Code not defined");
   }

   return (Msg);
}

// --------------------------------------------------------------------
// ---   function VerificationStep                                  ---
// --------------------------------------------------------------------

static void VerificationStep(interval& y, int& unique)
{
   const int kmax= 10;
   interval yIn, yOld, fc, dfy;
   double c,eps;
   int k;

   k = 0; 
   yIn = y; 
   eps = 0.25; 
   unique = 0;

   while (!unique && (k < kmax) ) {
      yOld = blow(y,eps);
      dfy = df(y); 
      if (in(0.0,dfy)) break;
      k++;
      c = mid(yOld);
      fc=f(_interval(c));
      y = c - fc / dfy;
      if (disjoint(y, yOld)) break;
      unique = in(y,yOld);
      y = y & yOld;
      if (y == yOld) eps=eps*8.0;
   }
   if (!unique) y = yIn;
}

// --------------------------------------------------------------------
// ---   function AllZeros (root finding and verification)          ---
// --------------------------------------------------------------------

void AllZeros(interval Start, double Epsilon, list* &zerolist, 
                   int& NumberOfZeros, int& Err, int MaxNumberOfZeros) {

   double MinEpsilon;
   list* pointer;

   if (1<=MaxNumberOfZeros && MaxNumberOfZeros <= MaxCount) {
      MaxZeroNo = MaxNumberOfZeros;
      Err = NoError; 
      NumberOfZeros = 0;
      MinEpsilon = q_succ(1.0)-1.0;			// 1ulp
      if (Epsilon < MinEpsilon)
         Epsilon = MinEpsilon;

      xinewton(f(Start), Start, Epsilon, 0 , zerolist, NumberOfZeros);

      if (NumberOfZeros > MaxNumberOfZeros) {	
         Err = NotAllZeros; 
         NumberOfZeros = MaxNumberOfZeros;
      }
   } else {
      Err = IllMaxZeroNo;
      NumberOfZeros = 0;
   }

   pointer = zerolist;
   while (pointer != NULL) {
      if (pointer->info==0) {
        VerificationStep(pointer->intval,pointer->info); 
      }
      pointer = pointer->next;
   }
}

// --------------------------------------------------------------------
// ---   function main                                              ---
// --------------------------------------------------------------------

int main ()
{
   interval SearchInterval;
   double Tolerance;
   int NumberOfZeros, Error;
   list* Zero;
   Zero = NULL;

   cout << "Extended-Interval-Newton method in C++ with fi_lib" << endl;
   cout << "==================================================" << endl << endl;
 
   cout << "Computing all zeros of the function f(x) = cosh(x) + 10 * x^2 * sin(x)^2 - 34 " 
        << endl;
   cout << "Search interval (e.g. '-10 10'): " ;
   cin >> SearchInterval;
   cout << endl << "Search interval = " << SearchInterval << endl;
   cout << "Tolerance (relative) (e.g. '1e-3' or '1e-15'): ";
   cin >> Tolerance;
   cout << endl;

   AllZeros(SearchInterval, Tolerance, Zero, NumberOfZeros, 
                                                        Error, MaxCount);


   if (Zero==NULL)
     cout << "Function contains no zeros in the search interval!" << endl;
   else {
     cout << "Ranges for zeros:" << endl;
     while (Zero != NULL) {
       cout << Zero->intval << endl;
       if (Zero->info)
          cout << " encloses a locally unique zero !" << endl;
       else
          cout << " may contain a zero (not verified unique)!"<< endl;
       Zero = Zero->next;
     }
   }
   cout << endl << NumberOfZeros << " interval enclosure(s)" << endl;
   if (Error) cout << endl << AllZerosErrMsg(Error) << endl;

   return 0;
}
