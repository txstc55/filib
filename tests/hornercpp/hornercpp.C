/*********************************************************************/
/*                                                                   */
/*   Example: hornercpp.C  (Interval Horner's scheme in C++)         */
/*          (For copyright and info's see file "fi_lib.h")           */
/*                                                                   */
/*********************************************************************/

#include<iostream>
#include<string>
#include"filib/interval.hpp"    /* use library fi_lib with C++ - interface */

using namespace std;

/* --- main program ----------------------------------------------------- */

int main()
{
  interval coeff[16];
  interval p, x, res;
  int i;

  /* --- Computation of the coefficients -------------------------------- */
  coeff[0] = _interval(1.0);
  coeff[1] = _interval(1.0);
  p = _interval(1.0);
  for (i=2; i<=15; i++){
    p = p*(double)i;
    coeff[i] = 1.0/p;
  }
    
  cout << "Interval Horner's scheme in C++ with fi_lib" << endl;
  cout << "===========================================" << endl << endl;

  cout << "Computation of the polynom (sum_i=0^15 1/i! x^i)," << endl;
  cout << "that means the first terms of the taylor series" << endl;
  cout << "of the exponential function." << endl << endl;
  cout << "Enclosures for the polynom coefficients:" << endl;

  for (i=0; i<=15; i++)
    cout << "   coeff[" << setw(2) << i << "] = " << coeff[i] << endl;
  cout << endl;

  cout << "Now, you can choose an interval argument (e.g. 'x = 1 1', " << endl;
  cout << "'x = 1.01 1.02', 'x= -1 -1' or 'x= -2.0 -1.99' ...): " << endl;
  cout << endl;
  cout << "x = ";
  cin >> x;

  /* --- interval Horner's scheme ---------------------------------------- */
  res = coeff[15];
  for (i=14; i>=0; i--) 
    res = res*x + coeff[i];

  cout << "Result: " << res << endl; 

  return 0;
}

