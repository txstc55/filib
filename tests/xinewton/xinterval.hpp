/*********************************************************************/
/*                                                                   */
/*   Example: xinterval.hpp  (Extended interval arithmetic)          */
/*       (For copyright and info's see file "fi_lib.h")              */
/*                                                                   */
/*********************************************************************/

#include"filib/interval.hpp"

typedef enum { Finite, PlusInfty, MinusInfty, Double, Empty } KindType;

//----------------------------------------------------------------------------

class xinterval {                  // Extended intervals according
  public:                          // to the above definition
    KindType  kind;                //-----------------------------
    double    inf, sup;

    xinterval ( );
    xinterval ( const KindType&, const double&, const double& );
    xinterval ( const xinterval& );

    xinterval& operator= ( const xinterval& );

    friend xinterval operator% (const interval& A, const interval& B );
    friend xinterval operator- ( const double a, const xinterval B );
    friend interval* operator& ( interval X, const xinterval Y );
};

//----------------------------------------------------------------------------
 
interval EmptyIntval ( );          // Irregular (empty) interval

//----------------------------------------------------------------------------

interval EmptyIntval ( )           // Irregular (empty) interval
{                                  //---------------------------
  interval x;
  x.INF = 999999999.0;
  x.SUP = -999999999.0;
  return x;
}

xinterval::xinterval ( )
{
  kind = Finite;
  inf  = 0.0;
  sup  = 0.0;
}

xinterval::xinterval ( const KindType& k, const double& i, const double& s )
{
  kind = k;
  inf  = i;
  sup  = s;
}

xinterval::xinterval ( const xinterval& a )
{
  kind = a.kind;
  inf  = a.inf;
  sup  = a.sup;
}

xinterval& xinterval::operator= ( const xinterval& a )
{
  kind = a.kind;
  inf  = a.inf;
  sup  = a.sup;
  return *this;
}

//----------------------------------------------------------------------------
// Extended interval division 'A / B' where 0 in 'B' is allowed.
//----------------------------------------------------------------------------

xinterval operator% (const interval& A, const interval& B )
{
  interval  c;
  xinterval Q;

  if ( in(0.0, B) ) {
    if ( in(0.0, A) ) {
      Q.kind = Double;                    // Q = [-oo,+oo] = [-oo,0] v [0,+oo]
      Q.sup  = 0.0;                       //----------------------------------
      Q.inf  = 0.0;
    }
    else if ( B == 0.0 ) {                                          // Q = [/]
      Q.kind = PlusInfty;                                           //--------
      Q.inf  = q_pred(sup(A)/inf(B));
    }
    else if ( (sup(A) < 0.0) && (sup(B) == 0.0) ) {         // Q = [Q.inf,+oo]
      Q.kind = PlusInfty;                                   //----------------
      Q.inf  = q_pred(sup(A)/inf(B));
    }
    else if ( (sup(A) < 0.0) && (inf(B) < 0.0) && (sup(B) > 0.0) ) {
      Q.kind = Double;                        // Q = [-oo,Q.sup] v [Q.inf,+oo]
      Q.sup  = q_succ(sup(A)/sup(B));         //------------------------------
      Q.inf  = q_pred(sup(A)/inf(B));
    }
    else if ( (sup(A) < 0.0) && (inf(B) == 0.0) ) {         // Q = [-oo,Q.sup]
      Q.kind = MinusInfty;                                  //----------------
      Q.sup  = q_succ(sup(A)/sup(B));
    }
    else if ( (inf(A) > 0.0) && (sup(B) == 0.0) ) {         // Q = [-oo,Q.sup]
      Q.kind = MinusInfty;                                  //----------------
      Q.sup  = q_succ(inf(A)/inf(B));
    }
    else if ( (inf(A) > 0.0) && (inf(B) < 0.0) && (sup(B) > 0.0) ) {
      Q.kind = Double;                        // Q = [-oo,Q.sup] v [Q.inf,+oo]
      Q.sup  = q_succ(inf(A)/inf(B));         //------------------------------
      Q.inf  = q_pred(inf(A)/sup(B));
    }
    else { // if ( (Inf(A) > 0.0) && (Inf(B) == 0.0) )
      Q.kind = PlusInfty;                                   // Q = [Q.inf,+oo]
      Q.inf  = q_pred(inf(A)/sup(B));                       //----------------
    }
  } // in(0.0,B)
  else {  // !in(0.0,B)
    c = A / B;                                            // Q = [C.inf,C.sup]
    Q.kind = Finite;                                      //------------------
    Q.inf  = inf(c);
    Q.sup  = sup(c);
  }

  return Q;
} 

//----------------------------------------------------------------------------
// Subtraction of an extended interval 'B' from a double value 'a'.
//----------------------------------------------------------------------------

xinterval operator- ( const double a, const xinterval B )
{
  xinterval D;

  switch (B.kind) {
    case Finite     : D.kind = Finite;                    // D = [D.inf,D.sup]
                      D.inf  = q_pred(a-B.sup);           //------------------
                      D.sup  = q_succ(a-B.inf);
                      break;
    case PlusInfty  : D.kind = MinusInfty;                    // D = [inf,+oo]
                      D.sup  = q_succ(a-B.inf);               //--------------
                      break;
    case MinusInfty : D.kind = PlusInfty;                     // D = [-oo,sup]
                      D.inf  = q_pred(a-B.sup);               //--------------
                      break;
    case Double     : D.kind = Double;        // D = [-oo,D.sup] v [D.inf,+oo]
                      D.inf  = q_pred(a-B.sup); //----------------------------
                      D.sup  = q_succ(a-B.inf);
                      if (D.inf < D.sup) D.inf = D.sup;
                      break;
    case Empty      : D.kind = Empty;                               // D = [/]
                      D.inf  = q_pred(a-B.sup);                     //--------
                      break;
  } // switch
  return D;
}

//----------------------------------------------------------------------------
// Intersection of an interval 'X' and an extended interval 'Y'. The result
// is given as a pair (vector) of intervals, where one or both of them can
// be empty intervals.
//----------------------------------------------------------------------------

interval* operator& ( interval X, const xinterval Y )
{
  interval H;
  interval* IS = new interval[3];

  IS[1] = EmptyIntval();
  IS[2] = EmptyIntval();

  switch (Y.kind) {
    case Finite     : // [X.inf,X.sup] & [Y.inf,Y.sup]
                      //------------------------------
                      H = _interval(Y.inf,Y.sup);
                      if ( !disjoint(X,H) ) IS[1] = X & H;

                      break;
    case PlusInfty  : // [X.inf,X.sup] & [Y.inf,+oo]
                      //----------------------------
                      if (sup(X) >= Y.inf)
                        if (inf(X) > Y.inf)
                          IS[1] = X;
                        else
                          IS[1] = _interval(Y.inf,sup(X));

                      break;
    case MinusInfty : // [X.inf,X.sup] & [-oo,Y.sup]
                      //----------------------------
                      if (Y.sup >= inf(X))
                        if (sup(X)<Y.sup)
                          IS[1] = X;
                        else
                          IS[1] = _interval(inf(X),Y.sup);

                      break;
    case Double     : if ( (inf(X) <= Y.sup) && (Y.inf <= sup(X)) ) {
                        IS[1] = _interval(inf(X),Y.sup);    // X & [-oo,Y.sup]
                        IS[2] = _interval(Y.inf,sup(X));    // X & [Y.inf,+oo]
                      }
                      else if (Y.inf <= sup(X)) // [X.inf,X.sup] & [Y.inf,+oo]
                        if (inf(X) >= Y.inf)    //----------------------------
                          IS[1] = X;
                        else
                          IS[1] = _interval(Y.inf,sup(X));
                      else if (inf(X) <= Y.sup) // [X.inf,X.sup] & [-oo,Y.sup]
                        if (sup(X) <= Y.sup)    //----------------------------
                          IS[1] = X;
                        else
                          IS[1] = _interval(inf(X),Y.sup);

                      break;
    case Empty      : break;                           // [X.inf,X.sup] ** [/]
  } // switch                                          //---------------------

  return IS;
} // operator&
