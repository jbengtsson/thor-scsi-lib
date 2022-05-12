/* Author:        Johan Bengtsson

    Definitions:  Interface to Fortran library for Truncated Power
		  Series Algebra.

    Note, linear case is a special case, see e.g. daexp

*/

#ifndef TPSA_FOR_PM_H
#define TPSA_FOR_PM_H

#include <armadillo>
#include <tps/ss_vect.h>
#include <tps/tps_type.h>
#include <cassert>

extern int  bufsize;  // Note, max no of monomials is (no+nv)!/(nv!*no!)


long int fact(long int n);

long int nok(long int n, long int k);

double getmat(const ss_vect<tps> &map, const int i, const int j);

void putmat(ss_vect<tps> &map, const int i, const int j, const double r);

/**
 * tpsa_for_pm.cc does not yet compile ...
 *
 * thus now implemented as an inline aborting program.
 * make it only a definition as soon as tpsa_for_pm.cc has been ported
 */
inline void getlinmat(const int nv, const ss_vect<tps> &map, arma::mat &mat){
  assert(0);
}

/**
 * tpsa_for_pm.cc does not yet compile ...
 *
 * thus now implemented as an inline aborting program.
 * make it only a definition as soon as tpsa_for_pm.cc has been ported
 */
inline ss_vect<tps> putlinmat(const int nv, const arma::mat &mat){
  assert(0);
}

void idprset(const int level);

tps atan2(const tps &b,const tps &a);

#endif /* TPSA_FOR_PM_H */
