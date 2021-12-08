/* Tracy-2

   J. Bengtsson  CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -

*/

#ifndef _THOR_SCSI_MATH_ETY_H_
#define _THOR_SCSI_MATH_ETY_H_

#include <armadillo>

void ETY(int n, int low, int high, arma::mat &a, arma::vec &ort);

void ETYT(int n, int low, int high, arma::mat &a, arma::vec &ort, arma::mat &z);

void ety2(int n, int low, int high, arma::mat &h, arma::vec &wr, arma::vec &wi,
	  arma::mat &z, int &ierr);

#endif /* _THOR_SCSI_MATH_ETY_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
