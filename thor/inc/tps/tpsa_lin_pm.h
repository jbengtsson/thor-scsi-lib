 /* Author:	 Johan Bengtsson

   Definitions:  Polymorphic number class.              */

#ifndef TPSA_LIN_PM_H
#define TPSA_LIN_PM_H
// long int fact(long int n);

// long int nok(long int n, long int k);

void idprset(const int level);

void TPSA_Ini(void);

tps atan2(const tps &b, const tps &a);


void prt_coeff(std::ostringstream &s, const tps &a, const long int jj[],
	       const int ord, int &n_coeff);

void prt_header(std::ostream &s, const bool res_basis);
#endif
