#include <thor_scsi/elements/element_helpers.h>
#include <tps/tps_type.h>

namespace tse = thor_scsi::elements;

// template tps sqr(const tps &);

void tse::get_twoJ(const int n_DOF, const ss_vect<double> &ps, const ss_vect<tps> &A,
	      double twoJ[])
{
  int             j, no;
  long int        jj[ps_dim];
  ss_vect<double> z;

  no = no_tps;
  danot_(1);

  for (j = 0; j < ps_dim; j++)
    jj[j] = (j < 2*n_DOF)? 1 : 0;

  z = (PInv(A, jj)*ps).cst();

  for (j = 0; j < n_DOF; j++)
    twoJ[j] = sqr(z[2*j]) + sqr(z[2*j+1]);

  danot_(no);
}
