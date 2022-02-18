#include "drift.h"
#include <thor_scsi/core/config.h>
#include <thor_scsi/core/element_helpers.h>
#include <tps/ss_vect.h>
#include <tps/tps.h>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

namespace thor_scsi::elements {
template<typename T>
void Drift(tsc::ConfigType &conf, const double L, ss_vect<T> &ps)
{
  T u;

  if (!conf.H_exact) {
    // Small angle axproximation.
    u = L/(1e0+ps[delta_]);
    ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
    ps[ct_] += u*(sqr(ps[px_])+sqr(ps[py_]))/(2e0*(1e0+ps[delta_]));
  } else {
	  u = L/tse::get_p_s(conf, ps);
    ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
    ps[ct_] += u*(1e0+ps[delta_]) - L;
  }
  if (conf.pathlength) ps[ct_] += L;
}
}

template<typename T>
void tse::DriftType::_advance(tsc::ConfigType &conf, ss_vect<T> &ps)
{
  Drift(conf, PL, ps);

  if (conf.emittance && !conf.Cavity_on)
    // Needs A^-1.
    curly_dH_x = is_tps<tps>::get_curly_H(ps);
}

template void tse::DriftType::_advance(tsc::ConfigType &conf, ss_vect<double> &ps);
template void tse::DriftType::_advance(tsc::ConfigType &conf, ss_vect<tps> &ps);


/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
