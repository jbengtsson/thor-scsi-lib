#ifndef _TPS_ENUMS_H_
#define _TPS_ENUMS_H_ 1

#include <tps/config.h>
const int
  max_str   = 132,
  n_m2      = 21,       ///< No of 1st & 2nd order moments: (6 over 5) + 6.
  ps_tr_dim = 4,        ///< Transverse phase-space dim.
/*
 *
 */
  ps_dim    = 6,        ///< 6D phase-space dim.
#ifdef NO_TPSA
  ss_dim    = ps_dim,
#else
#warning "unchecked code"
  ss_dim    = ps_dim +1,
#endif
  tps_n     = ps_dim+1; ///< 6D phase-space linear terms & constant

// Spatial components.
enum spatial_ind { X_ = 0, Y_ = 1, Z_ = 2 };

// Phase-space components.
// (Note, e.g. spin components should be added here)
enum phase_space_ind { x_ = 0, px_ = 1, y_ = 2, py_ = 3, delta_ = 4, ct_ = 5 };

// Truncated Power Series Algebra (TPSA)
const int
  nv_tps   = ps_dim,    ///< No of variables.
  nd_tps   = 3,         ///< No of degrees of freedom.
  iref_tps = 0;         /* File with resonances to be excluded from the map
			   normal form: fort.7. */


// Truncated Power Series Algebra (TPSA)
const  int  no_tps = NO_TPSA, ndpt_tps = 5;



#endif
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
