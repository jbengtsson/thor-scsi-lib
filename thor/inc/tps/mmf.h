#ifndef _TPS_MMF_H_
#define _TPS_MMF_H_ 1

#include <tps/forward_decl.h>
#include <tps/ss_vect.h>

typedef struct MNF_struct
{
  tps
    K,       // New effective Hamiltonian.
    g;       // Generator for nonlinear transformation to Floquet space.
  ss_vect<tps>
    A0,      // Transformation to fixed point.
    A1,      // Transformation (linear) to Floquet space.
    map_res; // Residual map.
} MNF_struct;
#endif /* _TPS_MMF_H_ */
