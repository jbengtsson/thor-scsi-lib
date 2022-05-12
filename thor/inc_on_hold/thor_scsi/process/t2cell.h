/* Tracy-3

   J. Bengtsson, BNL 2007

*/

#ifndef THOR_SCSI_T2CELL_H
#define THOR_SCSI_T2CELL_H

#include <tps/ss_vect.h>
#include <tps/tps_type.h>

extern tps  sigma_;

bool GetCOD(long imax, double eps, double dP, long &lastpos);

// template<typename T>
// void Elem_Pass(ss_vect<T> &x);

template<typename T>
void Cell_Pass(const long i0, const long i1, ss_vect<T> &x, long &lastpos);

void Cell_Pass(const long i0, const long i1, tps &sigma, long &lastpos);

#endif /* THOR_SCSI_T2CELL_H */
