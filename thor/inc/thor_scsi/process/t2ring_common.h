#ifndef _THOR_SCSI_T2RING_COMMON_H_
#define _THOR_SCSI_T2RING_COMMON_H_
/**

   Todo:
       better naming?
 */
#include <tps/ss_vect.h>
#include <tps/tps_type.h>

void get_twoJ(const int n_DOF, const ss_vect<double> &ps,
	      const ss_vect<tps> &A, double twoJ[]);

#endif /*_THOR_SCSI_T2RING_COMMON_H_ */
