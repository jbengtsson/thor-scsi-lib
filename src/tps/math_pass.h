#ifndef _TPS_MATH_PATH_H_
#define _TPS_MATH_PATH_H_ 1

/**
   Todo:
       I assume this is only valid for the linear configuration
       NO_TPSA 1
*/
#include <tps/ss_vect_utils.h>
#include <vector>

inline ss_vect<double> mat_pass(std::vector< std::vector<double> > &M,
				ss_vect<double> &ps)
{ return vectops(stlmattomat(M)*pstovec(ps)); }

inline ss_vect<tps> mat_pass(std::vector< std::vector<double> > &M,
			     ss_vect<tps> &ps)
{ return mattomap(stlmattomat(M)*maptomat(ps)); }

#endif
