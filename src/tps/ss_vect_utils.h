#ifndef _TPS_SS_VECT_UTILS_H_
#define _TPS_SS_VECT_UTILS_H_ 1

#include <armadillo>
#include <tps/ss_vect.h>

inline arma::vec pstovec(const ss_vect<double> &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_], 1e0}; }

inline ss_vect<double> vectops(const arma::vec vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_]}; }

inline std::vector<double> vectostl(const arma::vec &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_], 1e0}; }

inline std::vector<double> stltovec(const arma::vec &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_], 1e0}; }

inline std::vector<double> pstostl(const ss_vect<double> &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_], 1e0}; }

inline ss_vect<double> stltops(const std::vector<double> &vec)
{ return {vec[x_], vec[px_], vec[y_], vec[py_], vec[delta_], vec[ct_]}; }

#endif /* _TPS_SS_VECT_UTILS_H_ */
