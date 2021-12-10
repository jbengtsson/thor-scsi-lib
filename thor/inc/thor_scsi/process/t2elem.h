/* Tracy-3

   J. Bengtsson, BNL 2007

*/

#ifndef THOR_SCSI_PROCESS_T2ELEM_H
#define THOR_SCSI_PROCESS_T2ELEM_H 1

#include <tps/ss_vect.h>
#include <tps/tps_type.h>
#include <vector>

extern double q_fluct;

namespace thor_scsi {

extern bool   sympl;
extern int    FieldMap_filetype;
extern double C_u, C_gamma, C_q, cl_rad,

template<typename T>
T get_p_s(ConfigType &conf, const ss_vect<T> &);

void getelem(long i, CellType *cellrec);
void putelem(long i, CellType *cellrec);

template<typename T>
void GtoL(ss_vect<T> &X, std::vector<double> &S, std::vector<double> &R,
	  const double c0, const double c1, const double s1);
template<typename T>
void LtoG(ss_vect<T> &X, std::vector<double> &S, std::vector<double> &R,
	  double c0, double c1, double s1);
template<typename T>
void p_rot(ConfigType &conf, double phi, ss_vect<T> &x);


template<typename T>
void get_B2(const double h_ref, const T B[], const ss_vect<T> &xp,
	    T &B2_perp, T &B2_par);

template<typename T>
void radiate(ConfigType &conf, ss_vect<T> &x, const double L,
	     const double h_ref, const T B[]);

template<typename T>
void Drift(ConfigType &conf, const double L, ss_vect<T> &x);
template<typename T>
void bend_fringe(ConfigType &conf, const double hb, ss_vect<T> &x);
template<typename T>
void EdgeFocus(ConfigType &conf, const double irho, double phi,
	       double gap, ss_vect<T> &x);
template<typename T>
void quad_fringe(ConfigType &conf, const double b2, ss_vect<T> &x);
template<typename T>
void thin_kick(ConfigType &conf, const int Order, const MpoleArray &MB,
	       const double L, const double h_bend, const double h_ref,
	       ss_vect<T> &x);
template<typename T>
void sol_pass(ConfigType &conf, const ElemType *elem, ss_vect<T> &x);


arma::mat get_edge_lin_map(const double h, const double phi,
			   const double gap, const double delta);
arma::mat get_sbend_lin_map(const double L, const double h, const double b2,
			    const double delta);
arma::mat get_thin_kick_lin_map(const double b2L, const double delta);
arma::mat get_lin_mat(ElemType *Elem, const double delta);
void get_lin_mats(const double delta);

double Elem_GetKval(int Fnum1, int Knum1, int Order);

};

#endif /*  THOR_SCSI_PROCESS_T2ELEM_H*/
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
