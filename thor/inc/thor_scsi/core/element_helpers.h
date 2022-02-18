#ifndef _THOR_SCSI_CORE_ELEMENTS_HELPERS_H_
#define _THOR_SCSI_CORE_ELEMENTS_HELPERS_H_ 1

/**
   Contains functions not part of the global API but rather
   code common to the elemetns.

 */
#include <vector>
#include <string>
#include <tps/ss_vect.h>
#include <tps/ss_vect_utils.h>
#include <tps/tps.h>
#include <thor_scsi/core/cells.h>
#include <thor_scsi/core/elements_enums.h>
#include <thor_scsi/core/config.h>

#include <exception>
#include <iostream>
#include <string>

#include <tps/tpsa_lin.h>
// #include <thor_scsi/process/t2ring_common.h>


extern double q_fluct; /// < Track down if it is a constant or a global variable reflecting current status

namespace thor_scsi {
	namespace elements {

void get_twoJ(const int n_DOF, const ss_vect<double> &ps,
	      const ss_vect<tps> &A, double twoJ[]);

/**
 *
 * Todo: check naming:
 *
 *  Function exist that are named
 *
 *
 *  get_p_s
 *  get_ps
 */

template<typename T>
inline T get_p_s(thor_scsi::core::ConfigType &conf, const ss_vect<T> &ps)
{
  T p_s, p_s2;

  if (!conf.H_exact)
    // Small angle axproximation.
    p_s = 1e0 + ps[delta_];
  else {
    p_s2 = sqr(1e0+ps[delta_]) - sqr(ps[px_]) - sqr(ps[py_]);
    if (p_s2 >= 0e0)
      p_s = sqrt(p_s2);
    else {
//      printf("get_p_s: *** Speed of light exceeded!\n");
      p_s = NAN;
    }
  }
  return(p_s);
}


// partial template-class specialization
// primary version
template<typename T>
class is_tps { };

// partial specialization
template<>
class is_tps<double> {
public:
	static inline void get_ps(const ss_vect<double> &x, thor_scsi::elements::CellType *Cell)
		{ Cell->BeamPos = pstostl(x); }

	static inline double set_prm(const int k) { return 1e0; }

	static inline double get_curly_H(const ss_vect<tps> &x){
		std::cout << "get_curly_H: operation not defined for double" << std::endl;
		throw std::domain_error("get_curly_H: operation not defined for double");
		return 0e0;
	}

	static inline double get_dI_eta(const ss_vect<tps> &A){
		std::cout << "get_dI_eta: operation not defined for double" << std::endl;
		throw std::domain_error("get_dI_eta: operation not defined for double");
		return 0e0;
	}

	static inline void emittance(thor_scsi::core::ConfigType &conf, const double B2,
				     const double u, const double ps0,
				     const ss_vect<double> &xp) { }

	static inline void diff_mat(const double B2, const double u,
				    const double ps0, const ss_vect<double> &xp) { }

};

// partial specialization
template<>
class is_tps<tps> {
public:
	static inline void get_ps(const ss_vect<tps> &x, thor_scsi::elements::CellType *Cell)
		{ Cell->BeamPos = pstostl(x.cst()); Cell->A = maptostlmat(x); }

	static inline tps set_prm(const int k) { return tps(0e0, k); }

	static inline double get_curly_H(const ss_vect<tps> &A){
		int             j;
		double          curly_H[2];
		ss_vect<double> eta;

		eta.zero();
		for (j = 0; j < 4; j++)
			eta[j] = A[j][delta_];

		get_twoJ(2, eta, A, curly_H);

		return curly_H[X_];
	}

	static inline double get_dI_eta(const ss_vect<tps> &A){
		return A[x_][delta_];
	}

	/**

	  see
	  M. Sands "The Physics of Electron Storage Rings" SLAC-121, Eq. (5.20),
	  p. 118:
	    dN<u^2>/E^2 =
	      3*C_U*C_gamma*h_bar*c*E_0^5*(1+delta)^4*(B_perp/(Brho))^3/(4*pi*m_e^3 [eV/c^2])

	  A contains the eigenvectors.
	 */
	static inline void emittance(thor_scsi::core::ConfigType &conf, const tps &B2_perp,
				     const tps &ds, const tps &p_s0,
				     const ss_vect<tps> &A){

		int          j;
		double       B_66;
		ss_vect<tps> A_inv;

		if (B2_perp > 0e0) {
			B_66 = (q_fluct*pow(B2_perp.cst(), 1.5)*pow(p_s0, 4)*ds).cst();
			A_inv = Inv(A);
			// D_11 = D_22 = curly_H_x,y * B_66 / 2,
			// curly_H_x,y = eta_Fl^2 + etap_Fl^2
			for (j = 0; j < 3; j++)
				conf.D_rad[j] +=
					(sqr(A_inv[j*2][delta_])+sqr(A_inv[j*2+1][delta_]))*B_66/2e0;
		}
	}

	static inline void diff_mat(const tps &B2_perp, const tps &ds,
				    const tps &p_s0, ss_vect<tps> &x)
		{ }

};

	}
}
#endif /*  _THOR_SCSI_CORE_ELEMENTS_HELPERS_H_  */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
