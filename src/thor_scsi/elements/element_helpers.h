#ifndef _THOR_SCSI_CORE_ELEMENTS_HELPERS_H_
#define _THOR_SCSI_CORE_ELEMENTS_HELPERS_H_ 1

/**
   Contains functions not part of the global API but rather
   code common to the elemetns.

 */
#include <vector>
#include <string>
// #include <tps/ss_vect.h>
// #include <tps/ss_vect_utils.h>
#include <tps/tps_type.h>
#include <gtpsa/ss_vect.h>
#include <thor_scsi/core/config.h>
#include <thor_scsi/core/exceptions.h>

#include <exception>
#include <iostream>
#include <string>

// #include <tps/tpsa_lin.h>
// #include <thor_scsi/process/t2ring_common.h>


extern double q_fluct; /// < Track down if it is a constant or a global variable reflecting current status

namespace thor_scsi::elements {

/*
 *
 *
 *                          2
 *         K1*gap*h*(1 + sin phi)
 *  psi = ----------------------- * (1 - K2*g*gap*tan phi)
 *              cos phi
 *
 *
 */
/**
 * @brief Correction for magnet gap (longitudinal fringe field)
 *
 *
 * irho = h,
 * @f[
 *  h = \frac{1}{\rho} [1/m]
 * @f]
 * @f[\phi@f]  edge angle
 * gap  full gap between poles
 *
 *  @f[K_1 @f] is usually 1/2,
 *  @f[K_2 @f] is zero here
 *
 *
 * @f[
 *     \psi = \frac{K_1 \ gap \  h (1 + \left[sin(phi)\right]^2)}{\cos{(\phi)}}
 *             * (1 - K_2 \ g  \ gap \tan(\phi))
 * @f]
 *
 * \verbatim embed:rst:leading-asterisk
 * .. Warning::
 *
 *    phi: value expected in degree!
 * \endverbatim
 */
	double get_psi(const double irho, const double phi, const double gap);

/**
 * @brief compute the linear action J
 *
 * \verbatim embed:rst:leading-asterisk
 *
 * .. Todo::
 *
 *      Check interface. See if phi and J should be computed at the same time
 *      Check return double (e.g. std::array)
 *
 * \endverbatim
 *
 */
template<typename T>
void get_twoJ(const int n_DOF, const gtpsa::ss_vect<double> &ps,
	      const gtpsa::ss_vect<T> &A, double twoJ[]);

/**
 *
 * Todo: check naming:
 *
 *  Function exist that are named
 *
 *  Compute longitudinal momentum
 *  get_p_s
 *  get_ps
 */

/**
 *  @brief Compute longitudinal momentum
 */
template<typename T>
inline T get_p_s(const thor_scsi::core::ConfigType &conf, const gtpsa::ss_vect<T> &ps)
{
	T p_s(ps[0]), p_s2(ps[0]);

	if (!conf.H_exact) {
		// Small angle axproximation.
		p_s = 1e0 + ps[delta_];
	} else {
		p_s2 = sqr(1e0+ps[delta_]) - sqr(ps[px_]) - sqr(ps[py_]);
		if (p_s2 >= 0e0){
			p_s = sqrt(p_s2);
		}  else {
			throw PhysicsViolation("Speed of light exceeded");
			p_s = NAN;
		}
	}
	return(p_s);
}


	template<typename T>
	double get_curly_H(const gtpsa::ss_vect<T>      &A);
	double get_curly_H(const gtpsa::ss_vect<tps>    &A);
	double get_curly_H(const gtpsa::ss_vect<double> &A);


	//template<>
	template<typename T>
	static inline double get_dI_eta(const gtpsa::ss_vect<T> &A){
#warning "optimise dI eta"
		arma::mat jac = A.jacobian();
		return jac(x_, delta_);
	}

	static inline double get_dI_eta(const gtpsa::ss_vect<tps> &A){
		return A[x_][delta_];
	}

	static inline double get_dI_eta(const gtpsa::ss_vect<double> &A){
		std::cout << "get_dI_eta: operation not defined for double" << std::endl;
		throw std::domain_error("get_dI_eta: operation not defined for tps");
		return 0e0;
	}

// partial template-class specialization
// primary version
template<typename T>
class is_tps { };

/**
 * partial specialization
 *
 * Todo: check if moved to oberver
 */
template<>
class is_tps<double> {
public:
	/**
	 *  @brief Compute phase space
	 */
	/*
	static inline void get_ps(const ss_vect<double> &x, thor_scsi::core::CellType *Cell)
		{ Cell->BeamPos = pstostl(x); }
	*/
	static inline double set_prm(const int k) { return 1e0; }

	/**
	 * @brief Compute (linear) dispersion action
	 */
	/**
	 * @brief Synchrotron integral
	 *
	 * Todo: Check which one
	 */
	static inline double get_dI_eta(const gtpsa::ss_vect<tps> &A){
		std::cout << "get_dI_eta: operation not defined for double" << std::endl;
		throw std::domain_error("get_dI_eta: operation not defined for tps");
		return 0e0;
	}
	static inline double get_dI_eta(const gtpsa::ss_vect<gtpsa::tpsa> &A){
		std::cout << "get_dI_eta: operation not defined for gtpsa::tpsa" << std::endl;
		throw std::domain_error("get_dI_eta: operation not defined for gtpsa::tpsa");
		return 0e0;
	}


	static inline void diff_mat(const double B2, const double u,
				    const double ps0, const gtpsa::ss_vect<double> &xp) { }

};

// partial specialization
template<>
class is_tps<tps> {
public:
	/*
	static inline void get_ps(const ss_vect<tps> &x, thor_scsi::core::CellType *Cell)
		{ Cell->BeamPos = pstostl(x.cst()); Cell->A = maptostlmat(x); }
	*/

	static inline tps set_prm(const int k) { return tps(0e0, k); }

	/*
	 *  dN<u^2>/E^2 =
	 * 3*C_U*C_gamma*h_bar*c*E_0^5*(1+delta)^4*(B_perp/(Brho))^3/(4*pi*m_e^3 [eV/c^2])
	 */
	/**
	 *
	 *
	 * @f[
	 *  \frac{dN<u^2>}{E^2} =
	 *    3 \ C_U \ C_\gamma \hbar \ c \sqrt{E_0} (1+\delta)^4
	 *       \frac{\left(\frac{B_{perp}}{B\rho}\right)^3}{(4*pi*m_e^3 [eV/c^2])}
	 * @f]
	 * A contains the eigenvectors.
	 *
	 * \verbatim embed:rst:leading-asterisk
	 * see [Sands]_
	 *
	 * .. [Sands] M. Sands "The Physics of Electron Storage Rings" SLAC-121, Eq. (5.20),
	 *            p. 118
	 *
	 * \endverbatim
	 */

	static inline void diff_mat(const tps &B2_perp, const tps &ds,
				    const tps &p_s0, gtpsa::ss_vect<tps> &x)
		{ }


};
		/**
		 * @brief forward a phase space with a drift
		 *
		 *
		 * Required as helper for drift -- kick -- drift implementationa
		 * See drift.h for DriftType, which is an implementation as separate lattice lement
		 */
		template<typename T, typename T2>
		void drift_propagate(const thor_scsi::core::ConfigType &conf, const T2& L, gtpsa::ss_vect<T> &ps);

		/**
		 * @brief implementation of the thin kick (thin lens approximation)
		 *
		 *   @param h_bend: 1/rho_bend (radius of the the bend, magnetic rigidty of the dipoles)
		 *   @param h_ref:  1/rho which is the curvature of the design / reference orbit
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * .. Todo::
		 *
		 *     * Split up function in different parts or functions:
		 *       E.g.: one for dipoles and one for anything else.
		 *     * Consider to combine interpolation with length L and thus not using field interpolation
		 *       but field integral interpolation
		 *     * consider renaming function to thin lens
		 *
		 *
		 * .. Warning::
		 *
		 *    * cartesion bend with evaluation of gradients not supported curently (code not
		 *      active
		 *
		 * The vector potential for the combined-function sector bend is from: [Iselin:85]_
		 *
		 * .. [Iselin:85]  C. Iselin "Lie Transformations and Transport Equations for Combined-
		 *                 Function Dipoles" Part. Accel. 17, 143-155 (1985).
		 * \endverbatim
		 */
		template<typename T>
		void thin_kick(const thor_scsi::core::ConfigType &conf,
			       // const thor_scsi::core::Field2DInterpolation& intp,
			       const T BxoBrho, const T ByoBrho,
			       const double L,
			       const double h_bend, const double h_ref,
			       const gtpsa::ss_vect<T> &ps0, gtpsa::ss_vect<T> &ps);

}// namespace thor_scsi::elements
#endif /*  _THOR_SCSI_CORE_ELEMENTS_HELPERS_H_  */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
