#ifndef _THOR_SCSI_ELEMENTS_MPOLE_H_
#define _THOR_SCSI_ELEMENTS_MPOLE_H_

#include <thor_scsi/elements/element_local_coordinates.h>
#include <thor_scsi/elements/enums.h>
#include <thor_scsi/elements/utils.h>
// move to API
#include <thor_scsi/core/multipoles.h>

/* Calculate multipole kick. The kick is given by

                   e L      L delta      L x
	theta  = - --- B  + -------  -  -----  ,
             x     p    y     rho           2
                    0                    rho
                 e L
	theta  = --- B
             y   p    x
                  0
    where
                           ====
                           \                       n-1
	(B + iB  ) = B rho  >   (ia  + b ) (x + iy)
	  y    x           /       n    n
	                   ====
    where
			e      1
			-- = -----
			p    B rho
			 0
*/
namespace thor_scsi::elements {
	/**
	 * Calculate multipole kick. The kick is given by
	 * @f[
	 *     \vartheta_x =  \frac{e}{p_\mathrm{0}} L B_y + \frac{L}{\rho}
	 *                  \left(\delta - \frac{x}{\rho} \right)
	 * @f]
	 * @f[
	 \vartheta_y =  \frac{e}{p_\mathrm{0}} L B_x
	 * @f]
	 * with
	 * @f[
	 * \frac{e}{p_\mathrm{0}} = \frac{1}{B\, \rho}
	 * @f]
	 *
	 * \verbatim embed:rst:leading-asterisk
	 *   The magnetic field expansion is represented by a
	 *   :any:`Field`Field2DInterpolation object:
	 * \endverbatim
	 */
	const int HOMmax = 21;
	typedef std::vector<double> MpoleArray;
	class MpoleType : public LocalGalileanPRot {
	public:

		inline MpoleType(const Config &config) : LocalGalileanPRot(config){
			// Field interpolation type
			intp = nullptr;
		}
		// friend MpoleType* Mpole_Alloc(void);
		//	ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
		// void print(const std::string &);
		/// std::string repr(void);
		///< additional information of mpole for mpole
		///< more than just solely type name
		///< used by meth:`repr`
		// std::string repr_add(void);

		const char* type_name(void) const override final { return "mpole"; };

		inline double GetKpar(int order){
			throw std::logic_error("Implement harmonic access ");
			// return PBpar[order+HOMmax];
		}

		inline void SetKpar(int order, double val){
			throw std::logic_error("Implement harmonic access ");
			// PBpar[order+HOMmax] = val;
		}


		virtual void localPass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps) override final
		{ _localPass(conf, ps);}
		// virtual void pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps) override final
		// { _pass(conf, ps);}

		///< Make smart pointer
		thor_scsi::core::Field2DInterpolation* intp;

	  private:
		template<typename T>
			void _localPass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
#if 0
		template<typename T>
			void Mpole_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
		void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps)
		{ Mpole_Pass(conf, ps); };
		void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps)
		{ Mpole_Pass(conf, ps); };
#endif
	  public:
#if 0

		/**
		 *  sets up the  calculation of the Euclidian group: translation
		 */
		void SetdS(void);
		void SetdT(void);
		/**
		 * triggers calculation of member PB from members PBrms, PBrnd, PBsys
		 */
		void SetPB(const int n);
		double GetdT(void);
		double GetPB(const int n);
		double
			// should be covered by transform
			PdTpar,                    ///< Roll angle [deg]: design
			PdTsys,                    //<                    systematic
			PdTrms,                    ///<                    rms
			PdTrnd;                    ///<                    random number.
		std::vector<double>
			PdSsys{0e0, 0e0},          ///< Displacement errors [m]: systematic
			PdSrms{0e0, 0e0},          ///<                          rms
			PdSrnd{0e0, 0e0};          ///<                          random number.
#endif
		int
			Pmethod,                   ///< Integration Method.
			PN,                        ///< Number of integration steps.
			Porder,                    ///< The highest order in PB.
			n_design;                  ///< multipole order (design).
		double
			PTx1,                      ///<  Bend angle [deg]:  hor. entrance angle
			PTx2,                      ///<  Bend angle [deg]: hor. exit angle.
			Pgap,                      ///< Total magnet gap [m].
			Pirho,                     ///< 1/rho [1/m].
			Pc0,                       ///< Corrections for roll error of bend.
			Pc1,
			Ps1;
#if 0
		/* todo: review to implement is a complex arrays */
		MpoleArray
			PBpar,                             ///< Multipole strengths: design
			PBsys,                     ///< Multipole strengths: systematic
			PBrms,                     ///< Multipole strengths: random variation \f$ <\sigma> \f$
			PBrnd,                     ///< random number used for computing PBrms.
			PB;                        ///< Multipole strengths: total used by integrator
#endif
		pthicktype
			Pthick;
		std::vector< std::vector<double> >
			M_elem                     ///< Transport matrix & orbit.
			{{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}};
	private:
		// SI units are used internally
		// apart from energy [GeV]
		/*  c_1 = 1/(2*(2-2^(1/3))),    c_2 = (1-2^(1/3))/(2*(2-2^(1/3)))
		    d_1 = 1/(2-2^(1/3)),        d_2 = -2^(1/3)/(2-2^(1/3))                 */

		/**
		 * Symplectic integrator
		 * 2nd order
		 *
		 *  @f[ \frac{L}{2} \to bnl \to  \frac{L}{2} ]@
		 *
		 * Here 4th order: calculating coefficients
		 *
		 * 4th order:
		 *
		 *  @f[ d_1 L \to c_1  bnL \to d_1 L \to c_2 bnl \to d_1 L ]@
		 *
		 * with
		 *
		 * @f[ d_1 + d_2 + d_1 = L @f]
		 *
		 *  @f[c_1 + c_2 = 1 @f]
		 *
		 * d_2:  negative thus creates a negative drift
		 */
		const double c_1 = 1e0/(2e0*(2e0-thirdroot(2e0)));
		const double c_2 = 0.5e0 - c_1;
		const double d_1 = 2e0*c_1;
		const double d_2 = 1e0 - 2e0*d_1;

	};
} // Name space

#endif // _THOR_SCSI_ELEMENTS_MPOLE_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
