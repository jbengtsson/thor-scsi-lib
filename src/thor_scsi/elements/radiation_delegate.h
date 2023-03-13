#ifndef _THOR_SCSI_RADIATION_DELEGATE_H_
#define _THOR_SCSI_RADIATION_DELEGATE_H_ 1
#include <thor_scsi/elements/radiation_delegate_api.h>
#include <thor_scsi/elements/constants.h>
#include <tps/tpsa_lin.h>
#include <array>

namespace thor_scsi::elements {
	using thor_scsi::core::ElemType;

	template<class EC>
	class RadiationDelegateKnobbed: public  RadiationDelegateInterfaceKnobbed<EC>{
	public:
		~RadiationDelegateKnobbed(void){};
		inline RadiationDelegateKnobbed(void)
			: curly_dH_x(0e0)
			, delegator_name("")
			, delegator_index(-1)
			{ this->reset(); }

		inline void reset(void) {
			this->curly_dH_x = 0e0;
		}

		inline double getCurlydHx(void) const {
			return this->curly_dH_x;
		}
		/*
		 * Used for computing curly_dHx
		 */
		virtual void view(const EC& kick, const gtpsa::ss_vect<double>      &ps, const enum ObservedState state, const int cnt) override;
		// virtual void view(const ElemType& kick, const gtpsa::ss_vect<tps>         &ps, const enum ObservedState state, const int cnt) override;
		virtual void view(const EC& kick, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum ObservedState state, const int cnt) override;

		virtual void show(std::ostream& strm, int level) const override final;

		std::string getDelegatorName(void){
			return this->delegator_name;
		}
		int getDelegatorIndex(void){
			return this->delegator_index;
		}
	private:
		template<typename T>
		inline void _view(const EC&, const gtpsa::ss_vect<T> &ps, const enum ObservedState state, const int cnt);

		template<typename T>
		inline void computeAndStoreCurlyH(const gtpsa::ss_vect<T> &ps);

		double curly_dH_x;
		std::string delegator_name;
		int delegator_index;
	};

	typedef RadiationDelegateKickInterfaceKnobbed<ElemType> RadiationDelegateKickInterface;
	/**
	 *
	 * @note synchrotron integrals are not used nor has their
	 *       functionality been checked
	 */

	/*
	 * EC typically is FieldKickAPIKnobbed<C>
	 */
	template<class FC>
	class RadiationDelegateKickKnobbed: public RadiationDelegateKickInterfaceKnobbed<FC> {
	public:
		RadiationDelegateKickKnobbed(void)
			: curly_dH_x(0e0)
			, index(0e0)
			, dI( {0e0, 0e0, 0e0, 0e0, 0e0, 0e0} )
			, D_rad( {0e0, 0e0, 0e0} )
			{
			this->reset();
		}
		~RadiationDelegateKickKnobbed(void){};

		/*
		 * @brief: reset parameters for radiation
		 *
		 * passing ps to template the function. currently these functions are void
		 * (perhaps required later on)
		 */
		inline void reset(void) {
			this->curly_dH_x = 0e0;
			this->dI.fill({0e0});
			this->D_rad.fill({0e0});
			// this->dEnergy = 0e0;

		}

		// should be renamed to diffusion coefficient ...
		// add in coment that it is a prerequisite for emittance calculations
		/*
		 * @brief computing diffusion coefficients
		 * @todo review if it would be a better approach to store B66
		 */
		inline void computeDiffusion(const bool flag){
			this->compute_diffusion = flag;
		}
		inline bool isComputingDiffusion() const {
			return this->compute_diffusion;
		}

		// See if not link to a global machine setting property
		void setEnergy(const double val);

		inline double getEnergy(void) const {
			return this->energy;
		}

		/*
		 * Used for computing synchrotron integrals
		 */
		virtual void view(const FC& kick, const gtpsa::ss_vect<double>      &ps, const enum ObservedState state, const int cnt) override;
		// virtual void view(const FieldKickAPI& kick, const gtpsa::ss_vect<tps>         &ps, const enum ObservedState state, const int cnt) override;
		virtual void view(const FC& kick, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum ObservedState state, const int cnt) override;

		virtual void show(std::ostream& strm, int level) const override final;
		//virtual void view(const ElemType& kick, const gtpsa::ss_vect<double> &ps, const enum ObservedState state, const int cnt) override final;
		//virtual void view(const ElemType& kick, const gtpsa::ss_vect<tps> &ps, const enum ObservedState state, const int cnt) override final;
		/**
		 * @brief Radiation effect due to local field
		 *
		 * @f[
		 *     \frac{\mathrm{d}\delta}{\mathrm{d}(ds)} =
		 *           -C_{\gamma} \, E_0^3 \, (1 + \delta)^2 \, \left( \frac{B_{perp}}{B \rho}\right)^2
		 *                \frac{1}{2\pi}
		 * @f]
		 *
		 * @todo: depends on energy of ring .... currently taken from config ...
		 *
		 * M. Sands "The hysics of Electron Storage Rings" SLAC-121, p. 98.
		 */
		template<typename T>
		void radiate(const thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<T> &x, const double L, const double h_ref, const std::array<T, 3> B);

		inline auto getSynchrotronIntegralsIncrement(void) const {
			return this->dI;
		}

		inline auto getDiffusionCoefficientsIncrement(void) const {
			return this->D_rad;
		}

		inline double getCurlydHx(void) const {
			return this->curly_dH_x;
		}

		std::string getDelegatorName(void) const {
			return this->delegator_name;
		}

		int getDelegatorIndex(void) const {
			return this->delegator_index;
		}

	private:
		//using RadiationDelegateKickInterface::view;

		/*
		 * @brief: finalise calculation of radiation integrals
		 *
		 * Applies appropriate balance to curly_H and integral 4.
		 * Calculates all other integrals
		 */
		//
		template<typename T>
		inline void synchrotronIntegralsFinish(const FC &kick, const gtpsa::ss_vect<T> &ps);

		// calculate the effect of radiation
		template<typename T>
		inline void //radiate
		synchrotronIntegralsStep(const gtpsa::ss_vect<T> &ps);

		template<typename T>
		inline void _view(const FC&, const gtpsa::ss_vect<T> &ps, const enum ObservedState state, const int cnt);

		template<typename T>
		void diffusion(const T &B2_perp, const T &ds, const T &p_s0,  const gtpsa::ss_vect<T> &A);

		double curly_dH_x;
		int index;
		std::array<double, 6> dI;           ///< Local contributions to the synchrotron integrals
		std::array<double, 3> D_rad;        //< Diffusion coefficients (Floquet space).
		bool compute_diffusion = false;
		double energy = NAN;
		double q_fluct = NAN;

		std::string delegator_name = "";
		int delegator_index = -1;

	};

    typedef RadiationDelegateKnobbed<thor_scsi::core::ElemTypeKnobbed> RadiationDelegate;
    typedef RadiationDelegateKickKnobbed<thor_scsi::elements::FieldKickAPIKnobbed<thor_scsi::core::StandardDoubleType>> RadiationDelegateKick;
} // namespace thor_scsi::elements

#endif /* _THOR_SCSI_RADIATION_DELEGATOR_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
