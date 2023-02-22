#ifndef _THOR_SCSI_ELEMENTS_CLASSICAL_MAGNET_H_
#define _THOR_SCSI_ELEMENTS_CLASSICAL_MAGNET_H_

#include <thor_scsi/elements/mpole.h>
#include <thor_scsi/core/multipoles.h>
#include <thor_scsi/core/multipole_types.h>

namespace thor_scsi::elements {
	/**
	 * @brief a magnet with a single strong multipole
	 *
	 * typically not directly used but using derived classes as Quadrupole, Sextupole or Bending
	 */
    template<class C>
	class ClassicalMagnetWithKnob : public MpoleTypeWithKnob<C> {
    protected:
        using complex_type = typename C::complex_type;
        using double_type = typename C::double_type;

        /*
         * no type name as not directly used
         */
	public:
		inline ClassicalMagnetWithKnob(const Config &config) : MpoleTypeWithKnob<C>(config){
			/*
			 * don't now how to access the multipole number derived by
			 *  a subclass
			 */
		}
		/**
		 * get the major harmonic number
		 */
		virtual int getMainMultipoleNumber(void) const = 0;
		/**
		 * is it a skew multipole ?
		 */
		virtual bool isSkew(void) const = 0;

		/**
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * .. Todo::
		 *    Check if double function is appropriate ...
		 *
		 *
		 * \endverbatim
		 */
		inline const std::shared_ptr<thor_scsi::core::TwoDimensionalMultipolesKnobbed<C>> getMultipoles(void) const {
			auto tmp = this->getFieldInterpolator();
			if(!tmp){
				std::stringstream strm;
				strm << __FILE__ << " ClassicalMAgnet.getMultipoles: this->getFieldInterpolator = NULL ";
				std::cerr << strm.str() << std::endl;
				abort();
				throw std::logic_error(strm.str());
			}
			return tmp;
		}

		inline void setMainMultipoleStrength(const Config &config){
			const double_type K = config.get<double>("K");
			// Watch the apersand ...
			this->setMainMultipoleStrength(K);
		}

		inline void setMainMultipoleStrength(const complex_type mul){
			const int n = this->getMainMultipoleNumber();
			this->getMultipoles()->setMultipole(n, mul);
		}

		inline void setMainMultipoleStrength(const Config &config, const int n){
			double K = config.get<double>("K");
			this->setMainMultipoleStrength(K);
		}

		inline void setMainMultipoleStrength(const double_type part){
			double_type re=0e0, im=0e0;
			if(!this->isSkew()){re = part; } else {	im = part;}
			const complex_type Cn(re, im);
			this->setMainMultipoleStrength(Cn);
		}

		inline complex_type getMainMultipoleStrength(void) const {
			auto n = this->getMainMultipoleNumber();
			return this->getMultipoles()->getMultipole(n);
		};

		inline double_type getMainMultipoleStrengthComponent(void) const {
			auto cm = this->getMainMultipoleStrength();
			if(this->isSkew()){
				return cm.imag();
			}
			return cm.real();
		}
	};

    typedef ClassicalMagnetWithKnob<thor_scsi::core::StandardDoubleType> ClassicalMagnet;
    typedef ClassicalMagnetWithKnob<thor_scsi::core::TpsaVariantType> ClassicalMagnetTpsa;
} // Name space

#endif // _THOR_SCSI_ELEMENTS_CLASSICAL_MAGNET_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
