#ifndef _THOR_SCSI_ELEMENTS_CLASSICAL_MAGNET_H_
#define _THOR_SCSI_ELEMENTS_CLASSICAL_MAGNET_H_

#include <thor_scsi/elements/mpole.h>
#include <thor_scsi/core/multipoles.h>

namespace thor_scsi::elements {
	/**
	 * @brief a magnet with a single strong multipole
	 *
	 * typically not directly used
	 */
	class ClassicalMagnet : public MpoleType {
		/*
		 * no type name as not directly used
		 */
	public:
		inline ClassicalMagnet(const Config &config) : MpoleType(config){
			/*
			 * don't now how to access the multipole number derived by
			 *  a subclass
			 */
		}


		inline void setMainMultipoleStrength(const Config &config, const int n){
			auto muls = dynamic_cast<thor_scsi::core::PlanarMultipoles &> (*this->intp);
			muls.setMultipole(n, config.get<double>("K"));
		}

		inline void setMainMultipoleStrength(const Config &config){
			const int n = this->getMainMultipoleNumber();
			const double K = config.get<double>("K");
			// Watch the apersand ...
			auto& muls = dynamic_cast<thor_scsi::core::PlanarMultipoles &> (*this->intp);
			muls.setMultipole(n, K);
		}

		inline void setMainMultipoleStrength(const thor_scsi::core::cdbl mul){
			const int n = this->getMainMultipoleNumber();
			// Watch the apersand ...
			auto& muls = dynamic_cast<thor_scsi::core::PlanarMultipoles &> (*this->intp);
			muls.setMultipole(n, mul);
		}

		inline void setMainMultipoleStrength(const double normal){
			this->setMainMultipoleStrength(thor_scsi::core::cdbl(normal, 0));
		}
		/**
		 * get the major harmonic number
		 */
		inline thor_scsi::core::cdbl getMainMultipoleStrength(void) const {
			auto& muls = dynamic_cast<thor_scsi::core::PlanarMultipoles &> (*this->intp);
			return muls.getMultipole(this->getMainMultipoleNumber());
		};

		virtual int getMainMultipoleNumber(void) const = 0;
		//thor_scsi::core::PlanarMultipoles* intp;

		/**
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * .. Warning::
		 *
		 *     fix memory handling of interpolation object
		 *     add move constructor
		 *
		 * \endverbatim
		 */
		inline const thor_scsi::core::PlanarMultipoles& getMultipoles(void) const {
			return static_cast<thor_scsi::core::PlanarMultipoles&>(*this->intp);
		}

	};


} // Name space

#endif // _THOR_SCSI_ELEMENTS_CLASSICAL_MAGNET_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
