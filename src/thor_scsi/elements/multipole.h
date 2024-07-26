#ifndef _THOR_SCSI_ELEMENTS_MULTIPOLE_H_
#define _THOR_SCSI_ELEMENTS_MULTIPOLE_H_

#include <thor_scsi/elements/classical_magnet.h>

namespace thor_scsi::elements {

	template<class C>
	class MultipoleTypeWithKnob : public ClassicalMagnetWithKnob<C> {
	public:
		inline MultipoleTypeWithKnob(const Config &config)
			: ClassicalMagnetWithKnob<C> (config)
			{
				this->setMainMultipoleStrength(config);
				const double b_3 = config.get<double>("B_3");
				this->getMultipoles()->setMultipole(3, b_3);
			}

		inline int getMainMultipoleNumber(void) const override final {
			return 2;
		};
		inline bool isSkew(void) const override final {
			return false;
		};

		const char* type_name(void) const override final
			{ return "Multipole"; };
	};
	typedef MultipoleTypeWithKnob<thor_scsi::core::StandardDoubleType>
        MultipoleType;
	typedef MultipoleTypeWithKnob<thor_scsi::core::TpsaVariantType>
	MultipoleTypeTpsa;
} // Name space

#endif // _THOR_SCSI_ELEMENTS_MULTIPOLE_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
