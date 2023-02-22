#ifndef _THOR_SCSI_ELEMENTS_SEXTUPOLE_H_
#define _THOR_SCSI_ELEMENTS_SEXTUPOLE_H_

#include <thor_scsi/elements/classical_magnet.h>

namespace thor_scsi::elements {

	template<class C>
	class SextupoleTypeWithKnob : public ClassicalMagnetWithKnob<C> {
	public:
		inline SextupoleTypeWithKnob(const Config &config)
			: ClassicalMagnetWithKnob<C> (config)
			{ this->setMainMultipoleStrength(config); }

		inline int getMainMultipoleNumber(void) const override final {
			return 3;
		};
		inline bool isSkew(void) const override final {
			return false;
		};

		const char* type_name(void) const override final { return "Sextupole"; };
	};
	typedef SextupoleTypeWithKnob<thor_scsi::core::StandardDoubleType> SextupoleType;
	typedef SextupoleTypeWithKnob<thor_scsi::core::TpsaVariantType> SextupoleTypeTpsa;
} // Name space

#endif // _THOR_SCSI_ELEMENTS_SEXTUPOLE_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
