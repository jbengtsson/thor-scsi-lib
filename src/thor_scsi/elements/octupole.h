#ifndef _THOR_SCSI_ELEMENTS_OCTUPOLE_H_
#define _THOR_SCSI_ELEMENTS_OCTUPOLE_H_

#include <thor_scsi/elements/classical_magnet.h>

namespace thor_scsi::elements {

	template<class C>
	class OctupoleTypeWithKnob : public ClassicalMagnetWithKnob<C> {
	public:
		inline OctupoleTypeWithKnob(const Config &config)
			: ClassicalMagnetWithKnob<C> (config)
			{
				this->setMainMultipoleStrength(config);
				const double
					b_4 = config.get<double>("B_4", 0.0);
				this->getMultipoles()->setMultipole(4, b_4);
			}

		inline int getMainMultipoleNumber(void) const override final {
			return 2;
		};
		inline bool isSkew(void) const override final {
			return false;
		};

		const char* type_name(void) const override final
			{ return "Octupole"; };

	};

	typedef OctupoleTypeWithKnob<thor_scsi::core::StandardDoubleType>
	OctupoleType;
	typedef OctupoleTypeWithKnob<thor_scsi::core::TpsaVariantType>
	OctupoleTypeTpsa;
} // Name space

#endif // _THOR_SCSI_ELEMENTS_OCTUPOLE_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
