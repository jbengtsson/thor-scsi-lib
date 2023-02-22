#ifndef _THOR_SCSI_ELEMENTS_QUADRUPOLE_H_
#define _THOR_SCSI_ELEMENTS_QUADRUPOLE_H_

#include <thor_scsi/elements/classical_magnet.h>

namespace thor_scsi::elements {
    template<class C>
	class QuadrupoleTypeWithKnob : public ClassicalMagnetWithKnob<C> {
	public:
		inline QuadrupoleTypeWithKnob(const Config &config)
           : ClassicalMagnetWithKnob<C>(config)
           {
			    this->setMainMultipoleStrength(config);
           }

		inline int getMainMultipoleNumber(void) const override final {
			return 2;
		};
		inline bool isSkew(void) const override final {
			return false;
		};
		const char* type_name(void) const override final { return "Quadrupole"; };
	};

    typedef QuadrupoleTypeWithKnob<thor_scsi::core::StandardDoubleType> QuadrupoleType;
} // Name space

#endif // _THOR_SCSI_ELEMENTS_QUADRUPOLE_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
