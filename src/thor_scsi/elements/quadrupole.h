#ifndef _THOR_SCSI_ELEMENTS_QUADRUPOLE_H_
#define _THOR_SCSI_ELEMENTS_QUADRUPOLE_H_

#include <thor_scsi/elements/classical_magnet.h>

namespace thor_scsi::elements {
	class QuadrupoleType : public ClassicalMagnet {
	public:
		inline QuadrupoleType(const Config &config) : ClassicalMagnet(config){
			this->setMainMultipoleStrength(config);
		}

		inline int getMainMultipoleNumber(void) const override final {
			return 2;
		};
		const char* type_name(void) const override final { return "Quadrupole"; };
	};

} // Name space

#endif // _THOR_SCSI_ELEMENTS_QUADRUPOLE_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
