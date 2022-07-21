#ifndef _THOR_SCSI_ELEMENTS_OCTUPOLE_H_
#define _THOR_SCSI_ELEMENTS_OCTUPOLE_H_

#include <thor_scsi/elements/classical_magnet.h>

namespace thor_scsi::elements {
	class OctupoleType : public ClassicalMagnet {
	public:
		inline OctupoleType(const Config &config) : ClassicalMagnet(config){
			this->setMainMultipoleStrength(config);
		}

		inline int getMainMultipoleNumber(void) const override final {
			return 4;
		};
		inline bool isSkew(void) const override final {
			return false;
		};
		const char* type_name(void) const override final { return "Octupole"; };
	};
} // Name space

#endif // _THOR_SCSI_ELEMENTS_OCTUPOLE_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
