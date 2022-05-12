#ifndef _THOR_SCSI_ELEMENTS_SEXTUPOLE_H_
#define _THOR_SCSI_ELEMENTS_SEXTUPOLE_H_

#include <thor_scsi/elements/classical_magnet.h>

namespace thor_scsi::elements {
	class SextupoleType : public ClassicalMagnet {
	public:
		inline SextupoleType(const Config &config) : ClassicalMagnet(config){
			this->setMainMultipoleStrength(config);
		}

		inline int getMainMultipoleNumber(void) const override final {
			return 3;
		};
		const char* type_name(void) const override final { return "Sextupole"; };
	};
} // Name space

#endif // _THOR_SCSI_ELEMENTS_SEXTUPOLE_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
