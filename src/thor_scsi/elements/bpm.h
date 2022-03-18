#ifndef _THOR_SCSI_ELEMENTS_BPM_H_
#define _THOR_SCSI_ELEMENTS_BPM_H_

#include <thor_scsi/elements/marker.h>

namespace thor_scsi::elements {
	/**
	 * @brief: BPM: derived from marker as lattice file contains two names
	 *
	 * No differnce from marker currently
	 */
	class BPMType : public MarkerType {
	public:
		inline BPMType(const Config &config) : MarkerType(config) {
		}

		const char* type_name(void) const override final { return "BPM"; };

	};

} // Name space

#endif // _THOR_SCSI_ELEMENTS_BPM_H_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
