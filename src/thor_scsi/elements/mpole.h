#ifndef _THOR_SCSI_ELEMENTS_MPOLE_H_
#define _THOR_SCSI_ELEMENTS_MPOLE_H_

#include <thor_scsi/elements/field_kick.h>
#include <thor_scsi/core/multipoles.h>

namespace thor_scsi::elements {
	class MpoleType : public FieldKick {
	public:
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
		inline MpoleType(const Config &config) : FieldKick(config){
			intp = new thor_scsi::core::PlanarMultipoles;
		}

		inline ~MpoleType(void){
			if(intp){
				delete intp;
			}
		}
		const char* type_name(void) const override { return "mpole"; };

		//thor_scsi::core::PlanarMultipoles* intp;
	};

} // Name space

#endif // _THOR_SCSI_ELEMENTS_MPOLE_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
