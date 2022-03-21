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
			auto* muls = new thor_scsi::core::PlanarMultipoles;
			intp.reset(muls);
		}

		const char* type_name(void) const override { return "mpole"; };

		inline const auto getFieldInterpolator(void) const {
			return std::static_pointer_cast<thor_scsi::core::PlanarMultipoles>(this->intp);
		}


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
