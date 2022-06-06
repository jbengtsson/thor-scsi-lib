#ifndef _THOR_SCSI_ELEMENTS_MPOLE_H_
#define _THOR_SCSI_ELEMENTS_MPOLE_H_

#include <thor_scsi/elements/field_kick.h>
#include <thor_scsi/core/multipoles.h>
#include <iomanip>

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
			std::shared_ptr<thor_scsi::core::TwoDimensionalMultipoles> tmp(new thor_scsi::core::TwoDimensionalMultipoles);
			this->intp = std::move(std::dynamic_pointer_cast<thor_scsi::core::Field2DInterpolation>(tmp));
			/*
			  std::cerr << "Mpole Type: " << std::setw(10) << name << " interpolation object " << this->intp  << std::endl;
			  auto parent = dynamic_cast<FieldKick*>(this);
			  std::cerr << "FieldKick:             interpolation object " << parent->getFieldInterpolator()  << std::endl;
			*/
		}

		const char* type_name(void) const override { return "mpole"; };

		inline auto getFieldInterpolator(void) const {
			auto p = std::dynamic_pointer_cast<thor_scsi::core::TwoDimensionalMultipoles>(this->intp);
			return std::const_pointer_cast<thor_scsi::core::TwoDimensionalMultipoles>(p);
		}

		inline void setFieldInterpolator(std::shared_ptr<thor_scsi::core::TwoDimensionalMultipoles> intp) {
			auto p = std::dynamic_pointer_cast<thor_scsi::core::Field2DInterpolation>(intp);
			if(!p){
				throw std::runtime_error("Could not cast multipole to Field2DInterpolation");
			}
			this->intp = p;
		}

		//thor_scsi::core::TwoDimensionalMultipoles* intp;
	};

} // Name space

#endif // _THOR_SCSI_ELEMENTS_MPOLE_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
