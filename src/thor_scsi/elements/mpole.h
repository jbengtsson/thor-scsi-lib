#ifndef _THOR_SCSI_ELEMENTS_MPOLE_H_
#define _THOR_SCSI_ELEMENTS_MPOLE_H_

#include <thor_scsi/elements/field_kick.h>
#include <thor_scsi/core/multipoles.h>
#include <iomanip>

namespace thor_scsi::elements {
    template<class C, typename = typename C::complex_type, typename = typename C::double_type> // CK for complex knob
	class MpoleTypeWithKnob : public FieldKickKnobbed<C> {
    protected:
        using complex_type = typename C::complex_type;
        using double_type = typename C::double_type;
        using multipoles_knobbed = thor_scsi::core::TwoDimensionalMultipolesKnobbed<C>;
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
		inline MpoleTypeWithKnob(const Config &config) : FieldKickKnobbed<C>(config){
			auto tmp = std::make_shared<multipoles_knobbed>(0e0);
			this->intp = std::move(std::dynamic_pointer_cast<thor_scsi::core::Field2DInterpolationKnobbed<C>>(tmp));
			/*
			  std::cerr << "Mpole Type: " << std::setw(10) << name << " interpolation object " << this->intp  << std::endl;
			  auto parent = dynamic_cast<FieldKick*>(this);
			  std::cerr << "FieldKick:             interpolation object " << parent->getFieldInterpolator()  << std::endl;
			*/
		}

		const char* type_name(void) const override { return "mpole"; };

		inline auto getFieldInterpolator(void) const {
			std::stringstream strm;
			strm << __FILE__ << ":mpole.getFieldInterpolator: this->intp " << static_cast<void *>(this->intp.get());
			auto p = std::dynamic_pointer_cast<multipoles_knobbed>(this->intp);
			strm << " p " << static_cast<void *>(p.get());
			if(!p){
				strm<<" cast failed!";
				throw std::runtime_error(strm.str());
			}
			auto cp =  std::const_pointer_cast<multipoles_knobbed>(p);
			strm << " cp " << static_cast<void *>(cp.get());
			if(!cp){
				strm<<" cast failed!";
				throw std::runtime_error(strm.str());
			}
					    // std::cerr << strm.str() << std::endl;
			return cp;
		}

		inline void setFieldInterpolator(std::shared_ptr<thor_scsi::core::TwoDimensionalMultipolesKnobbed<C>> intp) {
			auto p = std::dynamic_pointer_cast<thor_scsi::core::Field2DInterpolationKnobbed<C>>(intp);
			if(!p){
				throw std::runtime_error("Could not cast multipole to Field2DInterpolation");
			}
			this->intp = p;
		}

		//thor_scsi::core::TwoDimensionalMultipoles* intp;
	};


    typedef class MpoleTypeWithKnob<core::StandardDoubleType> MpoleType;
    typedef class MpoleTypeWithKnob<core::TpsaVariantType> MpoleTypeTpsa;

} // Name space

#endif // _THOR_SCSI_ELEMENTS_MPOLE_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
