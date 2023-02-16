#ifndef _THOR_SCSI_ELEMENTS_FIELD_KICK_API_H_
#define _THOR_SCSI_ELEMENTS_FIELD_KICK_API_H_
#include <thor_scsi/elements/element_local_coordinates.h>
#include <thor_scsi/core/field_interpolation.h>
#include <thor_scsi/core/multipole_types.h>

namespace thor_scsi::elements {
	// declaration of API ... to be used by radiation observer
    template<class C>
	class FieldKickAPIKnobbed : public LocalGalileanPRotKnobbed<C> {

	public:
		inline FieldKickAPIKnobbed(const Config &conf)
			: LocalGalileanPRotKnobbed<C>(conf)
			, intp(nullptr)
			{}
		inline FieldKickAPIKnobbed(FieldKickAPIKnobbed&& O)
			: LocalGalileanPRotKnobbed<C>(std::move(O)),
			  intp(std::move(O.intp))
			{
			}
		/*
		 * @brief:
		 *
		 * returns
		 * @f[ \frac{1}{\rho} = irho == 0 @f]
		 */

		inline double getCurvature(void) const {
			return this->Pirho;
		}

		inline void setCurvature(const double val) {
			this->Pirho = val;
		}
		inline int getNumberOfIntegrationSteps(void) const {
			return this->integration_steps;
		}

		inline void setFieldInterpolator(std::shared_ptr<thor_scsi::core::Field2DInterpolationKnobbed<C>> a_intp){
			this->intp = a_intp;
		}
		inline auto getFieldInterpolator(void) const {
			/*
			  std::cerr << "Getting field interpolator " <<  std::endl;
			  std::cerr.flush();
			  std::cerr << " address " << this->intp << std::endl;
			*/
			if(!this->intp){
				throw std::logic_error("interpolator pointer null");
			}
			/*
			  std::cerr << "use count " << this->intp.use_count() << std::endl;
			  std::cerr.flush();
			*/
			auto tmp =  std::const_pointer_cast<thor_scsi::core::Field2DInterpolationKnobbed<C>>(this->intp);
			return tmp;
			// return ;
		}
	protected:
		double Pirho = 0;             ///< 1/rho [1/m].
		std::shared_ptr<thor_scsi::core::Field2DInterpolationKnobbed<C>> intp;
		int integration_steps = 1;
	};

    typedef  FieldKickAPIKnobbed<thor_scsi::core::StandardDoubleType> FieldKickAPI;
} // Name space

#endif // _THOR_SCSI_ELEMENTS_FIELD_KICK_API_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
