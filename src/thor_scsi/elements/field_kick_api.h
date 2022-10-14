#ifndef _THOR_SCSI_ELEMENTS_FIELD_KICK_API_H_
#define _THOR_SCSI_ELEMENTS_FIELD_KICK_API_H_
#include <thor_scsi/elements/element_local_coordinates.h>
#include <thor_scsi/core/field_interpolation.h>

namespace thor_scsi::elements {
	// declaration of API ... to be used by radiation observer
	class FieldKickAPI : public LocalGalileanPRot {

	public:
		inline FieldKickAPI(const Config &conf)
			: LocalGalileanPRot(conf)
			, intp(nullptr)
			{}
		inline FieldKickAPI(FieldKickAPI&& O)
			: LocalGalileanPRot(std::move(O)),
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

		inline void setFieldInterpolator(std::shared_ptr<thor_scsi::core::Field2DInterpolation> a_intp){
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
			auto tmp =  std::const_pointer_cast<thor_scsi::core::Field2DInterpolation>(this->intp);
			return tmp;
			// return ;
		}
	protected:
		double Pirho = 0;             ///< 1/rho [1/m].
		std::shared_ptr<thor_scsi::core::Field2DInterpolation> intp;
		int integration_steps = 1;
	};

} // Name space

#endif // _THOR_SCSI_ELEMENTS_FIELD_KICK_API_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
