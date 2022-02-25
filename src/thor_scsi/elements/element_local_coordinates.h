#ifndef _THOR_SCSI_ELEMENTS_LOCAL_COORDINATES_
#define _THOR_SCSI_ELEMENTS_LOCAL_COORDINATES_

#include <thor_scsi/core/elements_basis.h>
#include <thor_scsi/core/transform_phase_space.h>

namespace thor_scsi::elements {
	using thor_scsi::core::ElemType;
	/**
	 * @ brief: provide a local_pass method which will be executed in local coordiantes
	 *
	 * This class is never expected to be instanciated by it self so it gets no type name
	 *
	 * Todo:
	 *      transformation should contain a constant part and a random part
	 */
	class LocalCoordinates : public ElemType {

	public:
		inline LocalCoordinates(const Config &config) : ElemType(config){}

		inline virtual void global2Local(ss_vect<double> &ps) = 0;
		inline virtual void local2Global(ss_vect<double> &ps) = 0;

		//inline virtual void global2Local(ss_vect<tps> &ps) = 0;
		//inline virtual void local2Global(ss_vect<tps> &ps) = 0;

		virtual  void localPass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps) = 0;

		inline void pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps) override final
			{ _pass(conf, ps); };
		// inline void pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps) override final
		// { _pass(conf, ps); };

	private:
		template<typename T>
		void _pass(const thor_scsi::core::ConfigType &conf, ss_vect<T> &ps){
			this->global2Local(ps);
			this->localPass(conf, ps);
			this->local2Global(ps);
		}
	};

	/*
	 * @brief: Device in local coordinates translated and rotated
	 *
	 * see thor_scsi::core::PhaseSpaceGalilean2DTransform for implementation
	 */
	class LocalGalilean : public LocalCoordinates {

	public:
		inline LocalGalilean(const Config &config) : LocalCoordinates(config) {}
		inline virtual void global2Local(ss_vect<double> &ps) override final {
			this->transform.forward(ps);
		}

		inline virtual void local2Global(ss_vect<double> &ps) override final {
			this->transform.backward(ps);
		}

		thor_scsi::core::PhaseSpaceGalilean2DTransform transform;
	};

	/*
	 * @brief: Device in local coordinates translated and rotated
	 *
	 * see thor_scsi::core::PhaseSpaceGalileanPRot2DTransform for implementation
	 */
	class LocalGalileanPRot  : public LocalCoordinates {

	public:
		inline LocalGalileanPRot(const Config &config) : LocalCoordinates(config) {}
		inline virtual void global2Local(ss_vect<double> &ps) override final {
			this->transform.forward(ps);
		}
		inline virtual void local2Global(ss_vect<double> &ps) override final {
			this->transform.backward(ps);
		}

		thor_scsi::core::PhaseSpaceGalileanPRot2DTransform transform;
	};


} //namespace thor_scsi::elements

#endif // _THOR_SCSI_ELEMENTS_LOCAL_COORDINATES_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
