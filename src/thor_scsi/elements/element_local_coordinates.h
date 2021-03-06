#ifndef _THOR_SCSI_ELEMENTS_LOCAL_COORDINATES_
#define _THOR_SCSI_ELEMENTS_LOCAL_COORDINATES_

#include <thor_scsi/core/elements_basis.h>
#include <thor_scsi/core/transform_phase_space.h>
#include <tps/tps_type.h>

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
		virtual ~LocalCoordinates(){}
		LocalCoordinates(LocalCoordinates&& o) : ElemType(std::move(o)) {
			/*
			std::cerr << __FILE__ << "::" << __FUNCTION__ << " ctor @ " << __LINE__
				  << " name " << this->name << std::endl;
			*/
		}

		inline virtual void global2Local(ss_vect<double> &ps) = 0;
		inline virtual void local2Global(ss_vect<double> &ps) = 0;

		inline virtual void global2Local(ss_vect<tps> &ps) = 0;
		inline virtual void local2Global(ss_vect<tps> &ps) = 0;

		virtual void localPass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps) = 0;
		virtual void localPass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps) = 0;

		inline void pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps) override final
			{ _pass(conf, ps); };
		inline void pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps) override final
			{ _pass(conf, ps); };

	private:
		template<typename T>
		void _pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps){
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
		virtual ~LocalGalilean(){}
		inline LocalGalilean(LocalGalilean&& o) :
			LocalCoordinates(std::move(o)),
			transform(std::move(o.transform))
			{

			}

		inline virtual void global2Local(ss_vect<double> &ps) override final {
			this->_global2Local(ps);
		}

		inline virtual void local2Global(ss_vect<double> &ps) override final {
			this->_local2Global(ps);
		}
		inline virtual void global2Local(ss_vect<tps> &ps) override final {
			this->_global2Local(ps);
		}

		inline virtual void local2Global(ss_vect<tps> &ps) override final {
			this->_local2Global(ps);
		}

		inline auto* getTransform(void){
			return &this->transform;
		}
		thor_scsi::core::PhaseSpaceGalilean2DTransform transform;

	private:
		template<typename T>
		void _global2Local(ss_vect<T> &ps){
			this->transform.forward(ps);
		}
		template<typename T>
		void _local2Global(ss_vect<T> &ps){
			this->transform.backward(ps);
		}
	};

	/*
	 * @brief: Device in local coordinates translated and rotated
	 * @todo: could the template functions be defind by deriving ...
	 * see thor_scsi::core::PhaseSpaceGalileanPRot2DTransform for implementation
	 */
	class LocalGalileanPRot  : public LocalCoordinates {

	public:
		inline LocalGalileanPRot(const Config &config) : LocalCoordinates(config) {}
		virtual ~LocalGalileanPRot(){}
		inline LocalGalileanPRot(LocalGalileanPRot&& o) :
			LocalCoordinates(std::move(o))
			{
				this->transform = o.transform;
				//transform(std::move(o.transform));
			}

		inline virtual void global2Local(ss_vect<double> &ps) override final {
			this->transform.forward(ps);
		}
		inline virtual void local2Global(ss_vect<double> &ps) override final {
			this->transform.backward(ps);
		}
		inline virtual void global2Local(ss_vect<tps> &ps) override final {
			this->_global2Local(ps);
		}

		inline virtual void local2Global(ss_vect<tps> &ps) override final {
			this->_local2Global(ps);
		}

		inline auto* getTransform(void){
			return &this->transform;
		}

		thor_scsi::core::PhaseSpaceGalileanPRot2DTransform transform;

	private:
		template<typename T>
		void _global2Local(ss_vect<T> &ps){
			this->transform.forward(ps);
		}
		template<typename T>
		void _local2Global(ss_vect<T> &ps){
			this->transform.backward(ps);
		}
	};


} //namespace thor_scsi::elements

#endif // _THOR_SCSI_ELEMENTS_LOCAL_COORDINATES_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
