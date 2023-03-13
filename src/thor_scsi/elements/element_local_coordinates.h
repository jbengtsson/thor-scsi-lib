#ifndef _THOR_SCSI_ELEMENTS_LOCAL_COORDINATES_
#define _THOR_SCSI_ELEMENTS_LOCAL_COORDINATES_

#include <thor_scsi/core/elements_basis.h>
#include <thor_scsi/core/config.h>
#include <thor_scsi/core/transform_phase_space.h>
#include <tps/tps_type.h>

namespace thor_scsi::elements {
	using thor_scsi::core::ElemTypeKnobbed;
	using thor_scsi::core::ConfigType;
	/**
	 * @ brief: provide a local_pass method which will be executed in local coordiantes
	 *
	 * This class is never expected to be instantiated by it self so it gets no type name
	 *
	 * Todo:
	 *      transformation should contain a constant part and a random part
	 */
	 template<class C>
	 class LocalCoordinatesKnobbed : public ElemTypeKnobbed /* <C> */ {

	public:
		 inline LocalCoordinatesKnobbed(const Config &config) : ElemTypeKnobbed /* <C> */ (config){}
		virtual ~LocalCoordinatesKnobbed(){}
		 LocalCoordinatesKnobbed(LocalCoordinatesKnobbed&& o) : ElemTypeKnobbed /* <C> */ (std::move(o)) {
			/*
			std::cerr << __FILE__ << "::" << __FUNCTION__ << " ctor @ " << __LINE__
				  << " name " << this->name << std::endl;
			*/
		}

		// inline virtual void global2Local(ss_vect<double>             &ps) = 0;
		// inline virtual void local2Global(ss_vect<double>             &ps) = 0;
		// inline virtual void global2Local(ss_vect<tps>                &ps) = 0;
		// inline virtual void local2Global(ss_vect<tps>                &ps) = 0;

		inline virtual void global2Local(gtpsa::ss_vect<double>      &ps) = 0;
		inline virtual void global2Local(gtpsa::ss_vect<gtpsa::tpsa> &ps) = 0;
		 // inline virtual void global2Local(gtpsa::ss_vect<tps>         &ps) = 0;

		inline virtual void local2Global(gtpsa::ss_vect<double>      &ps) = 0;
		inline virtual void local2Global(gtpsa::ss_vect<gtpsa::tpsa> &ps) = 0;
		 // inline virtual void local2Global(gtpsa::ss_vect<tps>         &ps) = 0;

		// virtual void localPropagate(ConfigType &conf, ss_vect<double>             &ps)  = 0;
		// virtual void localPropagate(ConfigType &conf, ss_vect<tps>                &ps)  = 0;

		virtual void localPropagate(ConfigType &conf, gtpsa::ss_vect<double>      &ps)  = 0;
		virtual void localPropagate(ConfigType &conf, gtpsa::ss_vect<gtpsa::tpsa> &ps)  = 0;
		 // virtual void localPropagate(ConfigType &conf, gtpsa::ss_vect<tps>         &ps)  = 0;

		// inline void propagate(ConfigType &conf, ss_vect<double>             &ps) override final { _propagate(conf, ps); };
		// inline void propagate(ConfigType &conf, ss_vect<tps>                &ps) override final { _propagate(conf, ps); };
		virtual inline void propagate(ConfigType &conf, gtpsa::ss_vect<double>      &ps) override final { _propagate(conf, ps); };
		virtual inline void propagate(ConfigType &conf, gtpsa::ss_vect<gtpsa::tpsa> &ps) override final { _propagate(conf, ps); };
		 // virtual inline void propagate(ConfigType &conf, gtpsa::ss_vect<tps>         &ps) override final { _propagate(conf, ps); };

	private:
		// template<typename T>
		// void _propagate(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps){
		//	this->global2Local(ps);
		//	this->localPropagate(conf, ps);
		//	this->local2Global(ps);
		// }

		template<typename T>
		void _propagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<T> &ps){
			this->global2Local(ps);
			this->localPropagate(conf, ps);
			this->local2Global(ps);
		}
	};

	/*
	 * @brief: Device in local coordinates translated and rotated
	 *
	 * see thor_scsi::core::PhaseSpaceGalilean2DTransform for implementation
	 */
    template<class C>
	class LocalGalileanKnobbed : public LocalCoordinatesKnobbed<C> {

	public:
		inline LocalGalileanKnobbed(const Config &config)
			: LocalCoordinatesKnobbed<C>(config)
			, transform()
			{}
		virtual ~LocalGalileanKnobbed(){}
		inline LocalGalileanKnobbed(LocalGalileanKnobbed&& o) :
			LocalCoordinatesKnobbed<C>(std::move(o)),
			transform(std::move(o.transform))
			{

			}

		// inline virtual void global2Local(ss_vect<double>             &ps) override { this->_global2Local(ps); }
		// inline virtual void global2Local(ss_vect<tps>                &ps) override { this->_global2Local(ps); }
		// inline virtual void local2Global(ss_vect<tps>                &ps) override { this->_local2Global(ps); }
		// inline virtual void local2Global(ss_vect<double>             &ps) override { this->_local2Global(ps); }

		inline virtual void global2Local(gtpsa::ss_vect<double>      &ps) override { this->_global2Local(ps); }
	        inline virtual void global2Local(gtpsa::ss_vect<gtpsa::tpsa> &ps) override { this->_global2Local(ps); }
	        // inline virtual void global2Local(gtpsa::ss_vect<tps>         &ps) override { this->_global2Local(ps); }

		inline virtual void local2Global(gtpsa::ss_vect<double>      &ps) override { this->_local2Global(ps); }
		inline virtual void local2Global(gtpsa::ss_vect<gtpsa::tpsa> &ps) override { this->_local2Global(ps); }
	        // inline virtual void local2Global(gtpsa::ss_vect<tps>         &ps) override { this->_local2Global(ps); }


		inline auto* getTransform(void){
			return &this->transform;
		}
		thor_scsi::core::PhaseSpaceGalilean2DTransformKnobbed<C> transform;

	private:
		// template<typename T> void _global2Local(ss_vect<T>        &ps){ this->transform.forward(ps);	}
		// template<typename T> void _local2Global(ss_vect<T>        &ps){ this->transform.backward(ps);	}
		template<typename T> void _global2Local(gtpsa::ss_vect<T> &ps){	this->transform.forward(ps);	}
		template<typename T> void _local2Global(gtpsa::ss_vect<T> &ps){	this->transform.backward(ps);	}
	};

	/*
	 * @brief: Device in local coordinates translated and rotated
	 * @todo: could the template functions be defind by deriving ...
	 * see thor_scsi::core::PhaseSpaceGalileanPRot2DTransform for implementation
	 */
    template<class C>
	class LocalGalileanPRotKnobbed  : public LocalCoordinatesKnobbed<C> {

	public:
		inline LocalGalileanPRotKnobbed(const Config &config)
			: LocalCoordinatesKnobbed<C>(config)
			, transform()
			{}

		virtual ~LocalGalileanPRotKnobbed(){}
		inline LocalGalileanPRotKnobbed(LocalGalileanPRotKnobbed&& o)
			: LocalCoordinatesKnobbed<C>(std::move(o) )
			, transform(std::move(o.transform))
			{
				// this->transform = o.transform;
				//transform(std::move(o.transform));
			}

		// inline virtual void global2Local(ss_vect<double> &ps) override final { this->_global2Local(ps);  }
		// inline virtual void global2Local(ss_vect<tps>    &ps) override final { this->_global2Local(ps);  }
		// inline virtual void local2Global(ss_vect<double> &ps) override final { this->_local2Global(ps);  }
		// inline virtual void local2Global(ss_vect<tps>    &ps) override final { this->_local2Global(ps);  }

		inline virtual void global2Local(gtpsa::ss_vect<double>      &ps) override final { this->_global2Local(ps);  }
	    // inline virtual void global2Local(gtpsa::ss_vect<tps>         &ps) override final { this->_global2Local(ps);  }
		inline virtual void global2Local(gtpsa::ss_vect<gtpsa::tpsa> &ps) override final { this->_global2Local(ps);  }
		inline virtual void local2Global(gtpsa::ss_vect<double>      &ps) override final { this->_local2Global(ps);  }
	    // inline virtual void local2Global(gtpsa::ss_vect<tps>         &ps) override final { this->_local2Global(ps);  }
		inline virtual void local2Global(gtpsa::ss_vect<gtpsa::tpsa> &ps) override final { this->_local2Global(ps);  }


		inline auto* getTransform(void){return &this->transform;		}

		thor_scsi::core::PhaseSpaceGalileanPRot2DTransformKnobbed<C> transform;

	private:
		// template<typename T> void _global2Local(ss_vect<T> &ps){ this->transform.forward(ps);  }
		// template<typename T> void _local2Global(ss_vect<T> &ps){ this->transform.backward(ps); }
		template<typename T> void _global2Local(gtpsa::ss_vect<T> &ps){ this->transform.forward(ps);  }
		template<typename T> void _local2Global(gtpsa::ss_vect<T> &ps){ this->transform.backward(ps); }
	};

    typedef LocalGalileanKnobbed<thor_scsi::core::StandardDoubleType> LocalGalilean;
    typedef LocalGalileanPRotKnobbed<thor_scsi::core::StandardDoubleType> LocalGalileanPRot;

} //namespace thor_scsi::elements

#endif // _THOR_SCSI_ELEMENTS_LOCAL_COORDINATES_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
