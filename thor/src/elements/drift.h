#ifndef _THOR_SCSI_CORE_ELEMENTS_DRIFT_H_
#define _THOR_SCSI_CORE_ELEMENTS_DRIFT_H_ 1

#include <thor_scsi/core/elements_basis.h>


namespace thor_scsi {
	namespace elements {
		/**

		   Empty space between two "typical accelerator components"
		 */
		class DriftType : public ElemType {
		public:
			inline DriftType(const Config &config) : ElemType(config){
				PL = 0;
				// ... done by Elemtype initialisation
				// ... pleonamsmus
				this->transform.setDx(0.0);
				this->transform.setDy(0.0);
				this->transform.setRoll(0.0);
			}


			const char* type_name() const override final { return "drift"; };
			virtual std::string repr(void) override final{ return "Dift repr: implement me!";} ;

			/**
			 *
			 * Todo:
			 *     replace function with mv operator
			 */
			virtual void assign(const ElementVoidBase *other) override{
				const DriftType *O = static_cast<const DriftType*>(other);
				// transform = O->transform
				ElementVoidBase::assign(other);
			}

#if 0
			void print(const std::string &);
			std::string repr(void);
			void SetdS(void) {};
			void SetdT(void) {};
			void SetPB(const int n) {};
			double GetdT(void) { return 0e0; };
#endif
			// double GetPB(const int n) { return 0e0; };

			inline void pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps) override final
				{ //_advance(conf, ps);
				};
			inline void pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps) override final
				{ //_advance(conf, ps);
				};

		private:
			template<typename T>
			void _pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
		};
	}
}

#endif // _THOR_SCSI_CORE_ELEMENTS_DRIFT_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
