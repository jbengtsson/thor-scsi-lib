#ifndef _THOR_SCSI_CORE_ELEMENTS_DRIFT_H_
#define _THOR_SCSI_CORE_ELEMENTS_DRIFT_H_ 1

#include <thor_scsi/core/elements_basis.h>
#include <thor_scsi/core/multipole_types.h>


namespace thor_scsi::elements {
		/**

		   Empty space between two "typical accelerator components"
		 */
		using thor_scsi::core::ElemTypeKnobbed;
		using thor_scsi::core::ConfigType;
		template<class C>
		class DriftTypeWithKnob : public ElemTypeKnobbed<C> {
		public:
			using base = ElemTypeKnobbed<C>;
			inline DriftTypeWithKnob(const Config &config) :
				base(config)
				{
				// transformation done by transfrom
				// ... done by Elemtype initialisation
				// ... pleonamsmus
				}


			const char* type_name() const override final { return "Drift"; };
			/**
			 *
			 * Todo:
			 *     replace function with mv operator
			 */
			//virtual void assign(const ElementVoidBase *other) override{
			//	const DriftType *O = static_cast<const DriftType*>(other);
			//	// transform = O->transform
			//	ElementVoidBase::assign(other);
			//}
			// double GetPB(const int n) { return 0e0; };

		        // inline void propagate(ConfigType &conf, ss_vect<double>             &ps) override final { _propagate(conf, ps); };
			// inline void propagate(ConfigType &conf, ss_vect<tps>                &ps) override final { _propagate(conf, ps); };
			inline virtual void propagate(ConfigType &conf, gtpsa::ss_vect<double>      &ps) override final { _propagate(conf, ps); };
			inline virtual void propagate(ConfigType &conf, gtpsa::ss_vect<gtpsa::tpsa> &ps) override final { _propagate(conf, ps); };
			// inline virtual void propagate(ConfigType &conf, gtpsa::ss_vect<tps>         &ps) override final { _propagate(conf, ps); };

		private:
			// template<typename T> void _propagate(const ConfigType &conf, ss_vect<T>        &ps);
			template<typename T> void _propagate(const ConfigType &conf, gtpsa::ss_vect<T> &ps);
		};

		typedef DriftTypeWithKnob<thor_scsi::core::StandardDoubleType> DriftType;
	        typedef DriftTypeWithKnob<thor_scsi::core::TpsaVariantType> DriftTypeTpsa;

} // namespace thor_scsi::elements

#endif // _THOR_SCSI_CORE_ELEMENTS_DRIFT_H_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
