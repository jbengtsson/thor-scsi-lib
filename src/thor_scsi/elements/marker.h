#ifndef _THOR_SCSI_ELEMENTS_MARKER_H_
#define _THOR_SCSI_ELEMENTS_MARKER_H_

#include <thor_scsi/elements/element_local_coordinates.h>
#include <thor_scsi/elements/element_helpers.h>
#include <thor_scsi/elements/radiation_delegate_api.h>

namespace thor_scsi::elements {
	/**
	 * @brief: Marker: just collect current position
	 *
	 * Doing nothing yet as observers are not yet implemented
	 * @todo: implement radiation observer in a base class?
	 */
	class MarkerType : public LocalGalilean {
	public:
		inline MarkerType(const Config &config)
		    : LocalGalilean(config)
		    , rad_del(nullptr)
		{}

		const char* type_name(void) const override { return "Marker"; };

		inline void setRadiationDelegate(std::shared_ptr<thor_scsi::elements::RadiationDelegateInterface> rad_del){
			this->rad_del = rad_del;
		}
		inline auto getRadiationDelegate(void) {
			return this->rad_del;
		}

	    virtual void localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<double>      &ps) override final { _localPropagate(conf, ps);}
	    // virtual void localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<tps>         &ps) override final { _localPropagate(conf, ps);}
	    virtual void localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<gtpsa::tpsa> &ps) override final { _localPropagate(conf, ps);}

	private:
		template<typename T>
		inline void _localPropagate(const thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<T> &ps){
			if (conf.emittance && !conf.Cavity_on && this->rad_del){
				rad_del->view(*this, ps, thor_scsi::core::ObservedState::end, 0);
				// // Needs A^-1.
				// curly_dH_x = is_tps<tps>::get_curly_H(ps);
			}
		}

		std::shared_ptr<thor_scsi::elements::RadiationDelegateInterface> rad_del;

	};

} // Name space

#endif // _THOR_SCSI_ELEMENTS_MARKER_H_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
