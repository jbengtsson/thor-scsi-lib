#ifndef _THOR_SCSI_ELEMENTS_MARKER_H_
#define _THOR_SCSI_ELEMENTS_MARKER_H_

#include <thor_scsi/elements/element_local_coordinates.h>

namespace thor_scsi::elements {
	/**
	 * @brief: Marker: just collect current position
	 *
	 * Doing nothing yet as observers are not yet implemented
	 */
	class MarkerType : public LocalGalilean {
	public:
		inline MarkerType(const Config &config) : LocalGalilean(config){
		}

		const char* type_name(void) const override { return "Marker"; };

		virtual void localPass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps) override final
		{ _localPass(conf, ps);}
		virtual void localPass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps) override final
		{ _localPass(conf, ps);}

	private:
		template<typename T>
		inline void _localPass(const thor_scsi::core::ConfigType &conf, ss_vect<T> &ps){
#ifdef THOR_SCSI_USE_RADIATION
			if (conf.emittance && !conf.Cavity_on){
				// Needs A^-1.
				curly_dH_x = is_tps<tps>::get_curly_H(ps);
			}
#endif /* THOR_SCSI_USE_RADIATION */
		}
	};

} // Name space

#endif // _THOR_SCSI_ELEMENTS_MARKER_H_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
