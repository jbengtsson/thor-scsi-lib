#ifndef _THOR_SCSI_ELEMENTS_CAVITY_H_
#define _THOR_SCSI_ELEMENTS_CAVITY_H_

#include <thor_scsi/elements/element_local_coordinates.h>
#include <tps/tps_type.h>

namespace thor_scsi::elements {
	/**
	 * @brief: Marker: just collect current position
	 *
	 * Doing nothing yet as observers are not yet implemented
	 */
	class CavityType : public LocalGalilean {
	public:
		CavityType(const Config &config);

		const char* type_name(void) const override final { return "Cavity"; };

		// virtual void localPropagate(thor_scsi::core::ConfigType &conf, ss_vect<double>             &ps) override final { _localPropagate(conf, ps); }
		// virtual void localPropagate(thor_scsi::core::ConfigType &conf, ss_vect<tps>                &ps) override final { _localPropagate(conf, ps); }
		virtual void localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<double>      &ps) override final { _localPropagate(conf, ps); }
	    // virtual void localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<tps>         &ps) override final { _localPropagate(conf, ps); }
		virtual void localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<gtpsa::tpsa> &ps) override final { _localPropagate(conf, ps); }

		inline void setVoltage(const double val){
			this->Pvolt = val;
		}

		inline double getVoltage(void) const {
			return this->Pvolt;
		}

		inline void setFrequency(const double val){
			this->Pfreq = val;
		}

		inline double getFrequency(void) const {
			return this->Pfreq;
		}

		inline void setPhase(const double val){
			this->phi = val;
		}
		inline double getPhase(void) const {
			return this->phi;
		}

		inline void setHarmonicNumber(const int n){
			this->Ph = n;
		}

		inline int getHarmonicNumber(void) const {
			return this->Ph;
		}

	private:
		/**
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * .. Todo::
		 *
		 *    Energy deviation stored in configuration: review where to store it
		 *
		 * \endverbatim
		 */
		template<typename T> void _localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<T> &ps);
		// template<typename T>
		// void _localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<tpsa> &ps);

		double
		Pvolt = 0.0,                     ///< Vrf [V].
			Pfreq = 0.0,             ///< Vrf [Hz].
			phi = 0.0;               ///< RF phase.

		bool
		entry_focus = false,             ///< Edge focusing at entry.
			exit_focus = false;      ///< Edge focusing at exit.
		int
		// PN = 1,                          ///< Number of integration step (currently unused).
			Ph = 1;                  ///< Harmonic number.

	};
} // Name space

#endif // _THOR_SCSI_ELEMENTS_CAVITY_H_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
