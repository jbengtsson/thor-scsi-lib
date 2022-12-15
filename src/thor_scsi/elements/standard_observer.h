#ifndef _THOR_SCSI_STD_MACHINE_OBSERVER_
#define _THOR_SCSI_STD_MACHINE_OBSERVER_
#include <thor_scsi/core/cell_void.h>

namespace thor_scsi::elements{
	/**
	 * @brief observer making copies of the data received at end state
	 */

	class StandardObserver : public thor_scsi::core::Observer{
	public:

		void view(std::shared_ptr<const thor_scsi::core::CellVoid> elem, const gtpsa::ss_vect<double>      &ps, const enum thor_scsi::core::ObservedState, const int cnt) override final;
		void view(std::shared_ptr<const thor_scsi::core::CellVoid> elem, const gtpsa::ss_vect<tps>         &ps, const enum thor_scsi::core::ObservedState, const int cnt) override final;
		void view(std::shared_ptr<const thor_scsi::core::CellVoid> elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum thor_scsi::core::ObservedState, const int cnt) override final;

		inline std::string getObservedName(void) const {
			return this->m_observed_name;
		}

		inline int getObservedIndex(void) const {
			return this->m_observed_index;
		}

		inline bool hasPhaseSpace(void) const {
			return this->m_has_ps;
		}

		inline  gtpsa::ss_vect<double>& getPhaseSpace(void){
			return this->m_ps;
		}

		/*
		 * @brief
		 *
		 * @todo  decide if it should return tps
		 */
		inline bool hasTruncatedPowerSeries(void) const {
			return this->m_has_tps;
		}
		inline bool hasTruncatedPowerSeriesA(void) const {
			return this->m_has_tpsa;
		}

		inline  gtpsa::ss_vect<tps>& getTruncatedPowerSeries(void){
			return this->m_tps;
		}
		inline auto getTruncatedPowerSeriesA(void){
			return this->m_ptpsa;
		}

		inline void reset(void){
			this->m_has_ps = this->m_has_tps = this->m_has_tpsa = false;
		}

		void show(std::ostream& strm, const int level) const override final;

		// python support
		std::string repr(void) const;
		// python support
		std::string pstr(void) const;

	private:
		std::string _repr(const int level) const;

		inline void store(const gtpsa::ss_vect<tps>         &a_tps) {
			this->m_tps = a_tps.clone();
			this->m_has_tps = true;
		}
		inline void store(const gtpsa::ss_vect<gtpsa::tpsa> &a_tps)  {
			this->m_ptpsa = std::make_shared<gtpsa::ss_vect<gtpsa::tpsa>>(a_tps[0]);
			this->m_ptpsa->_copyInPlace(a_tps);
			this->m_has_tpsa = true;
		}

		inline void store(const gtpsa::ss_vect<double>       &a_ps) {
			this->m_ps = a_ps.clone();
			this->m_has_ps = true;
		}

		template<typename T>
		void _view(std::shared_ptr<const thor_scsi::core::CellVoid> elem, const gtpsa::ss_vect<T> &ps, const enum thor_scsi::core::ObservedState, const int cnt);

		gtpsa::ss_vect<double> m_ps = {0, 0, 0, 0, 0, 0};
		gtpsa::ss_vect<tps> m_tps = {tps(), tps(), tps(), tps(), tps(), tps()};
		std::shared_ptr<gtpsa::ss_vect<gtpsa::tpsa>> m_ptpsa = nullptr;
		bool m_has_ps = false, m_has_tps = false, m_has_tpsa = false;
		std::string m_observed_name = "";
		int m_observed_index=-1;
	};



}
#endif /* _THOR_SCSI_STD_MACHINE_OBSERVER_ */
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
