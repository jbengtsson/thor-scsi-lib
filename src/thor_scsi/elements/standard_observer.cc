#include <thor_scsi/elements/standard_observer.h>

namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;

template<typename T>
inline void tse::StandardObserver::_view(std::shared_ptr<const tsc::CellVoid> elem, const ss_vect<T> &ps, const enum tsc::ObservedState state, const int cnt)
{
	switch(state){
	case tsc::ObservedState::start:
		this->reset();
		this->m_delegator_name = elem->name;
		this->m_delegator_index = elem->index;
		return;
		break;
	case tsc::ObservedState::end:
		this->store(ps);
		return;
		break;
	default:
		return;
	}
}

void tse::StandardObserver::view(std::shared_ptr<const tsc::CellVoid> elem, const ss_vect<double> &ps, const enum tsc::ObservedState state, const int cnt)
{
	_view(elem, ps, state, cnt);
}

void tse::StandardObserver::view(std::shared_ptr<const tsc::CellVoid> elem, const ss_vect<tps> &ps, const enum tsc::ObservedState state, const int cnt)
{
	_view(elem, ps, state, cnt);
}


void tse::StandardObserver::show(std::ostream& strm, int level) const
{
	strm << "StandardObserver("
	     << "delegator_name=\"" << this->m_delegator_name
	     << "\", delegator_index=" << this->m_delegator_index
	     << ", has_ps=" << this->m_has_ps <<", ps=";
	this->m_ps.show(strm, 4, false);
		strm << ", has_tps=" << this->m_has_tps <<", tps=";
	this->m_tps.show(strm, 4, false);
	strm << ")";
}
