#include <thor_scsi/elements/standard_observer.h>

namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;

template<typename T>
inline void tse::StandardObserver::_view(std::shared_ptr<const tsc::CellVoid> elem, const gtpsa::ss_vect<T> &ps, const enum tsc::ObservedState state, const int cnt)
{
	switch(state){
	case tsc::ObservedState::start:
		this->reset();
		this->m_observed_name = elem->name;
		this->m_observed_index = elem->index;
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

void tse::StandardObserver::view(std::shared_ptr<const tsc::CellVoid> elem, const gtpsa::ss_vect<double> &ps, const enum tsc::ObservedState state, const int cnt)
{
	_view(elem, ps, state, cnt);
}

void tse::StandardObserver::view(std::shared_ptr<const tsc::CellVoid> elem, const gtpsa::ss_vect<tps> &ps, const enum tsc::ObservedState state, const int cnt)
{
	_view(elem, ps, state, cnt);
}
void tse::StandardObserver::view(std::shared_ptr<const tsc::CellVoid> elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum tsc::ObservedState state, const int cnt)
{
	_view(elem, ps, state, cnt);
}


void tse::StandardObserver::show(std::ostream& strm, const int level) const
{
	strm << "StandardObserver("
	     << "observed_name=\"" << this->m_observed_name
	     << "\", observed_index=" << this->m_observed_index
	     << ", has_ps=" << this->m_has_ps <<", ps=";
	this->m_ps.show(strm, 4, false);
		strm << ", has_tps=" << this->m_has_tps
		     << ", has_tpsa=" << this->m_has_tpsa;
	this->m_tps.show(strm, 4, false);
	strm << ")";
}

std::string tse::StandardObserver::_repr(const int level) const
{
	std::ostringstream strm;
	this->show(strm, level);
	return strm.str();
}

std::string tse::StandardObserver::repr(void) const
{
	return this->_repr(10);
}

std::string tse::StandardObserver::pstr(void) const
{
	return this->_repr(0);
}
