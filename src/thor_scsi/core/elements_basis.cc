#include <thor_scsi/core/elements_basis.h>

namespace tsc = thor_scsi::core;

void tsc::ElemType::show(std::ostream& strm, int level) const
{
	CellVoid::show(strm, level);
	if(level >= 1){
		strm << " L="<<this->PL<<"";
	}
	if(!this->m_aperture){
		strm << " aperture=None";
	} else {
		strm << " aperture=";
		this->m_aperture->show(strm, level);
	}
	if(!(this->observer())){
		strm << " observer=None";
	} else {
		strm << " observer=";
		this->observer()->show(strm, level);
	}

}
