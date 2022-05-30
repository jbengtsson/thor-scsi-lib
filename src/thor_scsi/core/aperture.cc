#include <thor_scsi/core/aperture.h>
#include <sstream>

namespace tsc = thor_scsi::core;

void tsc::TwoDimensionalAperture::show(std::ostream& strm, int level) const
{
	strm << this->type_name();
}

std::string tsc::TwoDimensionalAperture::repr(void) const
{
	std::ostringstream strm;

	strm << "< some aperture"
	     << " @ " << (const void *)(this) << ">(";
	this->show(strm, 10);
	strm << ")";
	return strm.str();
}
