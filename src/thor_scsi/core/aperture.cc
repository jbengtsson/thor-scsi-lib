#include <thor_scsi/core/aperture.h>

namespace tsc = thor_scsi::core;

void tsc::TwoDimensionalAperture::show(std::ostream& strm, int level) const
{
	strm << this->type_name();
}
