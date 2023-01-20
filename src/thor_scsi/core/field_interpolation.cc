#include <thor_scsi/core/field_interpolation.h>
#include <boost/type_index.hpp>
#include <sstream>

namespace tsc = thor_scsi::core;
namespace bti = boost::typeindex;

template<>
std::string tsc::Field2DInterpolation::prettyClassname(void) const
{
	return bti::type_id_with_cvr<decltype(this)>().pretty_name();
}

template<>
std::string tsc::Field2DInterpolation::pstr(void) const
{
	std::stringstream strm;

	strm << "<" << this->prettyClassname()
	     << " @ " << (const void *)(this) << ">(";
	this->show(strm, 0);
	strm << ")";
	return strm.str();
}

template<>
std::string tsc::Field2DInterpolation::repr(void) const
{
	std::stringstream strm;

	strm << "<" << this->prettyClassname()
	     << " @ " << (const void *)(this) << ">(";
	this->show(strm, 10);
	strm << ")";
	return strm.str();
}
