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
std::string tsc::Field2DInterpolationDependence::prettyClassname(void) const
{
    return bti::type_id_with_cvr<decltype(this)>().pretty_name();
}

template<typename T>
static std::string show_helper(T& o, int level)
{
	std::stringstream strm;

	strm << "<" << o.prettyClassname()
	     << " @ " << (const void *)(&o) << ">(";
	o.show(strm, level);
	strm << ")";
	return strm.str();
}

template<>
std::string tsc::Field2DInterpolation::pstr(void) const { return show_helper(*this,  0); }
template<>
std::string tsc::Field2DInterpolation::repr(void) const { return show_helper(*this, 10); }
template<>
std::string tsc::Field2DInterpolationDependence::pstr(void) const { return show_helper(*this,  0); }
template<>
std::string tsc::Field2DInterpolationDependence::repr(void) const { return show_helper(*this, 10); }
