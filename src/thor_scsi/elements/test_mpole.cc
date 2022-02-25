#define BOOST_TEST_MODULE mpole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <thor_scsi/elements/mpole.h>

namespace tse = thor_scsi::elements;

BOOST_AUTO_TEST_CASE(test01_drift_print)
{
	Config C;
	C.set<std::string>("name", "d");
	tse::MpoleType mpole(C);
}
