#define BOOST_TEST_MODULE mpole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <thor_scsi/elements/mpole.h>

namespace tse = thor_scsi::elements;

BOOST_AUTO_TEST_CASE(test01_kick_print)
{
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("Method", 4.0);
	tse::FieldKick kick(C);

	std::cout<< "Field kick " << kick <<std::endl;
	std::cout<< "     ";
	kick.show(std::cout, 4);
	std::cout<<std::endl;
}

BOOST_AUTO_TEST_CASE(test02_mpole_print)
{
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("Method", 4.0);
	tse::MpoleType mpole(C);

	std::cout<< "mpole " << mpole <<std::endl;
	std::cout<< "     ";
	mpole.show(std::cout, 4);
	std::cout<<std::endl;
}
