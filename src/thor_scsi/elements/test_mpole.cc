#define BOOST_TEST_MODULE mpole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <thor_scsi/elements/mpole.h>
#include <ostream>

namespace tsc = thor_scsi::core;
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


BOOST_AUTO_TEST_CASE(test03_mpole_wrong_method)
{
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("Method", 5.0);
	//tse::MpoleType mpole(C);

	BOOST_CHECK_THROW(((tse::MpoleType)(C)), thor_scsi::NotImplemented);

}

BOOST_AUTO_TEST_CASE(test10_mpole_kick_zero)
{
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");

	tse::MpoleType mpole(C);

	const ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		ss_vect<double> ps = ps_orig;
		mpole.pass(calc_config, ps);

		BOOST_CHECK_SMALL(ps[x_],  1e-14);
		BOOST_CHECK_SMALL(ps[y_],  1e-14);
		BOOST_CHECK_SMALL(ps[px_], 1e-14);
		BOOST_CHECK_SMALL(ps[py_], 1e-14);
		BOOST_CHECK_SMALL(ps[ct_], 1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

	}
}

BOOST_AUTO_TEST_CASE(test11_mpole_kick_longitudinal_one)
{
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("L", 1.0);
	tse::MpoleType mpole(C);


	BOOST_CHECK_CLOSE(mpole.getLength(), 1.0, 1e-14);

	std::cout << "mpole "; mpole.show(std::cout, 4); std::cout<<std::endl;
	// is this correct for a thin kick
	const ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		ss_vect<double> ps = ps_orig;
		mpole.pass(calc_config, ps);

		BOOST_CHECK_SMALL(ps[x_],  1e-14);
		BOOST_CHECK_SMALL(ps[y_],  1e-14);
		BOOST_CHECK_SMALL(ps[px_], 1e-14);
		BOOST_CHECK_SMALL(ps[py_], 1e-14);
		BOOST_CHECK_SMALL(ps[ct_], 1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

		std::cout << "thin kick ps " << ps << std::endl;
	}

	{
		ss_vect<double> ps = ps_orig;
		mpole.asThick(true);
		mpole.pass(calc_config, ps);

		BOOST_CHECK_SMALL(ps[x_],  1e-14);
		BOOST_CHECK_SMALL(ps[y_],  1e-14);
		BOOST_CHECK_SMALL(ps[px_], 1e-14);
		BOOST_CHECK_SMALL(ps[py_], 1e-14);
		BOOST_CHECK_SMALL(ps[ct_], 1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

		std::cout << "thick kick ps " << ps << std::endl;
	}
}
