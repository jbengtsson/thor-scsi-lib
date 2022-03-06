#define BOOST_TEST_MODULE mpole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
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

	{
		boost::test_tools::output_test_stream output;
		output << kick;
		BOOST_CHECK( !output.is_empty( false ) );
	}
	{
		boost::test_tools::output_test_stream output;
		kick.show(output, 4);
		BOOST_CHECK( !output.is_empty( false ) );
	}
}

BOOST_AUTO_TEST_CASE(test01_kick_integral_zero_length)
{
	Config C;
	C.set<std::string>("name", "test");

	{
		C.set<double>("L", 0.0);
		tse::FieldKick kick(C);
		BOOST_CHECK_EQUAL(kick.isIntegral(), true);
	}
	{
		C.set<double>("L", 1.0);
		tse::FieldKick kick(C);
		BOOST_CHECK_EQUAL(int(kick.isIntegral()), 0);
		BOOST_CHECK_EQUAL(kick.isIntegral(), false);
	}
	{
		C.set<double>("L", 0.0);
		tse::FieldKick kick(C);
		BOOST_CHECK_EQUAL(kick.isIntegral(), true);
	}
}

BOOST_AUTO_TEST_CASE(test02_mpole_print)
{
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("Method", 4.0);
	tse::MpoleType mpole(C);

	{
		boost::test_tools::output_test_stream output;
		output << mpole;
		BOOST_CHECK( !output.is_empty( false ) );
	}
	{
		boost::test_tools::output_test_stream output;
		mpole.show(output, 4);
		BOOST_CHECK( !output.is_empty( false ) );
	}
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

	// initialised to 0 length by default
	BOOST_CHECK_EQUAL(mpole.isIntegral(), true);

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

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	// is this correct for a thin kick
	const ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		ss_vect<double> ps = ps_orig;
		mpole.pass(calc_config, ps);

		BOOST_CHECK_SMALL(ps[x_],     1e-14);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[px_],    1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

		// std::cout << "thin kick ps " << ps << std::endl;
	}

	{
		ss_vect<double> ps = ps_orig;
		mpole.asThick(true);
		mpole.pass(calc_config, ps);

		BOOST_CHECK_SMALL(ps[x_],     1e-14);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[px_],    1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

		//std::cout << "thick kick ps " << ps << std::endl;
	}

}

BOOST_AUTO_TEST_CASE(test11_mpole_kick_dipole_component_thin_kick)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("L", 0.0);
	tse::MpoleType mpole(C);

	/* */
	(dynamic_cast<tsc::PlanarMultipoles*>(mpole.intp))->setMultipole(1, tsc::cdbl(1e0,0e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	const ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		ss_vect<double> ps = ps_orig;
		mpole.pass(calc_config, ps);

		/* Length 0 -> harmonics turn integral ? */
		BOOST_CHECK_CLOSE(ps[px_],    1, 1e-14);

		BOOST_CHECK_SMALL(ps[x_],     1e-14);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

	}
}

BOOST_AUTO_TEST_CASE(test11_mpole_kick_dipole_component_thin_kick_l1)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("L", 1.0);
	tse::MpoleType mpole(C);

	/* */
	(dynamic_cast<tsc::PlanarMultipoles*>(mpole.intp))->setMultipole(1, tsc::cdbl(1e-3, 0e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	const ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		ss_vect<double> ps = ps_orig;
		mpole.pass(calc_config, ps);

		std::cout << "Ps after kick" << ps << std::endl;
		/* Length 0 -> harmonics turn integral ? */
		BOOST_CHECK_CLOSE(ps[x_],   -.5e-3, 1e-12);
		BOOST_CHECK_CLOSE(ps[px_],   -1e-3, 1e-12);

		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    2e-7);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

	}
}
