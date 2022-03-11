#define BOOST_TEST_MODULE mpole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/elements/mpole.h>
#include <ostream>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

#if 0
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

	mpole.show(std::cout, 4);

	const ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		ss_vect<double> ps = ps_orig;
		mpole.pass(calc_config, ps);

		/* Length 0 -> harmonics turn integral ? */
		/*
		 * todo: cross check the sign with coordinate system conventions
		 *       compare results to tracy
		 */
		BOOST_CHECK_CLOSE(ps[px_],    -1, 1e-14);

		BOOST_CHECK_SMALL(ps[x_],     1e-14);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

	}
}
#endif

/*
 * review if this test is sensible
 */
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
		BOOST_CHECK_CLOSE(ps[px_],   -1e-3, 1e-12);

		BOOST_CHECK_SMALL(ps[y_],     1e-12);
		BOOST_CHECK_SMALL(ps[x_],     1e-12);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    2e-7);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);
	}
}


BOOST_AUTO_TEST_CASE(test21_mpole_kick_dipole_component_thick_kick_polar_ideal)
{


	/*
	 * todo: check that different length give the same result
	 */
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");

	const double phi = 1e-3; // 1 mrad
	const double length = 1.0;
	const double irho = phi / length;

	const double px_expected = phi;
	const double x_expected = 1e0 / 2e0 * phi * length;

	C.set<double>("L", length);
	tse::MpoleType mpole(C);

        mpole.asThick(true);
	mpole.Pirho = irho;
	mpole.PN = 1;


	/* */
	(dynamic_cast<tsc::PlanarMultipoles*>(mpole.intp))->setMultipole(1, tsc::cdbl(0e0,0e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	std::cerr << "M pole "; mpole.show(std::cerr, 4); std::cerr<<std::endl;
	const ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		ss_vect<double> ps = ps_orig;
		mpole.pass(calc_config, ps);

		std::cerr << "ps out " << ps << std::endl;
		/* Length 0 -> harmonics turn integral ? */
		BOOST_CHECK_SMALL(ps[x_],     1e-14);
		BOOST_CHECK_SMALL(ps[px_],    1e-14);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

	}
}


BOOST_AUTO_TEST_CASE(test21_mpole_kick_dipole_component_thick_kick_off_momentum)
{
	/*
	 * todo: check that different length give the same result
	 */
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");

	// const double phi = 1e-3; // 1 mrad
	const double phi = 1/180e0 * M_PI; // 1 deg
	const double length = 1.0;
	const double rho = length / phi;

	// relative momentum deviation
	const double delta = 1e-3;

	const double x_expected = rho * (1 - cos(phi)) * 1e0/(1 + delta) * delta;
	const double px_expected = sin(phi) * delta;

	std::cerr  << "length " << length << " rho " << rho
		   <<  "phi " << phi << " compare " << rho * length << std::endl;
	C.set<double>("L", length);
	tse::MpoleType mpole(C);

        mpole.asThick(true);
	mpole.Pirho = 1/rho;
	mpole.PN = 10;

	calc_config.Cart_Bend = false;

	/* */
	(dynamic_cast<tsc::PlanarMultipoles*>(mpole.intp))->setMultipole(1, tsc::cdbl(0e0,0e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	std::cerr << "M pole "; mpole.show(std::cerr, 4); std::cerr<<std::endl;
	ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	ps_orig[delta_] = delta;
	{
		ss_vect<double> ps = ps_orig;
		std::cerr << "ps in " << ps << ", ps orig " << ps_orig << " delta " << delta << std::endl;
		mpole.pass(calc_config, ps);

		std::cerr << "ps out " << ps << std::endl;
		/* Length 0 -> harmonics turn integral ? */
		BOOST_CHECK_CLOSE(ps[x_],     x_expected,  1e-14);
		BOOST_CHECK_CLOSE(ps[px_],    px_expected, 1e-14);
		BOOST_CHECK_CLOSE(ps[delta_], delta,       1e-14);

		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);

		// BOOST_WARN_CLOSE(ps[x_],     x_expected, 1e-14);
		// BOOST_WARN_CLOSE(ps[px_],    px_expected, 1e-14);
	}
}
