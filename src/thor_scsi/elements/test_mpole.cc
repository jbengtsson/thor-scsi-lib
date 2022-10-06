#define BOOST_TEST_MODULE mpole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/elements/mpole.h>
#include <ostream>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

auto a_desc = gtpsa::desc(1, 0);
auto tpsa_ref = gtpsa::tpsa(a_desc, mad_tpsa_default);


BOOST_AUTO_TEST_CASE(test02_mpole_print)
{
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("Method", 4.0);
	C.set<double>("N", 1);
	C.set<double>("L", .2);
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
	C.set<double>("N", 1);
	C.set<double>("L", .2);
	//tse::MpoleType mpole(C);

	BOOST_CHECK_THROW(((tse::MpoleType)(C)), thor_scsi::NotImplemented);

}

BOOST_AUTO_TEST_CASE(test10_mpole_kick_zero)
{
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("N", 1);
	C.set<double>("L", 0e0);

	tse::MpoleType mpole(C);

	// initialised to 0 length by default
	BOOST_CHECK_EQUAL(mpole.isThick(), false);

	const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		gtpsa::ss_vect<double> ps = ps_orig;
		mpole.propagate(calc_config, ps);

		BOOST_CHECK_SMALL(ps[x_],     1e-14);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[px_],    1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

	}
}

BOOST_AUTO_TEST_CASE(test11_mpole_kick_longitudinal_one)
{
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("N", 1);
	C.set<double>("L", 1.0);
	tse::MpoleType mpole(C);


	BOOST_CHECK_CLOSE(mpole.getLength(), 1.0, 1e-14);

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	// is this correct for a thin kick
	const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		gtpsa::ss_vect<double> ps = ps_orig;
		mpole.propagate(calc_config, ps);

		BOOST_CHECK_SMALL(ps[x_],     1e-14);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[px_],    1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

		// std::cout << "thin kick ps " << ps << std::endl;
	}

	{
		gtpsa::ss_vect<double> ps = ps_orig;
		mpole.asThick(true);
		mpole.propagate(calc_config, ps);

		BOOST_CHECK_SMALL(ps[x_],     1e-14);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[px_],    1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

		//std::cout << "thick kick ps " << ps << std::endl;
	}
}


/*
 * test a orbit trim: horizontal
 */
BOOST_AUTO_TEST_CASE(test12_orbit_trim_horizontal)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("N", 1);
	C.set<double>("L", 0.0);
	tse::MpoleType mpole(C);

	/* */
	mpole.getFieldInterpolator()->setMultipole(1, tsc::cdbl(1e0,0e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	// mpole.show(std::cout, 4);

	const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		gtpsa::ss_vect<double> ps = ps_orig;
		mpole.propagate(calc_config, ps);

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

/*
 * test a orbit trim: vertical
 */
BOOST_AUTO_TEST_CASE(test13_orbit_trim_vertical)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("N", 1);
	C.set<double>("L", 0.0);
	tse::MpoleType mpole(C);

	/* */
	mpole.getFieldInterpolator()->setMultipole(1, tsc::cdbl(0, 1e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	// mpole.show(std::cout, 4);

	const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		gtpsa::ss_vect<double> ps = ps_orig;
		mpole.propagate(calc_config, ps);

		/* Length 0 -> harmonics turn integral ? */
		/*
		 * todo: cross check the sign with coordinate system conventions
		 *       compare results to tracy
		 */
		BOOST_CHECK_CLOSE(ps[py_],    1, 1e-14);

		BOOST_CHECK_SMALL(ps[x_],     1e-14);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[px_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);

	}
}

BOOST_AUTO_TEST_CASE(test14_higher_orders_normal_multipole)
{
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("N", 1);
	C.set<double>("L", 0.0);

	for (int n=2; n<=4; ++n){
		tse::MpoleType mpole(C);
		mpole.getFieldInterpolator()->setMultipole(n, tsc::cdbl(1, 0));

		BOOST_CHECK(mpole.isThick() == false);
		/* on axis */
		{
			gtpsa::ss_vect<double> ps = {0, 0, 0, 0, 0, 0};
			mpole.propagate(calc_config, ps);

			BOOST_CHECK_SMALL(ps[x_],     1e-14);
			BOOST_CHECK_SMALL(ps[px_],    1e-14);
			BOOST_CHECK_SMALL(ps[y_],     1e-14);
			BOOST_CHECK_SMALL(ps[py_],    1e-14);
			BOOST_CHECK_SMALL(ps[ct_],    1e-14);
			BOOST_CHECK_SMALL(ps[delta_], 1e-14);
		}

		/* along x */
		for (int xi=-10; xi<10; ++xi){

			if(xi == 0){
				/* special comparisons at zero */
				continue;
			}
			const double x = xi * 1e-3;
			/* Todo: check sign */
			const double px_expected = -pow(x, (n-1));

			gtpsa::ss_vect<double> ps = {x, 0, 0, 0, 0, 0};
			mpole.propagate(calc_config, ps);
			BOOST_CHECK_CLOSE(ps[x_],     x,           1e-14);
			BOOST_CHECK_CLOSE(ps[px_],    px_expected, 1e-14);

			BOOST_CHECK_SMALL(ps[y_],     1e-14);
			BOOST_CHECK_SMALL(ps[py_],    1e-14);
			BOOST_CHECK_SMALL(ps[ct_],    1e-14);
			BOOST_CHECK_SMALL(ps[delta_], 1e-14);
		}

		/* along y */
		for (int yi=-10; yi<10; ++yi){
			if(yi == 0){
				continue;
			}
			const double y = yi * 1e-3;

			/* Todo: check sign */
			const tsc::cdbl p_expected = pow(tsc::cdbl(0e0, y), (n-1));

			gtpsa::ss_vect<double> ps = {0, 0, y, 0, 0, 0};
			mpole.propagate(calc_config, ps);
			BOOST_CHECK_CLOSE(ps[y_],     y,                 1e-14);
			BOOST_CHECK_CLOSE(ps[px_],   -p_expected.real(), 1e-14);
			BOOST_CHECK_CLOSE(ps[py_],    p_expected.imag(), 1e-14);

			BOOST_CHECK_SMALL(ps[x_],     1e-14);
			BOOST_CHECK_SMALL(ps[ct_],    1e-14);
			BOOST_CHECK_SMALL(ps[delta_], 1e-14);
		}

	}
}

BOOST_AUTO_TEST_CASE(test15_higher_orders_skew_multipole)
{
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("N", 1);
	C.set<double>("L", 0.0);

	for (int n=2; n<=4; ++n){
		tse::MpoleType mpole(C);
		tsc::cdbl t_mul = tsc::cdbl(0, 355e0/113e0/double(n));
		mpole.getFieldInterpolator()->setMultipole(n, t_mul);

		BOOST_CHECK(mpole.isThick() == false);
		/* on axis */
		{
			gtpsa::ss_vect<double> ps = {0, 0, 0, 0, 0, 0};
			mpole.propagate(calc_config, ps);

			BOOST_CHECK_SMALL(ps[x_],     1e-14);
			BOOST_CHECK_SMALL(ps[px_],    1e-14);
			BOOST_CHECK_SMALL(ps[y_],     1e-14);
			BOOST_CHECK_SMALL(ps[py_],    1e-14);
			BOOST_CHECK_SMALL(ps[ct_],    1e-14);
			BOOST_CHECK_SMALL(ps[delta_], 1e-14);
		}

		/* along x */
		for (int xi=-10; xi<10; ++xi){

			if(xi == 0){
				/* special comparisons at zero */
				continue;
			}
			const double x = xi * 1e-3;
			/* Todo: check sign */
			const tsc::cdbl p_expected = t_mul *  pow(tsc::cdbl(x, 0), (n-1));

			gtpsa::ss_vect<double> ps = {x, 0, 0, 0, 0, 0};
			mpole.propagate(calc_config, ps);
			BOOST_CHECK_CLOSE(ps[x_],     x,                 1e-14);
			BOOST_CHECK_CLOSE(ps[px_],    p_expected.real(), 2e-12);
			BOOST_CHECK_CLOSE(ps[py_],    p_expected.imag(), 2e-12);

			BOOST_CHECK_SMALL(ps[y_],     1e-14);
			BOOST_CHECK_SMALL(ps[ct_],    1e-14);
			BOOST_CHECK_SMALL(ps[delta_], 1e-14);

			// accuracy found was 1.5e-16
			BOOST_WARN_CLOSE(ps[px_],     p_expected.real(), 2e-14);
			BOOST_WARN_CLOSE(ps[py_],     p_expected.imag(), 3e-14);
		}

		/* along y */
		for (int yi=-10; yi<10; ++yi){
			if(yi == 0){
				continue;
			}
			const double y = yi * 1e-3;

			/* Todo: check sign */
			const tsc::cdbl p_expected = t_mul *  pow(tsc::cdbl(0, y), (n-1));

			gtpsa::ss_vect<double> ps = {0, 0, y, 0, 0, 0};
			mpole.propagate(calc_config, ps);
			BOOST_CHECK_CLOSE(ps[y_],     y,                 1e-14);
			BOOST_CHECK_CLOSE(ps[px_],   -p_expected.real(), 2e-12);
			BOOST_CHECK_CLOSE(ps[py_],    p_expected.imag(), 2e-12);

			BOOST_CHECK_SMALL(ps[x_],     1e-14);
			BOOST_CHECK_SMALL(ps[ct_],    1e-14);
			BOOST_CHECK_SMALL(ps[delta_], 1e-14);

			// accuracy found was 1.5e-16
			BOOST_WARN_CLOSE(ps[px_],    -p_expected.real(), 2e-14);
			BOOST_WARN_CLOSE(ps[py_],     p_expected.imag(), 3e-14);
		}

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
	C.set<double>("N", 1);
	tse::MpoleType mpole(C);

        mpole.asThick(true);
	mpole.setCurvature(irho);
	mpole.setNumberOfIntegrationSteps(1);


	/* */
	mpole.getFieldInterpolator()->setMultipole(1, tsc::cdbl(0e0,0e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	std::cerr << "M pole "; mpole.show(std::cerr, 4); std::cerr<<std::endl;
	const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		gtpsa::ss_vect<double> ps = ps_orig;
		mpole.propagate(calc_config, ps);

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


/*
 * test thick kick for 1 mrad
 */
BOOST_AUTO_TEST_CASE(test21_mpole_kick_dipole_component_thick_kick_off_momentum)
{
	/*
	 * todo: check that different length give the same result
	 */
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("N", 1);

	const double phi = 1e-3; // 1 mrad
	// const double phi = 1/180e0 * M_PI; // 1 deg
	const double length = 1.0;
	const double rho = length / phi;

	// relative momentum deviation
	const double delta = 1e-3;

	const double x_expected = rho * (1 - cos(phi)) * 1e0/(1 + delta) * delta;
	const double px_expected = sin(phi) * delta;

	C.set<double>("L", length);
	tse::MpoleType mpole(C);

        mpole.asThick(true);
	mpole.setCurvature(1/rho);
	mpole.setNumberOfIntegrationSteps(1);

	calc_config.Cart_Bend = false;

	/* */
	mpole.getFieldInterpolator()->setMultipole(1, tsc::cdbl(0e0,0e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	std::cerr << "M pole "; mpole.show(std::cerr, 4); std::cerr<<std::endl;
	gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	ps_orig[delta_] = delta;
	{
		gtpsa::ss_vect<double> ps = ps_orig;
		mpole.propagate(calc_config, ps);

		/* Length 0 -> harmonics turn integral ? */
		BOOST_CHECK_CLOSE(ps[x_],     x_expected,  1e-6);
		BOOST_CHECK_CLOSE(ps[px_],    px_expected, 1e-6);
		BOOST_CHECK_CLOSE(ps[delta_], delta,       1e-14);
		BOOST_CHECK_SMALL(ps[ct_],                 1e-7);

		// values for 1 mrad
		BOOST_WARN_CLOSE(ps[x_],     x_expected,  1e-8);
		BOOST_WARN_CLOSE(ps[px_],    px_expected, 2e-8);

		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);

		// BOOST_WARN_CLOSE(ps[x_],     x_expected, 1e-14);
		// BOOST_WARN_CLOSE(ps[px_],    px_expected, 1e-14);
	}
}
