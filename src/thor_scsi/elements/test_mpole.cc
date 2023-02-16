#define BOOST_TEST_MODULE mpole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/elements/mpole.h>
#include <ostream>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

typedef typename tsc::StandardDoubleType::complex_type cdbl;


auto a_desc = std::make_shared<gtpsa::desc>(1, 1);
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

#if 0
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
		gtpsa::ss_vect<double> ps = ps_orig.clone();
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
		gtpsa::ss_vect<double> ps = ps_orig.clone();
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
		gtpsa::ss_vect<double> ps = ps_orig.clone();
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
	mpole.getFieldInterpolator()->setMultipole(1, cdbl(1e0,0e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	// mpole.show(std::cout, 4);

	const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		gtpsa::ss_vect<double> ps = ps_orig.clone();
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
	mpole.getFieldInterpolator()->setMultipole(1, cdbl(0, 1e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	// mpole.show(std::cout, 4);

	const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		gtpsa::ss_vect<double> ps = ps_orig.clone();
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
		mpole.getFieldInterpolator()->setMultipole(n, cdbl(1, 0));

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
			const cdbl p_expected = pow(cdbl(0e0, y), (n-1));

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
		cdbl t_mul = cdbl(0, 355e0/113e0/double(n));
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
			const cdbl p_expected = t_mul *  pow(cdbl(x, 0), (n-1));

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
			const cdbl p_expected = t_mul *  pow(cdbl(0, y), (n-1));

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
	mpole.getFieldInterpolator()->setMultipole(1, cdbl(0e0,0e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	std::cerr << "M pole "; mpole.show(std::cerr, 4); std::cerr<<std::endl;
	const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		gtpsa::ss_vect<double> ps = ps_orig.clone();
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
	mpole.getFieldInterpolator()->setMultipole(1, cdbl(0e0,0e0));

	boost::test_tools::output_test_stream output;
	mpole.show(output, 4);
	BOOST_CHECK( !output.is_empty( false ) );

	std::cerr << "M pole "; mpole.show(std::cerr, 4); std::cerr<<std::endl;
	gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	ps_orig[delta_] = delta;
	{
		gtpsa::ss_vect<double> ps = ps_orig.clone();
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
#endif
/*
 * test mpole for parameters
 */
BOOST_AUTO_TEST_CASE(test30_mpole_engineering)
{
    tsc::ConfigType calc_config;
    Config C;
    C.set<std::string>("name", "test");
    C.set<double>("N", 1);

    C.set<double>("L", 0e0);
    tse::MpoleTypeTpsa mpole(C);

    // not yet using knobs ... to be explored
    int nv = 8;
    auto a_desc = std::make_shared<gtpsa::desc>(nv, 5);
    // auto L = gtpsa::tpsa(a_desc, 1);
    //mpole.setLength(L);
    auto t = gtpsa::CTpsaOrComplex(0e0);
    //tsc::TwoDimensionalMultipolesKnobbed<tsc::TpsaVariantType> h(t);
    auto h = std::make_shared<tsc::TwoDimensionalMultipolesTpsa>(t, 20);

    auto c3 = gtpsa::ctpsa(a_desc, 3);
    c3.set({0,0},{1e-4, 0});
    // a decapole
    auto c5 = gtpsa::ctpsa(a_desc, 3);
    c5.set({0,0},{3e-4, 0});

    // gradient for these variables
    c3.setv(1, {0,0, 0,0, 0,0, 1,0});
    c5.setv(1, {0,0, 0,0, 0,0, 0,1});

    std::complex<double> c1 (1e0, 0e0);
    h->setMultipole(1, c1);
    h->setMultipole(3, c3);
    h->setMultipole(5, c5);

    mpole.setFieldInterpolator(h);

    gtpsa::ss_vect<gtpsa::tpsa> ss_vect(a_desc, 2, 9);
    ss_vect.set_identity();
    ss_vect[0].setName("x");
    ss_vect[1].setName("px");
    ss_vect[2].setName("y");
    ss_vect[3].setName("py");
    ss_vect[4].setName("delta");
    ss_vect[5].setName("ct");
    ss_vect[6].setName("c3");
    ss_vect[7].setName("c5");

    // offset and angle in the horizontal plane
    ss_vect[0].set(0, 1e-3);
    ss_vect[1].set(0, 1e-3);

    calc_config = tsc::ConfigType();
    mpole.propagate(calc_config, ss_vect);
    double eps = 1e-12;
    for(size_t i =0; i<ss_vect.size(); ++i){
        auto& t = ss_vect[i];
        t.print(0, eps);
    }
}
