#define BOOST_TEST_MODULE quadrupole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/core/field_interpolation.h>
#include <thor_scsi/core/multipoles.h>
#include <thor_scsi/elements/quadrupole.h>
#include "check_multipole.h"
#include <ostream>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;
typedef typename tsc::StandardDoubleType::complex_type cdbl;

auto desc = std::make_shared<gtpsa::desc>(6, 6);
const auto t_ref = gtpsa::tpsa(desc, 1);

template<class C>
static void check_only_quad_set(std::shared_ptr<tsc::TwoDimensionalMultipolesKnobbed<C>> muls, const cdbl ref)
{
	auto check = CheckMultipoles<C>();
    check.only_major_multipole_set(muls, ref, 2);
}

BOOST_AUTO_TEST_CASE(test10_quadrupole_ostream)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("K", 0.0);
	C.set<double>("L", 0e0);
	C.set<double>("N", 1);

	auto quad = tse::QuadrupoleType(C);

	boost::test_tools::output_test_stream output;
	output << quad;
	BOOST_CHECK( !output.is_empty( false ) );
}

BOOST_AUTO_TEST_CASE(test10_quadrupole_K)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	const double grad = 27e0;
	C.set<double>("K", grad);
	C.set<double>("L", 0e0);
	C.set<double>("N", 1);

	const cdbl ref(grad, 0e0);
	auto quad = tse::QuadrupoleType(C);
	BOOST_CHECK(quad.getMainMultipoleNumber() == 2);
	// quadrupoles are always thin elements
	BOOST_CHECK(quad.isThick() == false);
	check_only_quad_set(quad.getMultipoles(), ref);

	auto c = quad.getMainMultipoleStrength();
	BOOST_CHECK_CLOSE(c.real(), grad, 1e-12);
	BOOST_CHECK_SMALL(c.imag(),       1e-12);

}


BOOST_AUTO_TEST_CASE(test10_quadrupole_setMul)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	const double grad = 27e0;
	C.set<double>("K", grad);
	C.set<double>("L", 0e0);
	C.set<double>("N", 1);

	/* check that setting 0 gets set */
	{
		auto quad = tse::QuadrupoleType(C);
		{
			auto val = quad.getMainMultipoleStrength();
			BOOST_CHECK_CLOSE(val.real(), grad, 1e-12);
			BOOST_CHECK_SMALL(val.imag(),       1e-12);
		}

		quad.setMainMultipoleStrength(0);
		{
			auto val = quad.getMainMultipoleStrength();
			BOOST_CHECK_SMALL(val.real(),  1e-12);
			BOOST_CHECK_SMALL(val.imag(),  1e-12);
		}

		const cdbl ref(0e0, 0e0);
		check_only_quad_set(quad.getMultipoles(), ref);
	}
	{
		auto quad = tse::QuadrupoleType(C);
		quad.setMainMultipoleStrength(cdbl(0, 0));
		auto val = quad.getMainMultipoleStrength();
		BOOST_CHECK_SMALL(val.real(),  1e-12);
		BOOST_CHECK_SMALL(val.imag(),  1e-12);
	}

	/* check that setting to a double value  gets set */
	{
		auto quad = tse::QuadrupoleType(C);
		{
			auto val = quad.getMainMultipoleStrength();
			BOOST_CHECK_CLOSE(val.real(), grad, 1e-12);
			BOOST_CHECK_SMALL(val.imag(),       1e-12);
		}

		const double val = -7e0;
		quad.setMainMultipoleStrength(val);
		{
			auto c = quad.getMainMultipoleStrength();
			BOOST_CHECK_CLOSE(c.real(),  val,  1e-12);
			BOOST_CHECK_SMALL(c.imag(),        1e-12);
		}
		const cdbl ref(val, 0e0);
		check_only_quad_set(quad.getMultipoles(), ref);
	}

        /* check that setting double works */
	{
		auto quad = tse::QuadrupoleType(C);
		const double val = 3e0;
		const cdbl ref(val, 0e0);
		quad.setMainMultipoleStrength(ref);
		auto c = quad.getMainMultipoleStrength();
		BOOST_CHECK_CLOSE(c.real(),  val, 1e-12);
		BOOST_CHECK_SMALL(c.imag(),       1e-12);

		check_only_quad_set(quad.getMultipoles(), ref);
	}
	{
		auto quad = tse::QuadrupoleType(C);
		const double val = 42e0;
		const cdbl ref(0e0, val);
		quad.setMainMultipoleStrength(-ref);
		auto c = quad.getMainMultipoleStrength();
		BOOST_CHECK_SMALL(c.real(),      1e-12);
		BOOST_CHECK_CLOSE(c.imag(), -val, 1e-12);

		check_only_quad_set(quad.getMultipoles(), -ref);
	}

}


BOOST_AUTO_TEST_CASE(test20_quadrupole_thin_eval)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("K", 0e0);
	C.set<double>("L", 0e0);
	C.set<double>("N", 1);

	{
		auto quad = tse::QuadrupoleType(C);
		quad.setMainMultipoleStrength(0.0);

		/* no multipole */
		for(int i = -1; i <= 1; ++i){
			gtpsa::ss_vect<double> ps_ref{0,0,0,0,0,0};
			const double x = i;

			gtpsa::ss_vect<double> ps{0,0,0,0,0,0};
			ps[x_] = x;

			quad.propagate(calc_config, ps);

			BOOST_CHECK_CLOSE(ps[x_],      x, 1e-14);
			BOOST_CHECK_SMALL(ps[y_],         1e-14);
			BOOST_CHECK_SMALL(ps[px_],        1e-14);
			BOOST_CHECK_SMALL(ps[py_],        1e-14);
			BOOST_CHECK_SMALL(ps[ct_],        1e-14);
			BOOST_CHECK_SMALL(ps[delta_],     1e-14);


			gtpsa::ss_vect<gtpsa::tpsa> ss_vect(t_ref);
			ss_vect[x_] = x;

			quad.propagate(calc_config, ss_vect);

			auto cst = ss_vect.cst();
			BOOST_CHECK_CLOSE(cst[x_],     x, 1e-14);
			BOOST_CHECK_SMALL(cst[y_],        1e-14);
			BOOST_CHECK_SMALL(cst[px_],       1e-14);
			BOOST_CHECK_SMALL(cst[py_],       1e-14);
			BOOST_CHECK_SMALL(cst[ct_],       1e-14);
			BOOST_CHECK_SMALL(cst[delta_],    1e-14);

		}
	}

	for(int i = -1; i <= 1; ++i){
		auto quad = tse::QuadrupoleType(C);
		if (i == 0){
			/* checked above */
			continue;
		}

		quad.setMainMultipoleStrength(cdbl(1e0/i, 0e0));
		gtpsa::ss_vect<double> ps{0,0,0,0,0,0};;
		quad.propagate(calc_config, ps);

		BOOST_CHECK_SMALL(ps[x_],      1e-14);
		BOOST_CHECK_SMALL(ps[y_],      1e-14);
		BOOST_CHECK_SMALL(ps[px_],     1e-14);
		BOOST_CHECK_SMALL(ps[py_],     1e-14);
		BOOST_CHECK_SMALL(ps[ct_],     1e-14);
		BOOST_CHECK_SMALL(ps[delta_],  1e-14);

		gtpsa::ss_vect<gtpsa::tpsa> ss_vect(t_ref);
		auto cst = ss_vect.cst();
		quad.propagate(calc_config, ss_vect);
		BOOST_CHECK_SMALL(cst[x_],     1e-14);
		BOOST_CHECK_SMALL(cst[y_],     1e-14);
		BOOST_CHECK_SMALL(cst[px_],    1e-14);
		BOOST_CHECK_SMALL(cst[py_],    1e-14);
		BOOST_CHECK_SMALL(cst[ct_],    1e-14);
		BOOST_CHECK_SMALL(cst[delta_], 1e-14);

	}

	{
		/* normal quadrupole */
		const double grad = 355e0 / 113e0;
		auto quad = tse::QuadrupoleType(C);
		quad.setMainMultipoleStrength(cdbl(grad, 0));
		for(int i = -1; i <= 1; ++i){
			if (i == 0){
				/* checked above */
				continue;
			}

			const double x = i;
			std::vector<double> vec2 {x, 0e0, 0e0, 0e0, 0e0, 0e0};
			gtpsa::ss_vect<double> ps(vec2);

			quad.propagate(calc_config, ps);

			BOOST_CHECK_CLOSE(ps [px_],  - x* grad, 1e-12);
			BOOST_CHECK_CLOSE(ps [x_],           x, 1e-12);
			BOOST_CHECK_SMALL(ps [y_],              1e-12);
			BOOST_CHECK_SMALL(ps [py_],             1e-12);
			BOOST_CHECK_SMALL(ps [ct_],             1e-12);
			BOOST_CHECK_SMALL(ps [delta_],          1e-12);


			gtpsa::ss_vect<gtpsa::tpsa> ss_vect(t_ref);
			ss_vect.set_identity();
			ss_vect[x_] = x;
			BOOST_CHECK_CLOSE(ss_vect[x_].cst(), x, 1e-12);

			quad.propagate(calc_config, ss_vect);

			auto cst = ss_vect.cst();

			BOOST_CHECK_CLOSE(cst[px_],  - x* grad, 1e-12);
			BOOST_CHECK_CLOSE(cst[x_],           x, 1e-12);
			BOOST_CHECK_SMALL(cst[y_],              1e-12);
			BOOST_CHECK_SMALL(cst[py_],             1e-12);
			BOOST_CHECK_SMALL(cst[ct_],             1e-12);
			BOOST_CHECK_SMALL(cst[delta_],          1e-12);

			// test that the gradient is found in the first derivative
			std::vector<num_t> vec(ss_vect.size());
			// first order
			ss_vect[py_].getv(1, &vec);
			BOOST_CHECK_CLOSE(vec[y_],      grad, 1e-12);
			BOOST_CHECK_SMALL(vec[x_],            1e-12);
			BOOST_CHECK_SMALL(vec[px_],           1e-12);
			BOOST_CHECK_CLOSE(vec[py_],      1e0, 1e-12);
			BOOST_CHECK_SMALL(vec[ct_],           1e-12);
			BOOST_CHECK_SMALL(vec[delta_],        1e-12);
		}
	}

	{
		/* skew quadrupole */
		const double grad = 1/28.0;
		auto quad = tse::QuadrupoleType(C);
		quad.setMainMultipoleStrength(cdbl(0, grad));

		for(int i = -1; i <= 1; ++i){
			if (i == 0){
				/* checked above */
				continue;
			}

			const double x = i;
			std::vector<double> vec2 {x, 0e0, 0e0, 0e0, 0e0, 0e0};
			gtpsa::ss_vect<double> ps(vec2);

			quad.propagate(calc_config, ps);

			BOOST_CHECK_CLOSE(ps [py_], x * grad, 1e-12);
			BOOST_CHECK_CLOSE(ps [x_],         x, 1e-12);
			BOOST_CHECK_SMALL(ps [y_],            1e-12);
			BOOST_CHECK_SMALL(ps [px_],           1e-12);
			BOOST_CHECK_SMALL(ps [ct_],           1e-12);
			BOOST_CHECK_SMALL(ps [delta_],        1e-12);


			gtpsa::ss_vect<gtpsa::tpsa> ss_vect(t_ref);
			ss_vect.set_identity();
			ss_vect[x_] = x;

			quad.propagate(calc_config, ss_vect);
			auto cst = ss_vect.cst();

			BOOST_CHECK_CLOSE(cst[py_], x * grad, 1e-12);
			BOOST_CHECK_CLOSE(cst[x_],         x, 1e-12);
			BOOST_CHECK_SMALL(cst[y_],            1e-12);
			BOOST_CHECK_SMALL(cst[px_],           1e-12);
			BOOST_CHECK_SMALL(cst[ct_],           1e-12);
			BOOST_CHECK_SMALL(cst[delta_],        1e-12);

			// test that the gradient is found in the first derivative
			std::vector<num_t> vec(ss_vect.size());
			// first order
			ss_vect[py_].getv(1, &vec);
			BOOST_CHECK_CLOSE(vec[x_],      grad, 1e-12);
			BOOST_CHECK_SMALL(vec[px_],           1e-12);
			BOOST_CHECK_SMALL(vec[y_],            1e-12);
			BOOST_CHECK_CLOSE(vec[py_],        1, 1e-12);
			BOOST_CHECK_SMALL(vec[ct_],           1e-12);
			BOOST_CHECK_SMALL(vec[delta_],        1e-12);
		}
	}
}

BOOST_AUTO_TEST_CASE(test20_quadrupole_typical_length_eval)
{

	const double gdl = 5e0, length=1e-3;
	const double gradient = gdl / length;
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("K", gradient);
	C.set<double>("L", length);
	C.set<double>("N", 1);

	{
		/* normal quadrupole */
		auto quad = tse::QuadrupoleType(C);
		quad.setNumberOfIntegrationSteps(1);
		// quadrupoles are always thin elements
		BOOST_CHECK(quad.isThick());
		BOOST_CHECK_SMALL(quad.getCurvature(), 1e-12);
		BOOST_CHECK(quad.assumingCurvedTrajectory() == false);

		BOOST_CHECK_CLOSE(quad.getLength(), length,   1e-12);
		auto c = quad.getMainMultipoleStrength();
		BOOST_CHECK_CLOSE(c.real(), gradient, 1e-12);
		BOOST_CHECK_SMALL(c.imag(),           1e-12);

		for(int i = -5; i <= 5; i+=5){
			if (i == 0){
				/* checked above */
				continue;
			}

			const double xs = i * (1e-3), l2 = length / 2e0;
			const double By = -xs * (gdl), xe = xs + By * l2;
			gtpsa::ss_vect<double> ps{0,0,0,0,0,0};


			ps[x_] = xs;
			quad.propagate(calc_config, ps);

			BOOST_CHECK_CLOSE(ps[px_],   By, 2);
			BOOST_CHECK_CLOSE(ps[x_],    xe, 0.5);
			BOOST_WARN_CLOSE( ps[px_],   By, 1.8);
			BOOST_WARN_CLOSE( ps[x_],    xe, 0.06);
			BOOST_CHECK_SMALL(ps[ct_],       4e-5);
			BOOST_CHECK_SMALL(ps[y_],        1e-14);
			BOOST_CHECK_SMALL(ps[py_],       1e-14);
			BOOST_CHECK_SMALL(ps[delta_],    1e-14);

			gtpsa::ss_vect<gtpsa::tpsa> ss_vect(t_ref);
			ss_vect[x_] = xs;
			quad.propagate(calc_config, ss_vect);

			auto cst = ss_vect.cst();

			BOOST_CHECK_CLOSE(cst[px_],   By, 2);
			BOOST_CHECK_CLOSE(cst[x_],    xe, 0.5);
			BOOST_WARN_CLOSE( cst[px_],   By, 1.8);
			BOOST_WARN_CLOSE( cst[x_],    xe, 0.06);
			BOOST_CHECK_SMALL(cst[ct_],       4e-5);
			BOOST_CHECK_SMALL(cst[y_],        1e-14);
			BOOST_CHECK_SMALL(cst[py_],       1e-14);
			BOOST_CHECK_SMALL(cst[delta_],    1e-14);

		}
	}
}
