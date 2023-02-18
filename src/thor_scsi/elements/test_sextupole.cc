#define BOOST_TEST_MODULE sextupole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/core/field_interpolation.h>
#include <thor_scsi/core/multipoles.h>
#include <thor_scsi/elements/sextupole.h>
#include "check_multipole.h"
#include <ostream>
#include <chrono>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

using cdbl = std::complex<double>;

template<class C>
static void check_only_sext_set(std::shared_ptr<tsc::TwoDimensionalMultipolesKnobbed<C>> muls, const cdbl ref)
{
    auto check = CheckMultipoles<C>();
	check.only_major_multipole_set(muls, ref, 3);
}


BOOST_AUTO_TEST_CASE(test10_sextupole_ostream)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("K", 0.0);
	C.set<double>("L", 0e0);
	C.set<double>("N", 1);

	auto sext = tse::SextupoleType(C);

	boost::test_tools::output_test_stream output;
	output << sext;
	BOOST_CHECK( !output.is_empty( false ) );
}

BOOST_AUTO_TEST_CASE(test10_sextupole_K)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	const double chroma = 27e0;
	C.set<double>("K", chroma);
	C.set<double>("L", 0e0);
	C.set<double>("N", 1);
	const cdbl ref(chroma, 0e0);

	auto sext = tse::SextupoleType(C);
	BOOST_CHECK(sext.getMainMultipoleNumber() == 3);
	// sextrupoles are always thin elements
	BOOST_CHECK(sext.isThick() == false);
	check_only_sext_set(sext.getMultipoles(), ref);

	auto c = sext.getMainMultipoleStrength();
	BOOST_CHECK_CLOSE(c.real(), chroma, 1e-12);
	BOOST_CHECK_SMALL(c.imag(),       1e-12);

}


BOOST_AUTO_TEST_CASE(test10_sextupole_setMul)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	const double chroma = 36.96e0;
	C.set<double>("K", chroma);
	C.set<double>("L", 0e0);
	C.set<double>("N", 1);

	/* check that setting 0 gets set */
	{
		auto sext = tse::SextupoleType(C);
		{
			auto val = sext.getMainMultipoleStrength();
			BOOST_CHECK_CLOSE(val.real(), chroma, 1e-12);
			BOOST_CHECK_SMALL(val.imag(),       1e-12);
		}

		sext.setMainMultipoleStrength(0);
		{
			auto val = sext.getMainMultipoleStrength();
			BOOST_CHECK_SMALL(val.real(),  1e-12);
			BOOST_CHECK_SMALL(val.imag(),  1e-12);
		}

		const cdbl ref(0e0, 0e0);
		check_only_sext_set(sext.getMultipoles(), ref);
	}
	{
		auto sext = tse::SextupoleType(C);
		sext.setMainMultipoleStrength(cdbl(0, 0));
		auto val = sext.getMainMultipoleStrength();
		BOOST_CHECK_SMALL(val.real(),  1e-12);
		BOOST_CHECK_SMALL(val.imag(),  1e-12);
	}

	/* check that setting to a double value  gets set */
	{
		auto sext = tse::SextupoleType(C);
		{
			auto val = sext.getMainMultipoleStrength();
			BOOST_CHECK_CLOSE(val.real(), chroma, 1e-12);
			BOOST_CHECK_SMALL(val.imag(),       1e-12);
		}

		const double val = -7e0;
		sext.setMainMultipoleStrength(val);
		{
			auto c = sext.getMainMultipoleStrength();
			BOOST_CHECK_CLOSE(c.real(),  val,  1e-12);
			BOOST_CHECK_SMALL(c.imag(),        1e-12);
		}
		const cdbl ref(val, 0e0);
		check_only_sext_set(sext.getMultipoles(), ref);
	}

        /* check that setting double works */
	{
		auto sext = tse::SextupoleType(C);
		const double val = 3e0;
		const cdbl ref(val, 0e0);
		sext.setMainMultipoleStrength(ref);
		auto c = sext.getMainMultipoleStrength();
		BOOST_CHECK_CLOSE(c.real(),  val, 1e-12);
		BOOST_CHECK_SMALL(c.imag(),       1e-12);

		check_only_sext_set(sext.getMultipoles(), ref);
	}
	{
		auto sext = tse::SextupoleType(C);
		const double val = 42e0;
		const cdbl ref(0e0, val);
		sext.setMainMultipoleStrength(-ref);
		auto c = sext.getMainMultipoleStrength();
		BOOST_CHECK_SMALL(c.real(),      1e-12);
		BOOST_CHECK_CLOSE(c.imag(), -val, 1e-12);

		check_only_sext_set(sext.getMultipoles(), -ref);
	}

}


BOOST_AUTO_TEST_CASE(test20_sextupole_thin_eval)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("K", 0e0);
	C.set<double>("L", 0e0);
	C.set<double>("N", 1);

	{
		auto sext = tse::SextupoleType(C);
		sext.setMainMultipoleStrength(0.0);

		/* no multipole */
		for(int i = -1; i <= 1; ++i){
			const double x = i;
			gtpsa::ss_vect<double> ps(x);
			ps[x_] = x;


			sext.propagate(calc_config, ps);

			BOOST_CHECK_CLOSE(ps[x_],     x, 1e-14);
			BOOST_CHECK_SMALL(ps[y_],        1e-14);
			BOOST_CHECK_SMALL(ps[px_],       1e-14);
			BOOST_CHECK_SMALL(ps[py_],       1e-14);
			BOOST_CHECK_SMALL(ps[ct_],       1e-14);
			BOOST_CHECK_SMALL(ps[delta_],    1e-14);
		}
	}

	for(int i = -1; i <= 1; ++i){
		auto sext = tse::SextupoleType(C);
		if (i == 0){
			/* checked above */
			continue;
		}

		sext.setMainMultipoleStrength(cdbl(1e0/i, 0e0));
		gtpsa::ss_vect<double> ps(0e0);
		sext.propagate(calc_config, ps);

		BOOST_CHECK_SMALL(ps[x_],     1e-14);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[px_],    1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);
	}

	{
		/* normal sextrupole */
		const double chroma = 355e0 / 113e0;
		auto sext = tse::SextupoleType(C);
		sext.setMainMultipoleStrength(cdbl(chroma, 0));
		for(int i = -1; i <= 1; ++i){
			if (i == 0){
				/* checked above */
				continue;
			}

			const double x = i;
			gtpsa::ss_vect<double> ps(x);
			ps[x_] = x;

			sext.propagate(calc_config, ps);

			// Todo: check sign!!!!
			const double xc =  - (x* chroma) * i;

			// check for signAPIAPI
			BOOST_CHECK_CLOSE(ps[x_],           x, 1e-14);
			BOOST_CHECK_CLOSE(ps[px_],         xc, 1e-14);
			BOOST_CHECK_SMALL(ps[y_],              1e-14);
			BOOST_CHECK_SMALL(ps[py_],             1e-14);
			BOOST_CHECK_SMALL(ps[ct_],             1e-14);
			BOOST_CHECK_SMALL(ps[delta_],          1e-14);
		}
	}

	{
		/* skew sextrupole */
		const double chroma = 1/28.0;
		auto sext = tse::SextupoleType(C);
		sext.setMainMultipoleStrength(cdbl(0, chroma));

		for(int i = -1; i <= 1; ++i){
			if (i == 0){
				/* checked above */
				continue;
			}

			const double x = i;
			gtpsa::ss_vect<double> ps(x);
			ps[x_] = x;

			sext.propagate(calc_config, ps);

			// Todo: check sign!!!!
			const double yc = x * chroma * i, pspy = ps[py_];
			std::cerr << "ps[px_] "  << ps[px_] << " x  " << x  << std::endl;
			std::cerr << "ps[py_] "  << ps[py_] << " yc " << yc << std::endl;

			BOOST_CHECK_CLOSE(ps[py_], x * chroma * i, 1e-12);
			BOOST_CHECK_CLOSE(ps[x_],              x, 1e-14);
			BOOST_CHECK_SMALL(ps[y_],              1e-14);
			BOOST_CHECK_SMALL(ps[px_],             1e-14);
			BOOST_CHECK_SMALL(ps[ct_],             1e-14);
			BOOST_CHECK_SMALL(ps[delta_],          1e-14);

			BOOST_CHECK_CLOSE(pspy, yc, 1e-12);

		}
	}
}

BOOST_AUTO_TEST_CASE(test20_sextupole_typical_length_eval)
{

	const double cdl = 6.336, length=0.16;
	const double chroma = cdl / length;
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("K", chroma);
	C.set<double>("L", length);
	C.set<double>("N", 1);

	{
		/* normal sextrupole */
		auto sext = tse::SextupoleType(C);
		sext.setNumberOfIntegrationSteps(1);
		// sextrupoles are always thin elements
		BOOST_CHECK(sext.isThick());
		BOOST_CHECK_SMALL(sext.getCurvature(), 1e-12);
		BOOST_CHECK(sext.assumingCurvedTrajectory() == false);

		BOOST_CHECK_CLOSE(sext.getLength(), length,   1e-12);
		auto c = sext.getMainMultipoleStrength();
		BOOST_CHECK_CLOSE(c.real(), chroma, 1e-12);
		BOOST_CHECK_SMALL(c.imag(),         1e-12);

		for(int i = -5; i <= 5; i+=5){
			if (i == 0){
				/* checked above */
				continue;
			}

			const double xs = i * (1e-3), l2 = length / 2e0;
			const double By = -(xs * xs) * (cdl), xe = xs + By * l2;
			gtpsa::ss_vect<double> ps(xs);

			ps[x_] = xs;
			sext.propagate(calc_config, ps);

			BOOST_CHECK_CLOSE(ps[px_],  By, 2);
			BOOST_CHECK_CLOSE(ps[x_],   xe, 0.5);
			BOOST_WARN_CLOSE(ps[px_],   By, 1.8);
			BOOST_WARN_CLOSE(ps[x_],    xe, 0.06);
			BOOST_CHECK_SMALL(ps[ct_],       4e-5);
			BOOST_CHECK_SMALL(ps[y_],        1e-14);
			BOOST_CHECK_SMALL(ps[py_],       1e-14);
			BOOST_CHECK_SMALL(ps[delta_],    1e-14);
		}
	}
}

BOOST_AUTO_TEST_CASE(test30_sextupole_tpsa)
{

	const double cdl = 6.336, length=0.16;
	const double chroma = cdl / length;
	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("K", chroma);
	C.set<double>("L", length);
	C.set<double>("N", 1);

	const int mo = 9;
	auto desc = std::make_shared<gtpsa::desc>(9, mo);

	/* normal sextrupole */
	auto sext = tse::SextupoleTypeTpsa(C);
	auto muls = std::make_shared<tsc::TwoDimensionalMultipolesTpsa>(0e0, 20);
	sext.setFieldInterpolator(muls);
	sext.setNumberOfIntegrationSteps(1);
	auto K2 = gtpsa::ctpsa(desc, mo);
	K2.set(std::complex(0e0, 0e0), std::complex(6.3, 0e0));
	/* gradient on K2 */
	std::vector<cdbl> K2_grad({0, 0, 0,  0, 0, 0,  1, 0, 0});
	K2.setv(1, K2_grad);
	// sext.getMultipoles()->setMultipole(3, K2);
	sext.setMainMultipoleStrength(K2);

	auto c9 = gtpsa::ctpsa(desc, mo);
	c9.set(std::complex(1e0, 0e0), std::complex(1e-4, 0e0));
	/* gradient on K2 */
	std::vector<cdbl> c9_grad({0, 0, 0,  0, 0, 0,  0, 0, 1});

	c9.setv(1, c9_grad);
	sext.getMultipoles()->setMultipole(9, c9);

	auto dx = gtpsa::tpsa(desc, mo);
	dx.set(0, 1e-3);
	/* gradient on dx */
	dx.setv(1, {0, 0, 0,  0, 0, 0,  0, 1, 0});
	sext.getTransform()->setDx(dx);

	std::cout << sext.repr() << std::endl;

	gtpsa::ss_vect<gtpsa::tpsa> ssv(desc, mo);
	ssv.set_identity();

	std::vector<std::string> names = {"x", "px", "y", "py", "delta", "ct", "K2", "dx"};

	{
		for(size_t i=0; i<ssv.size(); ++i){
			ssv[i].print(names[i].c_str(), 1e-6);
		}
	}


	auto t_start = std::chrono::system_clock::now();
	sext.propagate(calc_config, ssv);
	auto t_end = std::chrono::system_clock::now();
	std::chrono::duration<double> t_diff = t_end - t_start;
	std::cout <<" Propgagte required " << t_diff.count() << " seconds" << std::endl;

	{
		for(size_t i=0; i<ssv.size(); ++i){
			ssv[i].print(names[i].c_str(), 1e-6);
		}
	}

}
