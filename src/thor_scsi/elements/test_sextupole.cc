#define BOOST_TEST_MODULE sextupole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/core/field_interpolation.h>
#include <thor_scsi/core/multipoles.h>
#include <thor_scsi/elements/sextupole.h>
#include "check_multipole.h"
#include <tps/ss_vect.h>
#include <ostream>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;


static void check_only_sext_set(std::shared_ptr<tsc::PlanarMultipoles> muls, const tsc::cdbl ref)
{
	check_only_major_multipole_set(muls, ref, 3);
}


BOOST_AUTO_TEST_CASE(test10_sextupole_ostream)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("K", 0.0);

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
	const tsc::cdbl ref(chroma, 0e0);

	auto sext = tse::SextupoleType(C);
	BOOST_CHECK(sext.getMainMultipoleNumber() == 2);
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

		const tsc::cdbl ref(0e0, 0e0);
		check_only_sext_set(sext.getMultipoles(), ref);
	}
	{
		auto sext = tse::SextupoleType(C);
		sext.setMainMultipoleStrength(tsc::cdbl(0, 0));
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
		const tsc::cdbl ref(val, 0e0);
		check_only_sext_set(sext.getMultipoles(), ref);
	}

        /* check that setting double works */
	{
		auto sext = tse::SextupoleType(C);
		const double val = 3e0;
		const tsc::cdbl ref(val, 0e0);
		sext.setMainMultipoleStrength(ref);
		auto c = sext.getMainMultipoleStrength();
		BOOST_CHECK_CLOSE(c.real(),  val, 1e-12);
		BOOST_CHECK_SMALL(c.imag(),       1e-12);

		check_only_sext_set(sext.getMultipoles(), ref);
	}
	{
		auto sext = tse::SextupoleType(C);
		const double val = 42e0;
		const tsc::cdbl ref(0e0, val);
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

	{
		auto sext = tse::SextupoleType(C);
		sext.setMainMultipoleStrength(0.0);

		/* no multipole */
		for(int i = -1; i <= 1; ++i){
			ss_vect<double> ps;
			const double x = i;
			ps[x_] = x;


			sext.pass(calc_config, ps);

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

		sext.setMainMultipoleStrength(tsc::cdbl(1e0/i, 0e0));
		ss_vect<double> ps;
		sext.pass(calc_config, ps);

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
		sext.setMainMultipoleStrength(tsc::cdbl(chroma, 0));
		for(int i = -1; i <= 1; ++i){
			if (i == 0){
				/* checked above */
				continue;
			}

			ss_vect<double> ps;
			const double x = i;
			ps[x_] = x;

			sext.pass(calc_config, ps);

			BOOST_CHECK_CLOSE(ps[px_],  - x* chroma, 1e-12);
			BOOST_CHECK_CLOSE(ps[x_],           x, 1e-14);
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
		sext.setMainMultipoleStrength(tsc::cdbl(0, chroma));

		for(int i = -1; i <= 1; ++i){
			if (i == 0){
				/* checked above */
				continue;
			}

			ss_vect<double> ps;
			const double x = i;
			ps[x_] = x;

			sext.pass(calc_config, ps);

			BOOST_CHECK_CLOSE(ps[py_], x * chroma, 1e-12);
			BOOST_CHECK_CLOSE(ps[x_],         x, 1e-14);
			BOOST_CHECK_SMALL(ps[y_],            1e-14);
			BOOST_CHECK_SMALL(ps[px_],           1e-14);
			BOOST_CHECK_SMALL(ps[ct_],           1e-14);
			BOOST_CHECK_SMALL(ps[delta_],        1e-14);
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

	{
		/* normal sextrupole */
		auto sext = tse::SextupoleType(C);
		sext.setNumberOfIntegrationSteps(1);
		// sextrupoles are always thin elements
		BOOST_CHECK(sext.isThick());
		BOOST_CHECK_SMALL(sext.getInverseRigidity(), 1e-12);
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
			ss_vect<double> ps;

			ps[x_] = xs;
			sext.pass(calc_config, ps);

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
