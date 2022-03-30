#define BOOST_TEST_MODULE thin_kick
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <flame/core/config.h>
#include <thor_scsi/elements/element_helpers.h>
#include <thor_scsi/core/multipoles.h>
#include <ostream>
#include <vector>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

typedef std::array<bool, 6> flags_t;
typedef std::array<double, 6> eps_rel_t;

/* formerly a function of element helpers */
template<typename T>
void thin_kick(const thor_scsi::core::ConfigType &conf,
	 const thor_scsi::core::Field2DInterpolation& intp,
	 const double L, const double h_bend, const double h_ref,
	 ss_vect<T> &ps)
{
 	T BxoBrho, ByoBrho;

	intp.field(ps[x_], ps[y_], &BxoBrho, &ByoBrho);
	tse::thin_kick(conf, BxoBrho, ByoBrho, L, h_bend, h_ref, ps);
}

static void inline test_ps_small(const ss_vect<double>& ps)
{

	static const double eps_percent = 1e-14;

	BOOST_CHECK_SMALL(ps[x_],     eps_percent);
	BOOST_CHECK_SMALL(ps[y_],     eps_percent);
	BOOST_CHECK_SMALL(ps[px_],    eps_percent);
	BOOST_CHECK_SMALL(ps[py_],    eps_percent);
	BOOST_CHECK_SMALL(ps[ct_],    eps_percent);
	BOOST_CHECK_SMALL(ps[delta_], eps_percent);
}


BOOST_AUTO_TEST_CASE(test01_wrong_thin_kick_L0)
{
	const ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	const double h_ref = 0e0, h_bend = 0e0;
	const double length = 0.0;

	tsc::ConfigType calc_config;
	const tsc::PlanarMultipoles muls;

	{
		ss_vect<double> ps = ps_orig;
		thin_kick(calc_config, muls, length, h_bend, h_ref, ps);
		test_ps_small(ps);

	}

	for(int i=1; i <= 4; ++i)
	{
		ss_vect<double> ps = ps_orig;
		tsc::PlanarMultipoles muls2 = muls.clone();
		muls2.setMultipole(i, tsc::cdbl(1.0, 0.0));

		thin_kick(calc_config, muls2, length, h_bend, h_ref, ps);
		test_ps_small(ps);
	}

	{
		// Leave an extra instance here ... in case an extra test is
		// required
		ss_vect<double> ps = ps_orig;
		tsc::PlanarMultipoles muls2 = muls.clone();
		muls2.setMultipole(5, tsc::cdbl(1.0, 0.0));

		thin_kick(calc_config, muls2, length, h_bend, h_ref, ps);
		test_ps_small(ps);
	}
}

BOOST_AUTO_TEST_CASE(test02_wrong_thin_kick_L0_ps_off)
{
	const double x = 1, y = 1;
	const ss_vect<double> ps_orig = {x, 0, y, 0, 0, 0};
	const double h_ref = 0e0, h_bend = 0e0;
	const double length = 0.0;

	tsc::ConfigType calc_config;
	const tsc::PlanarMultipoles muls;

	{
		ss_vect<double> ps = ps_orig;
		thin_kick(calc_config, muls, length, h_bend, h_ref, ps);

		BOOST_CHECK_CLOSE(ps[x_],     x, 1e-14);
		BOOST_CHECK_CLOSE(ps[y_],     y, 1e-14);
		BOOST_CHECK_SMALL(ps[px_],    1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);
	}

	for(int i=1; i <= 4; ++i)
	{
		ss_vect<double> ps = ps_orig;
		tsc::PlanarMultipoles muls2 = muls.clone();
		muls2.setMultipole(i, tsc::cdbl(1.0, 0.0));

		thin_kick(calc_config, muls2, length, h_bend, h_ref, ps);


		BOOST_CHECK_CLOSE(ps[x_],     x, 1e-14);
		BOOST_CHECK_CLOSE(ps[y_],     y, 1e-14);
		BOOST_CHECK_SMALL(ps[px_],    1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);
	}

}

BOOST_AUTO_TEST_CASE(test03_wrong_thin_kick_L0_ps_off)
{
	const double x = -3e-3, px = 7e-3, y = -5e-3, py = 11e-4, delta=13e-4, ct=17e-4;
	const ss_vect<double> ps_orig = {x, px, y, py, delta, ct};
	const double h_ref = 0e0, h_bend = 0e0;
	const double length = 0;

	tsc::ConfigType calc_config;
	const tsc::PlanarMultipoles muls;

	{
		ss_vect<double> ps = ps_orig;
		thin_kick(calc_config, muls, length, h_bend, h_ref, ps);

		BOOST_CHECK_CLOSE(ps[x_],     x,     1e-14);
		BOOST_CHECK_CLOSE(ps[px_],    px,    1e-14);
		BOOST_CHECK_CLOSE(ps[y_],     y,     1e-14);
		BOOST_CHECK_CLOSE(ps[py_],    py,    1e-14);
		BOOST_CHECK_CLOSE(ps[ct_],    ct,    1e-14);
		BOOST_CHECK_CLOSE(ps[delta_], delta, 1e-14);
	}

	for(int i=1; i <= 4; ++i)
	{
		ss_vect<double> ps = ps_orig;
		tsc::PlanarMultipoles muls2 = muls.clone();
		muls2.setMultipole(i, tsc::cdbl(1.0, 0.0));

		thin_kick(calc_config,  muls2, length, h_bend, h_ref, ps);
		BOOST_CHECK_CLOSE(ps[x_],     x,     1e-14);
		BOOST_CHECK_CLOSE(ps[px_],    px,    1e-14);
		BOOST_CHECK_CLOSE(ps[y_],     y,     1e-14);
		BOOST_CHECK_CLOSE(ps[py_],    py,    1e-14);
		BOOST_CHECK_CLOSE(ps[ct_],    ct,    1e-14);
		BOOST_CHECK_CLOSE(ps[delta_], delta, 1e-14);

	}

}


BOOST_AUTO_TEST_CASE(test10_thin_kick_L1_no_field)
{
	const ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	const double h_ref = 0e0, h_bend = 0e0;
	const double length = 1.0;

	tsc::ConfigType calc_config;
	const tsc::PlanarMultipoles muls;

	{
		ss_vect<double> ps = ps_orig;
		thin_kick(calc_config, muls, length, h_bend, h_ref, ps);
		test_ps_small(ps);

	}

	for(int i=1; i <= 4; ++i)
	{
		ss_vect<double> ps = ps_orig;
		tsc::PlanarMultipoles muls2 = muls.clone();

		thin_kick(calc_config, muls2, length / (double(i)), h_bend, h_ref, ps);
		test_ps_small(ps);
	}

}

BOOST_AUTO_TEST_CASE(test11_thin_kick_L1_no_ps_off)
{
	const double x = -28e-4, px = -464e-5, y = +37e-4, py = 31e-4, delta=29e-4, ct=23e-4;
	const ss_vect<double> ps_orig = {x, px, y, py, delta, ct};
	const double h_ref = 0e0, h_bend = 0e0;
	const double length = 355/113;

	tsc::ConfigType calc_config;
	const tsc::PlanarMultipoles muls;

	for(int i=1; i <= 4; ++i)
	{
		ss_vect<double> ps = ps_orig;
		tsc::PlanarMultipoles muls2 = muls.clone();

		thin_kick(calc_config, muls2, length * i, h_bend, h_ref, ps);

		BOOST_CHECK_CLOSE(ps[x_],     x,     1e-14);
		BOOST_CHECK_CLOSE(ps[px_],    px,    1e-14);
		BOOST_CHECK_CLOSE(ps[y_],     y,     1e-14);
		BOOST_CHECK_CLOSE(ps[py_],    py,    1e-14);
		BOOST_CHECK_CLOSE(ps[ct_],    ct,    1e-14);
		BOOST_CHECK_CLOSE(ps[delta_], delta, 1e-14);

	}

}

BOOST_AUTO_TEST_CASE(test20_thin_kick_L1_dipole)
{
	tsc::ConfigType calc_config;
	tsc::PlanarMultipoles muls;

	muls.setMultipole(1, tsc::cdbl(1e-3, 0e0));
	const double h_ref = 0e0, h_bend = 0e0;

	// length is set to zero so that the multipoles are treated as integrals ....
	const double length = 1.0;

	const ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
	{
		ss_vect<double> ps = ps_orig;

		thin_kick(calc_config, muls, length, h_bend, h_ref, ps);

		BOOST_CHECK_CLOSE(ps[px_],   -1e-3, 1e-12);

		// Small kick should x deviate or not ?
		BOOST_CHECK_SMALL(ps[x_],     1e-12);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    2e-7);
		BOOST_CHECK_SMALL(ps[delta_], 1e-14);
	}
}


static void compute_irho_px(const double phi, const double length, const double delta,
			    double *irho, double *px)
{
	*irho = phi/length;
	*px = *irho * length * delta;
}

BOOST_AUTO_TEST_CASE(test51_sector_bend_delta)
{
	tsc::ConfigType calc_config;
	const tsc::PlanarMultipoles muls;


	const double phi = 1e-3; // 1 mrad
	const double length = 1.0;
	// relative momentum deviation
	const double delta = 1e-3;

	double irho, px_expected;
	compute_irho_px(phi, length, delta, &irho, &px_expected);
	const double h_bend = irho, h_ref = irho;


	const ss_vect<double> ps_orig = {0, 0, 0, 0, delta, 0};
	BOOST_CHECK_CLOSE(ps_orig[delta_],  delta, 1e-14);

	{
		ss_vect<double> ps = ps_orig;
		// std::cerr << "ps in " << ps << ", ps orig " << ps_orig << " delta " << delta << std::endl;
		thin_kick(calc_config,  muls, length, h_bend, h_ref, ps);

		// std::cerr << "ps out " << ps << std::endl;
		/* Length 0 -> harmonics turn integral ? */
		BOOST_CHECK_CLOSE(ps[px_],    px_expected, 1e-23);
		BOOST_CHECK_CLOSE(ps[delta_], delta,       1e-14);

		BOOST_CHECK_SMALL(ps[x_],     1e-24);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
		BOOST_CHECK_SMALL(ps[ct_],    1e-14);

	}
}


BOOST_AUTO_TEST_CASE(test52_sector_bend_delta_pars)
{
	tsc::ConfigType calc_config;
	const tsc::PlanarMultipoles muls;

	/* test different length */
	{
		// relative momentum deviation
		const double delta = 1e-3;
		const ss_vect<double> ps_orig = {0, 0, 0, 0, delta, 0};
		BOOST_CHECK_CLOSE(ps_orig[delta_],  delta, 1e-14);

		const double phi = 1e-3; // 1 mrad

		for(int i = 1; i < 10; ++i)
		{
			const double length = 355e0/113e0 /(double(i));

			ss_vect<double> ps = ps_orig;

			double irho, px_expected;
			compute_irho_px(phi, length, delta, &irho, &px_expected);
			const double h_bend = irho, h_ref = irho;

			thin_kick(calc_config,  muls, length, h_bend, h_ref, ps);

			// std::cerr << "ps out " << ps << std::endl;
			/* Length 0 -> harmonics turn integral ? */
			BOOST_CHECK_CLOSE(ps[px_],    px_expected, 1e-12);
			BOOST_CHECK_CLOSE(ps[delta_], delta,       1e-12);

			BOOST_CHECK_SMALL(ps[x_],     1e-12);
			BOOST_CHECK_SMALL(ps[y_],     1e-14);
			BOOST_CHECK_SMALL(ps[py_],    1e-14);
			BOOST_CHECK_SMALL(ps[ct_],    1e-14);

		}
	}

	/* test different momenta */
	{

		const double length = 7e0/22e0;
		const double phi = 1e-3; // 1 mrad

		int last_i = 1;
		for(int i = 2; i < 10; ++i)
		{

			// relative momentum deviation
			const double delta = 1e-3 * double(last_i)/double(i);
			last_i = i;

			const ss_vect<double> ps_orig = {0, 0, 0, 0, delta, 0};
			BOOST_CHECK_CLOSE(ps_orig[delta_],  delta, 1e-14);
			ss_vect<double> ps = ps_orig;

			double irho, px_expected;
			compute_irho_px(phi, length, delta, &irho, &px_expected);
			const double h_bend = irho, h_ref = irho;

			thin_kick(calc_config,  muls, length, h_bend, h_ref, ps);

			// std::cerr << "ps out " << ps << std::endl;
			/* Length 0 -> harmonics turn integral ? */
			BOOST_CHECK_CLOSE(ps[px_],    px_expected, 1e-12);
			BOOST_CHECK_CLOSE(ps[delta_], delta,       1e-12);

			BOOST_CHECK_SMALL(ps[x_],     1e-12);
			BOOST_CHECK_SMALL(ps[y_],     1e-14);
			BOOST_CHECK_SMALL(ps[py_],    1e-14);
			BOOST_CHECK_SMALL(ps[ct_],    1e-14);

		}
	}

	/* test different angles */
	{
		const double length = 1.98;
		const double delta = 41e-4;

		const ss_vect<double> ps_orig = {0, 0, 0, 0, delta, 0};
		BOOST_CHECK_CLOSE(ps_orig[delta_],  delta, 1e-14);

		for(int i = 2; i < 10; ++i)
		{

			const double phi = 5e-3 / double(i);
			double irho, px_expected;
			compute_irho_px(phi, length, delta, &irho, &px_expected);
			const double h_bend = irho, h_ref = irho;

			ss_vect<double> ps = ps_orig;

			thin_kick(calc_config,  muls, length, h_bend, h_ref, ps);

			// std::cerr << "ps out " << ps << std::endl;
			/* Length 0 -> harmonics turn integral ? */
			BOOST_CHECK_CLOSE(ps[px_],    px_expected, 1e-12);
			BOOST_CHECK_CLOSE(ps[delta_], delta,       1e-12);

			BOOST_CHECK_SMALL(ps[x_],     1e-12);
			BOOST_CHECK_SMALL(ps[y_],     1e-14);
			BOOST_CHECK_SMALL(ps[py_],    1e-14);
			BOOST_CHECK_SMALL(ps[ct_],    1e-14);

		}

	}
}


BOOST_AUTO_TEST_CASE(test60_higher_orders)
{
	tsc::ConfigType calc_config;
	Config C;
	const double length = 1, h_bend = 0.0, h_ref = 0.0;

	for (int n=2; n<=4; ++n){
		tsc::PlanarMultipoles muls;
		muls.setMultipole(n, tsc::cdbl(1, 0));

		/* on axis */
		{
			ss_vect<double> ps = {0, 0, 0, 0, 0, 0};
			thin_kick(calc_config, muls, length, h_bend, h_ref, ps);
			BOOST_CHECK_SMALL(ps[x_],     1e-14);
			BOOST_CHECK_SMALL(ps[px_],    1e-14);
			BOOST_CHECK_SMALL(ps[y_],     1e-14);
			BOOST_CHECK_SMALL(ps[py_],    1e-14);
			BOOST_CHECK_SMALL(ps[ct_],    1e-14);
			BOOST_CHECK_SMALL(ps[delta_], 1e-14);
		}

		for (int xi=1; xi<10; ++xi){
			const double x = xi * 1e-3;
			/* Todo: check sign */
			const double px_expected = -pow(x, (n-1));

			ss_vect<double> ps = {x, 0, 0, 0, 0, 0};
			thin_kick(calc_config, muls, length, h_bend, h_ref, ps);
			BOOST_CHECK_CLOSE(ps[x_],     x,           1e-14);
			BOOST_CHECK_CLOSE(ps[px_],    px_expected, 1e-14);

			BOOST_CHECK_SMALL(ps[y_],     1e-14);
			BOOST_CHECK_SMALL(ps[py_],    1e-14);
			BOOST_CHECK_SMALL(ps[ct_],    1e-14);
			BOOST_CHECK_SMALL(ps[delta_], 1e-14);
		}

	}
}

/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
