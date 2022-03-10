#define BOOST_TEST_MODULE thin_kick
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/elements/element_helpers.h>
#include <thor_scsi/core/multipoles.h>
#include <ostream>
#include <vector>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

typedef std::array<bool, 6> flags_t; 
typedef std::array<double, 6> eps_rel_t; 

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
	const int order = 4;
	const double length = 0.0;

	tsc::ConfigType calc_config;	
	const tsc::PlanarMultipoles muls;
       
	{
		ss_vect<double> ps = ps_orig;
		tse::thin_kick(calc_config, order, muls, length, h_bend, h_ref, ps);
		test_ps_small(ps);
		
	}

	for(int i=1; i <= 4; ++i)
	{
		ss_vect<double> ps = ps_orig;
		tsc::PlanarMultipoles muls2 = muls.clone();
		muls2.setMultipole(i, tsc::cdbl(1.0, 0.0));
		
		tse::thin_kick(calc_config, order, muls2, length, h_bend, h_ref, ps);
		test_ps_small(ps);
	}
	
	{
		// Leave an extra instance here ... in case an extra test is
		// required
		ss_vect<double> ps = ps_orig;
		tsc::PlanarMultipoles muls2 = muls.clone();
		muls2.setMultipole(5, tsc::cdbl(1.0, 0.0));
		
		tse::thin_kick(calc_config, order, muls2, length, h_bend, h_ref, ps);
		test_ps_small(ps);
	}
}

BOOST_AUTO_TEST_CASE(test02_wrong_thin_kick_L0_ps_off)
{
	const double x = 1, y = 1;
	const ss_vect<double> ps_orig = {x, 0, y, 0, 0, 0};
	const double h_ref = 0e0, h_bend = 0e0;
	const int order = 4;
	const double length = 0.0;

	tsc::ConfigType calc_config;	
	const tsc::PlanarMultipoles muls;
       
	{
		ss_vect<double> ps = ps_orig;
		tse::thin_kick(calc_config, order, muls, length, h_bend, h_ref, ps);
		
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
		
		tse::thin_kick(calc_config, order, muls2, length, h_bend, h_ref, ps);


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
	const int order = 4;
	const double length = 0;

	tsc::ConfigType calc_config;	
	const tsc::PlanarMultipoles muls;
       
	{
		ss_vect<double> ps = ps_orig;
		tse::thin_kick(calc_config, order, muls, length, h_bend, h_ref, ps);
		
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
		
		tse::thin_kick(calc_config, order, muls2, length, h_bend, h_ref, ps);
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
	const int order = 4;
	const double length = 1.0;

	tsc::ConfigType calc_config;	
	const tsc::PlanarMultipoles muls;
       
	{
		ss_vect<double> ps = ps_orig;
		tse::thin_kick(calc_config, order, muls, length, h_bend, h_ref, ps);
		test_ps_small(ps);
		
	}

	for(int i=1; i <= 4; ++i)
	{
		ss_vect<double> ps = ps_orig;
		tsc::PlanarMultipoles muls2 = muls.clone();
		
		tse::thin_kick(calc_config, order, muls2, length / (double(i)), h_bend, h_ref, ps);
		test_ps_small(ps);
	}
	
}

BOOST_AUTO_TEST_CASE(test11_thin_kick_L1_no_ps_off)
{
	const double x = -28e-4, px = -464e-5, y = +37e-4, py = 31e-4, delta=29e-4, ct=23e-4;
	const ss_vect<double> ps_orig = {x, px, y, py, delta, ct};
	const double h_ref = 0e0, h_bend = 0e0;
	const int order = 4;
	const double length = 355/113;

	tsc::ConfigType calc_config;	
	const tsc::PlanarMultipoles muls;
       
	for(int i=1; i <= 4; ++i)
	{
		ss_vect<double> ps = ps_orig;
		tsc::PlanarMultipoles muls2 = muls.clone();
		
		tse::thin_kick(calc_config, order, muls2, length * i, h_bend, h_ref, ps);
		
		BOOST_CHECK_CLOSE(ps[x_],     x,     1e-14);
		BOOST_CHECK_CLOSE(ps[px_],    px,    1e-14);
		BOOST_CHECK_CLOSE(ps[y_],     y,     1e-14);
		BOOST_CHECK_CLOSE(ps[py_],    py,    1e-14);
		BOOST_CHECK_CLOSE(ps[ct_],    ct,    1e-14);
		BOOST_CHECK_CLOSE(ps[delta_], delta, 1e-14);

	}
	
}




/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
