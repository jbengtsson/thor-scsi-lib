#define BOOST_TEST_MODULE combination
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <thor_scsi/elements/drift.h>
#include <gtpsa/ss_vect.h>
#include <cmath>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

/* angle indices not defined, lets compiler check if parameters are consistent. */
enum {
	xi = X_, yi = Y_, zi = Z_, pxi = 3, pyi = 4, pzi=5
};

static void
diff_position_direction(const gtpsa::ss_vect<double>& start, const gtpsa::ss_vect<double>& end,
			double *dpos, double *ddir)
{
	const gtpsa::ss_vect<double> d = end - start;

	*dpos = sqrt(sqr(d[xi]) + sqr(d[yi]) + sqr(d[zi]));
	*ddir = sqrt(sqr(d[pxi]) + sqr(d[pyi]) + sqr(d[pzi]));
}


static void
test_zero_movement(tse::DriftType& drift, tsc::ConfigType& calc_config,
		   const gtpsa::ss_vect<double>& ps_start,
		   const double pos_diff=1e-15,
		   const double dir_diff=1e-15)
{

        gtpsa::ss_vect<double> ps = ps_start.clone();
	drift.propagate(calc_config, ps);

	double dpos, ddir;
	diff_position_direction(ps_start, ps, &dpos, &ddir);

	BOOST_CHECK_CLOSE(dpos, 0.0, pos_diff);
	BOOST_CHECK_CLOSE(ddir, 0.0, dir_diff);
}

BOOST_AUTO_TEST_CASE(test01_drift_print)
{
	Config C;
	C.set<std::string>("name", "d");
	tse::DriftType drift(C);

	BOOST_CHECK_CLOSE(drift.getLength(), 0.0, 1e-15);
	gtpsa::ss_vect<double> ps = {0, 0, 0, 0, 0, 0};

	// just check show works
	std::cout << "state: ";
	ps.show(std::cout, 3);

	// just check it is printing
	std::cout << "state: " << ps << std::endl;
}

BOOST_AUTO_TEST_CASE(test02_drift_zero_length)
{

	  Config C;
	  C.set<std::string>("name", "d");
	  tse::DriftType drift(C);

	  BOOST_CHECK_CLOSE(drift.getLength(), 0.0, 1e-15);

	  tsc::ConfigType calc_config;

	  {
		  const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
		  test_zero_movement(drift, calc_config, ps_orig);
	  }

	  {
		  const gtpsa::ss_vect<double> ps_orig = {0, 0, 254, 0, 0, 0};
		  test_zero_movement(drift, calc_config, ps_orig);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {0, 0, 254, 0, 0, 1};
		  test_zero_movement(drift, calc_config, ps_orig, 1e-300, 1e-300);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {1e-3,   1e-3,  254,  1e-3,  1e-3, 1};
		  test_zero_movement(drift, calc_config, ps_orig, 1e-300, 1e-300);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {-1e-3, -1e-3,  254,  1e-3,  1e-3, 1};
		  test_zero_movement(drift, calc_config, ps_orig, 1e-300, 1e-300);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {1e-3,   1e-3,  254, -1e-3, -1e-3, 1};
		  test_zero_movement(drift, calc_config, ps_orig, 1e-300, 1e-300);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {-1e-3, -1e-3,  254, -1e-3, -1e-3, 1};
		  test_zero_movement(drift, calc_config, ps_orig, 1e-300, 1e-300);
	  }
}

static void
test_length_drift(tse::DriftType& drift, tsc::ConfigType& calc_config,
		  const gtpsa::ss_vect<double>& ps_start, const double length = 1.0,
		  const double pos_diff=1e-15, const double dir_diff=1e-15)
{
	gtpsa::ss_vect<double> ps = ps_start.clone(), ps_ref = ps_start.clone();
	drift.propagate(calc_config, ps);

	/* just propagate*/
	const double
		dx = length * ps_start[pxi],
		dy = length * ps_start[pyi],
		dz = length * ps_start[pzi];

	ps_ref[xi] = ps_start[xi] + dx;
	ps_ref[yi] = ps_start[yi] + dy;
	ps_ref[zi] = ps_start[zi] + dz;

	double dpos, ddir;
	diff_position_direction(ps_start, ps_ref, &dpos, &ddir);

	const double expected_pos_diff = sqrt(dx*dx + dy*dy + dz*dz);
	BOOST_CHECK_CLOSE(dpos, expected_pos_diff, pos_diff);
	BOOST_CHECK_CLOSE(ddir, 0.0, dir_diff);
}

static void
test_unit_length_drift(tse::DriftType& drift, tsc::ConfigType& calc_config,
		       const gtpsa::ss_vect<double>& ps_start,
		       const double pos_diff=1e-15, const double dir_diff=1e-15)
{
	test_length_drift(drift, calc_config, ps_start, 1.0, pos_diff, dir_diff);
}

BOOST_AUTO_TEST_CASE(test03_drift_one_length_test)
{

	  Config C;
	  C.set<double>("L", 1);
	  C.set<std::string>("name", "d");

	  tsc::ConfigType calc_config;

	  tse::DriftType drift(C);
	  BOOST_CHECK_CLOSE(drift.getLength(), 1, 1e-6);

	  {
		  const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
		  test_unit_length_drift(drift, calc_config, ps_orig);
	  }

	  {
		  const gtpsa::ss_vect<double> ps_orig = {0, 0, 254, 0, 0, 0};
		  test_unit_length_drift(drift, calc_config, ps_orig);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {0, 0, 254, 0, 0, 1};
		  test_unit_length_drift(drift, calc_config, ps_orig, 1e-300, 1e-300);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {1e-3,   1e-3,  254,  1e-3,  1e-3, 1};
		  test_unit_length_drift(drift, calc_config, ps_orig, 1e-300, 1e-300);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {-1e-3, -1e-3,  254,  1e-3,  1e-3, 1};
		  test_unit_length_drift(drift, calc_config, ps_orig, 1e-300, 1e-300);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {1e-3,   1e-3,  254, -1e-3, -1e-3, 1};
		  test_unit_length_drift(drift, calc_config, ps_orig, 1e-300, 1e-300);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {-1e-3, -1e-3,  254, -1e-3, -1e-3, 1};
		  test_unit_length_drift(drift, calc_config, ps_orig, 1e-300, 1e-300);
	  }
}

/**
 * typical lengthes
 * short quadrupole
 */
BOOST_AUTO_TEST_CASE(test04_drift_short_quad_length_test)
{

	  Config C;
	  /* long quadrupole */
	  const double length=0.2;
	  C.set<double>("L", length);
	  C.set<std::string>("name", "d");

	  tsc::ConfigType calc_config;

	  tse::DriftType drift(C);
	  BOOST_CHECK_CLOSE(drift.getLength(), length, 1e-15);

	  {
		  const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
		  test_length_drift(drift, calc_config, ps_orig, length);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {-1e-3, -1e-3,  254,  1e-3,  1e-3, 1};
		  test_length_drift(drift, calc_config, ps_orig, length, 1e-11, 1e-300);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {-1e-3,   1e-3,  254, -1e-3,  1e-3, 1};
		  test_length_drift(drift, calc_config, ps_orig, length, 1e-10, 1e-300);
	  }
}
/**
 * typical lengthes
 * long quadrupole
 */
BOOST_AUTO_TEST_CASE(test05_drift_long_quad_length_test)
{

	  Config C;
	  /* long quadrupole */
	  const double length=0.5;
	  C.set<double>("L", length);
	  C.set<std::string>("name", "d");

	  tsc::ConfigType calc_config;

	  tse::DriftType drift(C);
	  BOOST_CHECK_CLOSE(drift.getLength(), length, 1e-15);

	  {
		  const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, 0};
		  test_length_drift(drift, calc_config, ps_orig, length);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {-1e-3, -1e-3,  254,  1e-3,  1e-3, 1};
		  test_length_drift(drift, calc_config, ps_orig, length, 1e-19, 1e-300);
	  }
	  {
		  const gtpsa::ss_vect<double> ps_orig = {-1e-3,   1e-3,  254, -1e-3,  1e-3, 1};
		  test_length_drift(drift, calc_config, ps_orig, length, 1e-18, 1e-300);
	  }
}

/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
