#define BOOST_TEST_MODULE cavity
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/elements/cavity.h>
#include <thor_scsi/elements/constants.h>
#include <ostream>

namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;


static double compute_delta(const double volt,  const double frequency, const double phase,
		     const double energy, const double ct)
{
	auto c0 = tse::speed_of_light;

	const double delta = - volt / (energy)
			*sin(2e0*M_PI*frequency/c0*ct+phase);
	return delta;
}


BOOST_AUTO_TEST_CASE(test20_cavity)
{

	const double frequency=500e6, voltage=0.500e5;
	const double ct = 1e-3;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("Frequency", frequency);
	C.set<double>("Voltage", voltage);

	tsc::ConfigType calc_config;
	tse::CavityType cavity(C);

	calc_config.Cavity_on = true;

	std::cout << "energy " << calc_config.Energy << std::endl;
	const double delta = compute_delta(voltage, frequency, 0.0, calc_config.Energy, ct);

	BOOST_CHECK_CLOSE(cavity.getVoltage(),   voltage,   1e-12);
	BOOST_CHECK_CLOSE(cavity.getFrequency(), frequency, 1e-12);

	const gtpsa::ss_vect<double> ps_orig = {0, 0, 0, 0, 0, ct};
	BOOST_CHECK_CLOSE(ps_orig[ct_],   ct,   1e-12);


	{
		gtpsa::ss_vect<double> ps = ps_orig;
		cavity.propagate(calc_config, ps);

		BOOST_CHECK_CLOSE(ps[ct_],    ct, 1e-14);
		BOOST_CHECK_CLOSE(ps[delta_], delta, 1e-14);

		BOOST_CHECK_SMALL(ps[x_],     1e-14);
		BOOST_CHECK_SMALL(ps[y_],     1e-14);
		BOOST_CHECK_SMALL(ps[px_],    1e-14);
		BOOST_CHECK_SMALL(ps[py_],    1e-14);
	}

}
