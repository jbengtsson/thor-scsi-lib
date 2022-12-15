#define BOOST_TEST_MODULE sector_bend
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/elements/bending.h>
#include <thor_scsi/elements/radiation_delegate.h>
#include <thor_scsi/std_machine/accelerator.h>
#include <thor_scsi/std_machine/std_machine.h>
#include <tps/enums.h>

namespace ts = thor_scsi;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;


int reg_done = ts::register_elements();

const double energy = 2.5e9;

BOOST_AUTO_TEST_CASE(test01_radiation_delegate)
{
	tsc::ConfigType calc_config;
	Config C;

	const double length = 1e0, b2 = 1e0, phi = 5e0, delta = 0e0;

	C.set<std::string>("name", "test");
	C.set<double>("K", b2);
	C.set<double>("L", length);
	C.set<double>("T", phi);
	C.set<double>("N", 1);

	const std::string drift_txt(
		"d05l2t8r: Drift, L = 0.42;\n"
		"mini_cell : LINE = (d05l2t8r);\n");

	GLPSParser parse;
	Config *C2 = parse.parse_byte(drift_txt);
	auto machine = ts::Accelerator(*C2);
	machine.set_log_level(0);



	calc_config.radiation = true;

	auto bend = tse::BendingType(C);
	auto p_del = std::make_shared<tse::RadiationDelegateKick>();
	p_del->setEnergy(energy);
	bend.setRadiationDelegate(p_del);


	const double unused = 0e0;
	gtpsa::ss_vect<double> ps_ref(unused);
	ps_ref.set_zero();

	try
	{
	        gtpsa::ss_vect<double> ps = ps_ref.clone();
		bend.propagate(calc_config, ps);
		std::cout << ps;

	}catch(std::exception& e){
		abort();
	}


}
