#define BOOST_TEST_MODULE radiation_delegate
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <thor_scsi/elements/radiation_delegate.h>
#include <thor_scsi/elements/field_kick.h>
#include <thor_scsi/elements/marker.h>
#include <cmath>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

BOOST_AUTO_TEST_CASE(test10_rad_del_strm)
{
	auto rad_del = tse::RadiationDelegate();

	boost::test_tools::output_test_stream output;
	output << rad_del;
	BOOST_CHECK( !output.is_empty( false ) );

	std::cout << rad_del << std::endl;
}

BOOST_AUTO_TEST_CASE(test10_rad_del_kick_strm)
{
	auto rad_del = tse::RadiationDelegateKick();

	boost::test_tools::output_test_stream output;
	output << rad_del;
	BOOST_CHECK( !output.is_empty( false ) );

	std::cout << rad_del << std::endl;
}



BOOST_AUTO_TEST_CASE(test20_rad_del_ps0)
{
	Config marker_C;
	marker_C.set<std::string>("name", "marker");

	auto marker = tse::MarkerType(marker_C);
	auto rad_del = tse::RadiationDelegate();
	ss_vect<double> ps;

	ps.zero();

	rad_del.view(marker, ps, tsc::ObservedState::start, 0);

	BOOST_CHECK_THROW(rad_del.view(marker, ps, tsc::ObservedState::end, 0), std::domain_error);

}

BOOST_AUTO_TEST_CASE(test21_rad_del_ps0)
{
	Config marker_C;
	marker_C.set<std::string>("name", "marker");

	tsc::ConfigType calc_config;
	calc_config.emittance = true;

	auto marker = tse::MarkerType(marker_C);
	auto tmp = new tse::RadiationDelegate();
	std::shared_ptr<tse::RadiationDelegate> rad_del(tmp);
	marker.setRadiationDelegate(rad_del);

	ss_vect<double> ps;
	ps.zero();

	BOOST_CHECK_THROW(marker.pass(calc_config, ps), std::domain_error);

}

BOOST_AUTO_TEST_CASE(test30_rad_del_tps)
{
	auto rad_del = tse::RadiationDelegate();
	Config marker_C;
	marker_C.set<std::string>("name", "marker_test_tps_i");

	auto marker = tse::MarkerType(marker_C);
	ss_vect<tps> tps;

	tps.set_identity();

	rad_del.view(marker, tps, tsc::ObservedState::start, 0);
	rad_del.view(marker, tps, tsc::ObservedState::end, 0);

	std::cout<< rad_del << std::endl;
	std::cout<< tps << std::endl;
}

BOOST_AUTO_TEST_CASE(test31_rad_del_tps_identity)
{
	Config marker_C;
	marker_C.set<std::string>("name", "marker_test_tps_i_prop");

	tsc::ConfigType calc_config;
	calc_config.emittance = true;

	auto marker = tse::MarkerType(marker_C);
	auto tmp = new tse::RadiationDelegate();
	std::shared_ptr<tse::RadiationDelegate> rad_del(tmp);
	marker.setRadiationDelegate(rad_del);

	ss_vect<tps> tps;

	tps.set_identity();

	rad_del->view(marker, tps, tsc::ObservedState::start, 0);
	rad_del->view(marker, tps, tsc::ObservedState::end, 0);

	std::cout<< rad_del << std::endl;
	std::cout<< tps << std::endl;

}
