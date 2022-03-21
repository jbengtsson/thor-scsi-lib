#define BOOST_TEST_MODULE quadrupole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/elements/field_kick.h>
#include <ostream>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

BOOST_AUTO_TEST_CASE(test01_kick_print)
{
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("Method", 4.0);
	tse::FieldKick kick(C);

	{
		boost::test_tools::output_test_stream output;
		output << kick;
		BOOST_CHECK( !output.is_empty( false ) );
	}
	{
		boost::test_tools::output_test_stream output;
		kick.show(output, 4);
		BOOST_CHECK( !output.is_empty( false ) );
	}
}

BOOST_AUTO_TEST_CASE(test02_kick_integral_zero_length)
{
	Config C;
	C.set<std::string>("name", "test");

	{
		C.set<double>("L", 0.0);
		tse::FieldKick kick(C);
		BOOST_CHECK_EQUAL(kick.isThick(), false);
	}
	{
		C.set<double>("L", 1.0);
		tse::FieldKick kick(C);
		BOOST_CHECK_EQUAL(kick.isThick(), true);
	}
	{
		C.set<double>("L", 0.0);
		tse::FieldKick kick(C);
		BOOST_CHECK_EQUAL(kick.isThick(), false);
	}
}

BOOST_AUTO_TEST_CASE(test03_default_values)
{
	Config C;
	const double length = 0.5;
	C.set<std::string>("name", "test");
	C.set<double>("L", length);

	auto kick = tse::FieldKick(C);

	BOOST_CHECK_CLOSE(kick.getLength(), length, 1e-12);
	BOOST_CHECK_CLOSE(kick.getInverseRigidity(), 0e0, 1e-12);

	BOOST_CHECK(kick.getNumberOfIntegrationSteps() ==  1);
	BOOST_CHECK(kick.assumingCurvedTrajectory() ==  false);

	BOOST_CHECK_CLOSE(kick.getBendingAngle(),  0e0, 1e-12);
	BOOST_CHECK_CLOSE(kick.getEntranceAngle(), 0e0, 1e-12);
	BOOST_CHECK_CLOSE(kick.getExitAngle(),     0e0, 1e-12);

}


BOOST_AUTO_TEST_CASE(test20_4Order_constants)
{
	Config C;
	const double length = 1.0;
	C.set<std::string>("name", "test");

	auto kick = tse::FieldKick(C);
	auto delegator = dynamic_cast<const tse::FieldKick::FieldKickForthOrder &>(kick.getFieldKickDelegator());

	// constants as computed by python
	const double c_1 =   0.6756035959798289;
	const double c_2 =  -0.1756035959798288;

	const double d_1 =   1.3512071919596578;
	const double d_2 =  -1.7024143839193153;


	double dL1, dL2, dkL1, dkL2;
	delegator.splitIntegrationStep(1.0, &dL1, &dL2, &dkL1, &dkL2);

	BOOST_CHECK_CLOSE(dL1,  c_1, 1e-12);
	BOOST_CHECK_CLOSE(dL2,  c_2, 1e-12);
	BOOST_CHECK_CLOSE(dkL1, d_1, 1e-12);
	BOOST_CHECK_CLOSE(dkL2, d_2, 1e-12);
}


BOOST_AUTO_TEST_CASE(test20_4Order_lengthes)
{
	Config C;
	const double length = 1.0;
	C.set<std::string>("name", "test");
	C.set<double>("L", length);

	auto kick = tse::FieldKick(C);
	auto delegator = dynamic_cast<const tse::FieldKick::FieldKickForthOrder &>(kick.getFieldKickDelegator());

}

/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
