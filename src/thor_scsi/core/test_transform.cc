#define BOOST_TEST_MODULE combination
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <thor_scsi/core/transform.h>
#include <cmath>

namespace tsc = thor_scsi::core;

BOOST_AUTO_TEST_CASE(test_init)
{

	auto tf = tsc::Euclidian2DTransform();

	BOOST_CHECK_CLOSE(tf.getDx(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getDy(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getRoll(), 0, 0);
}


BOOST_AUTO_TEST_CASE(test_roll)
{

	const double test_angle = 90.0/180.0*M_PI;

	auto tf = tsc::Euclidian2DTransform();
	tf.setRoll(test_angle);

	BOOST_CHECK_CLOSE(tf.getDx(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getDy(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getRoll(), test_angle, 1e-45);
}

BOOST_AUTO_TEST_CASE(test_roll_small_angle)
{

	const double test_angle = 1e-3 /* a milli radiant is rather more standard */;

	auto tf = tsc::Euclidian2DTransform();
	tf.setRoll(test_angle);

	BOOST_CHECK_CLOSE(tf.getDx(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getDy(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getRoll(), test_angle, 1e-45);
}


BOOST_AUTO_TEST_CASE(test_dpos)
{

	const double test_val = .3e-3 /* fraction of a millimeter offset  */;

	{
		auto tf = tsc::Euclidian2DTransform();
		tf.setDx(0);	    tf.setDy(test_val);

		BOOST_CHECK_CLOSE(tf.getDx(), 0, 1e-42);
		BOOST_CHECK_CLOSE(tf.getDy(), test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getRoll(), 0, 1e-45);

	}
	{
		auto tf = tsc::Euclidian2DTransform();
		tf.setDx(test_val); tf.setDy(0);

		BOOST_CHECK_CLOSE(tf.getDx(), test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getDy(), 0, 1e-42);
		BOOST_CHECK_CLOSE(tf.getRoll(), 0, 1e-45);

	}

}

BOOST_AUTO_TEST_CASE(test_dpos_dval)
{

	const double
		test_val = .3e-3 /* fraction of a millimeter offset  */,
		test_angle = 1e-3 /* a milli radiant is rather more standard */;

	{
		auto tf = tsc::Euclidian2DTransform();
		tf.setDx(test_val); tf.setDy(test_val); tf.setRoll(test_angle);

		BOOST_CHECK_CLOSE(tf.getDx(),   test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getDy(),   test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getRoll(), test_angle, 1e-45);

	}
	{
		auto tf = tsc::Euclidian2DTransform();
		tf.setDx(-test_val); tf.setDy(-test_val); tf.setRoll(-test_angle);

		BOOST_CHECK_CLOSE(tf.getDx(),   -test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getDy(),   -test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getRoll(), -test_angle, 1e-45);

	}
	{
		auto tf = tsc::Euclidian2DTransform();
		tf.setDx(-test_val); tf.setDy(+test_val); tf.setRoll(-test_angle);

		BOOST_CHECK_CLOSE(tf.getDx(),   -test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getDy(),    test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getRoll(), -test_angle, 1e-45);

	}

}
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
