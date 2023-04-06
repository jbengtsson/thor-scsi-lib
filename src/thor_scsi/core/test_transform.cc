#define BOOST_TEST_MODULE transform
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <thor_scsi/core/transform.h>
#include <cmath>
#include <iostream>

namespace tsc = thor_scsi::core;

BOOST_AUTO_TEST_CASE(test_init)
{

	auto tf = tsc::Galilean2DTransform();

	BOOST_CHECK_CLOSE(tf.getDx(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getDy(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getRoll(), 0, 0);
}


BOOST_AUTO_TEST_CASE(test_roll)
{

	const double test_angle = 90.0/180.0*M_PI;

	auto tf = tsc::Galilean2DTransform();
	tf.setRoll(test_angle);

	BOOST_CHECK_CLOSE(tf.getDx(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getDy(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getRoll(), test_angle, 1e-45);
}

BOOST_AUTO_TEST_CASE(test_roll_small_angle)
{

	const double test_angle = 1e-3 /* a milli radian is rather more standard */;

	auto tf = tsc::Galilean2DTransform();
	tf.setRoll(test_angle);

	BOOST_CHECK_CLOSE(tf.getDx(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getDy(), 0, 1e-42);
	BOOST_CHECK_CLOSE(tf.getRoll(), test_angle, 1e-45);
}


BOOST_AUTO_TEST_CASE(test_dpos)
{

	const double test_val = .3e-3 /* fraction of a millimeter offset  */;

	{
		auto tf = tsc::Galilean2DTransform();
		tf.setDx(0);	    tf.setDy(test_val);

		BOOST_CHECK_CLOSE(tf.getDx(), 0, 1e-42);
		BOOST_CHECK_CLOSE(tf.getDy(), test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getRoll(), 0, 1e-45);

	}
	{
		auto tf = tsc::Galilean2DTransform();
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
		test_angle = 1e-3 /* a milli radian is rather more standard */;

	{
		auto tf = tsc::Galilean2DTransform();
		tf.setDx(test_val); tf.setDy(test_val); tf.setRoll(test_angle);

		BOOST_CHECK_CLOSE(tf.getDx(),   test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getDy(),   test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getRoll(), test_angle, 1e-45);

	}
	{
		auto tf = tsc::Galilean2DTransform();
		tf.setDx(-test_val); tf.setDy(-test_val); tf.setRoll(-test_angle);

		BOOST_CHECK_CLOSE(tf.getDx(),   -test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getDy(),   -test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getRoll(), -test_angle, 1e-45);

	}
	{
		auto tf = tsc::Galilean2DTransform();
		tf.setDx(-test_val); tf.setDy(+test_val); tf.setRoll(-test_angle);

		BOOST_CHECK_CLOSE(tf.getDx(),   -test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getDy(),    test_val, 1e-42);
		BOOST_CHECK_CLOSE(tf.getRoll(), -test_angle, 1e-45);

	}

}

BOOST_AUTO_TEST_CASE(test_Prot_zero)
{

	auto tf = tsc::PRotTransform();

	BOOST_CHECK_SMALL(tf.getC0(), 1e-21);
	BOOST_CHECK_SMALL(tf.getC1(), 1e-21);
	BOOST_CHECK_SMALL(tf.getS1(), 1e-21);

	std::cout << "PRotTransform  "<< tf << std::endl;

}
BOOST_AUTO_TEST_CASE(test_Prot_set)
{
	auto tf = tsc::PRotTransform();

	const double angle = (355e0 / 113e0) * 1e-3;
	const double val0 = 1 - angle * angle;
	const double val1 = 1 - 100 * angle * angle;

	tf.setC0(val0);
	tf.setC1(val1);
	tf.setS1(angle);

	BOOST_CHECK_CLOSE(tf.getC0(), val0, 1e-21);
	BOOST_CHECK_CLOSE(tf.getC1(), val1, 1e-21);
	BOOST_CHECK_CLOSE(tf.getS1(), angle, 1e-21);

	auto tf2 = std::move(tf);

	BOOST_CHECK_CLOSE(tf2.getC0(), val0, 1e-21);
	BOOST_CHECK_CLOSE(tf2.getC1(), val1, 1e-21);
	BOOST_CHECK_CLOSE(tf2.getS1(), angle, 1e-21);
}
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
