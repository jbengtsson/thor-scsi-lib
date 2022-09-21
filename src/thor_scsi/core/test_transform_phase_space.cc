#define BOOST_TEST_MODULE transform
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <thor_scsi/core/transform_phase_space.h>
#include <gtpsa/ss_vect.h>
#include <cmath>

namespace tsc = thor_scsi::core;

#define CHECK_LONGITUDINAL_ZERO(ps) do{	BOOST_CHECK_SMALL(ps[ct_], 1e-21); BOOST_CHECK_SMALL(ps[delta_], 1e-21);}while(0)

BOOST_AUTO_TEST_CASE(test00_init_api)
{
	auto tf = tsc::PhaseSpaceGalilean2DTransform();
	std::cout << "Galilean phase space transform " << tf << std::endl;
}

BOOST_AUTO_TEST_CASE(test10_init_zero_ps)
{
	auto tf = tsc::PhaseSpaceGalilean2DTransform();
	double unused=0e0;
	gtpsa::ss_vect<double> ps(unused);
	ps.set_zero();
	// made code crash if assigning to variable with same name!
	const gtpsa::ss_vect<double> ps_ref = ps.clone();

	tf.forward(ps);
	BOOST_CHECK_SMALL(ps[x_], 1e-21);
	BOOST_CHECK_SMALL(ps[y_], 1e-21);
	BOOST_CHECK_SMALL(ps[px_], 1e-21);
	BOOST_CHECK_SMALL(ps[py_], 1e-21);
	CHECK_LONGITUDINAL_ZERO(ps);

	tf.backward(ps);
	BOOST_CHECK_SMALL(ps[x_], 1e-21);
	BOOST_CHECK_SMALL(ps[y_], 1e-21);
	BOOST_CHECK_SMALL(ps[px_], 1e-21);
	BOOST_CHECK_SMALL(ps[py_], 1e-21);
	CHECK_LONGITUDINAL_ZERO(ps);


}

BOOST_AUTO_TEST_CASE(test20_init_shift_x)
{
	auto tf = tsc::PhaseSpaceGalilean2DTransform();
	const gtpsa::ss_vect<double> ps_ref{0,0,0,0,0,0};
	gtpsa::ss_vect<double> ps = ps_ref.clone();
	const double val = 1;

	tf.setDx(val);
	tf.forward(ps);
	BOOST_CHECK_CLOSE(ps[x_], -val, 1e-14);
	BOOST_CHECK_SMALL(ps[y_], 1e-21);
	BOOST_CHECK_SMALL(ps[px_], 1e-21);
	BOOST_CHECK_SMALL(ps[py_], 1e-21);
	CHECK_LONGITUDINAL_ZERO(ps);

	tf.backward(ps);
	BOOST_CHECK_SMALL(ps[x_], 1e-21);
	BOOST_CHECK_SMALL(ps[y_], 1e-21);
	BOOST_CHECK_SMALL(ps[px_], 1e-21);
	BOOST_CHECK_SMALL(ps[py_], 1e-21);
	CHECK_LONGITUDINAL_ZERO(ps);

}

BOOST_AUTO_TEST_CASE(test21_init_shift_y)
{
	auto tf = tsc::PhaseSpaceGalilean2DTransform();
	const gtpsa::ss_vect<double> ps_ref{0,0,0,0,0,0};
	gtpsa::ss_vect<double> ps = ps_ref.clone();
	const double val = 1;

	tf.setDy(val);
	tf.forward(ps);
	BOOST_CHECK_CLOSE(ps[y_], -val, 1e-14);
	BOOST_CHECK_SMALL(ps[x_], 1e-21);
	BOOST_CHECK_SMALL(ps[px_], 1e-21);
	BOOST_CHECK_SMALL(ps[py_], 1e-21);
	CHECK_LONGITUDINAL_ZERO(ps);

	tf.backward(ps);
	BOOST_CHECK_SMALL(ps[y_], 1e-21);
	BOOST_CHECK_SMALL(ps[x_], 1e-21);
	BOOST_CHECK_SMALL(ps[px_], 1e-21);
	BOOST_CHECK_SMALL(ps[py_], 1e-21);
	CHECK_LONGITUDINAL_ZERO(ps);

}


BOOST_AUTO_TEST_CASE(test30_realistic_shift_x)
{
	auto tf = tsc::PhaseSpaceGalilean2DTransform();

	const double x = 2e-3, y=-3e-3, px=-5e-4, py=-5e-7;
	const double val = .5;
	const gtpsa::ss_vect<double> ps_ref{x,px,y,py,0,0};
	gtpsa::ss_vect<double> ps = ps_ref.clone();

	BOOST_CHECK_CLOSE(ps[x_],         x, 1e-14);
	BOOST_CHECK_CLOSE(ps[y_],         y, 1e-14);
	BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
	BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
	CHECK_LONGITUDINAL_ZERO(ps);

	tf.setDx(val);
	tf.forward(ps);
	BOOST_CHECK_CLOSE(ps[x_],  -val + x, 1e-12);
	BOOST_CHECK_CLOSE(ps[y_],         y, 1e-14);
	BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
	BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
	CHECK_LONGITUDINAL_ZERO(ps);

	tf.backward(ps);
	BOOST_CHECK_CLOSE(ps[x_],         x, 1e-12);
	BOOST_CHECK_CLOSE(ps[y_],         y, 1e-14);
	BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
	BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
	CHECK_LONGITUDINAL_ZERO(ps);

}


BOOST_AUTO_TEST_CASE(test30_realistic_shift_y)
{
	auto tf = tsc::PhaseSpaceGalilean2DTransform();

	const double x =-2e-3, y= 3e-3, px=-5e4, py=-5e-7;
	const double val = .5;
	const gtpsa::ss_vect<double> ps_ref{x,px,y,py,0,0};
	gtpsa::ss_vect<double> ps = ps_ref.clone();

	tf.setDy(val);
	tf.forward(ps);
	BOOST_CHECK_CLOSE(ps[x_],         x, 1e-14);
	BOOST_CHECK_CLOSE(ps[y_],  -val + y, 1e-12);
	BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
	BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
	CHECK_LONGITUDINAL_ZERO(ps);

	tf.backward(ps);
	BOOST_CHECK_CLOSE(ps[x_],         x, 1e-14);
	BOOST_CHECK_CLOSE(ps[y_],         y, 1e-12);
	BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
	BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
	CHECK_LONGITUDINAL_ZERO(ps);

}

BOOST_AUTO_TEST_CASE(test40_zero_ps)
{
	auto tf = tsc::PhaseSpaceGalilean2DTransform();

	const gtpsa::ss_vect<double> ps_ref{0,0,0,0,0,0};
	{
		gtpsa::ss_vect<double> ps = ps_ref.clone();

		tf.forward(ps);
		BOOST_CHECK_SMALL(ps[x_],   1e-21);
		BOOST_CHECK_SMALL(ps[x_],   1e-21);
		BOOST_CHECK_SMALL(ps[px_],  1e-21);
		BOOST_CHECK_SMALL(ps[py_],  1e-21);
		CHECK_LONGITUDINAL_ZERO(ps);

		tf.backward(ps);
		BOOST_CHECK_SMALL(ps[x_],   1e-21);
		BOOST_CHECK_SMALL(ps[y_],   1e-21);
		BOOST_CHECK_SMALL(ps[px_],  1e-21);
		BOOST_CHECK_SMALL(ps[py_],  1e-21);
		CHECK_LONGITUDINAL_ZERO(ps);
	}
	{
		gtpsa::ss_vect<double> ps = ps_ref.clone();
		tf.setRoll(M_PI);
		tf.forward(ps);
		BOOST_CHECK_SMALL(ps[x_],   1e-21);
		BOOST_CHECK_SMALL(ps[x_],   1e-21);
		BOOST_CHECK_SMALL(ps[px_],  1e-21);
		BOOST_CHECK_SMALL(ps[py_],  1e-21);
		CHECK_LONGITUDINAL_ZERO(ps);

		tf.backward(ps);
		BOOST_CHECK_SMALL(ps[x_],   1e-21);
		BOOST_CHECK_SMALL(ps[y_],   1e-21);
		BOOST_CHECK_SMALL(ps[px_],  1e-21);
		BOOST_CHECK_SMALL(ps[py_],  1e-21);
		CHECK_LONGITUDINAL_ZERO(ps);
	}
}

BOOST_AUTO_TEST_CASE(test41_zero_roll)
{
	auto tf = tsc::PhaseSpaceGalilean2DTransform();
	const double x = 2, y=-3, px=5, py=7;
	const gtpsa::ss_vect<double> ps_ref{x,px,y,py,0,0};
	gtpsa::ss_vect<double> ps = ps_ref.clone();

	BOOST_CHECK_CLOSE(ps[x_],         x, 1e-14);
	BOOST_CHECK_CLOSE(ps[y_],         y, 1e-14);
	BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
	BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
	CHECK_LONGITUDINAL_ZERO(ps);

	tf.forward(ps);
	BOOST_CHECK_CLOSE(ps[x_],         x, 1e-12);
	BOOST_CHECK_CLOSE(ps[y_],         y, 1e-14);
	BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
	BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
	CHECK_LONGITUDINAL_ZERO(ps);

	tf.backward(ps);
	BOOST_CHECK_CLOSE(ps[x_],         x, 1e-12);
	BOOST_CHECK_CLOSE(ps[y_],         y, 1e-14);
	BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
	BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
	CHECK_LONGITUDINAL_ZERO(ps);
}


BOOST_AUTO_TEST_CASE(test50_roll_half_quater)
{
	auto tf = tsc::PhaseSpaceGalilean2DTransform();
	const double x = 2, y=-3, px=5, py=7;
	const gtpsa::ss_vect<double> ps_ref{x,px,y,py,0,0};

	{
		gtpsa::ss_vect<double> ps = ps_ref.clone();
		BOOST_CHECK_CLOSE(ps[x_],         x, 1e-14);
		BOOST_CHECK_CLOSE(ps[y_],         y, 1e-14);
		BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
		BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
		CHECK_LONGITUDINAL_ZERO(ps);

		tf.setRoll(M_PI);
		tf.forward(ps);
		BOOST_CHECK_CLOSE(ps[x_],        -x, 1e-12);
		BOOST_CHECK_CLOSE(ps[y_],        -y, 1e-12);
		BOOST_CHECK_CLOSE(ps[px_],      -px, 1e-12);
		BOOST_CHECK_CLOSE(ps[py_],      -py, 1e-12);
		CHECK_LONGITUDINAL_ZERO(ps);

		tf.backward(ps);
		BOOST_CHECK_CLOSE(ps[x_],         x, 1e-14);
		BOOST_CHECK_CLOSE(ps[y_],         y, 1e-14);
		BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
		BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
		CHECK_LONGITUDINAL_ZERO(ps);
	}

	{
		gtpsa::ss_vect<double> ps = ps_ref.clone();
		BOOST_CHECK_CLOSE(ps[x_],         x, 1e-14);
		BOOST_CHECK_CLOSE(ps[y_],         y, 1e-14);
		BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
		BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
		CHECK_LONGITUDINAL_ZERO(ps);

		tf.setRoll(M_PI/2);
		tf.forward(ps);
		BOOST_CHECK_CLOSE(ps[x_],         y, 1e-12);
		BOOST_CHECK_CLOSE(ps[y_],        -x, 1e-12);
		BOOST_CHECK_CLOSE(ps[px_],       py, 1e-12);
		BOOST_CHECK_CLOSE(ps[py_],      -px, 1e-12);
		CHECK_LONGITUDINAL_ZERO(ps);

		tf.backward(ps);
		BOOST_CHECK_CLOSE(ps[x_],         x, 1e-12);
		BOOST_CHECK_CLOSE(ps[y_],         y, 1e-12);
		BOOST_CHECK_CLOSE(ps[px_],       px, 1e-12);
		BOOST_CHECK_CLOSE(ps[py_],       py, 1e-12);
		CHECK_LONGITUDINAL_ZERO(ps);
	}
	{
		gtpsa::ss_vect<double> ps = ps_ref.clone();
		BOOST_CHECK_CLOSE(ps[x_],         x, 1e-14);
		BOOST_CHECK_CLOSE(ps[y_],         y, 1e-14);
		BOOST_CHECK_CLOSE(ps[px_],       px, 1e-14);
		BOOST_CHECK_CLOSE(ps[py_],       py, 1e-14);
		CHECK_LONGITUDINAL_ZERO(ps);

		tf.setRoll(-M_PI/2);
		tf.forward(ps);
		BOOST_CHECK_CLOSE(ps[x_],        -y, 1e-12);
		BOOST_CHECK_CLOSE(ps[y_],         x, 1e-12);
		BOOST_CHECK_CLOSE(ps[px_],      -py, 1e-12);
		BOOST_CHECK_CLOSE(ps[py_],       px, 1e-12);
		CHECK_LONGITUDINAL_ZERO(ps);

		tf.backward(ps);
		BOOST_CHECK_CLOSE(ps[x_],         x, 1e-12);
		BOOST_CHECK_CLOSE(ps[y_],         y, 1e-12);
		BOOST_CHECK_CLOSE(ps[px_],       px, 1e-12);
		BOOST_CHECK_CLOSE(ps[py_],       py, 1e-12);
		CHECK_LONGITUDINAL_ZERO(ps);
	}
}

/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
