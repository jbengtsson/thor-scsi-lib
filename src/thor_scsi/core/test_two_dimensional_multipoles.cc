#define BOOST_TEST_MODULE two_dimensional_multipoles
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <thor_scsi/core/multipoles.h>
#include <gtpsa/utils_tps.hpp>
#include <tps/tps_type.h>
#include <cmath>
#include <array>

/* Second implemntation of .field method available? */
// #define THOR_SCSI_PLANAR_MULTIPOLES_FIELD2
namespace tsc = thor_scsi::core;

static const double rad2deg = 180.0/M_PI;
typedef typename tsc::StandardDoubleType::complex_type cdbl;
static const cdbl I(0, 1);

auto desc = std::make_shared<gtpsa::desc>(6, 1);
const auto t_ref = gtpsa::tpsa(desc, 1);

static auto create_ref_pos(void)
{
	auto x_ref =  gtpsa::tpsa(desc, 1);
	auto y_ref =  gtpsa::tpsa(desc, 1);
	x_ref.setv(1, {1,0, 0,0, 0,0});
	y_ref.setv(1, {0,0, 1,0, 0,0});

	std::array<gtpsa::tpsa, 2> pos{x_ref, y_ref};
	return pos;
}

auto ref_pos = create_ref_pos();

BOOST_AUTO_TEST_CASE(test01_complex_angle)
{

	{
		const double angle = M_PI/16.0;
		cdbl cang = exp(angle * I);
		BOOST_CHECK_CLOSE(std::arg(cang) * rad2deg, angle * rad2deg, 1e-4);
	}

	{
		const double angle = M_PI/4.0;
		cdbl cang = exp(angle * I);
		BOOST_CHECK_CLOSE(std::arg(cang) * rad2deg, angle * rad2deg, 1e-4);

	}
	{
		const double angle = M_PI/2.0;
		cdbl cang = exp(angle * I);
		BOOST_CHECK_CLOSE(std::arg(cang) * rad2deg, angle * rad2deg, 1e-4);

	}
	{
		const double angle = -M_PI/2.0;
		cdbl cang = exp(angle * I);
		BOOST_CHECK_CLOSE(std::arg(cang) * rad2deg, angle * rad2deg, 1e-4);
	}
	{
		const double angle = -M_PI*7/8.0;
		cdbl cang = exp(angle * I);
		BOOST_CHECK_CLOSE(std::arg(cang) * rad2deg, angle * rad2deg, 1e-4);
	}
}

BOOST_AUTO_TEST_CASE(test01_init_zero)
{
    // taking a reference here
	auto h = tsc::TwoDimensionalMultipoles(0e0);
	auto hmax = tsc::max_multipole;
	// coefficient by interface
	BOOST_CHECK_CLOSE(h.getMultipole(1).real(), 0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(1).imag(), 0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(hmax).real(), 0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(hmax).imag(), 0, 1e-12);

	// access to coefficients
	auto coeffs = h.getCoeffs();
	BOOST_CHECK_CLOSE(coeffs.at(0).real(), 0, 1e-12);
	BOOST_CHECK_CLOSE(coeffs.at(0).imag(), 0, 1e-12);
	BOOST_CHECK_CLOSE(coeffs.at(hmax - 1).real(), 0, 1e-12);
	BOOST_CHECK_CLOSE(coeffs.at(hmax - 1).imag(), 0, 1e-12);

}

BOOST_AUTO_TEST_CASE(test02_show)
{
	auto h = tsc::TwoDimensionalMultipoles(0e0);
	{
		boost::test_tools::output_test_stream output;
		output << h;
		BOOST_CHECK( !output.is_empty( false ) );
	}
	{
		boost::test_tools::output_test_stream output;
		h.show(output, 10);
		BOOST_CHECK( !output.is_empty( false ) );
	}
}

BOOST_AUTO_TEST_CASE(test03_set_hmax)
{
	// taking a reference here
	int h_max = 23;
	auto h = tsc::TwoDimensionalMultipoles(0e0, h_max);
    // more output ...
	BOOST_CHECK_CLOSE(double(h.getMultipoleMaxIndex()), double(h_max), 1e-12);
}

// check that it compiles and accepts the order
BOOST_AUTO_TEST_CASE(test04_add_test)
{
	// taking a reference here
	int h_max = 4;
	auto h = tsc::TwoDimensionalMultipoles(0e0, h_max);
	{
		std::vector<std::complex<double>> offset = {1.0, 1e-3, 1e-5, 1e-4};
		const auto h2 = h.clone();
		const auto h3 = h + offset;
	}
}

// check that it compiles and accepts the order
BOOST_AUTO_TEST_CASE(test04_mul_test)
{
	// taking a reference here
	int h_max = 4;
	auto h = tsc::TwoDimensionalMultipoles(0e0, h_max);
	{
		std::vector<std::complex<double>> scale = {1.0, 1e-3, 1e-5, 1e-4};
		const auto h2 = h.clone();
		const auto h3 = h * scale;
	}
}


BOOST_AUTO_TEST_CASE(test10_set_harmonic_dipole)
{
	auto h = tsc::TwoDimensionalMultipoles(0e0);

	// pure dipole
	h.setMultipole(1, cdbl(1, 0));
	BOOST_CHECK_CLOSE(h.getMultipole(1).real(), 1, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(1).imag(), 0, 1e-12);

	BOOST_CHECK_CLOSE(h.getMultipole(2).real(), 0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(2).imag(), 0, 1e-12);

	auto coeffs = h.getCoeffs();
	BOOST_CHECK_CLOSE(coeffs.at(0).real(), 1, 1e-12);
	BOOST_CHECK_CLOSE(coeffs.at(0).imag(), 0, 1e-12);
	BOOST_CHECK_CLOSE(coeffs.at(1).real(), 0, 1e-12);
	BOOST_CHECK_CLOSE(coeffs.at(1).imag(), 0, 1e-12);

	{
		auto field =  h.field(cdbl(0,0));
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-12);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-12);
	}
	{
		auto field =  h.field(1);
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-12);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-12);
	}

	{
		auto field =  h.field(cdbl(0,1));
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-12);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-12);
	}
	{
		auto field =  h.field(cdbl(3,.2));
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-12);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-12);
	}

	// tpsa ...
	{
		auto Bx = t_ref.clone(), By = t_ref.clone();
		h.field(gtpsa::tpsa(t_ref), gtpsa::tpsa(t_ref), &Bx, &By);
		BOOST_CHECK_CLOSE(By.cst(), 1, 1e-12);
		BOOST_CHECK_CLOSE(Bx.cst(), 0, 1e-12);
	}
	{
		auto x = t_ref.clone(), y = t_ref.clone(), Bx = t_ref.clone(), By = t_ref.clone();
		x  =  1;
		h.field(x, y, &Bx, &By);
		BOOST_CHECK_CLOSE(By.cst(), 1, 1e-12);
		BOOST_CHECK_CLOSE(Bx.cst(), 0, 1e-12);
	}

	{
		auto x = t_ref.clone(), y = t_ref.clone(), Bx = t_ref.clone(), By = t_ref.clone();
		y = 1;
		h.field(x, y, &Bx, &By);
		BOOST_CHECK_CLOSE(By.cst(), 1, 1e-12);
		BOOST_CHECK_CLOSE(Bx.cst(), 0, 1e-12);
	}
	{
		auto x = t_ref.clone(), y = t_ref.clone(), Bx = t_ref.clone(), By = t_ref.clone();
		y = .2;
		x = 3;
		h.field(x, y, &Bx, &By);
		BOOST_CHECK_CLOSE(By.cst(), 1, 1e-12);
		BOOST_CHECK_SMALL(Bx.cst(),    1e-12);
	}

	h.applyRollAngle(M_PI/2);
	BOOST_CHECK_CLOSE(h.getMultipole(1).imag(), 1, 1e-15);
	BOOST_CHECK_SMALL(h.getMultipole(2).real(), 1e-12);
	BOOST_CHECK_SMALL(h.getMultipole(2).imag(), 1e-12);
	BOOST_CHECK_SMALL(double(h.getMultipole(1).real()),  1e-13);

	// Notice By + i Bx
	auto field =  h.field(0);
	// should be better than that
	BOOST_CHECK_SMALL(double(field.real()), 1e-10);
	// tests should it could achieve this value
	BOOST_WARN_SMALL(double(field.real()), 1e-12);
	BOOST_CHECK_CLOSE(field.imag(),  1, 1e-14);

}


BOOST_AUTO_TEST_CASE(test20_set_harmonic_quadrupole)
{

	auto h = tsc::TwoDimensionalMultipoles(0e0);

	// pure quadrupole
	h.setMultipole(2, cdbl(1, 0));
	BOOST_CHECK_CLOSE(h.getMultipole(2).real(), 1, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(2).imag(), 0, 1e-12);

	{
		auto field =  h.field(cdbl(0,0));
		BOOST_CHECK_SMALL(double(field.real()), 1e-12);
		BOOST_CHECK_SMALL(double(field.imag()), 1e-12);
	}
	{
		auto field =  h.field(1);
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-12);
		BOOST_CHECK_SMALL(double(field.imag()), 1e-12);
	}

	{
		auto field =  h.field(cdbl(0,1));
		BOOST_CHECK_SMALL(double(field.real()), 1e-12);
		BOOST_CHECK_CLOSE(field.imag(), 1, 1e-12);
	}
	{
		auto field =  h.field(cdbl(3,.2));
		BOOST_CHECK_CLOSE(double(field.real()), 3, 1e-12);
		BOOST_CHECK_CLOSE(double(field.imag()), .2, 1e-12);
	}

	{
		auto x = ref_pos[0].clone(), y = ref_pos[1].clone(), Bx = t_ref.clone(), By = t_ref.clone();
		std::cout << "x\n " << x << "\n"
			  << "y\n " << y << "\n";
		h.field(x, y, &Bx, &By);
		BOOST_CHECK_SMALL(double(By.cst()), 1e-12);
		BOOST_CHECK_SMALL(double(Bx.cst()), 1e-12);
		std::cout << "Bx\n " << Bx << "\n"
			  << "By\n " << By << "\n";

	}
	{

		auto x = tps(0, 0+1);
		auto y = tps(0, 2+1);
		tps Bx, By;

		std::cout << "x\n " << x << "\n"
			  << "y\n " << y << "\n";
		h.field(x, y, &Bx, &By);
		std::cout << "Bx\n " << Bx << "\n"
			  << "By\n " << By << "\n";
	}
	{
		auto x = t_ref.clone(), y = t_ref.clone(), Bx = t_ref.clone(), By = t_ref.clone();
		x = 1;
		h.field(x, y, &Bx, &By);
		BOOST_CHECK_CLOSE(double(By.cst()), 1., 1e-12);
		BOOST_CHECK_SMALL(double(Bx.cst()), 1e-12);
	}

	{
		auto x = t_ref.clone(), y = t_ref.clone(), Bx = t_ref.clone(), By = t_ref.clone();
		y = 1;
		h.field(x, y, &Bx, &By);
		BOOST_CHECK_SMALL(double(By.cst()), 1e-12);
		BOOST_CHECK_CLOSE(double(Bx.cst()), 1., 1e-12);
	}
	{
		auto field =  h.field(cdbl(3,.2));
		auto x = t_ref.clone(), y = t_ref.clone(), Bx = t_ref.clone(), By = t_ref.clone();
		x = 3;
		y = .2;
		h.field(x, y, &Bx, &By);
		BOOST_CHECK_CLOSE(double(By.cst()), 3, 1e-12);
		BOOST_CHECK_CLOSE(double(Bx.cst()), .2, 1e-12);
	}

	const double angle = M_PI/2;
	h.applyRollAngle(angle);
	BOOST_CHECK_CLOSE(h.getMultipole(2).real(), -1, 1e-15);
	// to help boost on which accuracy is meant
	BOOST_CHECK_CLOSE(h.getMultipole(2).imag() + 1, 1, 1e-12);
}


BOOST_AUTO_TEST_CASE(test30_set_scale_sextupole)
{
	// pure quadrupole
	auto h = tsc::TwoDimensionalMultipoles(0e0);
	h.setMultipole(3, cdbl(.1, -.3));
	BOOST_CHECK_CLOSE(h.getMultipole(3).real(),  0.1, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(3).imag(), -0.3, 1e-12);

	h *= 3;
	BOOST_CHECK_CLOSE(h.getMultipole(3).real(),  0.3, 1e-13);
	BOOST_CHECK_CLOSE(h.getMultipole(3).imag(), -0.9, 1e-13);

	auto h2 = h * (1/3.0);
	BOOST_CHECK_CLOSE(h2.getMultipole(3).real(),  0.1, 1e-13);
	BOOST_CHECK_CLOSE(h2.getMultipole(3).imag(), -0.3, 1e-13);
}


BOOST_AUTO_TEST_CASE(test40_rotate_quadrupole)
{

	auto h = tsc::TwoDimensionalMultipoles(0e0);

	// pure quadrupole
	h.setMultipole(2, cdbl(1, 0));
	BOOST_CHECK_CLOSE(h.getMultipole(2).real(), 1, 1e-12);
	BOOST_CHECK_SMALL(h.getMultipole(2).imag(),    1e-12);


	h.applyRollAngle(M_PI/2);
	BOOST_CHECK_CLOSE(h.getMultipole(2).real(), -1, 1e-15);
	// to help boost on which accuracy is meant
	BOOST_CHECK_CLOSE(h.getMultipole(2).imag() + 1, 1, 1e-12);

	h.applyRollAngle(M_PI/4);
	// to help boost on which accuracy is meant
	BOOST_CHECK_CLOSE(h.getMultipole(2).real() + 1, 1, 1e-13);
	// todo: check with paper that it is correct ..
	BOOST_CHECK_CLOSE(h.getMultipole(2).imag(), -1, 1e-12);

}

// check that rotation of magnet and rotated position matches
BOOST_AUTO_TEST_CASE(test41_rotate_quadrupole_consistency)
{

	return;

	auto h = tsc::TwoDimensionalMultipoles(0e0);

	const int n = 2;
	// pure quadrupole
	h.setMultipole(n, cdbl(1, 0));

	BOOST_CHECK_CLOSE(h.getMultipole(n).real(), 1, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(n).imag(), 0, 1e-12);


	const double angle = M_PI/2.0, phase = angle * n;
	const cdbl cang = exp(angle * I);

	const cdbl
		pos0(  1,   0),
		pos1(  0,   1),
		pos2(  1,  -1),
		pos3(-0.2,  3);

	// Turn to position that will be found at the positions in question
	const cdbl
		field0 = h.field(pos0 - cang),
		field1 = h.field(pos1 - cang),
		field2 = h.field(pos2 - cang),
		field3 = h.field(pos3 - cang);

	h.applyRollAngle(angle);

	const cdbl
		check0 = h.field(pos0),
		check1 = h.field(pos1),
		check2 = h.field(pos2),
		check3 = h.field(pos3);

	BOOST_CHECK_SMALL(double(field0.real()), 1e-10);
	BOOST_CHECK_SMALL(double(check0.real()), 1e-10);
	// Should reach this perfomance
	BOOST_WARN_SMALL(double(field0.real()), 1e-14);
	BOOST_WARN_SMALL(double(check0.real()), 1e-14);
	BOOST_CHECK_CLOSE(field0.imag(), check0.imag(), 1e-16);

	// Assuming that the angle is specified to 90 degrees
	BOOST_CHECK_CLOSE(field1.real(), 1, 1e-16);
	BOOST_CHECK_CLOSE(field1.real(), check1.real(), 1e-16);
	BOOST_CHECK_SMALL((double)field1.imag(), 1e-16);
	BOOST_CHECK_SMALL((double)check1.imag(), 1e-10);
	BOOST_WARN_SMALL(double(field1.imag()), 1e-14);
	BOOST_WARN_SMALL(double(check1.imag()), 1e-14);

	BOOST_CHECK_CLOSE(field2.real(), check2.real(), 1e-10);
	BOOST_CHECK_CLOSE(field2.imag(), check2.imag(), 1e-10);
	BOOST_WARN_CLOSE(field2.real(), check2.real(), 1e-14);
	BOOST_WARN_CLOSE(field2.imag(), check2.imag(), 1e-14);

	BOOST_CHECK_CLOSE(field3.real(), check3.real(), 1e-10);
	BOOST_CHECK_CLOSE(field3.imag(), check3.imag(), 1e-10);

	// Close to achievable performance
	BOOST_WARN_CLOSE(field3.real(), check3.real(), 1e-14);
	BOOST_WARN_CLOSE(field3.imag(), check3.imag(), 1e-14);

}


BOOST_AUTO_TEST_CASE(test50_octupole_field)
{
	auto h = tsc::TwoDimensionalMultipoles(0e0);

	const int n = 4;

	// pure octupole
	h.setMultipole(n, cdbl(1, 0));

	BOOST_CHECK_CLOSE(h.getMultipole(n).real(), 1, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(n).imag(), 0, 1e-12);


	// cubic  interpolation
	{
		auto field =  h.field(0);
		BOOST_CHECK_CLOSE(double(field.real()), 0, 1e-12);
		BOOST_CHECK_CLOSE(double(field.imag()), 0, 1e-12);
	}
	{
		double val = -.5;
		auto field =  h.field(-0.5);
		BOOST_CHECK_CLOSE(double(field.real()), val * val * val, 1e-12);
		BOOST_CHECK_CLOSE(double(field.imag()), 0, 1e-12);
	}
	{
		auto field =  h.field(1);
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-12);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-12);
	}
	{
		double val = -.5;
		auto field =  h.field(-0.5 * I);
		BOOST_CHECK_CLOSE(field.real(), 0, 1e-12);
		BOOST_CHECK_CLOSE(field.imag(), -val * val * val, 1e-12);
	}
	{
		auto field =  h.field(I);
		BOOST_CHECK_CLOSE(field.real(), 0, 1e-12);
		// 90 degrees ...
		BOOST_CHECK_CLOSE(field.imag(), -1, 1e-12);
	}


}

BOOST_AUTO_TEST_CASE(test62_rotate_dipole)
{
	auto h = tsc::TwoDimensionalMultipoles(0e0);
	const int n = 1;
	// pure dipole
	h.setMultipole(n, cdbl(1, 0));

	const double angle = M_PI/4, phase = angle * n;
	const cdbl
		pos0(1, 0),
		pos1(1, 1);

	const cdbl
		field0 = h.field(pos0),
		field1 = h.field(pos1);

	h.applyRollAngle(angle);
	BOOST_CHECK_CLOSE(h.getMultipole(n).real(), cos(phase), 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(n).imag(), sin(phase), 1e-12);
	BOOST_CHECK_CLOSE(std::arg(h.getMultipole(n)), phase, 1e-12);
	BOOST_CHECK_CLOSE(std::arg(h.getMultipole(n))/angle, n, 1e-12);

	{
		// check that the phase advances as expected
		const cdbl
			check0 = h.field(pos0),
			check1 = h.field(pos1);

		BOOST_CHECK_CLOSE((std::arg(field0) + phase) * rad2deg, std::arg(check0) * rad2deg, 1e-10);
		BOOST_CHECK_CLOSE((std::arg(field1) + phase) * rad2deg, std::arg(check1) * rad2deg, 1e-10);
	}
	{
		// Same position in interial frame ... changed coefficients
		const cdbl cang = exp(-angle * I);
		const cdbl
			rpos0 = pos0 * cang,
			rpos1 = pos1 * cang;

		// Check that the positions are rotated backwards
		// when the coordinate system is rotated these will be the cofficients of the
		// original position
		BOOST_CHECK_CLOSE((std::arg(rpos0) + angle) * rad2deg, std::arg(pos0) * rad2deg, 1e-13);
		BOOST_CHECK_CLOSE((std::arg(rpos1) + angle) * rad2deg, std::arg(pos1) * rad2deg, 1e-13);

		const cdbl
			check0 = h.field(rpos0),
			check1 = h.field(rpos1);

		// same field in the sense as seens as picture
		// but its coordinates are turned ...
		BOOST_CHECK_CLOSE((std::arg(check0) - angle) * rad2deg, std::arg(field0) * rad2deg, 1e-10);
		BOOST_CHECK_CLOSE((std::arg(check1) - angle) * rad2deg, std::arg(field1) * rad2deg, 1e-10);

	}
}

BOOST_AUTO_TEST_CASE(test63_rotate_octupole)
{
	auto h = tsc::TwoDimensionalMultipoles(0e0);

	const int n = 4;
	// pure quadrupole
	h.setMultipole(n, cdbl(1, 0));

	const double angle = M_PI/45, phase = angle * n;

	const cdbl
		pos0(0, 0),
		pos1(1, 0),
		pos2 = I,
		pos3(-0.5, 0),
		pos4 = -0.5 * I;


	const cdbl
		field0 = h.field(pos0),
		field1 = h.field(pos1),
		field2 = h.field(pos2),
		field3 = h.field(pos3),
		field4 = h.field(pos4);


	h.applyRollAngle(angle);
	BOOST_CHECK_CLOSE(h.getMultipole(4).real(), cos(phase), 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(4).imag(), sin(phase), 1e-12);
	BOOST_CHECK_CLOSE(std::arg(h.getMultipole(4)), phase, 1e-12);
	BOOST_CHECK_CLOSE(std::arg(h.getMultipole(4))/angle, n, 1e-12);
	// BOOST_CHECK_SMALL(double(h.getMultipole(4).imag()),  1e-14);

	BOOST_CHECK_CLOSE(h.getMultipole(5).real(),  0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(5).imag(),  0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(3).real(),  0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(3).imag(),  0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(2).real(),  0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(2).imag(),  0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(1).real(),  0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(1).imag(),  0, 1e-12);

	{
		// check that the phase advances as expected
		const cdbl
			check0 = h.field(pos0),
			check1 = h.field(pos1),
			check2 = h.field(pos2),
			check3 = h.field(pos3),
			check4 = h.field(pos4);


		// No length no angle
		// BOOST_CHECK_CLOSE(std::arg(field0) * rad2deg, std::arg(check0) * rad2deg, 1e-12);
		BOOST_CHECK_CLOSE((std::arg(field1) + phase) * rad2deg, std::arg(check1) * rad2deg, 1e-10);
		BOOST_CHECK_CLOSE((std::arg(field2) + phase) * rad2deg, std::arg(check2) * rad2deg, 1e-10);
		BOOST_CHECK_CLOSE(-360 + (std::arg(field3) + phase) * rad2deg, std::arg(check3) * rad2deg, 1e-10);
		BOOST_CHECK_CLOSE((std::arg(field4) + phase) * rad2deg, std::arg(check4) * rad2deg, 1e-140);
	}

	{
		const cdbl cang = exp(-angle * I);
		const cdbl
			rpos0 = pos0 * cang,
			rpos1 = pos1 * cang,
			rpos2 = pos2 * cang,
			rpos3 = pos3 * cang,
			rpos4 = pos4 * cang;

		// Check that the positions are rotated backwards
		// when the coordinate system is rotated these will be the cofficients of the
		// original position
		BOOST_CHECK_CLOSE((std::arg(rpos1) + angle) * rad2deg, std::arg(pos1) * rad2deg, 1e-13);
		BOOST_CHECK_CLOSE((std::arg(rpos2) + angle) * rad2deg, std::arg(pos2) * rad2deg, 1e-13);
		BOOST_CHECK_CLOSE((std::arg(rpos3) + angle) * rad2deg, std::arg(pos3) * rad2deg, 1e-13);
		BOOST_CHECK_CLOSE((std::arg(rpos4) + angle) * rad2deg, std::arg(pos4) * rad2deg, 1e-13);

		// check that the field has the same angle as before if probed at the same
		// posititon ins pace is used (thus different coordinate coefficients)
		const cdbl
			check0 = h.field(rpos0),
			check1 = h.field(rpos1),
			check2 = h.field(rpos2),
			check3 = h.field(rpos3),
			check4 = h.field(rpos4);

		// No length no angle
		// BOOST_CHECK_CLOSE(std::arg(field0) * rad2deg, std::arg(check0) * rad2deg, 1e-12);
		// angle of field is now rotated
		BOOST_CHECK_SMALL(       (std::arg(check1) - angle) * rad2deg, 1e-10);
		BOOST_CHECK_SMALL(       (std::arg(field1)        ) * rad2deg, 1e-10);
		BOOST_CHECK_CLOSE(       (std::arg(check2) - angle) * rad2deg, std::arg(field2) * rad2deg, 1e-10);
		BOOST_CHECK_CLOSE(360 +  (std::arg(check3) - angle) * rad2deg, std::arg(field3) * rad2deg, 1e-10);
		BOOST_CHECK_CLOSE(       (std::arg(check4) - angle) * rad2deg, std::arg(field4) * rad2deg, 1e-10);

		// Magnitude should be still the same
		BOOST_CHECK_CLOSE(std::norm(check0), std::norm(field0), 1e-10);
		BOOST_CHECK_CLOSE(std::norm(check1), std::norm(field1), 1e-10);
		BOOST_CHECK_CLOSE(std::norm(check2), std::norm(field2), 1e-10);
		BOOST_CHECK_CLOSE(std::norm(check3), std::norm(field3), 1e-10);
		BOOST_CHECK_CLOSE(std::norm(check4), std::norm(field4), 1e-10);

	}
}


BOOST_AUTO_TEST_CASE(test60_rotate_octupole_small_angle)
{
	auto h = tsc::TwoDimensionalMultipoles(0e0);


	const int n = 4;
	// pure octupole
	h.setMultipole(n, cdbl(1, 0));

	const double angle = 1e-3, phase = angle * n;
	const cdbl cang = exp(angle * I);


	const cdbl
		pos0(0,0),
		pos1 = -0.5 * I,
		pos2(1,0);

	const cdbl
		field0 = h.field(pos0),
		field1 = h.field(pos1),
		field2 = h.field(pos2);

	h.applyRollAngle(angle);

	/* small angle: first check with approximate change */
	BOOST_CHECK_CLOSE(h.getMultipole(n).real(), 1, 1e-3);
	// 1 for calculating relative size
	BOOST_CHECK_CLOSE(h.getMultipole(n).imag()+1, 1 + phase, 1e-3);

	// calculate phase as given by equation
	BOOST_CHECK_CLOSE(h.getMultipole(n).real(), 1 * cos(phase), 1e-12);
	// 1 for calculating relative size
	BOOST_CHECK_CLOSE(h.getMultipole(n).imag(), sin(phase), 1e-13);

	{
		const cdbl
			check0 = h.field(pos0),
			check1 = h.field(pos1),
			check2 = h.field(pos2);

		// No directon no angle
		//BOOST_CHECK_CLOSE((std::arg(field0) + phase) * rad2deg, std::arg(check0) * rad2deg, 1e-10);
		BOOST_CHECK_CLOSE((std::arg(field1) + phase) * rad2deg, std::arg(check1) * rad2deg, 1e-10);
		BOOST_CHECK_CLOSE((std::arg(field2) + phase) * rad2deg, std::arg(check2) * rad2deg, 1e-10);

		// but strength should be the same
		BOOST_CHECK_CLOSE(std::norm(check0), std::norm(field0), 1e-10);
		BOOST_CHECK_CLOSE(std::norm(check1), std::norm(field1), 1e-10);
		BOOST_CHECK_CLOSE(std::norm(check2), std::norm(field2), 1e-10);
	}

	{
		const cdbl cang = exp(-angle * I);
		const cdbl
			rpos0 = pos0 * cang,
			rpos1 = pos1 * cang,
			rpos2 = pos2 * cang;

		const cdbl
			check0 = h.field(rpos0),
			check1 = h.field(rpos1),
			check2 = h.field(rpos2);

		// Field direction close to zero
		BOOST_CHECK_CLOSE(       (std::arg(check1) - angle) * rad2deg, std::arg(field1) * rad2deg, 1e-10);
		BOOST_CHECK_SMALL((std::arg(check2) - angle) * rad2deg, 1e-10);
		BOOST_CHECK_SMALL( std::arg(field2) * rad2deg, 1e-10);
	}
}

BOOST_AUTO_TEST_CASE(test70_translate_quadrupole_zero)
{

	auto h = tsc::TwoDimensionalMultipoles(0e0);

	// pure quadrupole .. assuming harmonics have been checked
	h.setMultipole(2, cdbl(1, 0));


	// checking interface
	h.applyTranslation(0);
	h.applyTranslation(0, 0);

	cdbl dz(0,0);
	h.applyTranslation(dz);


	BOOST_CHECK_CLOSE(h.getMultipole(2).real(), 1, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(2).imag(), 0, 1e-12);

	// check that no dipole component has been added
	BOOST_CHECK_SMALL(double(h.getMultipole(1).real()), 1e-13);
	BOOST_CHECK_SMALL(double(h.getMultipole(1).imag()), 1e-13);

}


BOOST_AUTO_TEST_CASE(test80_translate_quadrupole_one)
{

	tsc::TwoDimensionalMultipolesKnobbed h = tsc::TwoDimensionalMultipoles(0e0);

	// pure quadrupole .. assuming harmonics have been checked
	h.setMultipole(2, cdbl(1, 0));

	const cdbl dz(1, 0);
	const cdbl pos0(0, 0), pos1(0, .5), pos2(1., 0);

	// probe field where coordinates will be
	const cdbl
		field0 = h.field(pos0),
		field1 = h.field(pos1),
		field2 = h.field(pos2);

	BOOST_CHECK_SMALL(field0.real(), 1e-12);
	BOOST_CHECK_SMALL(field0.imag(), 1e-12);
	BOOST_CHECK_CLOSE(field1.imag(), 0.5, 1e-12);
	BOOST_CHECK_CLOSE(field2.real(),   1, 1e-12);

	// checking interface
	h.applyTranslation(dz);

	BOOST_CHECK_CLOSE(h.getMultipole(2).real(), 1, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(2).imag(), 0, 1e-12);

	// Dipole harmonic
	BOOST_CHECK_CLOSE(h.getMultipole(1).real(), 1, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(1).imag(), 0, 1e-12);

	// no sextupole
	BOOST_CHECK_CLOSE(h.getMultipole(3).real(), 0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(3).imag(), 0, 1e-12);

	// no octupole
	BOOST_CHECK_CLOSE(h.getMultipole(4).real(), 0, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(4).imag(), 0, 1e-12);

	// probe field where coordinates were
	const cdbl
		check0 = h.field(pos0 - dz),
		check1 = h.field(pos1 - dz),
		check2 = h.field(pos2 - dz);

	BOOST_CHECK_CLOSE(field0.real(), check0.real(), 1e-12);
	BOOST_CHECK_CLOSE(field0.imag(), check0.imag(), 1e-12);

	BOOST_CHECK_CLOSE(field1.real(), check1.real(), 1e-12);
	BOOST_CHECK_CLOSE(field1.imag(), check1.imag(), 1e-12);

	BOOST_CHECK_CLOSE(field2.real(), check2.real(), 1e-12);
	BOOST_CHECK_CLOSE(field2.imag(), check2.imag(), 1e-12);
}


// A typical quadrupole for a MBA machine fairly misplaced
BOOST_AUTO_TEST_CASE(test90_translate_quadrupole_test_field)
{
	tsc::TwoDimensionalMultipolesKnobbed h = tsc::TwoDimensionalMultipoles(0e0);

	// pure quadrupole .. assuming harmonics have been checked
	h.setMultipole(2, cdbl(1, 0) * 95e0);
	const cdbl pos0(0, 0), pos1(0, -5e-3), pos2(-0.3e-3, -1e-3);
	const cdbl dz(-.2e-3, 0);

	// probe field where coordinates will be
	const cdbl
		field0 = h.field(pos0),
		field1 = h.field(pos1),
		field2 = h.field(pos2);

	// checking interface
	h.applyTranslation(dz);


	const cdbl
		check0 = h.field(pos0 - dz),
		check1 = h.field(pos1 - dz),
		check2 = h.field(pos2 - dz);

	BOOST_CHECK_CLOSE(field0.real(), check0.real(), 1e-12);
	BOOST_CHECK_CLOSE(field0.imag(), check0.imag(), 1e-12);

	BOOST_CHECK_CLOSE(field1.real(), check1.real(), 1e-12);
	BOOST_CHECK_CLOSE(field1.imag(), check1.imag(), 1e-12);

	BOOST_CHECK_CLOSE(field2.real(), check2.real(), 1e-12);
	BOOST_CHECK_CLOSE(field2.imag(), check2.imag(), 1e-12);

}

BOOST_AUTO_TEST_CASE(test91_translate_sextupole)
{
	const double sextupole = 100;
	const cdbl dz(1, 0);
	tsc::TwoDimensionalMultipolesKnobbed h = tsc::TwoDimensionalMultipoles(0e0);
	h.setMultipole(3, cdbl(1, 0) * sextupole);
	const cdbl pos0(0, 0), pos1(0, .5), pos2(1., 0);

	BOOST_CHECK_CLOSE(h.getMultipole(3).real(), sextupole, 1e-12);
	const cdbl
		field0 = h.field(pos0),
		field1 = h.field(pos1),
		field2 = h.field(pos2);

	h.applyTranslation(dz);

	// sextupole ... same as sextupole
	BOOST_CHECK_CLOSE(h.getMultipole(3).real(), sextupole, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(3).imag(), 0, 1e-12);

	// quadrupole ... half the sextupole
	BOOST_CHECK_CLOSE(h.getMultipole(2).real(), 2 * sextupole, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(2).imag(), 0, 1e-12);

	// dipole
	BOOST_CHECK_CLOSE(h.getMultipole(1).real(), sextupole, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(1).imag(), 0, 1e-12);

	// probe field where coordinates were
	const cdbl
		check0 = h.field(pos0 - dz),
		check1 = h.field(pos1 - dz),
		check2 = h.field(pos2 - dz);

	BOOST_CHECK_CLOSE(field0.real(), check0.real(), 1e-12);
	BOOST_CHECK_CLOSE(field0.imag(), check0.imag(), 1e-12);

	BOOST_CHECK_CLOSE(field1.real(), check1.real(), 1e-12);
	BOOST_CHECK_CLOSE(field1.imag(), check1.imag(), 1e-12);

	BOOST_CHECK_CLOSE(field2.real(), check2.real(), 1e-12);
	BOOST_CHECK_CLOSE(field2.imag(), check2.imag(), 1e-12);
}

BOOST_AUTO_TEST_CASE(test92_translate_octupole)
{
	const double octupole = 1000;

	const cdbl dz(1, 0);
	tsc::TwoDimensionalMultipolesKnobbed h = tsc::TwoDimensionalMultipoles(0e0);
	h.setMultipole(4, cdbl(1, 0) * octupole);
	const cdbl pos0(0, 0), pos1(0, .5), pos2(1., 0);

	BOOST_CHECK_CLOSE(h.getMultipole(4).real(), octupole, 1e-12);


	const cdbl
		field0 = h.field(pos0),
		field1 = h.field(pos1),
		field2 = h.field(pos2);

	h.applyTranslation(dz);

	BOOST_CHECK_CLOSE(h.getMultipole(4).real(), octupole, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(4).imag(), 0, 1e-12);

	// sextupole same as octupole
	BOOST_CHECK_CLOSE(h.getMultipole(3).real(), 3 * octupole , 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(3).imag(), 0, 1e-12);

	// quadrupole ... a third of
	BOOST_CHECK_CLOSE(h.getMultipole(2).real(), 3 * octupole , 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(2).imag(), 0, 1e-12);

	// dipole
	BOOST_CHECK_CLOSE(h.getMultipole(1).real(), octupole, 1e-12);
	BOOST_CHECK_CLOSE(h.getMultipole(1).imag(), 0, 1e-12);

	const cdbl
		check0 = h.field(pos0 - dz),
		check1 = h.field(pos1 - dz),
		check2 = h.field(pos2 - dz);

	BOOST_CHECK_CLOSE(field0.real(), check0.real(), 1e-12);
	BOOST_CHECK_CLOSE(field0.imag(), check0.imag(), 1e-12);

	BOOST_CHECK_CLOSE(field1.real(), check1.real(), 1e-12);
	BOOST_CHECK_CLOSE(field1.imag(), check1.imag(), 1e-12);

	BOOST_CHECK_CLOSE(field2.real(), check2.real(), 1e-12);
	BOOST_CHECK_CLOSE(field2.imag(), check2.imag(), 1e-12);

}


BOOST_AUTO_TEST_CASE(test100_translate_quadrupole_small)
{

	tsc::TwoDimensionalMultipolesKnobbed h = tsc::TwoDimensionalMultipoles(0e0);
	const double unit = 1e-4, rref=40e-3;
	//  realistic example
	const cdbl dz(1e-4/rref, .3e-3/rref);

	// pure quadrupole .. assuming harmonics have been checked
	h.setMultipole(1, cdbl(1, 0));
	{
		const unsigned int n = 6;
		const cdbl	dodecapole = cdbl(.2, 3) * unit;
		if(n <= tsc::max_multipole){
			h.setMultipole(n, dodecapole);
		}
	}
	{
		const unsigned int n = 10;
		const cdbl icosapole = cdbl(.5, -.7) * unit;
		if(n <= tsc::max_multipole){
			h.setMultipole(n, icosapole);
		}
	}


	const cdbl pos0(0, 0e0/rref) , pos1(0, 17e-3/rref), pos2(23e-3/rref, 0);
	const cdbl
		field0 = h.field(pos0),
		field1 = h.field(pos1),
		field2 = h.field(pos2);


	h.applyTranslation(dz);

	const cdbl
		check0 = h.field(pos0 - dz),
		check1 = h.field(pos1 - dz),
		check2 = h.field(pos2 - dz);
#ifdef THOR_SCSI_PLANAR_MULTIPOLES_FIELD2
	const cdbl
		check0_2 = h.field2(pos0 - dz),
		check1_2 = h.field2(pos1 - dz),
		check2_2 = h.field2(pos2 - dz);
#endif // THOR_SCSI_PLANAR_MULTIPOLES_FIELD2


	// limits set for quadmath .. in percent
	BOOST_CHECK_CLOSE(field0.real(), check0.real(), 3e-4);
	// limits set for?  .. in percent
	BOOST_WARN_CLOSE(field0.real(), check0.real(), 1.6e-4);
	BOOST_WARN_CLOSE(field0.real(), check0.real(), 1e-15);
	BOOST_CHECK_SMALL(field0.imag(),  1e-16);
	BOOST_CHECK_SMALL(check0.imag(),  1.2e-6);
	BOOST_WARN_SMALL(check0.imag(),  1e-15);

	BOOST_CHECK_CLOSE(field1.real(), check1.real(), 1.1e-4);
	BOOST_CHECK_SMALL(field1.imag() - check1.imag(), 5e-4);
	BOOST_WARN_SMALL(field1.imag()  - check1.imag(),  1e-15);

	BOOST_CHECK_CLOSE(field2.real(), check2.real(), 2e-3);
	BOOST_WARN_CLOSE(field2.real(), check2.real(), 1e-15);
	BOOST_CHECK_SMALL(field2.imag() - check2.imag(), 2e-3);
	BOOST_WARN_SMALL(field2.imag() - check2.imag(), 1e-15);

}

// tests for supporting engineering tolerances
BOOST_AUTO_TEST_CASE(test200_add_multipoles_inplace)
{

	tsc::TwoDimensionalMultipoles h(0e0), h2(0e0), h3(0e0);

	h.setMultipole(1, cdbl(1, 0));
	h2.setMultipole(2, cdbl(2, 0));
	{
		// cross check harmonics properly set
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);

		BOOST_CHECK_CLOSE(c1.real(), 1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_SMALL(c2.real(), 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}

	{
		// cross check harmonics properly set
		const cdbl c1 = h2.getMultipole(1);
		const cdbl c2 = h2.getMultipole(2);
		const cdbl c3 = h2.getMultipole(3);

		BOOST_CHECK_SMALL(c1.real(), 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), 2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}

	h3 += h2;
	h3 += h;
	{
		// cross check harmonics properly set and maintained after math
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);

		BOOST_CHECK_CLOSE(c1.real(), 1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_SMALL(c2.real(), 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}

	{
		// cross check harmonics properly set
		const cdbl c1 = h2.getMultipole(1);
		const cdbl c2 = h2.getMultipole(2);
		const cdbl c3 = h2.getMultipole(3);

		BOOST_CHECK_SMALL(c1.real(), 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), 2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}
	{
		// properly assigned
		// cross check harmonics properly set
		const cdbl c1 = h3.getMultipole(1);
		const cdbl c2 = h3.getMultipole(2);
		const cdbl c3 = h3.getMultipole(3);

		BOOST_CHECK_CLOSE(c1.real(), 1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), 2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}


}

BOOST_AUTO_TEST_CASE(test201_add_multipoles)
{

	tsc::TwoDimensionalMultipoles h(0e0), h2(0e0), h3(0e0);

	h.setMultipole(1, cdbl(1, 0));
	h2.setMultipole(2, cdbl(2, 0));
	{
		// cross check harmonics properly set
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);

		BOOST_CHECK_CLOSE(c1.real(), 1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_SMALL(c2.real(), 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}

	{
		// cross check harmonics properly set
		const cdbl c1 = h2.getMultipole(1);
		const cdbl c2 = h2.getMultipole(2);
		const cdbl c3 = h2.getMultipole(3);

		BOOST_CHECK_SMALL(c1.real(), 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), 2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}

	h3 = (h2 + h);

	{
		// cross check harmonics properly set and maintained after math
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);

		BOOST_CHECK_CLOSE(c1.real(), 1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_SMALL(c2.real(), 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}

	{
		// cross check harmonics properly set
		const cdbl c1 = h2.getMultipole(1);
		const cdbl c2 = h2.getMultipole(2);
		const cdbl c3 = h2.getMultipole(3);

		BOOST_CHECK_SMALL(c1.real(), 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), 2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}
	{
		// properly assigned
		// cross check harmonics properly set
		const cdbl c1 = h3.getMultipole(1);
		const cdbl c2 = h3.getMultipole(2);
		const cdbl c3 = h3.getMultipole(3);

		BOOST_CHECK_CLOSE(c1.real(), 1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), 2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}


}

BOOST_AUTO_TEST_CASE(test210_mul_multipoles_inplace)
{

	tsc::TwoDimensionalMultipoles h(0e0);
	const double b1 = 5e0, b2 = 7e0, scale = 3e0;
	h.setMultipole(1, cdbl(b1, 0));
	h.setMultipole(2, cdbl(b2, 0));

	{
		// properly assigned
		// cross check harmonics properly set
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);
		BOOST_CHECK_CLOSE(c1.real(), b1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), b2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}

	h *= scale;
	{
		// properly assigned
		// cross check harmonics properly set
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);
		BOOST_CHECK_CLOSE(c1.real(), b1 * scale , 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), b2 * scale , 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}
}

BOOST_AUTO_TEST_CASE(test220_mul_multipoles)
{

	tsc::TwoDimensionalMultipoles h(0e0);
	const double b1 = 5e0, b2 = 7e0, scale = 3e0;
	h.setMultipole(1, cdbl(b1, 0));
	h.setMultipole(2, cdbl(b2, 0));

	{
		// properly assigned
		// cross check harmonics properly set
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);
		BOOST_CHECK_CLOSE(c1.real(), b1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), b2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}

	tsc::TwoDimensionalMultipolesKnobbed h2 = h * scale;
	{
		// properly assigned
		// cross check harmonics not modified
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);
		BOOST_CHECK_CLOSE(c1.real(), b1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), b2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}
	{
		// properly assigned
		// cross check harmonics properly set
		const cdbl c1 = h2.getMultipole(1);
		const cdbl c2 = h2.getMultipole(2);
		const cdbl c3 = h2.getMultipole(3);
		BOOST_CHECK_CLOSE(c1.real(), b1 * scale , 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), b2 * scale , 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}
}


static std::complex<double> ret_zero(const std::complex<double> val)
{
	return std::complex<double>(0e0, 0e0);
}

BOOST_AUTO_TEST_CASE(test240_mul_multipoles_vector_inplace)
{

	const double b1 = 5e0, b2 = 7e0, s1 = 11e0, s2=13e0;

	tsc::TwoDimensionalMultipoles h(0e0);
	std::vector<std::complex<double>> scale(h.getCoeffs().size());
	h.setMultipole(1, cdbl(b1, 0));
	h.setMultipole(2, cdbl(b2, 0));

	std::transform(scale.begin(), scale.end(), scale.begin(),  ret_zero);
	scale[0] = s1;
	scale[1] = s2;

	{
		// properly assigned
		// cross check harmonics properly set
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);
		BOOST_CHECK_CLOSE(c1.real(), b1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), b2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}

	h *= scale;
    {
		// properly assigned
		// cross check harmonics properly set
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);
		BOOST_CHECK_CLOSE(c1.real(), b1 * s1 , 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), b2 * s2 , 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}
}

BOOST_AUTO_TEST_CASE(test250_mul_multipoles_vector)
{

	const double b1 = 5e0, b2 = 7e0, s1 = 1./11e0, s2=1./13e0;

	tsc::TwoDimensionalMultipoles h(0e0);
	std::vector<std::complex<double>> scale(h.getCoeffs().size());
	h.setMultipole(1, cdbl(b1, 0));
	h.setMultipole(2, cdbl(b2, 0));

	std::transform(scale.begin(), scale.end(), scale.begin(),  ret_zero);
	scale[0] = s1;
	scale[1] = s2;

	{
		// properly assigned
		// cross check harmonics properly set
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);
		BOOST_CHECK_CLOSE(c1.real(), b1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), b2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}

	tsc::TwoDimensionalMultipoles h2 = h * scale;

	{
		// properly assigned
		// cross check h unchanged
		const cdbl c1 = h.getMultipole(1);
		const cdbl c2 = h.getMultipole(2);
		const cdbl c3 = h.getMultipole(3);
		BOOST_CHECK_CLOSE(c1.real(), b1, 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), b2, 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}

	{
		// properly assigned
		// cross check harmonics properly set
		const cdbl c1 = h2.getMultipole(1);
		const cdbl c2 = h2.getMultipole(2);
		const cdbl c3 = h2.getMultipole(3);
		BOOST_CHECK_CLOSE(c1.real(), b1 * s1 , 1e-12);
		BOOST_CHECK_SMALL(c1.imag(), 1e-12);
		BOOST_CHECK_CLOSE(c2.real(), b2 * s2 , 1e-12);
		BOOST_CHECK_SMALL(c2.imag(), 1e-12);
		BOOST_CHECK_SMALL(c3.real(), 1e-12);
		BOOST_CHECK_SMALL(c3.imag(), 1e-12);
	}
}

#if 1
BOOST_AUTO_TEST_CASE(test300_mul_tpsa_)
{
    const int nv = 4, no = 4;
    auto a_desc = std::make_shared<gtpsa::desc>(nv, no);

    auto t = gtpsa::CTpsaOrComplex(0e0);
    //tsc::TwoDimensionalMultipolesKnobbed<tsc::TpsaVariantType> h(t);
    tsc::TwoDimensionalMultipolesTpsa h(t);

    auto x = gtpsa::tpsa(a_desc, 3);
    auto y = gtpsa::tpsa(a_desc, 3);
    auto Bx = gtpsa::tpsa(a_desc, 3);
    auto By = gtpsa::tpsa(a_desc, 3);

    x.setName("x");
    y.setName("y");
    Bx.setName("Bx");
    By.setName("By");

    // a sextupole
    auto c3 = gtpsa::ctpsa(a_desc, 3);
    c3.set({0,0},{1e-4, 0});
    // a decapole
    auto c5 = gtpsa::ctpsa(a_desc, 3);
    c5.set({0,0},{3e-4, 0});

    // Set identity by hand not avialable her
    x.setv(1, {1, 0, 0, 0, 0, 0});
    y.setv(1, {0, 1, 0, 0, 0, 0});

    c3.setv(1, {0, 0, 1, 0});
    c5.setv(1, {0, 0, 0, 1});
    // a rather standard dipole
    auto c1 =  gtpsa::CTpsaOrComplex(std::complex<double>(1e0, 0e0));
    h.setMultipole(1, c1);
    h.setMultipole(3, c3);
    h.setMultipole(5, c5);

    h.field(x, y, &Bx, &By);

    {
	// x "curvature"
	double field_from_c3 = By.get("2");
	double check = c3.cst().real();
	BOOST_CHECK_CLOSE(field_from_c3, check, 1e-12);
    }
    {
	// y "curvature"
	double field_from_c3 = By.get("02");
	double check = c3.cst().real();
	BOOST_CHECK_CLOSE(field_from_c3, -check, 1e-12);
    }
    {
	// c3 on c3
	double c3_on_c3 = By.get("201");
	double check = c3.cst().real();
	BOOST_CHECK_CLOSE(c3_on_c3, 1e0, 1e-12);
    }
    {
	// x "curvature"
	double c3_on_c3 = By.get("021");
	double check = c3.cst().real();
	BOOST_CHECK_CLOSE(c3_on_c3, 1e0, 1e-12);
    }
    {
	// y "curvature"
	double field_from_c3 = Bx.get("11");
	double check = c3.cst().real();
	BOOST_CHECK_CLOSE(field_from_c3, -check,  1e-12);
    }

    {
	const double eps = 1e-6;
	By.print(0, eps, 0, 0);
	Bx.print(0, eps, 0, 0);
    }

}
#endif

// include dedicated honer complex tests
// check that the last coefficient is honoured correctly ...
BOOST_AUTO_TEST_CASE(test310_hohner_complex_last_coeff)
{
    int h_max = 4;
    tsc::TwoDimensionalMultipolesTpsa h(0e0, h_max);
    const std::complex<double> c4(23e-4, -355e-4/113e0);
    std::vector<std::complex<double>> coeffs(h_max);
    coeffs.at(h_max -1) = c4;

    // consider x / rref
    const double x = 1, y = .5;
    std::complex<double> z_prop (x, y);
    const auto h2 = h.clone();
    const auto f = thor_scsi::core::honer_complex(coeffs, z_prop);

    const auto ref_field = std::pow(z_prop, (4-1)) * c4;
    BOOST_CHECK_CLOSE(ref_field.real(), f.real(), 1e-12);
    BOOST_CHECK_CLOSE(ref_field.imag(), f.imag(), 1e-12);
}
