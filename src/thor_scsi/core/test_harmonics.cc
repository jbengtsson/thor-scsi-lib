#define BOOST_TEST_MODULE harmonics
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <thor_scsi/core/harmonics.h>
#include <cmath>

namespace tsc = thor_scsi::core;

BOOST_AUTO_TEST_CASE(test00_init_zero)
{

	auto h = tsc::PlanarHarmonics();

	// coefficient by interface
	BOOST_CHECK_CLOSE(h.getHarmonic(1).real(), 0, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(1).imag(), 0, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(21).real(), 0, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(21).imag(), 0, 1e-42);

	// access to coefficients
	auto coeffs = h.getCoeffs();
	BOOST_CHECK_CLOSE(coeffs.at(0).real(), 0, 1e-42);
	BOOST_CHECK_CLOSE(coeffs.at(0).imag(), 0, 1e-42);
	BOOST_CHECK_CLOSE(coeffs.at(20).real(), 0, 1e-42);
	BOOST_CHECK_CLOSE(coeffs.at(20).imag(), 0, 1e-42);

}

BOOST_AUTO_TEST_CASE(test01_set_harmonic_dipole)
{
	auto h = tsc::PlanarHarmonics();

	// pure dipole
	h.setHarmonic(1, tsc::cdbl(1, 0));
	BOOST_CHECK_CLOSE(h.getHarmonic(1).real(), 1, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(1).imag(), 0, 1e-42);

	BOOST_CHECK_CLOSE(h.getHarmonic(2).real(), 0, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(2).imag(), 0, 1e-42);

	auto coeffs = h.getCoeffs();
	BOOST_CHECK_CLOSE(coeffs.at(0).real(), 1, 1e-42);
	BOOST_CHECK_CLOSE(coeffs.at(0).imag(), 0, 1e-42);
	BOOST_CHECK_CLOSE(coeffs.at(1).real(), 0, 1e-42);
	BOOST_CHECK_CLOSE(coeffs.at(1).imag(), 0, 1e-42);

	{
		auto field =  h.field(tsc::cdbl(0,0));
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-22);
	}
	{
		auto field =  h.field(1);
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-22);
	}

	{
		auto field =  h.field(tsc::cdbl(0,1));
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-22);
	}

	h.applyRollAngle(M_PI/2);
	BOOST_CHECK_CLOSE(h.getHarmonic(1).imag(), 1, 1e-15);
	BOOST_CHECK_CLOSE(h.getHarmonic(2).real(), 0, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(2).imag(), 0, 1e-42);
	// to help boost on which accuracy is meant
	BOOST_CHECK_CLOSE(h.getHarmonic(1).real()+1, 1, 1e-42);

	// Notice by + i Bx
	auto field =  h.field(0);
	BOOST_CHECK_CLOSE(field.real() + 1, 0 + 1, 1e-22);
	BOOST_CHECK_CLOSE(field.imag(),  1, 1e-22);


}


BOOST_AUTO_TEST_CASE(test02_set_harmonic_quadrupole)
{

	auto h = tsc::PlanarHarmonics();

	// pure quadrupole
	h.setHarmonic(2, tsc::cdbl(1, 0));
	BOOST_CHECK_CLOSE(h.getHarmonic(2).real(), 1, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(2).imag(), 0, 1e-42);


	h.applyRollAngle(M_PI/2);
	BOOST_CHECK_CLOSE(h.getHarmonic(2).real(), -1, 1e-15);
	// to help boost on which accuracy is meant
	BOOST_CHECK_CLOSE(h.getHarmonic(2).imag() + 1, 1, 1e-12);


}


BOOST_AUTO_TEST_CASE(test03_set_scale_sextupole)
{
	// pure quadrupole
	auto h = tsc::PlanarHarmonics();
	h.setHarmonic(3, tsc::cdbl(.1, -.3));
	BOOST_CHECK_CLOSE(h.getHarmonic(3).real(),  0.1, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(3).imag(), -0.3, 1e-42);

	h *= 3;
	BOOST_CHECK_CLOSE(h.getHarmonic(3).real(),  0.3, 1e-13);
	BOOST_CHECK_CLOSE(h.getHarmonic(3).imag(), -0.9, 1e-13);

	auto h2 = h * (1/3.0);
	BOOST_CHECK_CLOSE(h2.getHarmonic(3).real(),  0.1, 1e-13);
	BOOST_CHECK_CLOSE(h2.getHarmonic(3).imag(), -0.3, 1e-13);
}


BOOST_AUTO_TEST_CASE(test04_rotate_quadrupole)
{
  {
	auto h = tsc::PlanarHarmonics();

	// pure quadrupole
	h.setHarmonic(2, tsc::cdbl(1, 0));
	BOOST_CHECK_CLOSE(h.getHarmonic(2).real(), 1, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(2).imag(), 0, 1e-42);


	tsc::cdbl I(0,1);

	// Gradient interpolation
	{
		auto field =  h.field(0);
		BOOST_CHECK_CLOSE(field.real(), 0, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-22);
	}
	{
		auto field =  h.field(-0.5);
		BOOST_CHECK_CLOSE(field.real(), -.5, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-22);
	}
	// check convienience function
	{
		auto field =  h.field(.5, -0.5);
		BOOST_CHECK_CLOSE(field.real(), .5, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), -.5, 1e-22);
	}
	{
		auto field =  h.field(1);
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-22);
	}
	{
		auto field =  h.field(-0.5 * I);
		BOOST_CHECK_CLOSE(field.real(), 0, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), -0.5, 1e-22);
	}
	{
		auto field =  h.field(I);
		BOOST_CHECK_CLOSE(field.real(), 0, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), 1, 1e-22);
	}


	h.applyRollAngle(M_PI/2);
	BOOST_CHECK_CLOSE(h.getHarmonic(2).real(), -1, 1e-15);
	// to help boost on which accuracy is meant
	BOOST_CHECK_CLOSE(h.getHarmonic(2).imag() + 1, 1, 1e-12);

	h.applyRollAngle(M_PI/4);
	// to help boost on which accuracy is meant
	BOOST_CHECK_CLOSE(h.getHarmonic(2).real() + 1, 1, 1e-13);
	// todo: check with paper that it is correct ..
	BOOST_CHECK_CLOSE(h.getHarmonic(2).imag(), -1, 1e-12);

	// Gradient interpolation ... matches values above if
	// position is rotated too
	{
		auto field =  h.field(0);
		BOOST_CHECK_CLOSE(field.real(), 0, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-22);
	}
	{
		auto field =  h.field(-0.5 * I);
		BOOST_CHECK_CLOSE(field.real(), -.5, 1e-22);
		BOOST_CHECK_CLOSE(field.imag() + 1, 1, 1e-22);
	}
	{
		auto field =  h.field(1.0 * I);
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-22);
		BOOST_CHECK_CLOSE(field.imag() +1, 1, 1e-13);
	}
	{
		auto field =  h.field(-0.5 * I * I);
		BOOST_CHECK_CLOSE(field.real()+1 , 1, 1e-13);
		BOOST_CHECK_CLOSE(field.imag(), -0.5, 1e-22);
	}
	{
		auto field =  h.field(I * I);
		BOOST_CHECK_CLOSE(field.real()+1, 1, 1e-13);
		BOOST_CHECK_CLOSE(field.imag(), 1, 1e-22);
	}

  }
}


BOOST_AUTO_TEST_CASE(test05_rotate_octupole)
{
	auto h = tsc::PlanarHarmonics();

	const int n = 4;
	// pure quadrupole
	h.setHarmonic(n, tsc::cdbl(1, 0));
	BOOST_CHECK_CLOSE(h.getHarmonic(n).real(), 1, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(n).imag(), 0, 1e-42);

	tsc::cdbl I(0,1);
	// cubic  interpolation
	{
		auto field =  h.field(0);
		BOOST_CHECK_CLOSE(field.real(), 0, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-22);
	}
	{
		double val = -.5;
		auto field =  h.field(-0.5);
		BOOST_CHECK_CLOSE(field.real(), val * val * val, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-22);
	}
	{
		auto field =  h.field(1);
		BOOST_CHECK_CLOSE(field.real(), 1, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), 0, 1e-22);
	}
	{
		double val = -.5;
		auto field =  h.field(-0.5 * I);
		BOOST_CHECK_CLOSE(field.real(), 0, 1e-22);
		BOOST_CHECK_CLOSE(field.imag(), -val * val * val, 1e-22);
	}
	{
		auto field =  h.field(I);
		BOOST_CHECK_CLOSE(field.real(), 0, 1e-22);
		// 90 degrees ...
		BOOST_CHECK_CLOSE(field.imag(), -1, 1e-22);
	}


}

BOOST_AUTO_TEST_CASE(test06_rotate_octupole_small_angle)
{
	auto h = tsc::PlanarHarmonics();

	const int n = 4;
	// pure quadrupole
	h.setHarmonic(n, tsc::cdbl(1, 0));

	tsc::cdbl I(0,1);
	const double small_angle = 1e-3, phase = small_angle * n;
	const tsc::cdbl cang = exp(small_angle * I);


	const tsc::cdbl
		pos0(0,0),
		pos1 = -0.5 * I,
		pos2(1,0);

	const tsc::cdbl
		field0 = h.field(pos0 * cang),
		field1 = h.field(pos1 * cang),
		field2 = h.field(pos2 * cang);

	h.applyRollAngle(small_angle);
	/* small angle: first check with approximate change */
	BOOST_CHECK_CLOSE(h.getHarmonic(n).real(), 1, 1e-3);
	// 1 for calculating relative size
	BOOST_CHECK_CLOSE(h.getHarmonic(n).imag()+1, 1 + phase, 1e-3);

	// calculate phase as given by equation
	BOOST_CHECK_CLOSE(h.getHarmonic(n).real(), 1 * cos(phase), 1e-42);
	// 1 for calculating relative size
	BOOST_CHECK_CLOSE(h.getHarmonic(n).imag()+1, 1 + 1 * sin(phase), 1e-13);

	const tsc::cdbl
		check0 = h.field(pos0),
		check1 = h.field(pos1),
		check2 = h.field(pos2);

	BOOST_CHECK_CLOSE(field0.real(), check0.real(), 1e-42);
	BOOST_CHECK_CLOSE(field0.imag(), check0.imag(), 1e-42);

	BOOST_CHECK_CLOSE(field1.real(), check1.real(), 1e-42);
	BOOST_CHECK_CLOSE(field1.imag(), check1.imag(), 1e-42);

	BOOST_CHECK_CLOSE(field2.real(), check2.real(), 1e-42);
	BOOST_CHECK_CLOSE(field2.imag(), check2.imag(), 1e-42);

}

BOOST_AUTO_TEST_CASE(test06_translate_quadrupole_zero)
{

	auto h = tsc::PlanarHarmonics();

	// pure quadrupole .. assuming harmonics have been checked
	h.setHarmonic(2, tsc::cdbl(1, 0));


	// checking interface
	h.applyTranslation(0);
	h.applyTranslation(0, 0);

	tsc::cdbl dz = (0,0);
	h.applyTranslation(dz);


	BOOST_CHECK_CLOSE(h.getHarmonic(2).real(), 1, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(2).imag(), 0, 1e-42);

}


BOOST_AUTO_TEST_CASE(test07_translate_quadrupole_one)
{

	tsc::PlanarHarmonics h = tsc::PlanarHarmonics(4);

	// pure quadrupole .. assuming harmonics have been checked
	h.setHarmonic(2, tsc::cdbl(1, 0));


	// checking interface
	h.applyTranslation(1, 0);

	BOOST_CHECK_CLOSE(h.getHarmonic(2).real(), 1, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(2).imag(), 0, 1e-42);

	// Dipole harmonic
	BOOST_CHECK_CLOSE(h.getHarmonic(1).real(), 1, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(1).imag(), 0, 1e-42);

	// no sextupole
	BOOST_CHECK_CLOSE(h.getHarmonic(3).real(), 0, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(3).imag(), 0, 1e-42);

	// no octupole
	BOOST_CHECK_CLOSE(h.getHarmonic(4).real(), 0, 1e-42);
	BOOST_CHECK_CLOSE(h.getHarmonic(4).imag(), 0, 1e-42);

}

BOOST_AUTO_TEST_CASE(test08_translate_quadrupole_test_field)
{
	tsc::PlanarHarmonics h = tsc::PlanarHarmonics(4);

	// pure quadrupole .. assuming harmonics have been checked
	h.setHarmonic(2, tsc::cdbl(1, 0));
	const tsc::cdbl pos0(0, 0), pos1(0, .5), pos2(1., 0);
	const tsc::cdbl dz(1, 0);

	// probe field where coordinates will be
	const tsc::cdbl
		field0 = h.field(pos0 + dz),
		field1 = h.field(pos1 + dz),
		field2 = h.field(pos2 + dz);


	// checking interface
	h.applyTranslation(dz);


	const tsc::cdbl
		check0 = h.field(pos0),
		check1 = h.field(pos1),
		check2 = h.field(pos2);

	BOOST_CHECK_CLOSE(field0.real(), check0.real(), 1e-42);
	BOOST_CHECK_CLOSE(field0.imag(), check0.imag(), 1e-42);

	BOOST_CHECK_CLOSE(field1.real(), check1.real(), 1e-42);
	BOOST_CHECK_CLOSE(field1.imag(), check1.imag(), 1e-42);

	BOOST_CHECK_CLOSE(field2.real(), check2.real(), 1e-42);
	BOOST_CHECK_CLOSE(field2.imag(), check2.imag(), 1e-42);

}

BOOST_AUTO_TEST_CASE(test09_translate_quadrupole_small)
{
	tsc::PlanarHarmonics h = tsc::PlanarHarmonics();
	const double unit = 1e-4;
	//  realistic example
	const tsc::cdbl dz(1e-4, .3e-3);
	const tsc::cdbl
		dodecapole = tsc::cdbl(.2, 3) * unit,
		icosapole = tsc::cdbl(.5, -.7) * unit;

	return;
	// pure quadrupole .. assuming harmonics have been checked
	h.setHarmonic(2, tsc::cdbl(1, 0));
	h.setHarmonic( 6, dodecapole);
	h.setHarmonic(10, icosapole);


	const tsc::cdbl pos0(0, 0), pos1(0, 17e-3), pos2(23e-3, 0);
	const tsc::cdbl
		field0 = h.field(pos0 + dz),
		field1 = h.field(pos1 + dz),
		field2 = h.field(pos2 + dz);


	h.applyTranslation(dz);

	const tsc::cdbl
		check0 = h.field(pos0),
		check1 = h.field(pos1),
		check2 = h.field(pos2);

	BOOST_CHECK_CLOSE(field0.real(), check0.real(), 1e-42);
	BOOST_CHECK_CLOSE(field0.imag(), check0.imag(), 1e-42);

	BOOST_CHECK_CLOSE(field1.real(), check1.real(), 1e-42);
	BOOST_CHECK_CLOSE(field1.imag(), check1.imag(), 1e-42);

	BOOST_CHECK_CLOSE(field2.real(), check2.real(), 1e-42);
	BOOST_CHECK_CLOSE(field2.imag(), check2.imag(), 1e-42);

}
