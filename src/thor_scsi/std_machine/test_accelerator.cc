#define BOOST_TEST_MODULE accelerator
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <thor_scsi/std_machine/accelerator.h>
#include <thor_scsi/std_machine/std_machine.h>
#include <thor_scsi/elements/drift.h>
#include <thor_scsi/elements/marker.h>
#include <thor_scsi/elements/cavity.h>
#include <thor_scsi/elements/quadrupole.h>
#include <thor_scsi/elements/sextupole.h>
#include <thor_scsi/elements/octupole.h>

namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;
namespace ts = thor_scsi;

int reg_done = ts::register_elements();


BOOST_AUTO_TEST_CASE(test10_drift)
{
	const std::string drift_txt(
		"d05l2t8r: Drift, L = 0.42;\n"
		"mini_cell : LINE = (d05l2t8r);\n");

	GLPSParser parse;
	Config *C = parse.parse_byte(drift_txt);
	auto machine = ts::Accelerator(*C);

	// proper checks required here
	auto cell = machine.find("d05l2t8r");
	BOOST_TEST(cell);

	auto drift = std::dynamic_pointer_cast<tse::DriftType>(cell);
	BOOST_TEST(drift);

	const double length_check = drift->getLength();

	BOOST_CHECK_CLOSE(length_check, 0.42, 1e-6);
	BOOST_CHECK_EQUAL(drift->name, "d05l2t8r");

}


BOOST_AUTO_TEST_CASE(test20_marker)
{
	const std::string marker_txt(
		"start : Marker;\n"
		"mini_cell : LINE = (start);\n");

	GLPSParser parse;
	Config *C = parse.parse_byte(marker_txt);
	auto machine = ts::Accelerator(*C);

	// proper checks required here
	auto cell = machine.find("start");
	BOOST_TEST(cell != nullptr);

	auto marker = std::dynamic_pointer_cast<tse::MarkerType>(cell);
	BOOST_TEST(marker !=  nullptr);

	const double length_check = marker->getLength();
	BOOST_CHECK_SMALL(length_check,  1e-6);
	BOOST_CHECK_EQUAL(marker->name, "start");

}

BOOST_AUTO_TEST_CASE(test30_cavity)
{

	const std::string cavity_txt(
		"cavh1t8r: Cavity, Frequency = 500e6, Voltage = 0.5e6;\n"
		"mini_cell : LINE = (cavh1t8r);\n"
		);

	GLPSParser parse;
	Config *C = parse.parse_byte(cavity_txt);
	auto machine = ts::Accelerator(*C);

	// proper checks required here
	auto cell = machine.find("cavh1t8r");
	BOOST_TEST(cell != nullptr);

	auto cavity = std::dynamic_pointer_cast<tse::CavityType>(cell);
	BOOST_TEST(cavity !=  nullptr);

	BOOST_CHECK_CLOSE(cavity->getFrequency(),  500e6, 1e-6);
	BOOST_CHECK_CLOSE(cavity->getVoltage(),    0.5e6, 1e-6);

	const double length_check = cavity->getLength();
	BOOST_CHECK_SMALL(length_check,  1e-6);
	BOOST_CHECK_EQUAL(cavity->name, "cavh1t8r");

}

BOOST_AUTO_TEST_CASE(test40_dipole)
{

}


BOOST_AUTO_TEST_CASE(test50_quadrupole)
{
	const std::string txt(
		"Nquad = 12;"
		"q4m2d1r: Quadrupole, L = 0.5, K = 1.4, N = Nquad, Method = 4;"
		"mini_cell : LINE = (q4m2d1r);\n"
		);

	GLPSParser parse;
	Config *C = parse.parse_byte(txt);
	auto machine = ts::Accelerator(*C);

	auto cell = machine[0];
	auto elem = std::dynamic_pointer_cast<tse::QuadrupoleType>(cell);
	BOOST_TEST(elem !=  nullptr);

	BOOST_CHECK(elem->getMainMultipoleNumber() == 2);
	BOOST_CHECK_CLOSE(elem->getLength(),                        0.5,  1e-12);
	BOOST_CHECK_CLOSE(elem->getMainMultipoleStrength().real(),  1.4,  1e-12);
	BOOST_CHECK_SMALL(elem->getMainMultipoleStrength().imag(),        1e-12);

}

BOOST_AUTO_TEST_CASE(test60_sextupole)
{
	const std::string txt(
		"Nsext = 10;"
		"s4m2d1rl: Sextupole, L = 0.08, K = 28, N = Nsext, Method = 4;\n"
		"mini_cell : LINE = (s4m2d1rl);\n"
		);


	GLPSParser parse;
	Config *C = parse.parse_byte(txt);
	auto machine = ts::Accelerator(*C);
	auto cell = machine[0];
	auto elem = std::dynamic_pointer_cast<tse::SextupoleType>(cell);
	BOOST_TEST(elem !=  nullptr);

	BOOST_CHECK_CLOSE(elem->getLength(),                        0.08, 1e-12);
	BOOST_CHECK_CLOSE(elem->getMainMultipoleStrength().real(), 28,    1e-12);
	BOOST_CHECK_SMALL(elem->getMainMultipoleStrength().imag(),        1e-12);
	BOOST_CHECK(elem->getMainMultipoleNumber() == 3);
}


BOOST_AUTO_TEST_CASE(test70_octupole)
{
	const std::string txt(
		"Nsext = 10;"
		"missing: Octupole, L = 0.281, K = 355, N=1;\n"
		"mini_cell : LINE = (missing);\n"
		);


	GLPSParser parse;
	Config *C = parse.parse_byte(txt);
	auto machine = ts::Accelerator(*C);
	auto cell = machine[0];
	auto elem = std::dynamic_pointer_cast<tse::OctupoleType>(cell);
	BOOST_TEST(elem !=  nullptr);

	BOOST_CHECK_CLOSE(elem->getLength(),                         0.281, 1e-12);
	BOOST_CHECK_CLOSE(elem->getMainMultipoleStrength().real(), 355,     1e-12);
	BOOST_CHECK_SMALL(elem->getMainMultipoleStrength().imag(),          1e-12);
	BOOST_CHECK(elem->getMainMultipoleNumber() == 4);
}



/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
