#define BOOST_TEST_MODULE accelerator
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#include <thor_scsi/std_machine/accelerator.h>
#include <thor_scsi/std_machine/std_machine.h>
#include <thor_scsi/elements/drift.h>
#include <thor_scsi/elements/marker.h>
#include <thor_scsi/elements/cavity.h>
#include <thor_scsi/elements/quadrupole.h>
#include <thor_scsi/elements/sextupole.h>
#include <thor_scsi/elements/octupole.h>
#include <thor_scsi/elements/bending.h>
#include <thor_scsi/elements/standard_observer.h>
#include <thor_scsi/elements/standard_aperture.h>
#include <thor_scsi/core/config.h>
#include <sstream>
#include <string>

namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;
namespace ts = thor_scsi;

int reg_done = ts::register_elements();

auto a_desc = std::make_shared<gtpsa::desc>(1, 6);
auto tpsa_ref = gtpsa::tpsa(a_desc, mad_tpsa_default);


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
		"cavh1t8r: Cavity, Frequency = 500e6, Voltage = 0.5e6, HarmonicNumber=538;\n"
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
	const std::string txt(
		"bend   : Bending, L = 1.1, T = 20, K =-1.2, T1 = 5, T2 = 7,"
		"         N=9, method=4;"
		"mini_cell : LINE = (bend);"
		);

	GLPSParser parse;
	Config *C = parse.parse_byte(txt);
	auto machine = ts::Accelerator(*C);

	auto cell = machine[0];
	auto elem = std::dynamic_pointer_cast<tse::BendingType>(cell);
	BOOST_TEST(elem !=  nullptr);

	BOOST_CHECK(elem->getMainMultipoleNumber() == 1);
	BOOST_CHECK_CLOSE(elem->getLength(),                        1.1,  1e-12);
	BOOST_CHECK_CLOSE(elem->getMainMultipoleStrength().real(),    0,  1e-12);
	BOOST_CHECK_SMALL(elem->getMainMultipoleStrength().imag(),        1e-12);

	// K correctly set
	auto c2 = elem->getMultipoles()->getMultipole(2);
	BOOST_CHECK_CLOSE(c2.real(),   -1.2,  1e-12);
	BOOST_CHECK_SMALL(c2.imag(),          1e-12);

	BOOST_CHECK_CLOSE(elem->getBendingAngle(),            20,  1e-12);
	BOOST_CHECK_CLOSE(elem->getEntranceAngle(),            5,  1e-12);
	BOOST_CHECK_CLOSE(elem->getExitAngle(),                7,  1e-12);

	BOOST_CHECK_EQUAL(elem->getNumberOfIntegrationSteps(), 9);

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

BOOST_AUTO_TEST_CASE(test80_standard_observer)
{
	const std::string txt(
		"Nquad = 12;"
		"q4m2d1r: Quadrupole, L = 0.5, K = 1.4, N = Nquad, Method = 4;"
		"mini_cell : LINE = (q4m2d1r);\n"
		);

	GLPSParser parse;
	Config *C = parse.parse_byte(txt);
	auto machine = ts::Accelerator(*C);

	for(auto& cv: machine){
		auto ob = std::make_shared<tse::StandardObserver>();
		cv->set_observer(std::dynamic_pointer_cast<tsc::Observer>(ob));
	}

	{
		boost::test_tools::output_test_stream output;
		for(auto& cv: machine){
			cv->show(output, 10);
		}
		BOOST_CHECK( !output.is_empty( false ) );

	}

	gtpsa::ss_vect<gtpsa::tpsa> ps(tpsa_ref);
	ps.set_identity();
	auto calc_config = tsc::ConfigType();

	int next_element = machine.propagate(calc_config, ps);
	BOOST_CHECK(next_element == 1);

	// Iteration still works after processing ...
	{
		boost::test_tools::output_test_stream output;
		for(auto& cv: machine){
			cv->show(output, 10);
		}
		BOOST_CHECK( !output.is_empty( false ) );

	}
	auto ob = machine[0]->observer();
	auto std_ob = std::dynamic_pointer_cast<tse::StandardObserver>(ob);
	BOOST_CHECK(std_ob);
	BOOST_CHECK(std_ob->getObservedIndex() == 0);
	BOOST_CHECK(std_ob->getObservedName() == "q4m2d1r");
	BOOST_CHECK(std_ob->hasTruncatedPowerSeriesA());
	BOOST_CHECK(!std_ob->hasPhaseSpace());

}


class Log2Stream : public tsc::Machine::Logger
{
public:
	Log2Stream(std::ostream *strm){this->m_strm = strm;}
	virtual ~Log2Stream() {}
	virtual void log(const tsc::Machine::LogRecord &r) override final {
		std::ostream *p_strm = nullptr;
		if(this->m_strm){
			p_strm = this->m_strm;
		} else {
			std::cerr << "Lost stream to log to" << std::endl;
			p_strm = &std::cerr;
		}
		std::ostream &strm = *p_strm;

		std::string msg(r.strm.str());
		strm << r.fname << ": " << r.lnum << " : " << msg;
		if(msg.empty() || msg[msg.size()-1] != '\n'){
			strm.put('\n');
		}
	}

private:
	std::ostream *m_strm = &std::cerr;

};

class Log2StringStream : public tsc::Machine::Logger
{
public:
	Log2StringStream(void) { return; this->reset();}
	virtual ~Log2StringStream() {}
	virtual void log(const tsc::Machine::LogRecord &r) override final {
		std::ostream &strm = this->m_strm;

		std::string msg(r.strm.str());
		strm << r.fname << ": " << r.lnum << " : " << msg;
		if(msg.empty() || msg[msg.size()-1] != '\n'){
			strm.put('\n');
		}
		// strm << "That's all folks" << std::endl;
	}

	void reset(void){
		const static std::stringstream initial;

		this->m_strm.str(std::string());
		this->m_strm.clear();
		this->m_strm.copyfmt(initial);
	}

	std::string get(void){
		return this->m_strm.str();
	}

private:
	std::ostringstream m_strm;

};

BOOST_AUTO_TEST_CASE(test100_logger_stream)
{

	return;
	const std::string txt(
		"start : Marker;"
		"mini_cell  : LINE = (start);"
		);

	GLPSParser parse;
	Config *C = parse.parse_byte(txt);
	auto machine = ts::Accelerator(*C);


	{
		boost::test_tools::output_test_stream output;
		auto logger = std::make_shared<Log2Stream>(Log2Stream(&output));
		machine.set_logger(logger);

		THOR_SCSI_LOG_ALWAYS(FINE)    << ("Fine test");
		BOOST_CHECK( !output.is_empty( false ) );

		// Should not be required
		// machine.set_logger(nullptr);
	}

	//auto logger = Log2Stream(&std::cout);
	// logger.set_trace(std::cout);

	THOR_SCSI_LOG_ALWAYS(FINE)    << ("Fine test");
	THOR_SCSI_LOG_ALWAYS(DEBUG)   << ("Debug test");
	THOR_SCSI_LOG_ALWAYS(INFO)    << ("Info test");
	THOR_SCSI_LOG_ALWAYS(WARN)    << ("Warn test");
	THOR_SCSI_LOG_ALWAYS(WARNING) << ("Warning test");
	THOR_SCSI_LOG_ALWAYS(ERROR)   << ("Error test");
}

BOOST_AUTO_TEST_CASE(test101_logger_stringstream)
{
	const std::string txt(
		"start : Marker;"
		"mini_cell  : LINE = (start);"
		);

	GLPSParser parse;
	Config *C = parse.parse_byte(txt);
	auto machine = ts::Accelerator(*C);

 	{
		auto p_log = std::make_shared<Log2StringStream>();
		machine.set_logger(p_log);
		THOR_SCSI_LOG_ALWAYS(FINE) << ("Fine test");
		std::cout << "logger output " << p_log->get() << std::endl;
	}


	// logger.set_trace(std::cout);
	{
	}

	THOR_SCSI_LOG_ALWAYS(FINE)    << ("Fine test 2");
	THOR_SCSI_LOG_ALWAYS(DEBUG)   << ("Debug test");
	THOR_SCSI_LOG_ALWAYS(INFO)    << ("Info test");
	THOR_SCSI_LOG_ALWAYS(WARN)    << ("Warn test");
	THOR_SCSI_LOG_ALWAYS(WARNING) << ("Warning test");
	THOR_SCSI_LOG_ALWAYS(ERROR)   << ("Error test");


	auto logger = machine.get_logger();
	auto strm_logger = std::dynamic_pointer_cast<Log2StringStream>(logger);

	std::cout << "Stream Logger\n" << strm_logger->get() << std::endl;

	strm_logger->reset();

	THOR_SCSI_LOG(FINE)    << ("Fine test 2");
	THOR_SCSI_LOG(DEBUG)   << ("Debug test");
	THOR_SCSI_LOG(INFO)    << ("Info test");
	THOR_SCSI_LOG(WARN)    << ("Warn test");
	THOR_SCSI_LOG(WARNING) << ("Warning test");
	THOR_SCSI_LOG(ERROR)   << ("Error test");

	strm_logger->reset();
	machine.set_log_level(THOR_SCSI_DEBUG);
	THOR_SCSI_LOG(DEBUG)   << ("Debug test");

	std::cout << "Stream Logger level test\n" << strm_logger->get() << std::endl;

}

BOOST_AUTO_TEST_CASE(test120_rectangular_aperture)
{
	const std::string txt(
		"Nquad = 12;"
		"q4m2d1r: Quadrupole, L = 0.5, K = 1.4, N = Nquad, Method = 4;"
		"mini_cell : LINE = (q4m2d1r);\n"
		);

	GLPSParser parse;
	Config *C = parse.parse_byte(txt);
	auto machine = ts::Accelerator(*C);
	const double width = 100e-3, height = 50e-3;

	// register the aperture
	for(auto& cv: machine){
		auto ap = std::make_shared<tse::RectangularAperture>(width, height);
		auto elem = std::dynamic_pointer_cast<tse::ElemType>(cv);
		if(!elem){
			throw std::runtime_error("Can not cast to elemtype");
		}
		elem->setAperture(std::dynamic_pointer_cast<tsc::TwoDimensionalAperture>(ap));
	}

	auto cv = machine.find("q4m2d1r");
	auto elem = std::dynamic_pointer_cast<tse::ElemType>(cv);

	// check that it gives at least some output
	auto ap = elem->getAperture();
	{
		boost::test_tools::output_test_stream output;
		output << ap;
		BOOST_CHECK( !output.is_empty( false ) );
	}

	// check phase spaces vectors double ...
	{
		gtpsa::ss_vect<double> ps(0e0);
		ps.set_zero();

		bool not_lost = elem->checkAmplitude(ps);
		BOOST_CHECK( not_lost);
	}
	{
		gtpsa::ss_vect<double> ps(0e0);
		ps.set_zero();
		// Slightly out is out
		ps[x_] = width + 1e-3;

		bool not_lost = elem->checkAmplitude(ps);
		BOOST_CHECK( !not_lost);
	}

	{
		gtpsa::ss_vect<gtpsa::tpsa> ps(tpsa_ref);
		ps.set_identity();
		ps[x_] = width + 1e-3;
		bool not_lost = elem->checkAmplitude(ps);
		BOOST_CHECK( !not_lost);
	}

	// check that propagate will identify it
	{
		auto calc_config = tsc::ConfigType();
		gtpsa::ss_vect<double> ps(0e0);
		ps.set_zero();
		ps[x_] = width / 2.0;

		calc_config.lossplane = 0;
		int last_element = machine.propagate(calc_config, ps);
		BOOST_CHECK(calc_config.lossplane == 0);
	}
	// check that propagate will identify it
	{
		auto calc_config = tsc::ConfigType();
		gtpsa::ss_vect<double> ps(0e0);
		ps.set_zero();
		ps[x_] = width + 1e-3;

		calc_config.lossplane = 0;
		int last_element = machine.propagate(calc_config, ps);
		BOOST_CHECK(calc_config.lossplane == tse::PlaneKind::Horizontal);
	}

}


BOOST_AUTO_TEST_CASE(test130_add_marker_start)
{
	const std::string txt(
		"Nquad = 12;"
		"q4m2d1r: Quadrupole, L = 0.5, K = 1.4, N = Nquad, Method = 4;"
		"mini_cell : LINE = (q4m2d1r);\n"
		);

	GLPSParser parse;
	Config *C = parse.parse_byte(txt);

	boost::test_tools::output_test_stream output;
	auto logger = std::make_shared<Log2Stream>(Log2Stream(&output));

	auto machine = ts::Accelerator(*C, true);
	auto cv = machine.at(0);

	auto marker = std::dynamic_pointer_cast<tse::MarkerType>(cv);
	BOOST_CHECK( (marker) );
	BOOST_CHECK( (marker->name == "Start") );

	auto cv2 = machine.at(1);
	BOOST_CHECK( (cv2->name == "q4m2d1r") );
	auto mpole = std::dynamic_pointer_cast<tse::MpoleType>(cv2);
	BOOST_CHECK( (mpole) );
	auto quad = std::dynamic_pointer_cast<tse::QuadrupoleType>(cv2);
	BOOST_CHECK( (quad) );


}

BOOST_AUTO_TEST_CASE(test131_no_marker_add_start_required)
{
	const std::string txt(
		"Nquad = 12;"
		"start_different_name: Marker;"
		"q4m2d1r: Quadrupole, L = 0.5, K = 1.4, N = Nquad, Method = 4;"
		"mini_cell : LINE = (start_different_name, q4m2d1r);\n"
		);

	GLPSParser parse;
	Config *C = parse.parse_byte(txt);

	boost::test_tools::output_test_stream output;
	auto logger = std::make_shared<Log2Stream>(Log2Stream(&output));

	auto machine = ts::Accelerator(*C, true);
	auto cv = machine.at(0);

	auto marker = std::dynamic_pointer_cast<tse::MarkerType>(cv);
	BOOST_CHECK((marker) );
	BOOST_CHECK((marker->name == "start_different_name") );

	auto cv2 = machine.at(1);
	BOOST_CHECK( (cv2->name == "q4m2d1r") );
	auto mpole = std::dynamic_pointer_cast<tse::MpoleType>(cv2);
	BOOST_CHECK( (mpole) );
	auto quad = std::dynamic_pointer_cast<tse::QuadrupoleType>(cv2);
	BOOST_CHECK( (quad) );
}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
