#define BOOST_TEST_MODULE accelerator
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#include <thor_scsi/std_machine/accelerator.h>
#include <thor_scsi/std_machine/std_machine.h>
#include <thor_scsi/core/config.h>
#include <thor_scsi/elements/drift.h>
#include <thor_scsi/elements/bending.h>

#include <sstream>
#include <string>
#include <iostream>

namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;
namespace ts = thor_scsi;

static int reg_done = ts::register_elements();

BOOST_AUTO_TEST_CASE(test10_two_elements_list)
{
#if 0
    const std::string drift_txt(
	"d05l2t8r: Drift, L = 0.42;\n"
	"mini_cell : LINE = (d05l2t8r);\n"
	);

    const std::string bending_txt(
	"bend   : Bending, L = 1.1, N=10, T = 20, K =-1.2, T1 = 5, T2 = 7, N=9, method=4;\n"
	"mini_cell : LINE = (bend);\n"
	);

    GLPSParser parse;
    Config *C_drift = parse.parse_byte(drift_txt);
    Config *C_bending = parse.parse_byte(bending_txt);


    // tse::DriftType drift(*C_drift);
    // tse::BendingType bend(*C_bending);
    // std::cout << *C_drift;
    // std::cout << *C_bending;

    std::vector<std::shared_ptr<thor_scsi::core::ElemType>> elems;
    elems.reserve(2);
    std::cout <<  "Building elems"  << std::endl;
    elems.push_back(std::make_shared<tse::DriftType>(*C_drift));
    std::cout <<  "Built drift"  << std::endl;
    elems.push_back(std::make_shared<tse::BendingType>(*C_bending));
     std::cout <<  "Built bend"  << std::endl;
    auto machine = ts::Accelerator(elems);

    BOOST_CHECK_EQUAL(machine.size(), 2);

    BOOST_CHECK_EQUAL(machine.at(0)->type_name(), "Drift"    );
    BOOST_CHECK_EQUAL(machine.at(0)->name,        "d05l2t8r" );
    BOOST_CHECK_EQUAL(machine.at(1)->type_name(), "Bending"  );
    BOOST_CHECK_EQUAL(machine.at(1)->name,        "bend"     );

#endif
}

BOOST_AUTO_TEST_CASE(test20_two_elements_confg)
{
    Config config_drift;
    Config config_bend;
    Config config;

    config_drift.set<std::string>( "type", "Drift"    );
    config_drift.set<std::string>( "name", "d05l2t8r" );


    config_bend .set<std::string>( "type", "Bending"  );
    config_bend .set<std::string>( "name", "bend"     );
    config_bend .set<double>     ( "L"   ,   0.42     );
    config_bend .set<double>     ( "N"   ,   9        );
    config_bend .set<double>     ( "T"   ,  20        );
    config_bend .set<double>     ( "T1"  ,   5        );
    config_bend .set<double>     ( "T1"  ,   7        );
    config_bend .set<double>     ( "K"   , - 1.2      );
    config_bend .set<double>     ( "K"   , - 1.2      );


    std::vector<Config> configs;
    configs.push_back(config_drift);
    configs.push_back(config_bend);
    config.setAny("elements", configs);

    std::cout << "Config " << config << std::endl;

    auto machine = ts::Accelerator(config);

    BOOST_CHECK_EQUAL(machine.size(), 2);

    BOOST_CHECK_EQUAL(machine.at(0)->type_name(), "Drift"    );
    BOOST_CHECK_EQUAL(machine.at(0)->name,        "d05l2t8r" );
    BOOST_CHECK_EQUAL(machine.at(1)->type_name(), "Bending"  );
    BOOST_CHECK_EQUAL(machine.at(1)->name,        "bend"     );

}
