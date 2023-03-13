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

int main(int argc, char * argv[])
{
#if 0
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

#endif
    std::vector<std::shared_ptr<thor_scsi::core::ElemType>> elems;
    elems.reserve(2);
    std::cout <<  "Building elems"  << std::endl;

#if 1
    {
	/*
	  const std::string bending_txt(
	    "bend: Bending, L = 1.1, N=10.0, T = 20, K =-1.2, T1 = 5, T2 = 7, method=4;\n"
	    "mini_cell : LINE = (bend);\n"
	    );
	*/
	const std::string bending_txt(
	    "bend   : Bending, L = 1.1, T = 20, K =-1.2, T1 = 5, T2 = 7,"
	    "         N=9, method=4;"
	    "mini_cell : LINE = (bend);"
	    );


	GLPSParser parse;
	Config *C_bending = parse.parse_byte(bending_txt);
	std::cout << "bending check " << std::endl;
	auto val = C_bending->getAny("L");
	auto dval = C_bending->get<double>("N");
	std::cout << dval  << std::endl;
	//C_bending.
	elems.push_back(std::make_shared<tse::BendingType>(*C_bending));
	std::cout <<  "Built bend"  << std::endl;
    }
#endif
    {
	GLPSParser parse;

	const std::string drift_txt(
	    "d05l2t8r: Drift, L = 0.42;\n"
	    "mini_cell : LINE = (d05l2t8r);\n"
	    );

	Config *C_drift = parse.parse_byte(drift_txt);
	auto configs = C_drift->get<std::vector<Config>>("elements");
	const Config c = configs.at(0);
	tse::DriftType *drift = new tse::DriftType(c);
	elems.push_back(std::shared_ptr<tse::DriftType>(drift));
	std::cout <<  "Built drift"  << std::endl;
    }
    auto machine = ts::Accelerator(elems);



    std::cout << "drift ? " << machine.at(0)->repr() << std::endl;
    //std::cout << "bend  ? " << *(machine.at(1)) << std::endl;
    return 0;
}
