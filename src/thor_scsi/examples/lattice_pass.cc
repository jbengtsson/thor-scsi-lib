#include <boost/program_options.hpp>
#include <thor_scsi/std_machine/std_machine.h>
#include <thor_scsi/core/elements_basis.h>
#include <tps/ss_vect.h>
#include <stdlib.h>
#include <ostream>
#include <memory>

namespace po = boost::program_options;
namespace tsc = thor_scsi::core;

po::variables_map vm;
static void process_cmd_line(int argc, char *argv[])
{

	// Declare the supported options.
	std::string progname(argv[0]);

	po::options_description desc(progname + ":allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("lattice_file", po::value<std::string>(),                 "lattice file name")
		("x_pos,x",      po::value<double>()->default_value(0e0), "horizontal position        (left hand coordinate system)")
		("px",           po::value<double>()->default_value(0e0), "horizontal direction       (left hand coordinate system)")
		("y_pos,y",      po::value<double>()->default_value(0e0), "vertical position          (left hand coordinate system)")
		("py",           po::value<double>()->default_value(0e0), "vertical direction         (left hand coordinate system)")
		("ct",           po::value<double>()->default_value(0e0), "relative phase deviation   (left hand coordinate system)")
		("delta",        po::value<double>()->default_value(0e0), "relative impulse deviation (left hand coordinate system)")
		("inspect,i",    po::value<std::string>()->default_value(""), "name of element to inspect")
		("number,n",     po::value<int>()->default_value(0), "number of element to inspect")
		("n_turns",      po::value<int>()->default_value(1),      "propagate n turns (set to zero for none)" )
		("dump_lattice", po::value<bool>()->default_value(false), "dump read in lattice (to stdout)" )
		("verbose,v",    po::value<bool>()->default_value(false), "verbose output" )
		;

	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << "\n";
		exit(1);
	}
}


static void process(tsc::Machine& machine)
{
	bool verbose = vm["verbose"].as<bool>();

	if(vm["dump_lattice"].as<bool>()){
		std::cout << "Machine configuration " << std::endl;
		for(auto cv : machine){
			auto& elem = dynamic_cast<tsc::ElemType&>(*cv);
			std::cout << elem << std::endl;
		}
	}

	const int n_turns = vm["n_turns"].as<int>();

	ss_vect<double> ps;
	ps[x_]     = vm["x_pos"].as<double>();
	ps[px_]    = vm["px"].as<double>();
	ps[y_]     = vm["y_pos"].as<double>();
	ps[py_]    = vm["py"].as<double>();
	ps[ct_]    = vm["ct"].as<double>();
	ps[delta_] = vm["delta"].as<double>();


	tsc::ConfigType calc_config;

	auto element_name = vm["inspect"].as<std::string>();
	if(element_name != ""){
		int number =  vm["number"].as<int>();
		auto elem = machine.find(element_name, number);
		std::cout << "Machine element " << element_name << " number " << number << ": ";
		if(elem){
			std::cout << *elem;
		} else {
			std::cout << "none found";
		}

		std::cout << std::endl;
	}
	if(n_turns <= 0){
		if(verbose){
			std::cout << "No propagation calculation requested" << std::endl;
		}
		return;
	}

	std::cout << "Start    ps " << ps << std::endl;
	for(int turn = 0; turn < n_turns; ++turn){
		if(verbose && (turn > 0)){
			std::cout << "Turn " << turn << std::endl
				  << "        ps " << ps << std::endl;
		}

		for(auto cv : machine){
			auto& elem = dynamic_cast<tsc::ElemType&>(*cv);
			elem.pass(calc_config, ps);
		}
	}
	std::cout << "End      ps " << ps << std::endl;

}

int main(int argc, char *argv[])
{

	process_cmd_line(argc, argv);

	GLPSParser parse;
	std::string fname = vm["lattice_file"].as<std::string>();
	bool verbose = vm["verbose"].as<bool>();
	std::unique_ptr<Config> config;

	// needs to be called once and only once
	register_elements();

	if(verbose){
		std::cerr << "Parsing file " << fname << std::endl;
		std::cerr.flush();
	}
	try{
		config.reset(parse.parse_file(fname.c_str()));
	}catch(std::exception& e){
		std::cerr<<"Parse error: "<<e.what()<<"\n";
		return 1;
	}

	try{
		auto machine = tsc::Machine(*config);
		process(machine);
	}catch(std::exception& e){
		std::cerr<<"Machine configuration error: "<<e.what()<<"\n";
		return 1;
	}

}

/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
