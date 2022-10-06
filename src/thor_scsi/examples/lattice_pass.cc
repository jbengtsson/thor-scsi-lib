#include <boost/program_options.hpp>
#include <thor_scsi/std_machine/std_machine.h>
#include <thor_scsi/std_machine/accelerator.h>
#include <thor_scsi/elements/radiation_delegate.h>
#include <thor_scsi/elements/radiation_delegate.h>
#include <thor_scsi/elements/bending.h>
#include <thor_scsi/core/elements_basis.h>
//#include <tps/ss_vect.h>
#include <stdlib.h>
#include <ostream>
#include <memory>

namespace po = boost::program_options;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;
namespace ts = thor_scsi;

static int log_level = THOR_SCSI_WARN;

po::variables_map vm;
static void process_cmd_line(int argc, char *argv[])
{

	// Declare the supported options.
	std::string progname(argv[0]);

	po::options_description desc(progname + ":allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("lattice_file,L", po::value<std::string>(),                 "lattice file name")
		("x_pos,x",      po::value<double>()->default_value(0e0), "horizontal position        (left hand coordinate system)")
		("px",           po::value<double>()->default_value(0e0), "horizontal direction       (left hand coordinate system)")
		("y_pos,y",      po::value<double>()->default_value(0e0), "vertical position          (left hand coordinate system)")
		("py",           po::value<double>()->default_value(0e0), "vertical direction         (left hand coordinate system)")
		("ct",           po::value<double>()->default_value(0e0), "relative phase deviation   (left hand coordinate system)")
		("delta",        po::value<double>()->default_value(0e0), "relative impulse deviation (left hand coordinate system)")
		("inspect,i",    po::value<std::string>()->default_value(""), "name of element to inspect: if option verbose is used extra information will be printed")
		("radiate,r",    po::value<bool>()->default_value(false), "enable radiation togther with calculation")
		("energy,E",     po::value<double>()->default_value(2.5e9), "energy")
		("number,n",     po::value<int>()->default_value(0), "number of element to inspect")
		("n_turns",      po::value<int>()->default_value(0),      "propagate n turns (set to zero for none)" )
		("transport_matrix", po::value<bool>()->default_value(false), "compute transport matrix")
		("start_element_number,s", po::value<int>()->default_value(0), "first element to use")
		("end_element_number,e",   po::value<int>()->default_value(-1), "last element to use, (-1) for last element of lattice")
		("dump_lattice",           po::value<bool>()->default_value(false), "dump read in lattice (to stdout)" )
		("log_level,l",            po::value<std::string>()->default_value("WARNING"), "Logging level: FINE|DEBUG|INFO|WARN|ERROR")
		("verbose,v",              po::value<bool>()->default_value(false), "verbose output" )
		("very_verbose,V",         po::value<bool>()->default_value(false), "even more verbose output" )
		;

	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	const bool verbose = vm["verbose"].as<bool>();
	if (!vm.count("lattice_file")) {
		std::cout << "A lattice file must be defined" << std::endl;
		std::cout << desc << std::endl;
		exit(1);
	} else 	if (vm.count("help")) {
		std::cout << desc << std::endl;
		exit(1);
	}
	if(verbose){
		std::cout << "Finished processing "  << std::endl;
	}
}



static void compute_transport_matrix(ts::Accelerator& accelerator, const int first_element, const int last_element)
{

	bool verbose = vm["verbose"].as<bool>();
	bool very_verbose = vm["very_verbose"].as<bool>();
	auto desc = gtpsa::desc(1, 6);
	auto a_tps = gtpsa::tpsa(desc, mad_tpsa_default);
	gtpsa::ss_vect<gtpsa::tpsa> ps(a_tps);
	ps.set_identity();

	tsc::ConfigType calc_config;

	std::cout << "Starting poincare map computation with " << std::endl
		  << ps << std::endl;
	if(verbose){
		std::cout << "Processing elements: ";
	}

	std::shared_ptr<tsc::ElemType> elem;
	for(int n_elem = first_element; n_elem <= last_element; ++n_elem){
		if(very_verbose){
			std::cout << "Processing element number " << n_elem << std::endl;
			std::cout.flush();
		}
		auto cv = accelerator.at(n_elem);
		// auto cv = accelerator[n_elem];
		if(!cv){
			throw std::runtime_error("no cv ...");
		}
		if(very_verbose){
			std::cout << "cv[" << cv->index <<"] '"<< cv->name << "', # counts " << cv.use_count() << std::endl;
		}
		if(!very_verbose && verbose){
			std::cout << cv->name << " ";
		}
		elem = std::dynamic_pointer_cast<tsc::ElemType>(cv);
		if(elem){
			if(very_verbose){
				std::cout << "elem '" << elem->name << "', # counts " << cv.use_count() << std::endl;
			}
			elem->propagate(calc_config, ps);
		} else {
			std::cerr << "Element " << n_elem << "could not be cast to Elemtype"
				  << " element: " << cv << std::endl;
			return;
		}

	}
	if(verbose){
		std::cout << std::endl;
	}
	if(elem && verbose){
		std::cout << "Last element ";
		elem->show(std::cout, 10);
		std::cout << std::endl;

	}
	std::cout << "Computed poincare map " << std::endl
		  << ps << std::endl;

}


static void add_radiation_delegates(ts::Accelerator& accelerator)
{

	double energy = vm["energy"].as<double>();

	for(auto cv : accelerator){
		auto fk = std::dynamic_pointer_cast<tse::FieldKick>(cv);
		if(fk){
			auto p_del = std::make_shared<tse::RadiationDelegateKick>();
			p_del->setEnergy(energy);
			fk->setRadiationDelegate(p_del);
		}
	}
}
static void compute_transport_matrix_prop(ts::Accelerator& accelerator, const int first_element, const int last_element)
{
	bool verbose = vm["verbose"].as<bool>();
	bool very_verbose = vm["very_verbose"].as<bool>();
	bool radiate = vm["radiate"].as<bool>();

	auto desc = gtpsa::desc(1, 6);
	auto a_tps = gtpsa::tpsa(desc, mad_tpsa_default);
	gtpsa::ss_vect<gtpsa::tpsa> ps(a_tps);
	ps.set_identity();

	tsc::ConfigType calc_config;

	std::cout << "Starting poincare map computation with " << std::endl
		  << ps << std::endl;
	if(verbose){
		std::cout << "Processing elements: ";
	}
	if(radiate){
		calc_config.radiation = true;
		add_radiation_delegates(accelerator);

	} else {
		calc_config.radiation = false;
	}

	accelerator.propagate(calc_config, ps, first_element, last_element);
	std::cout << "Computed poincare map " << std::endl
		  << ps << std::endl;

}

static void compute_transport_matrix_dbl(ts::Accelerator& accelerator, const int first_element, const int last_element)
{

	bool verbose = vm["verbose"].as<bool>();
	gtpsa::ss_vect<double> ps(0e0);
	ps.set_zero();
	ps[x_]  = 1e-3;
	ps[px_] = 2e-3;
	ps[y_]  = 3e-3;
	ps[py_] = 4e-3;

	tsc::ConfigType calc_config;

	std::cout << "Starting poincare map computation with " << std::endl
		  << ps << std::endl;
	for(int n_elem = first_element; n_elem <= last_element; ++n_elem){
		auto cv = accelerator[n_elem];
		if(verbose){
			std::cout << "Processing element number " << n_elem
				  << " element ";
			cv->show(std::cout, 10);
			std::cout << std::endl;
		}
		auto elem = std::dynamic_pointer_cast<tsc::ElemType>(cv);
		if(elem){
			elem->propagate(calc_config, ps);
		} else {
			std::cerr << "Element " << n_elem << "could not be cast to Elemtype"
				  << " element: " << cv << std::endl;
			return;
		}

	}
	std::cout << "Computed poincare map " << std::endl
		  << ps << std::endl;

}


static void user_compute_transport_matrix(ts::Accelerator& accelerator)

{
	bool verbose = vm["verbose"].as<bool>();
	if (!vm["transport_matrix"].as<bool>()) {
		if(verbose){
			std::cout << "no one requested to compute a transport matrix" << std::endl;
		}
		return;
	}
	if(verbose){
		std::cout << "computing transport matrix" << std::endl;
	}
	const int first_element = vm["start_element_number"].as<int>();
	int last_element = vm["end_element_number"].as<int>();

	if(verbose){
		std::cout << "user requested " << first_element
			  << " to "<< last_element  << std::endl;
	}
	if(last_element == -1){
		last_element = accelerator.size();
		if(last_element >0){
			last_element--;
		}
	}
	if(verbose){
		std::cout << "Computing transport matrix from " << first_element
			  << " to " << last_element << " number " << std::endl;
	}
	// compute_transport_matrix(accelerator, first_element, last_element);


	accelerator.set_trace(&std::cout);
	compute_transport_matrix_prop(accelerator, first_element, last_element);
}

static void track_n_turns(ts::Accelerator& accelerator)
{
	const int n_turns = vm["n_turns"].as<int>();
	bool verbose = vm["verbose"].as<bool>();
	bool radiate = vm["radiate"].as<bool>();

	if (n_turns <= 0) {
		if (verbose) {
			std::cout << "No propagation calculation requested" << std::endl;
		}
		return;
	}

	gtpsa::ss_vect<double> ps(0e0);
	ps.set_zero();
	ps[x_]     = vm["x_pos"].as<double>();
	ps[px_]    = vm["px"].as<double>();
	ps[y_]     = vm["y_pos"].as<double>();
	ps[py_]    = vm["py"].as<double>();
	ps[ct_]    = vm["ct"].as<double>();
	ps[delta_] = vm["delta"].as<double>();

	tsc::ConfigType calc_config;

	if(radiate){
		calc_config.radiation = true;
		add_radiation_delegates(accelerator);
	} else {
		calc_config.radiation = false;
	}
	std::cout << "Start    ps " << ps << std::endl;
	for(int turn = 0; turn < n_turns; ++turn){
		if(verbose && (turn > 0)){
			std::cout << "Turn " << turn << std::endl
				  << "        ps " << ps << std::endl;
		}

		for(auto cv : accelerator){
			auto& elem = dynamic_cast<tsc::ElemType&>(*cv);
			elem.propagate(calc_config, ps);
		}
	}
	std::cout << "End      ps " << ps << std::endl;
}

static void process(ts::Accelerator& accelerator)
{
	bool verbose = vm["verbose"].as<bool>();

	accelerator.set_log_level(log_level);

	if(vm["dump_lattice"].as<bool>()){
		std::cout << "Accelerator configuration " << std::endl;
		for(auto cv : accelerator){
			auto& elem = dynamic_cast<tsc::ElemType&>(*cv);
			std::cout << elem << std::endl;
		}
	}

	auto element_name = vm["inspect"].as<std::string>();
	if(element_name != ""){
		int number =  vm["number"].as<int>();
		auto elem = accelerator.find(element_name, number);
		std::cout << "Accelerator element " << element_name << " number " << number << ": ";
		if(elem){
                        if(verbose){
				elem->show(std::cout, 12);
                        }else{
				std::cout << *elem;
                        }
		} else {
			std::cout << "none found";
		}

		std::cout << std::endl;
	}
	track_n_turns(accelerator);
	user_compute_transport_matrix(accelerator);

}

static void check_log_level()
{

	struct levels
	{
		const char* name;
		int level;
	};


	std::vector<struct levels> t_levels   = {
		{"FINE",    THOR_SCSI_FINE},
		{"DEBUG",   THOR_SCSI_DEBUG},
		{"INFO",    THOR_SCSI_INFO},
		{"WARN",    THOR_SCSI_WARN},
		{"WARNING", THOR_SCSI_WARNING},
		{"ERROR",   THOR_SCSI_ERROR}
	};

	bool verbose = vm["verbose"].as<bool>();
	std::string log_level_str = vm["log_level"].as<std::string>();
	if(verbose){
		std::cerr << "Checking log level " << log_level_str << "\n";
	}

	for(auto lvl : t_levels){

		if (lvl.name == log_level_str){
			log_level = lvl.level;
			if(verbose){
				std::cerr << "Found log leven " << lvl.name << "\n" ;
			}
			return;
		}

	}
	std::ostringstream strm;
	strm << "Unknown log level " << log_level_str << "\n";
	throw std::runtime_error(strm.str());

}
int main(int argc, char *argv[])
{

	process_cmd_line(argc, argv);
	check_log_level();

	GLPSParser parse;
	std::string fname = vm["lattice_file"].as<std::string>();
	bool verbose = vm["verbose"].as<bool>();
	std::unique_ptr<Config> config;

	// needs to be called once and only once
	ts::register_elements();

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
		auto accelerator = ts::Accelerator(*config);
		process(accelerator);
	}catch(std::exception& e){
		std::cerr<<"Accelerator configuration error: "<<e.what()<<"\n";
		return 1;
	}

}

/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
