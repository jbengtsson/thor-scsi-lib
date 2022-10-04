#include <thor_scsi/core/machine.h>
#include <flame/core/util.h>
#include <boost/thread/mutex.hpp>



namespace {
// This mutex guards the global Machine::p_state_infos
typedef boost::mutex info_mutex_t;
info_mutex_t info_mutex;
}

namespace tsc=thor_scsi::core;


void tsc::Machine::updateElementList(p_elements_t& elements)
{

    p_lookup_t lookup_name, lookup_type;

    // keep sequence order ...
    for(size_t index = 0; index < elements.size(); ++index){
	auto& e_ptr = elements[index];
	e_ptr->index = index;
    }
    for(auto& E: elements){
	lookup_name.insert(std::make_pair(LookupKey(E->name, E->index), E));
	lookup_type.insert(std::make_pair(LookupKey(E->type_name(), E->index), E));
    }
    p_elements.swap(elements);
    p_lookup.swap(lookup_name);
    p_lookup_type.swap(lookup_type);

}
/*
 * Original code:
 *    1. create state (or sim type)
 *    2. add all elements to it ....
 *
 *
 *  Now switching over: only simulation state
 */
tsc::Machine::Machine(const p_elements_t& elements)
    :p_elements()
    ,p_trace(nullptr)
{
    p_elements_t a_cpy(elements.size());
    std::copy(elements.begin(), elements.end(), a_cpy.begin());
    this->updateElementList(a_cpy);
}

tsc::Machine::Machine(p_elements_t& elements)
    :p_elements()
    ,p_trace(nullptr)
{
    this->updateElementList(elements);
}

tsc::Machine::Machine(const Config& c)
    :p_elements()
    ,p_trace(nullptr)
    ,p_conf(c)
     //,p_info()
{

#if 0  /* NO state */
    /* originally: first finding sim state */
    std::string type(c.get<std::string>("sim_type"));
    p_simtype = type;
    THOR_SCSI_LOG(INFO)<<"Constructing Machine w/ sim_type='"<<type<<'\'';
#endif /* NO state */

    info_mutex_t::scoped_lock G(info_mutex);

#if 0  /* NO state */
    p_element_infos_t::iterator it = p_element_infos.find(type);
    if(it==p_element_infos.end()) {
        std::ostringstream msg;
        msg<<"Unsupport sim_type '"<<type<<"'";
        throw key_error(msg.str());
    }

    p_info = it->second;
#endif /* NO state */

    typedef Config::vector_t elements_t;
    // std::cerr << "Processing config " << c << std::endl;
    elements_t Es(c.get<elements_t>("elements"));

    p_elements_t result;
    result.reserve(Es.size());

    size_t idx=0;
    for(elements_t::iterator it=Es.begin(), end=Es.end(); it!=end; ++it)
    {
        const Config& EC = *it;

	// Todo: I guess that the lattice will return each different element with
	//        a name returned as type
        const std::string& etype(EC.get<std::string>("type"));

        p_element_infos_t::iterator eit = p_element_infos.find(etype);
        if(eit==p_element_infos.end()){
	    std::cerr << "Could not find an element for type " << etype << std::endl;
            throw key_error(etype);
	}

        element_builder_t* builder = eit->second.builder;

	std::shared_ptr<CellVoid> E;
        try{
            auto tmp  = (builder->build(EC));
	    E = tmp;
        }catch(key_error& e){
            std::ostringstream strm;
            strm<<"Error while initializing element "<<idx<<" '"<<EC.get<std::string>("name", "<invalid>")
               <<"' : missing required parameter '"<<e.what()<<"'";
            throw key_error(strm.str());

        }catch(std::exception& e){
            std::ostringstream strm;
            strm<<"Error while constructing element "<<idx<<" '"<<EC.get<std::string>("name", "<invalid>")
               <<"' : "<<typeid(e).name()<<" : "<<e.what();
            throw std::runtime_error(strm.str());
        }

        if(E->type_name()!=etype) {
            std::ostringstream strm;
            strm<<"Element type inconsistent "<<etype<<" "<<E->type_name();
            throw std::logic_error(strm.str());
        }

        *const_cast<size_t*>(&E->index) = idx++; // ugly

        //result.push_back(std::make_shared<tsc::CellVoid>(E));
	result.push_back(E);
    }

    G.unlock();

    this->updateElementList(result);
    THOR_SCSI_LOG(DEBUG) << "\nLattice:";
}

//! Elements with a given name
tsc::p_elements_t
tsc::Machine::elementsWithName(const std::string& name)
{
	auto iters = this->equal_range(name);
	tsc::p_elements_t sel_elems(std::get<0>(iters), std::get<1>(iters));
	return sel_elems;
}

//! Elements with a given name type
tsc::p_elements_t
tsc::Machine::elementsWithNameType(const std::string& name)
{
	auto iters = this->equal_range_type(name);
	tsc::p_elements_t sel_elems(std::get<0>(iters), std::get<1>(iters));
	return sel_elems;
}

tsc::Machine::~Machine()
{
    return;
    // Should now be cleared up automatically as shared pointers are used.
    /*
    for(p_elements_t::iterator it=p_elements.begin(), end=p_elements.end(); it!=end; ++it)
    {
        delete *it;
    }
    */
}


#if 0

StateBase*
tsc::Machine::allocState(const Config &c) const
{
    return (*p_info.builder)(c);
}
#endif

void tsc::Machine::reconfigure(size_t idx, const Config& c)
{
    if(idx>=p_elements.size())
        throw std::invalid_argument("element index out of range");

    const std::string& etype(c.get<std::string>("type"));

    auto eit = p_element_infos.find(etype);
    if(eit==p_element_infos.end())
        throw key_error(etype);

    element_builder_t *builder = eit->second.builder;

    builder->rebuild(p_elements[idx], c, idx);
}

tsc::Machine::p_element_infos_t tsc::Machine::p_element_infos;

#if 0  /* NO state */
void tsc::Machine::p_registerState(const char *name, state_builder_t b)
{
    info_mutex_t::scoped_lock G(info_mutex);
    if(p_element_infos.find(name)!=p_element_infos.end()) {
        std::ostringstream strm;
        strm<<"attempt to register already registered sim_type=\""<<name<<"\"";
        throw std::logic_error(strm.str());
    }
    element_info I;
    I.name = name;
    I.builder = b;
    p_element_infos[name] = I;
}
#endif /* NO state */

void tsc::Machine::p_registerElement(
				     //const std::string& sname,
				     const char *ename, element_builder_t *b)
{
    info_mutex_t::scoped_lock G(info_mutex);

    // originally iterating first over known states than over known elements
    // Interstates are less mobile than voters
    // thus if I add state I rather iterate first over the elements and then
    // over the state

    p_element_infos_t::iterator I = p_element_infos.find(ename);
#if 0 /* NO state */
    if(it==p_element_infos.end()) {
        std::ostringstream strm;
        strm<<"can't add element \""<<ename<<"\" for unknown sim_type=\""<<sname<<"\"";
        throw std::logic_error(strm.str());
    }
    element_info& I = it->second;
#endif

    if(I!=p_element_infos.end()) {
        std::ostringstream strm;
        strm<<"element type \""<<ename<<"\" has already been registered for "
	    <<"this machine";
	//  /* NO state */ "sim_type=\""<<sname<<"\"";
        throw std::logic_error(strm.str());
    }
    p_element_infos[ename].name = std::string(ename);
    p_element_infos[ename].builder = b;
    // p_element_infos[ename].namep_element_infos[ename].namea= b;
}

void tsc::Machine::registeryCleanup()
{
    info_mutex_t::scoped_lock G(info_mutex);

    for(p_element_infos_t::iterator it=p_element_infos.begin(), end=p_element_infos.end();
        it!=end; ++it)
    {
      /*
       * does not exist any more ... only single elemets
        element_info::elements_t::iterator it2, end2;
        for(it2=it->second.elements.begin(), end2=it->second.elements.end(); it2!=end2; ++it2)
        {
            delete it2->second;
        }
      */

    }
    p_element_infos.clear();
}

namespace thor_scsi::core {
    std::ostream& operator<<(std::ostream& strm, const thor_scsi::core::Machine& m)
    {
        // Not handling sim type explict any more
        // strm<<"sim_type: "<<m.p_info.name<<"\n
         strm<<"#Elements: "<<m.p_elements.size()<<"\n";
	 for(auto it=m.p_elements.begin(),
	      end=m.p_elements.end(); it!=end; ++it)
	{
	    (*it)->show(strm, 0);
	}
	return strm;
    }
}
namespace {
struct Logcerr : public tsc::Machine::Logger
{
    virtual ~Logcerr() {}
    virtual void log(const tsc::Machine::LogRecord &r) override final
    {
        std::string msg(r.strm.str());
        std::cerr<<r.fname<<':'<<r.lnum<<" : "<<msg;
        if(msg.empty() || msg[msg.size()-1]!='\n')
            std::cerr.put('\n');
    }
    static Logcerr singleton;
    static void noopdtor(Logcerr*) {}
};
Logcerr Logcerr::singleton;
}

#if 1
int tsc::Machine::log_detail = THOR_SCSI_DEBUG;
#else
int tsc::Machine::log_detail = THOR_SCSI_WARN;
#endif
std::shared_ptr<tsc::Machine::Logger> tsc::Machine::p_logger(&Logcerr::singleton, Logcerr::noopdtor);

void tsc::Machine::set_logger(const std::shared_ptr<Logger> &p)
{
    std::shared_ptr<Logger> temp(p);
    if(!temp)
        temp.reset(&Logcerr::singleton, Logcerr::noopdtor);
    {
        info_mutex_t::scoped_lock G(info_mutex);
        p_logger.swap(temp);
    }
}

void tsc::Machine::set_log_level(int lvl)
{
    switch(lvl) {
    case THOR_SCSI_FINE:
    case THOR_SCSI_DEBUG:
    case THOR_SCSI_INFO:
    case THOR_SCSI_WARN:
    case THOR_SCSI_ERROR:
	break;

    default:
	std::ostringstream strm;
	strm << "Log level " << lvl << " unkown log level: known levels "
	     << "fine = "   << THOR_SCSI_FINE
	     << " debug = " << THOR_SCSI_DEBUG
	     << " info = "  << THOR_SCSI_INFO
	     << " warn = "  << THOR_SCSI_WARN
	     << " error ="  << THOR_SCSI_ERROR
	     << ".";
	throw std::runtime_error(strm.str());
    }

    log_detail = lvl;
}

std::shared_ptr<tsc::Machine::Logger> tsc::Machine::get_logger(void)
{
    return p_logger;
}

tsc::Machine::LogRecord::~LogRecord()
{
    std::shared_ptr<Logger> logger;
    {
        info_mutex_t::scoped_lock G(info_mutex);
        logger = p_logger;
    }
    logger->log(*this);
}
