#include <thor_scsi/core/cell_void.h>
#include <boost/type_index.hpp>

namespace tsc = thor_scsi::core;
namespace bti = boost::typeindex;

tsc::CellVoid::CellVoid(const Config& conf)
	:name(conf.get<std::string>("name"))
	,index(0)
	,p_observe(nullptr)
	,p_conf(conf)
{}

tsc::CellVoid::CellVoid(CellVoid&& o):
	//p_conf(std::move(o.p_conf)),
	index(std::exchange(o.index, -1))
{
	// std::cerr << "Cell void move ctor on " <<  o << std::endl;
	this->p_conf = o.p_conf;
	this->p_observe = o.p_observe;
	o.p_observe = nullptr;
	this->name = o.name;
}

tsc::CellVoid::~CellVoid()
{
	if (p_observe){
		std::cerr << "Would delete p_observe" << std::endl;
	}
	// delete p_observe;
}


void tsc::CellVoid::show(std::ostream& strm, int level) const
{
	strm<<"cell "<<index<<": "<<this->type_name()<<"("<<this->name<<")";
}

std::string tsc::CellVoid::prettyClassname(void) const
{
	return bti::type_id_with_cvr<decltype(this)>().pretty_name();
}

std::string tsc::CellVoid::pstr(void) const
{
	std::stringstream strm;

	strm << "<" << this->prettyClassname()
	     << " @ " << (const void *)(this) << ">(";
	this->show(strm, 0);
	strm << ")";
	return strm.str();
}

std::string tsc::CellVoid::repr(void) const
{
	std::stringstream strm;

	strm << "<" << this->prettyClassname()
	     << " @ " << (const void *)(this) << ">(";
	this->show(strm, 10);
	strm << ")";
	return strm.str();
}
