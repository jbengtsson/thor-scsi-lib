#include <thor_scsi/core/cell_void.h>

namespace tsc=thor_scsi::core;

tsc::CellVoid::CellVoid(const Config& conf)
    :name(conf.get<std::string>("name"))
    ,index(0)
    ,p_observe(nullptr)
    ,p_conf(conf)
{}

tsc::CellVoid::~CellVoid()
{
    delete p_observe;
}

void tsc::CellVoid::show(std::ostream& strm, int level) const
{
    strm<<"cell "<<index<<": "<<name<<" ("<<type_name()<<")";
}
