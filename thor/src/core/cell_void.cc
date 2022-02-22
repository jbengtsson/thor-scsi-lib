#include <thor_scsi/core/cell_void.h>

namespace tsc=thor_scsi::core;

tsc::CellVoid::CellVoid(const Config& conf)
    :name(conf.get<std::string>("name"))
    ,index(0)
     /* ,length(conf.get<double>("L",0.0)) */
    ,p_observe(NULL)
    ,p_conf(conf)
{}

tsc::CellVoid::~CellVoid()
{
    delete p_observe;
}

void tsc::CellVoid::show(std::ostream& strm, int level) const
{
    strm<<"Cell "<<index<<": "<<name<<" ("<<type_name()<<")\n";
}
