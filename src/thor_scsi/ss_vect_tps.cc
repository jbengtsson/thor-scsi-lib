#include <gtpsa/ss_vect.h>
#include <tps/tps_type.h>

template<>
void gtpsa::ss_vect<tps>::show(std::ostream& strm, int level, bool with_endl) const
{
    for(size_t i= 0; i<this->state_space.size(); ++i){
	this->state_space[i].show(strm, level);
	strm << " ";
    }
    if(with_endl) {	strm << "\n";        }
}
