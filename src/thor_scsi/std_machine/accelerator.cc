#include <thor_scsi/std_machine/accelerator.h>

namespace ts = thor_scsi;
namespace tsc = thor_scsi::core;

ts::Accelerator::Accelerator(const Config & conf) :
  tsc::Machine(conf)
{}


template<typename T>
void
ts::Accelerator::_propagate(thor_scsi::core::ConfigType& conf, ss_vect<T> &ps, size_t start, int max)// const
{

	const int nelem = static_cast<int>(this->size());

	int next_elem = static_cast<int>(start);
	bool retreat = std::signbit(max);

	for(int i=start; next_elem >= 0 && next_elem<nelem && i<abs(max); i++)
	{
		size_t n = next_elem;
		std::shared_ptr<tsc::CellVoid> cv = this->at(n);

		if(retreat) {
			next_elem--;
		} else {
			next_elem++;
		}
		auto elem = std::dynamic_pointer_cast<tsc::ElemType>(cv);
		if(!elem){
			std::cerr << "Failed to cast to element to Elemtype " << cv->name << std::endl;
			return;
		}

		std::shared_ptr<tsc::Observer> observer = elem->observer();
		if(observer){
			observer->view(std::const_pointer_cast<tsc::ElemType>(elem), ps, tsc::ObservedState::start, 0);
		}
		elem->pass(conf, ps);
		if(observer){
			observer->view(elem, ps, tsc::ObservedState::end, 0);
		}
		auto trace = this->trace();
		if(trace)
			(*trace) << "After ["<< n<< "] " << cv->name << " " <<std::endl << ps << std::endl;
	}
}


void
ts::Accelerator::propagate(thor_scsi::core::ConfigType& conf, ss_vect_dbl &ps, size_t start, int max)// const
{
    _propagate(conf, ps, start, max);
}
void
ts::Accelerator::propagate(thor_scsi::core::ConfigType& conf, ss_vect_tps  &ps, size_t start, int max)// const
{
    _propagate(conf, ps, start, max);
}
