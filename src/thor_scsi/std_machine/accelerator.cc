#include <thor_scsi/std_machine/accelerator.h>
#include <thor_scsi/elements/standard_aperture.h>
#include <thor_scsi/elements/elements_enums.h>
#include <thor_scsi/core/exceptions.h>

namespace ts = thor_scsi;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;


static tsc::p_elements_t
vec_elem_type_to_cell_void(const std::vector<std::shared_ptr<thor_scsi::core::ElemType>>&  elements)
{
       tsc::p_elements_t cells;
       cells.reserve(elements.size());
       /* order needs to be preserved */
       for(size_t i=0; i < elements.size(); ++i){
	   cells.push_back(std::dynamic_pointer_cast<tsc::CellVoid>(elements[i]));
       }
       return cells;

}
ts::Accelerator::Accelerator(const Config & conf)
    :tsc::Machine(conf)
{}

ts::Accelerator::Accelerator(const std::vector<std::shared_ptr<thor_scsi::core::ElemType>>&  elements)
    :tsc::Machine(vec_elem_type_to_cell_void(elements))
{

#if 0
    // Waiting that machine is reogranised ...
    this->updateElementList(cells);
#endif
}


template<typename T>
int
ts::Accelerator::_propagate(thor_scsi::core::ConfigType& conf, ss_vect<T> &ps, size_t start, int max_elements, size_t n_turns)// const
{

	const int nelem = static_cast<int>(this->size());

	bool retreat = std::signbit(max_elements);
	int next_elem = static_cast<int>(start);

	for(size_t turn=0; turn<n_turns; ++turn) {
#if 1
	    next_elem = static_cast<int>(start-1);
	    for(int i=start-1; next_elem >= 0 && next_elem<nelem && i<std::abs(max_elements); i++)
#else
	    next_elem = 0;
	    for(int i=0; next_elem >= 0 && next_elem<nelem && i<std::abs(max_elements); i++)
#endif
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
			// Should raise an exception!
			// THOR_SCSI_LOG(ERROR) << "Failed to cast to element to Elemtype " << cv->name << "\n";
			std::runtime_error("Could not cast cell void to elemtype");
			return next_elem;
		}

		std::shared_ptr<tsc::Observer> observer = elem->observer();
		if(observer){
			observer->view(std::const_pointer_cast<tsc::ElemType>(elem), ps, tsc::ObservedState::start, 0);
		}
		elem->pass(conf, ps);
		if(observer){
			observer->view(elem, ps, tsc::ObservedState::end, 0);
		}
		auto aperture = elem->getAperture();
		if(aperture){
			/* could be any aperture ... */
			bool flag = elem->checkAmplitude(ps);
			if(not flag){
				THOR_SCSI_LOG(INFO) << "Element lost at " << elem.get()
						    <<" with aperture " << aperture.get();
				/* get a guess on the plane */
				auto rect_apt = std::dynamic_pointer_cast<tse::RectangularAperture>(aperture);
				if(rect_apt){
					double x=-1e0, y=-1e0;
					if(x < 0e0){
						conf.lossplane = tse::PlaneKind::Horizontal;
					} else if(y < 0e0){
						conf.lossplane = tse::PlaneKind::Vertical;
					} else {
						throw ts::SanityCheckError("Particle lost but no loss plane identified!");
					}
				}
				return next_elem;
			}
		}
		auto trace = this->trace();
		if(trace)
			(*trace) << "After ["<< n<< "] " << cv->name << " " <<std::endl << ps << std::endl;
	    }
	}
	return next_elem;
}


int
ts::Accelerator::propagate(thor_scsi::core::ConfigType& conf, ss_vect_dbl &ps, size_t start, int max_elements, size_t n_turns)// const
{
    return _propagate(conf, ps, start, max_elements, n_turns);
}

int
ts::Accelerator::propagate(thor_scsi::core::ConfigType& conf, ss_vect_tps  &ps, size_t start, int max_elements, size_t n_turns)// const
{
    return _propagate(conf, ps, start, max_elements, n_turns);
}
