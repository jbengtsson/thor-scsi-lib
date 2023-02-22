#include <thor_scsi/std_machine/accelerator.h>
#include <thor_scsi/elements/standard_aperture.h>
#include <thor_scsi/elements/elements_enums.h>
#include <thor_scsi/elements/marker.h>
#include <thor_scsi/core/exceptions.h>
#include <sstream>

namespace ts = thor_scsi;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;


static tsc::p_elements_t
vec_elem_type_to_cell_void
(const std::vector<std::shared_ptr<thor_scsi::core::ElemType>>& elements)
{
  tsc::p_elements_t cells;
  cells.reserve(elements.size());
  /* order needs to be preserved */
  std::transform(elements.begin(), elements.end(), std::back_inserter(cells),
		 [](std::shared_ptr<thor_scsi::core::ElemType> elem)
		     { return std::dynamic_pointer_cast<tsc::CellVoid>(elem) ; }
      );
  return cells;

}

ts::Accelerator::Accelerator(const Config & conf, bool marker_at_start)
  :tsc::Machine(conf)
{
    if(marker_at_start){
	this->addMarkerAtStartIfRequired();
    }

}

ts::Accelerator::Accelerator
(const std::vector<std::shared_ptr<thor_scsi::core::ElemType>>& elements, bool marker_at_start)
  :tsc::Machine(vec_elem_type_to_cell_void(elements))
{
    if(marker_at_start){
	this->addMarkerAtStartIfRequired();
    }
}

void ts::Accelerator::addMarkerAtStart(void)
{
    Config C;
    C.setAny("name", "Start");
    C.setAny("Length", 0e0);

    /* perhaps declare this>elements it protected ? */
    tsc::p_elements_t elements;
    elements.reserve(this->size() + 1);
    elements.insert(elements.begin(), std::make_shared<tse::MarkerType>(C));
    elements.insert(std::next(elements.begin(), 1), this->begin(), this->end());
    this->updateElementList(elements);
}

void ts::Accelerator::addMarkerAtStartIfRequired(void)
{
    auto cv = this->at(0);
    auto elem = std::dynamic_pointer_cast<tsc::ElemType>(cv);
    auto marker =  std::dynamic_pointer_cast<tse::MarkerType>(elem);

    if( (marker) && (elem->getLength() == 0e0) ){
	THOR_SCSI_LOG(INFO) << "Not inserting a 'Marker' at the beginning as already one is there.\n";
	return;
    }

    if( (!marker) && (elem->getLength() == 0e0) ) {
	THOR_SCSI_LOG(WARNING) << "Inserting a 'Marker' at the beginning even if"
			       << " first element is of length zero. \n";
	this->addMarkerAtStart();
	return;
    }
    THOR_SCSI_LOG(WARNING) << "Inserting a 'Marker' at the beginning as requested:"
			   << " first element is neither 'Marker' nor of zero length. \n";
    this->addMarkerAtStart();
}

template<typename T>
int
ts::Accelerator::_propagate(thor_scsi::core::ConfigType& conf, gtpsa::ss_vect<T> &ps, size_t start_elem, int max_elements, size_t n_turns,  bool tracy_compatible_indexing)// const
{

	/* I guess Tobin would complain about this extra complexity */
	int nelem = static_cast<int>(this->size());
	bool retreat = std::signbit(max_elements);

	if(tracy_compatible_indexing){
		if(start_elem>0) {
			// take the one off
			start_elem -= 1;
		} else {
			std:: stringstream strm;
			strm << "Requested tracy compatible indexing, but start element "
			     << start_elem << "was smaller or equal to zero";
			throw std::runtime_error(strm.str());
		}
	}

	int next_elem = static_cast<int>(start_elem);
	auto trace = this->trace();


	for(size_t turn=0; turn<n_turns; ++turn) {
	    if(trace)
		(*trace) << "turn " <<  turn << ", processing elements from " <<  start_elem
			 << " for a maximum of elements " << max_elements << std::endl;

	    next_elem = static_cast<int>(start_elem);
	    for(int i=0; next_elem >= 0 && next_elem<nelem && i<std::abs(max_elements); i++)
	    {
		size_t n = next_elem;

		std::shared_ptr<tsc::CellVoid> cv = this->at(n);
		auto elem = std::dynamic_pointer_cast<tsc::ElemType>(cv);
		if(!elem){
		    // Should raise an exception!
		    THOR_SCSI_LOG(ERROR)
			<< "Failed to cast to element to Elemtype " << cv->name << "\n";
		    std::runtime_error("Could not cast cell void to elemtype");
		    return next_elem;
		}
		if(retreat) {
		    next_elem--;
		} else {
		    next_elem++;
		}
		std::shared_ptr<tsc::Observer> observer = elem->observer();
		if(observer){
			observer->view(std::const_pointer_cast<tsc::ElemType>(elem), ps, tsc::ObservedState::start, 0);
		}
		elem->propagate(conf, ps);
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
		if(trace)
			(*trace) << "After ["<< n<< "] " << cv->name << " " <<std::endl << ps << std::endl;
	    }
	}
	return next_elem;
}


int
ts::Accelerator::
propagate(thor_scsi::core::ConfigType& conf, ss_vect_dbl  &ps, size_t start,
	  int max_elements, size_t n_turns, bool tracy_compatible_indexing)// const
{
    return _propagate(conf, ps, start, max_elements, n_turns, tracy_compatible_indexing);
}

/*
int
ts::Accelerator::
propagate(thor_scsi::core::ConfigType& conf, ss_vect_tps  &ps, size_t start,
	  int max_elements, size_t n_turns, bool tracy_compatible_indexing)// const
{
    return _propagate(conf, ps, start, max_elements, n_turns, tracy_compatible_indexing);
}
*/
int
ts::Accelerator::
propagate(thor_scsi::core::ConfigType& conf, ss_vect_tpsa &ps, size_t start,
	  int max_elements, size_t n_turns,  bool tracy_compatible_indexing)// const
{
    return _propagate(conf, ps, start, max_elements, n_turns, tracy_compatible_indexing);
}
