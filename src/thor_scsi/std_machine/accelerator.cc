#include <thor_scsi/std_machine/accelerator.h>
#include <thor_scsi/elements/standard_aperture.h>
#include <thor_scsi/elements/elements_enums.h>
#include <thor_scsi/elements/marker.h>
#include <thor_scsi/elements/quadrupole.h>
#include <thor_scsi/core/exceptions.h>
#include <sstream>
#include <thor_scsi/core/multipole_types.h>

namespace ts = thor_scsi;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;


//template<class C>
static tsc::p_elements_t
vec_elem_type_to_cell_void
(const std::vector<std::shared_ptr<thor_scsi::core::ElemTypeKnobbed/*<C>*/>>&
 elements)
{
  tsc::p_elements_t cells;
  cells.reserve(elements.size());
  /* order needs to be preserved */
  std::transform(elements.begin(), elements.end(), std::back_inserter(cells),
		 [](std::shared_ptr<thor_scsi::core::ElemTypeKnobbed/*<C>*/>
		    elem)
		 { return std::dynamic_pointer_cast<tsc::CellVoid>(elem) ; }
		 );
  return cells;

}

template<class C>
ts::AcceleratorKnobbable<C>::AcceleratorKnobbable
(const Config & conf, bool marker_at_start) :tsc::Machine(conf)
{
  if(marker_at_start){
    this->addMarkerAtStartIfRequired();
  }

}


template<class C>
ts::AcceleratorKnobbable<C>::AcceleratorKnobbable
(const std::vector<std::shared_ptr<thor_scsi::core::ElemTypeKnobbed>>& elements,
 bool marker_at_start)
  :tsc::Machine(vec_elem_type_to_cell_void(elements))
{
  if(marker_at_start){
    this->addMarkerAtStartIfRequired();
  }
}

template<class C>
ts::AcceleratorKnobbable<C>::AcceleratorKnobbable
(std::vector<std::shared_ptr<thor_scsi::core::CellVoid>>& elements,
 bool marker_at_start)
  : tsc::Machine(elements)
{
  if(marker_at_start){
    this->addMarkerAtStartIfRequired();
  }
}

template<class CTypes>
void ts::AcceleratorKnobbable<CTypes>::addMarkerAtStart(void)
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

template<class C>
void ts::AcceleratorKnobbable<C>::addMarkerAtStartIfRequired(void)
{
  auto cv = this->at(0);
  auto elem = std::dynamic_pointer_cast<tsc::ElemType>(cv);
  auto marker =  std::dynamic_pointer_cast<tse::MarkerType>(elem);

  if( (marker) && (elem->getLength() == 0e0) ){
    THOR_SCSI_LOG(INFO)
      << "Not inserting a 'Marker' at the beginning as already one is there.\n";
    return;
  }

  if( (!marker) && (elem->getLength() == 0e0) ) {
    THOR_SCSI_LOG(WARNING) << "Inserting a 'Marker' at the beginning even if"
			   << " first element is of length zero. \n";
    this->addMarkerAtStart();
    return;
  }
  THOR_SCSI_LOG(WARNING)
    << "Inserting a 'Marker' at the beginning as requested:"
    << " first element is neither 'Marker' nor of zero length. \n";
  this->addMarkerAtStart();
}

template<class C>
template<typename T>
int
ts::AcceleratorKnobbable<C>::_propagate
(thor_scsi::core::ConfigType& conf, gtpsa::ss_vect<T> &ps, size_t start_elem,
 int max_elements, size_t n_turns, bool tracy_compatible_indexing)// const
{
  /* I guess Tobin would complain about this extra complexity */
  if(tracy_compatible_indexing){
    if(start_elem<=0) {
      std:: stringstream strm;
      strm << "Requested tracy compatible indexing, but start element "
	   << start_elem << "was smaller or equal to zero";
      throw std::runtime_error(strm.str());
    }
    // take the one off
    start_elem -= 1;
  }

  bool retreat = std::signbit(max_elements);

  // default arguments ... do all
  if(start_elem == 0 && max_elements == ts::max_elements_default){
    max_elements = static_cast<int>(this->size());
  }
  int next_elem = static_cast<int>(start_elem);
  auto trace = this->trace();

  for(size_t turn=0; turn<n_turns; ++turn) {
    if(trace)
      (*trace) << "turn " << turn << ", processing elements from " << start_elem
	       << " for a maximum of elements " << max_elements << std::endl;

    next_elem = static_cast<int>(start_elem);
    // iterate over the number of requested elements ...
    for(int i=0; i<std::abs(max_elements); i++)
      {
	// check if next element is acceptable ...
	// currently not required as I access elements using .at which checks
	// that the access is within range
	// assert(next_elem >= 0 && next_elem<nelem);
	size_t n = next_elem;

	std::shared_ptr<tsc::CellVoid> cv = this->at(n);
	// auto elem = std::dynamic_pointer_cast<tsc::ElemTypeKnobbed/*<C>*/>
	//   (cv);

#if 0
	// J.B. 14/07/23: debugging.
	printf("_propagate: %-8s\n", cv->name.c_str());
	if (std::dynamic_pointer_cast<tsc::CellVoid>(cv) != NULL)
	  printf("  CellVoid\n", cv->name.c_str());
	if (std::dynamic_pointer_cast<tsc::ElemTypeKnobbed>(cv) != NULL)
	  printf("  ElemTypeKnobbed\n", cv->name.c_str());
	if (std::dynamic_pointer_cast<tse::QuadrupoleType>(cv) != NULL)
	  printf("  QuadrupoleType\n", cv->name.c_str());
	if (std::dynamic_pointer_cast<tse::QuadrupoleTypeTpsa>(cv) != NULL)
	  printf("  QuadrupoleTypeTpsa\n", cv->name.c_str());
#endif

	std::shared_ptr<tsc::ElemTypeKnobbed>
	  elem = std::dynamic_pointer_cast<tsc::ElemTypeKnobbed>(cv);
	if(!elem){
	  // Should raise an exception!
	  THOR_SCSI_LOG(ERROR)
	    << "Failed to cast to element to ElemtypeKnobbed " << cv->name
	    << "\n";
	  std::runtime_error("Could not cast cell void to elemtype");
	  return next_elem;
	}
	if(retreat) {
	  next_elem--;
	} else {
	  next_elem++;
	}
	std::shared_ptr<tsc::Observer> observer = elem->observer();
	const auto celem =
	  std::const_pointer_cast<tsc::ElemTypeKnobbed/*<C>*/>(elem);
	/* keep view on observer the same , end always in forward upstream
	   direction */
	const auto observed_first =
	  (retreat) ? tsc::ObservedState::end : tsc::ObservedState::start;
	const auto observed_last  =
	  (retreat) ? tsc::ObservedState::start : tsc::ObservedState::end;
	if(observer){
	  observer->view(celem, ps, observed_first, 0);
	}
	elem->propagate(conf, ps);
	if(observer){
	  observer->view(celem, ps, observed_last, 0);
	}
	auto aperture = elem->getAperture();
	if(aperture){
	  /* could be any aperture ... */
	  bool flag = elem->checkAmplitude(ps);
	  if(not flag){
	    THOR_SCSI_LOG(INFO) << "Element lost at " << elem.get()
				<<" with aperture " << aperture.get();
	    /* get a guess on the plane */
	    auto rect_apt = std::dynamic_pointer_cast<tse::RectangularAperture>
	      (aperture);
	    if(rect_apt){
	      double x=-1e0, y=-1e0;
	      if(x < 0e0){
		conf.lossplane = tse::PlaneKind::Horizontal;
	      } else if(y < 0e0){
		conf.lossplane = tse::PlaneKind::Vertical;
	      } else {
		throw ts::SanityCheckError
		  ("Particle lost but no loss plane identified!");
	      }
	    }
	    return next_elem;
	  }
	}
	if(trace)
	  (*trace) << "After ["<< n<< "] " << cv->name << " " << std::endl
		   << ps << std::endl;
      }
  }
  return next_elem;
}


template<class C>
int
ts::AcceleratorKnobbable<C>::propagate
(thor_scsi::core::ConfigType& conf, ss_vect_dbl  &ps, size_t start,
 int max_elements, size_t n_turns, bool tracy_compatible_indexing) // const
{
  return
    _propagate
    (conf, ps, start, max_elements, n_turns, tracy_compatible_indexing);
}

/*
  int
  ts::AcceleratorKnobbable::
  propagate(thor_scsi::core::ConfigType& conf, ss_vect_tps  &ps, size_t start,
  int max_elements, size_t n_turns, bool tracy_compatible_indexing) // const
  {
    return
      _propagate
      (conf, ps, start, max_elements, n_turns, tracy_compatible_indexing);
  }
*/

template<class C>
int
ts::AcceleratorKnobbable<C>::propagate
(thor_scsi::core::ConfigType& conf, ss_vect_tpsa &ps, size_t start,
 int max_elements, size_t n_turns,  bool tracy_compatible_indexing)// const
{
  return
    _propagate
    (conf, ps, start, max_elements, n_turns, tracy_compatible_indexing);
}

template ts::AcceleratorKnobbable<tsc::StandardDoubleType>::AcceleratorKnobbable
(const Config &conf, bool add_marker_at_start);
template ts::AcceleratorKnobbable<tsc::TpsaVariantType>::AcceleratorKnobbable
(const Config &conf, bool add_marker_at_start);
template ts::AcceleratorKnobbable<tsc::StandardDoubleType>::AcceleratorKnobbable
(std::vector<std::shared_ptr<thor_scsi::core::CellVoid>>& elements,
 bool add_marker_at_start);
template ts::AcceleratorKnobbable<tsc::TpsaVariantType>::AcceleratorKnobbable
(std::vector<std::shared_ptr<thor_scsi::core::CellVoid>>& elements,
 bool add_marker_at_start);
template ts::AcceleratorKnobbable<tsc::StandardDoubleType>::AcceleratorKnobbable
(const std::vector<std::shared_ptr<thor_scsi::core::ElemTypeKnobbed>>& elements,
 bool add_marker_at_start);
template ts::AcceleratorKnobbable<tsc::TpsaVariantType>::AcceleratorKnobbable
(const std::vector<std::shared_ptr<thor_scsi::core::ElemTypeKnobbed>>& elements,
 bool add_marker_at_start);
template void ts::AcceleratorKnobbable<tsc::TpsaVariantType>::addMarkerAtStart
(void);
template void ts::AcceleratorKnobbable<tsc::StandardDoubleType>::
addMarkerAtStart(void);

template
int ts::AcceleratorKnobbable<tsc::StandardDoubleType>::propagate
(thor_scsi::core::ConfigType&, ss_vect_tpsa &ps, size_t start, int max_elements,
 size_t n_turns, bool tracy_compatible_indexing);
template
int ts::AcceleratorKnobbable<tsc::StandardDoubleType>::propagate
(thor_scsi::core::ConfigType&, ss_vect_dbl &ps, size_t start, int max_elements,
 size_t n_turns, bool tracy_compatible_indexing);

template
int ts::AcceleratorKnobbable<tsc::TpsaVariantType>::propagate
(thor_scsi::core::ConfigType&, ss_vect_tpsa &ps, size_t start, int max_elements,
 size_t n_turns, bool tracy_compatible_indexing);
template
int ts::AcceleratorKnobbable<tsc::TpsaVariantType>::propagate
(thor_scsi::core::ConfigType&, ss_vect_dbl  &ps, size_t start, int max_elements,
 size_t n_turns, bool tracy_compatible_indexing);
