#ifndef _THOR_SCSI_LATTICE_FILTERS_H_
#define _THOR_SCSI_LATTICE_FILTERS_H_ 1

#include <thor_scsi/core/elements_filter.h>
#include <thor_scsi/core/lattice.h>

/*
 * Heavy lifting delegated to elements_filters
 */
namespace thor_scsi{
	inline auto filter_mpole_types(thor_scsi::core::LatticeType &lat){
		return thor_scsi::elements::filter_mpole_types(lat.elems);
	}
}


#endif /*  _THOR_SCSI_LATTICE_FILTERS_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
