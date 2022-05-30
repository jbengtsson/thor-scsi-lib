#ifndef _THOR_SCSI_ELEMENTS_CHECK_MULTIPOLE_H_
#define _THOR_SCSI_ELEMENTS_CHECK_MULTIPOLE_H_
#include <memory>
#include <thor_scsi/core/multipoles.h>

void
check_only_major_multipole_set(std::shared_ptr<thor_scsi::core::TwoDimensionalMultipoles> muls,
			       const thor_scsi::core::cdbl ref, const size_t n_major);


#endif
