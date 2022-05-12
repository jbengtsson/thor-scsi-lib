#ifndef _THOR_SCSI_COMPAT_NAMESPACES_H_
#define _THOR_SCSI_COMPAT_NAMESPACES_H_ 1

/**

     Include the whole API in a standard namespace.

     Target:          Have access to the new API in a flat format

     Suggestions for shortcuts for the different name spaces

     .. :

        namespace ts = thor_scsi;
	namespace tsc = thor_scsi::core;
	namespace tse = thor_scsi::elements;
	namespace tsl = thor_scsi::legacy;
	namespace tsx = thor_scsi::compat;


 */
#warning "Please consider to use names spaces!"
#include <thor_scsi/compat/compat.h>
#include <thor_scsi/legacy/legacy.h>

using namespace thor_scsi;
using namespace thor_scsi::core;
using namespace thor_scsi::elements;
using namespace thor_scsi::legacy;
using namespace thor_scsi::compat;

using namespace std;

#endif /* _THOR_SCSI_COMPAT_NAMESPACES_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
