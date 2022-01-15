#ifndef _THOR_SCSI_COMPAT_H_
#define _THOR_SCSI_COMPAT_H_ 1
/**

   Provide a compatability layer to the former TRACY API
   Compatability is only provided as so far as it does not touch dangerous code

   Warning:
        If writing new code, please use the thor_scsi API. This compatability
	will only be supported in a limited fashion. Please consider to uprade
	your code!


 */
#include <thor_scsi/compat/constants.h>
#include <thor_scsi/compat/families.h>
#include <thor_scsi/compat/globals.h>
#include <thor_scsi/compat/namespaces.h>
#include <thor_scsi/compat/std_headers.h>
#include <thor_scsi/compat/typedefs.h>
#include <thor_scsi/thor_scsi.h>
// #include <thor_scsi/obsolete/nsls-ii_lib.h>
// #include <thor_scsi/obsolete/orb_corr.h>

#endif /*  _THOR_SCSI_COMPAT_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
