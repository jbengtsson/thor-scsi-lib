#ifndef _THOR_SCSI_CORE_LEGACY_IO_CODES_H_
#define _THOR_SCSI_CORE_LEGACY_IO_CODES_H_ 1
/**
   numerical type codes for output

   Todo:
       compare these to thor_scsi/core/element_enums.h

*/

/*
 * Spot the difference!

Source prtmfile.cc    Source rdmfile.cc

#define marker_   -1  #define marker_   -1
#define drift_     0  #define drift_     0
#define mpole_     1  #define mpole_     1
#define cavity_    2  #define cavity_    2
#define thinkick_  3  #define thinkick_  3
#define wiggler_   4  #define wiggler_   4
#define kick_map_  6  #define kick_map   6
#define map_       7  #define map_       7
*/

namespace thor_scsi {
	namespace legacy {
		enum parts_kind_io{
			marker_   = -1,
			drift_    =  0,
			mpole_    =  1,
			cavity_   =  2,
			thinkick_ =  3,
			wiggler_  =  4,
			kick_map_ =  6,
			map_      =  7,
		};
	}
}
#endif /* _THOR_SCSI_CORE_LEGACY_IO_CODES_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
