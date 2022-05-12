#ifndef _THOR_SCSI_VERSION_H_
#define _THOR_SCSI_VERSION_H_ 1
#include <string>
namespace thor_scsi {
	static const int thor_version_major = 3;
	static const int thor_version_minor = 5;
	static const int thor_version_patch_level = 0;

	const char * compiled_version = "3.5.0";
	/**
	 * stores run version when compile process is running
	 */
	extern const char *run_version;
}
#endif /* _THOR_SCSI_VERSION_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
