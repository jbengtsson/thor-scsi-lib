#ifndef _THOR_SCSI_VERSION_H_
#define _THOR_SCSI_VERSION_H_ 1
#include <string>
namespace thor_scsi {
	static const int thor_version_major = 3;
	static const int thor_version_minor = 5;
	static const int thor_version_micro = 0;

	std::string get_compiled_version(void);
}
#endif /* _THOR_SCSI_VERSION_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
