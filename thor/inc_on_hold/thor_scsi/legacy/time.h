#ifndef _THOR_SCSI_LEGACY_TIME_H_
#define _THOR_SCSI_LEGACY_TIME_H_ 1

#include <time.h>
// Same C asctime \n.
char *asctime2(const struct tm *timeptr);
struct tm* GetTime();

#endif /* _THOR_SCSI_LEGACY_TIME_H_ */
