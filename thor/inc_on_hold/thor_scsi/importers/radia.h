#ifndef _THOR_SCSI_IMPORTERS_RADIA_H_
#define _THOR_SCSI_IMPORTERS_RADIA_H_ 1

#include <thor_scsi/core/config.h>
#include <thor_scsi/core/elements.h>

void Read_IDfile(char *fic_radia, const thor_scsi::core::ConfigType &conf,  thor_scsi::elements::InsertionType *ID);

/**

   Todo:

      Move it to interpolation
 */
void Matrices4Spline(thor_scsi::elements::InsertionType *);
#endif /* _THOR_SCSI_IMPORTERS_RADIA_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
