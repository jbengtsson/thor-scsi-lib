/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -

*/

#ifndef _THOR_SCSI_PROCESS_LSOC_H_
#define _THOR_SCSI_PROCESS_LSOC_H_ 1

#include <thor_scsi/core/lattice.h>
#include <vector>
namespace thor_scsi {
extern int              n_bpm_[2], n_corr_[2];
extern std::vector<int> bpms_[2], corrs_[2];

void zero_trims(thor_scsi::core::LatticeType &lat);

void prt_gcmat(const int plane);

void gcmat(thor_scsi::core::LatticeType &lat, const int plane);

void gcmat(thor_scsi::core::LatticeType &lat, const int n_bpm, const std::vector<int> bpms,
	   const int n_corr, const std::vector<int> corrs, const int plane,
	   const bool svd);

void gcmat(thor_scsi::core::LatticeType &lat, const int bpm, const int corr,
	   const int plane);

void lsoc(thor_scsi::core::LatticeType &lat, const int plane, const double scl);

void gtcmat(thor_scsi::core::LatticeType &lat, const int plane);

void gtcmat(thor_scsi::core::LatticeType &lat, const int n_bpm, const std::vector<int> bpms,
	    const int n_corr, const std::vector<int> corrs, const int plane,
	    const bool svd);

void lstc(thor_scsi::core::LatticeType &lat, const int plane, const double scl);
}
#endif /* _THOR_SCSI_PROCESS_LSOC_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
