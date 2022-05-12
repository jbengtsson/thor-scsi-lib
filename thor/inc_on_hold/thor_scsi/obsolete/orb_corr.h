#ifndef _THOR_SCSI_ORB_CORR_H_
#define _THOR_SCSI_ORB_CORR_H_ 1

#include <thor_scsi/core/lattice.h>
#include <string>
#include <vector>

namespace thor_scsi{
class orb_corr_type {
private:
  double **A, **Ai, *w, **U, **V, *b, *bb, *x, *xx, eps, hcut, vcut;

public:
  bool                  hor, periodic;
  int                   m, n;
  std::vector<long int> bpms, corrs;

  void alloc(thor_scsi::core::LatticeType &lat, const long int i0, const long int i1,
	     const long int i2, const std::vector<std::string> &bpm_Fam_names,
	     const std::vector<std::string> &corr_Fam_names, const bool hor,
	     const bool periodic, const double eps);
  void alloc(thor_scsi::core::LatticeType &lat, const std::vector<std::string> &bpm_Fam_names,
	     const std::vector<std::string> &corr_Fam_names, const bool hor,
	     const bool periodic, const double eps);
  void dealloc(void);
  void get_trm_mat(thor_scsi::core::LatticeType &lat);
  void get_orm_mat(thor_scsi::core::LatticeType &lat);
  void svd_decomp(void);
  void prt_svdmat(thor_scsi::core::LatticeType &lat);
  void solve(thor_scsi::core::LatticeType &lat, const double scl) const;
  void clr_trims(thor_scsi::core::LatticeType &lat);
};


void codstat(thor_scsi::core::LatticeType &lat, double mean[], double sigma[], double xmax[],
	     const long lastpos, const bool all,
	     const std::vector<long int> &bpms);

void thread_beam(thor_scsi::core::LatticeType &lat, const int n_cell, const std::string &Fam_name,
		 const std::vector<std::string> &bpm_Fam_names,
		 const std::vector<std::string> corr_Fam_names[],
		 const int n_thread, const double scl);

void cod_ini(thor_scsi::core::LatticeType &lat, const std::vector<std::string> &bpm_Fam_names,
	     const std::vector<std::string> corr_Fam_names[],
	     orb_corr_type orb_corr[]);

bool cod_correct(thor_scsi::core::LatticeType &lat, const int n_orbit, const double scl,
		 orb_corr_type orb_corr[]);

};
#endif /* _THOR_SCSI_ORB_CORR_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
