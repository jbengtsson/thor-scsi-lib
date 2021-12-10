#ifndef _THOR_SCSI_PROCESS_DYNAP_H_
#define _THOR_SCSI_PROCESS_DYNAP_H_ 1

#include <thor_scsi/core/lattice.h>

class DA_data_type {
 private:
 public:

  bool track(thor_scsi::core::LatticeType &lat, const param_data_type &params,
	     const double x, const double px, const double y, const double py,
	     const double delta, const double f_rf, const bool prt);
  void get_r_stable(LatticeType &lat, const param_data_type &params,
		    double &r, const double phi, const double delta,
		    const double eps);
  double get_dynap(LatticeType &lat, param_data_type &params,
		   FILE *fp, const double r, const double delta,
		   const double eps, double x_min[], double x_max[]);
  void get_DA_bare(LatticeType &lat, param_data_type &params);
  void get_DA_real(LatticeType &lat, param_data_type &params,
		   orb_corr_type orb_corr[]);
  void get_mean_sigma(const int n, double &m, double &s);
};

#endif /* _THOR_SCSI_PROCESS_DYNAP_H_ */
