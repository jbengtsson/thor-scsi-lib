#ifndef _THOR_SCSI_NSLS_II_LIB_H_
#define _THOR_SCSI_NSLS_II_LIB_H_ 1

#include <thor_scsi/core/lattice.h>
#include <thor_scsi/compat/typedefs.h>

namespace thor_scsi {
void lwr_case(char str[]);

void upr_case(char str[]);

//void prt_trace (void);

void chk_cod(const bool cod, const char *proc_name);

void no_sxt(thor_scsi::core::LatticeType &lat);

void get_map(thor_scsi::core::LatticeType &lat, const bool cod);

tps get_h(void);

void get_m2(const ss_vect<tps> &ps, tps m2[]);

double get_curly_H(const double alpha_x, const double beta_x,
		   const double eta_x, const double etap_x);

void printcod(thor_scsi::core::LatticeType &lat, const char *file_name);


//------------------------------------------------------------------------------

double get_Wiggler_BoBrho(thor_scsi::core::LatticeType &lat, const int Fnum, const int Knum);

void set_Wiggler_BoBrho(thor_scsi::core::LatticeType &lat, const int Fnum, const int Knum,
			const double BoBrhoV);
void set_Wiggler_BoBrho(thor_scsi::core::LatticeType &lat, const int Fnum, const double BoBrhoV);

void set_ID_scl(thor_scsi::core::LatticeType &lat, const int Fnum, const int Knum,
		const double scl);

void SetFieldValues_fam(thor_scsi::core::LatticeType &lat, const int Fnum, const bool rms,
			const double r0, const int n, const double Bn,
			const double An, const bool new_rnd);

void SetFieldValues_type(thor_scsi::core::LatticeType &lat, const int N, const bool rms,
			 const double r0, const int n, const double Bn,
			 const double An, const bool new_rnd);

void SetFieldErrors(thor_scsi::core::LatticeType &lat, const char *name, const bool rms,
		    const double r0, const int n, const double Bn,
		    const double An, const bool new_rnd);

bool CorrectCOD(thor_scsi::core::LatticeType &lat, const int n_orbit, const double scl);

double Touschek(thor_scsi::core::LatticeType &lat, const double Qb, const double delta_RF,
		const double eps_x, const double eps_y,
		const double sigma_delta, const double sigma_s);

double Touschek(thor_scsi::core::LatticeType &lat, const double Qb, const double delta_RF,
		const bool consistent, const double eps_x, const double eps_y,
		const double sigma_delta, double sigma_s, const int n_turn,
		const bool aper_on, double sum_delta[][2],
		double sum2_delta[][2]);

double f_IBS(const double chi_m);

double get_int_IBS(void);

void IBS(thor_scsi::core::LatticeType &lat, const double Qb, const double eps_SR[], double eps[],
	 const bool prt1, const bool prt2);

void IBS_BM(thor_scsi::core::LatticeType &lat, const double Qb, const double eps_SR[],
	    double eps[], const bool prt1, const bool prt2);

void rm_space(char *name);

void get_bn(thor_scsi::core::LatticeType &lat, const char file_name[], int n, const bool prt);

double get_chi2(long int n, double x[], double y[], long int m, thor_scsi::compat::psVector b);

void pol_fit(int n, double x[], double y[], int order, thor_scsi::compat::psVector &b,
	     double &sigma, const bool prt);

void get_ksi2(thor_scsi::core::LatticeType &lat, const double d_delta);

bool find_nu(const int n, const double nus[], const double eps, double &nu);

bool get_nu(thor_scsi::core::LatticeType &lat, const double Ax, const double Ay,
	    const double delta, double &nu_x, double &nu_y);

void dnu_dA(thor_scsi::core::LatticeType &lat, const double Ax_max, const double Ay_max,
	    const double delta, const int n_ampl);

bool orb_corr(thor_scsi::core::LatticeType &lat, const int n_orbit);

void get_alphac(thor_scsi::core::LatticeType &lat);

void get_alphac2(thor_scsi::core::LatticeType &lat);

/* void bend_cal_Fam(const int Fnum); */

/* void bend_cal(void); */

double h_ijklm(const tps &h, const int i, const int j, const int k,
	       const int l, const int m);

void set_tune(thor_scsi::core::LatticeType &lat, const char file_name1[],
	      const char file_name2[], const int n);

void prt_H_long(const int n, const double phi_max, const double delta_max,
		const double U0);

void get_map_twiss(const ss_vect<tps> &M,
		   double beta0[], double beta1[], double nu[], bool stable[]);

void set_map(thor_scsi::elements::MapType *Map);

void set_map(thor_scsi::core::LatticeType &lat, const int Fnum, const double dnu[]);

void set_map_per(thor_scsi::elements::MapType *Map,
		 const double alpha0[], const double beta0[],
		 const double eta0[], const double etap0[]);

void set_map_per(thor_scsi::core::LatticeType &lat, const int Fnum, const double alpha0[],
		 const double beta0[], const double eta0[],
		 const double etap0[]);

void set_map_reversal(thor_scsi::core::LatticeType &lat, thor_scsi::elements::CellType &Cell);

void set_map_reversal(thor_scsi::core::LatticeType &lat, const long int Fnum);

void setmp(long ilat, long m, long n, double rr, double bnoff, double cmn);

void setmpall(thor_scsi::core::LatticeType &lat, double rref);
};
#endif /* _THOR_SCSI_NSLS_II_LIB_H_ */
