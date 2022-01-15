/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -                                         */

#ifndef THOR_SCSI_PROCESS_T2RING_H
#define THOR_SCSI_PROCESS_T2RING_H

#include <thor_scsi/core/elements_basis.h>
#include <thor_scsi/core/lattice.h>
#include <thor_scsi/process/t2ring_common.h>
#include <tps/ss_vect.h>
#include <tps/tps_type.h>
#include <vector>

namespace thor_scsi {
// Computation result files
const char beam_envelope_file[] = "beam_envelope";

/**
 *
 * returns true if stable
 *
 * Not assuming mid-plane symmetry, the charachteristic polynomial for a
 * symplectic periodic matrix is given by
 *
 * @f[
 *    P(lambda) = det(M-lambda*I)
 *    = (lambda-lambda0)(lambda-1/lambda0)
 *	   (lambda-lambda1)(lambda-1/lambda1)
 * @f]
 *
 * It follows that
 * @f[
 *      P(1) = (2-x)(2-y),     P(-1) = (2+x)(2+y)
 * @f]
 *
 * where
 * @f[
 *   x = (lambda0+1/lambda0)/2 = cos(2 pi nu_x)
 * @f]
 *
 * and similarly for y. Eliminating y
 * @f[
 *    x^2 + 4 b x + 4 c = 0
 * @f]
 *
 * where
 * @f[
 *      b = (P(1)-P(-1))/16,    c =(P(1)+P(-1))/8 - 1
 * @f]
 *
 * Solving for x
 * @f[
 *    x,y = -b +/- sqrt(b^2-c)
 * @f]
 *
 * where the sign is given by
 * @f[
 *     trace(hor) > trace(ver)
 *  @f]
 *
 * gives
 * @f[
 *     nu_x,y = arccos(x,y)/(2 pi)
 * @f]
 *
 * For mid-plane symmetry it simplies to
 * @f[
 *    nu_x = arccos((m11+m22)/2)/(2 pi)
 * @f]
 *
 */
bool GetNu(std::vector<double> &nu, std::vector< std::vector<double> > &M);

/**
 * compute twiss parameters
 *
 * Get Alpha beta gamma nu
 *
 */
bool Cell_GetABGN(std::vector< std::vector<double> > &M,
		  std::vector<double> &alpha, std::vector<double> &beta,
		  std::vector<double> &gamma, std::vector<double> &nu);

/**
 * Get the dispersion
 *
 * Todo:
 *   Rename in computeDispersion ?
 */
void Cell_Geteta(long i0, long i1, bool ring, double dP);

void Cell_Twiss(long i0, long i1, ss_vect<tps> &Ascr, bool chroma, bool ring,
		double dP);

void Ring_Getchrom(double dP);

// Not implemented code?
// void Ring_GetTwiss(bool chroma, double dP);

#if 0
void Ring_Fittune(std::vector<double> &nu, double eps, std::vector<int> &nq,
		  long qf[], long qd[], double dkL, long imax);

void Ring_Fitchrom(std::vector<double> &ksi, double eps,
		   std::vector<> &ns, long sf[], long sd[], double dkpL,
		   long imax);

void Ring_FitDisp(long pos, double eta, double eps, long nq, long q[],
		  double dkL, long imax);
#endif

void get_dI_eta_5(const int k, std::vector<thor_scsi::elements::ElemType*> Elem);

double get_code(const thor_scsi::core::ConfigType &conf, const thor_scsi::elements::ElemType &Cell);

void Cell_Twiss(const long int i0, const long int i1);

void findcod(thor_scsi::core::LatticeType &lat, double dP);

void prt_lin_map(const int n_DOF, const ss_vect<tps> &map);

ss_vect<tps> get_A(const std::vector<double> &alpha,
		   const std::vector<double> &beta,
		   const std::vector<double> &eta,
		   const std::vector<double> &etap);

ss_vect<tps> get_A_CS(const int n, const ss_vect<tps> &A,
		      std::vector<double> &dnu);

void get_ab(const ss_vect<tps> &A, std::vector<double> &alpha,
	    std::vector<double> &beta, std::vector<double> &dnu,
	    std::vector<double> &eta, std::vector<double> &etap);


void Trac(thor_scsi::core::LatticeType &lat, double x, double px, double y, double py, double dp,
	  double ctau, long nmax, long pos, long &lastn, long &lastpos,
	  FILE *outf1);

void SetKLpar(thor_scsi::core::LatticeType &lat, int Fnum, int Knum, int Order, double kL);

double GetKpar(thor_scsi::core::LatticeType &lat, int Fnum, int Knum, int Order);

double get_dynap(thor_scsi::core::LatticeType &lat, const double delta, const int n_aper,
		 const int n_track, const bool cod);


void computeFandJ(thor_scsi::core::LatticeType &lat, int n, double *x, ss_vect<double> *fjac,
		  double *fvect);
int Newton_Raphson(thor_scsi::core::LatticeType &lat, int n, ss_vect<double> &x, int ntrial,
		   double tolx);

void get_dnu(const int n, const ss_vect<tps> &A, std::vector<double> &dnu);

ss_vect<tps> tp_S(const int n_DOF, const ss_vect<tps> &A);

void dynap(FILE *fp, thor_scsi::core::LatticeType &lat, double r, const double delta,
	   const double eps, const int npoint, const int nturn,double x[],
	   double y[], const bool floqs, const bool cod, const bool print);

double get_aper(int n, double x[], double y[]);

/**
 * Todo: n_DOF: could it be derived from ss_vect thats allocated internally ?
 */
ss_vect<tps> get_S(const int n_DOF);

void getdynap(thor_scsi::core::LatticeType &lat, double &r, double phi, double delta, double eps,
	      int nturn, bool floqs);
}
#endif /* THOR_SCSI_PROCESS_T2RING_H */
