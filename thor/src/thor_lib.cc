/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -

*/
#include <thor_scsi/core/lattice.h>
#include <thor_scsi/core/elements_basis.h>
#include <thor_scsi/math/interpolation.h>
namespace tsc = thor_scsi::core;
namespace tsm = thor_scsi::math;
namespace tse = thor_scsi::elements;

#include "interpolation.cc"
#include "interpolation_elements.cc"
#include "t2elem.cc"
// #include "t2ring.cc"
// #include "sigma_track.cc"



// #include "orb_corr.cc"
// #include "param.cc"
// #include "dynap.cc"

// ALS.
// #include "physlib.cc"
// #include "fft.cc"

// Soleil.
// #include "naffutils.cc"
// #include "modnaff.cc"
// #include "soleillib.cc"

// NSLS-II.
// #include "nsls-ii_lib.cc"


//#include "at_thor.cc"

int
  Fnum_Cart,
  n_iter_Cart;

double
  u_Touschek,           // argument for Touschek D(ksi)
  chi_m;                // argument for IBS D(ksi)

// IBS (Bjorken-Mtingwa)
double a_IBS, b_IBS, c_IBS, a_k_IBS, b_k_IBS;

// for IBS
int    i_, j_;
double **C_;

// ss_vect<tps> map;
// MNF_struct MNF;

/**
 *
 * Todo:
 *    Check if not part of TPSA
 */

double eps_tps = 1e-25; // Floating point truncation.


// Instantiate templates.
#include <thor_scsi/core/elements.h>
#include <thor_scsi/core/config.h>
#include <tps/ss_vect.h>
#include <tps/tps_type.h>
#include <iostream>
#include <string>
#include <math.h>

namespace tse=thor_scsi::elements;
namespace tsc=thor_scsi::core;



template void GtoL(ss_vect<double> &, std::vector<double> &,
		   std::vector<double> &, const double, const double,
		   const double);
template void GtoL(ss_vect<tps> &, std::vector<double> &, std::vector<double> &,
		   const double, const double, const double);

template void LtoG(ss_vect<tps> &, std::vector<double> &, std::vector<double> &,
		   double, double, double);
template void LtoG(ss_vect<double> &, std::vector<double> &,
		   std::vector<double> &, double, double, double);

template void p_rot(tsc::ConfigType &conf, double, ss_vect<double> &);
template void p_rot(tsc::ConfigType &conf, double, ss_vect<tps> &);


template void get_B2(const double, const double [], const ss_vect<double> &,
		     double &, double &);
template void get_B2(const double, const tps [], const ss_vect<tps> &,
		     tps &, tps &);

template void radiate(tsc::ConfigType &conf, ss_vect<double> &, const double,
		      const double, const double []);
template void radiate(tsc::ConfigType &conf, ss_vect<tps> &, const double,
		      const double, const tps []);

template void radiate_ID(tsc::ConfigType &conf, ss_vect<double> &,
			 const double, const double &);
template void radiate_ID(tsc::ConfigType &conf, ss_vect<tps> &,
			 const double, const tps &);

template void Drift(tsc::ConfigType &conf, const double, ss_vect<double> &);
template void Drift(tsc::ConfigType &conf, const double, ss_vect<tps> &);

template void bend_fringe(tsc::ConfigType &conf, const double, ss_vect<double> &);
template void bend_fringe(tsc::ConfigType &conf, const double, ss_vect<tps> &);

template void EdgeFocus(tsc::ConfigType &conf, const double, const double,
			const double, ss_vect<double> &);
template void EdgeFocus(tsc::ConfigType &conf, const double, const double,
			const double, ss_vect<tps> &);

template void quad_fringe(tsc::ConfigType &conf, const double, ss_vect<double> &);
template void quad_fringe(tsc::ConfigType &conf, const double, ss_vect<tps> &);


template void thin_kick(tsc::ConfigType &conf, const int, const tse::MpoleArray &,
			const double, const double, const double,
			ss_vect<double> &);
template void thin_kick(tsc::ConfigType &conf, const int, const tse::MpoleArray &,
			const double, const double, const double,
			ss_vect<tps> &);

template void Cav_Focus(const double L, const double delta, const bool entrance,
			ss_vect<double> &ps);
template void Cav_Focus(const double L, const double delta, const bool entrance,
			ss_vect<tps> &ps);
template void Cav_Focus(const double L, const tps delta, const bool entrance,
			ss_vect<tps> &ps);

template void Wiggler_pass_EF(tsc::ConfigType &conf, const tse::ElemType *elem,
			      ss_vect<double> &x);
template void Wiggler_pass_EF(tsc::ConfigType &conf, const tse::ElemType *elem,
			      ss_vect<tps> &x);

template void Wiggler_pass_EF2(tsc::ConfigType &conf, int nstep, double L,
			       double kxV, double kxH, double kz,
			       double BoBrhoV, double BoBrhoH, double phi,
			       ss_vect<double> &x);
template void Wiggler_pass_EF2(tsc::ConfigType &conf, int nstep, double L,
			       double kxV, double kxH, double kz,
			       double BoBrhoV, double BoBrhoH, double phi,
			       ss_vect<tps> &x);

template void Wiggler_pass_EF3(tsc::ConfigType &conf, tse::ElemType *Cell,
			       ss_vect<double> &x);
template void Wiggler_pass_EF3(tsc::ConfigType &conf, tse::ElemType *Cell,
			       ss_vect<tps> &x);

template void sol_pass(tsc::ConfigType &conf, const tse::ElemType *, ss_vect<double> &);
template void sol_pass(tsc::ConfigType &conf, const tse::ElemType *, ss_vect<tps> &);

namespace thor_scsi{
  namespace math {
    template void LinearInterpolation2(double &x, double &z, double &tx, double &tz, double &B2,
				       thor_scsi::elements::ElemType *,
				       bool &, int);
    template void LinearInterpolation2(tps &, tps &, tps &, tps &, tps &,
				       thor_scsi::elements::ElemType *, bool &, int);
  }
}


template void tsm::SplineInterpolation2(double &, double &, double &, double &,
				   tse::ElemType *, bool &);
template void tsm::SplineInterpolation2(tps &, tps &, tps &, tps &,
				   tse::ElemType *, bool &);

template void spline(const double [], const double [], int const,
		     double const, const double, double []);
template void spline(const double [], const tps [], int const,
		     double const, const double, tps []);

template void splint(const double[], const double [], const double [],
		     const int, const double &, double &);
template void splint(const double[], const double [], const double [],
		     const int, const tps &, tps &);

template void splint(const double[], const tps [], const tps [],
		     const int, const tps &, tps &);

template void splin2(const double [], const double [],
		     double **, double **, const int, const int,
		     const double &, const double &, double &);
template void splin2(const double [], const double [],
		     double **, double **, const int, const int,
		     const tps &, const tps &, tps &);



#if 0
double d_sign(double a, double b)
{
  double x;

  x = (a >= 0 ? a : - a);
  return( b >= 0 ? x : -x);
}

void t2init(void)
{

  std::cerr << "Why is this code still be used" << std::endl;
  throw ts::SanityCheckError();

//  iniranf(0); /* initialise le generateur aleatoire: graine 0 */

//  fprintf(stdout,"pi = %20.16e \n",pi);

//  daini((long)no_, (long)nv_, 0);

//  lieini((long)no_, (long)nv_, (long)nd2_);
}
#endif

/*
// Matlab BS
void exit_(int exit_code)
{

  printf("fatal error, <ret> to continue "); std::cin.ignore(1, '\n');

  exit(exit_code);
}
*/


void prt_name_ascii(std::string &name);

void prt_name_ascii(std::string &name)
{
  int i;

  printf("  %-8s (", name.c_str());
  for (i = 0; i < (int)name.length(); i++)
    printf(" %3d", (int)name[i]);
  printf(" )\n");
}
#if 0
#endif

long int tsc::LatticeType::ElemIndex(const std::string &name)
{
  long        i;
  std::string name1 = name;

  const bool prt = false;

  for (i = 0; i < (int)name1.length(); i++)
    name1[i] = tolower(name1[i]);

  if (prt) {
    printf("\nElemIndex:\n     ");
    prt_name_ascii(name1);
    printf("\n");
  }

  for (i = 1; i <= (int)elemf.size(); i++) {
    if (prt) {
      printf("  %3ld", i);
      prt_name_ascii(elemf[i-1].ElemF->Name);
    }

    if (name1 == elemf[i-1].ElemF->Name) break;
  }

  if (name1 != elemf[i-1].ElemF->Name) {
    std::cerr << "ElemIndex: undefined element >" << name << "<" << std::endl;
    throw std::invalid_argument("No element with given name");
      // exit_(1);
  }

  return i;
}
