#include <thor_scsi/elements/element_helpers.h>
#include <thor_scsi/elements/elements_enums.h>
#include <thor_scsi/elements/utils.h>

#include <tps/tps_type.h>

namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;

// template tps sqr(const tps &);

double tse::get_psi(const double irho, const double phi, const double gap)
{

  double psi;
  const double phir = degtorad(phi);

  const double k1 = 0.5e0, k2 = 0e0;

  if (phi == 0e0)
    psi = 0e0;
  else
    psi = k1*gap*irho*(1e0+sqr(sin(phir)))/cos(phir)
          *(1e0 - k2*gap*irho*tan(phir));

  return psi;
}


void tse::get_twoJ(const int n_DOF, const ss_vect<double> &ps, const ss_vect<tps> &A,
	      double twoJ[])
{
  int             j, no;
  long int        jj[ps_dim];
  ss_vect<double> z;

  no = no_tps;
  danot_(1);

  for (j = 0; j < ps_dim; j++)
    jj[j] = (j < 2*n_DOF)? 1 : 0;

  z = (PInv(A, jj)*ps).cst();

  for (j = 0; j < n_DOF; j++)
    twoJ[j] = sqr(z[2*j]) + sqr(z[2*j+1]);

  danot_(no);
}


namespace thor_scsi::elements{
	template<typename T>
	void tse::drift_pass(const tsc::ConfigType &conf, const double L, ss_vect<T> &ps)
	{
		T u;

		if (!conf.H_exact) {
			// Small angle axproximation.
			u = L/(1e0+ps[delta_]);
			ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
			ps[ct_] += u*(sqr(ps[px_])+sqr(ps[py_]))/(2e0*(1e0+ps[delta_]));
		} else {
			u = L/tse::get_p_s(conf, ps);
			ps[x_]  += u*ps[px_]; ps[y_] += u*ps[py_];
			ps[ct_] += u*(1e0+ps[delta_]) - L;
		}
		if (conf.pathlength) ps[ct_] += L;
	}

	/**
	 *
	 * The vector potential for the combined-function sector bend is from:
	 *  C. Iselin "Lie Transformations and Transport Equations for Combined-
	 * Function Dipoles" Part. Accel. 17, 143-155 (1985).
	 *
	 * Todo:
	 *     Split up function in different parts or functions
	 *     E.g.: one for dipoles and one for anything else...
	 */
	template<typename T>
	void tse::thin_kick(tsc::ConfigType &conf, const int Order,
			    const thor_scsi::core::Field2DInterpolation& intp,
			    // const MpoleArray &MB,
			    const double L, const double h_bend, const double h_ref,
			    ss_vect<T> &ps)
	{
		int        j;
		T          BxoBrho, ByoBrho, ByoBrho1, B[3], u, p_s;
		ss_vect<T> ps0;

		// Todo: what kind of bend is that?
		if ((h_bend != 0e0) || ((1 <= Order)
#if 0
					//  why this switch
					&& (Order <= HOMmax)
#endif
			    )) {
			ps0 = ps;
#warning "feld interpolation for bends missing"
#if 0
			/*
			 * Compute magnetic field with Horner's rule.
			 *
			 *  Todo: replace it with a call to field ...
			 *        but requires a dedicated enty
			 */
			ByoBrho = MB[Order+HOMmax]; BxoBrho = MB[HOMmax-Order];
			for (j = Order-1; j >= 1; j--) {
				ByoBrho1 = ps0[x_]*ByoBrho - ps0[y_]*BxoBrho + MB[j+HOMmax];
				BxoBrho  = ps0[y_]*ByoBrho + ps0[x_]*BxoBrho + MB[HOMmax-j];
				ByoBrho  = ByoBrho1;
			}
#endif
			{
				double Bx, By;


			}

			if (conf.radiation || conf.emittance) {
				B[X_] = BxoBrho; B[Y_] = ByoBrho + h_bend; B[Z_] = 0e0;
				radiate(conf, ps, L, h_ref, B);
			}

			if (h_ref != 0e0) {
				// Sector bend.
				if (true) {
					ps[px_] -= L*(ByoBrho+(h_bend-h_ref)/2e0+h_ref*h_bend*ps0[x_]
						      -h_ref*ps0[delta_]);
					ps[ct_] += L*h_ref*ps0[x_];
				} else {
					// The Hamiltonian is split into: H_d + H_k; with [H_d, H_d] = 0.
					p_s = get_p_s(conf, ps0); u = L*h_ref*ps0[x_]/p_s;
					ps[x_]  += u*ps0[px_];
					ps[y_]  += u*ps0[py_];
					ps[ct_] += u*(1e0+ps0[delta_]);
					// ps[px_] -= L*(h_bend*(1e0+h_ref*ps0[x_])-h_ref*p_s);

#warning "feld interpolation for sector bends missing"
#if 0
					// Field expansion up to sextupole like terms.
					ByoBrho += h_bend - MB[Quad+HOMmax]*h_ref*sqr(ps0[y_])/2e0;
					ps[px_] -= L*((1e0+h_ref*ps0[x_])*ByoBrho-h_ref*p_s);
					ps[py_] += L*(1e0+h_ref*ps0[x_])*BxoBrho;
#endif
				}
			} else {
				// Cartesian bend.
				ps[px_] -= L*(h_bend+ByoBrho);
			}
			ps[py_] += L*BxoBrho;
		}
	}

}
//template void tse::drift_pass(tsc::ConfigType &conf, const double, ss_vect<tps> &);
template void tse::drift_pass(const tsc::ConfigType &conf, const double, ss_vect<double> &);
