#include <thor_scsi/elements/mpole.h>
#include <thor_scsi/elements/element_helpers.h>
#include <thor_scsi/elements/elements_enums.h>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

template<typename T>
void tse::MpoleType::_localPass(tsc::ConfigType &conf, ss_vect<T> &ps)
{
  int          seg = 0, i;
  double       dL = 0e0, dL1 = 0e0, dL2 = 0e0,
               dkL1 = 0e0, dkL2 = 0e0, h_ref = 0e0;

  // now handled by LocalCoordinateElement
  //GtoL(ps, dS, dT, Pc0, Pc1, Ps1);

  if (conf.emittance && !conf.Cavity_on) {
    // Needs A^-1.
    curly_dH_x = 0e0;
    for (i = 0; i <= 5; i++)
      /* Synchrotron integrals */ dI[i] = 0e0;
  }

  switch (Pmethod) {

  case Meth_Fourth:
    if (conf.mat_meth && (Porder <= Quad)) {
      ps = mat_pass(M_elem, ps);

      // if (conf.emittance && !conf.Cavity_on) &&
      // 	(PL != 0e0) && (Pirho != 0e0)) get_dI_eta_5(this);
    } else {
      // Fringe fields.
      if (conf.quad_fringe && (PB[Quad+HOMmax] != 0e0))
	quad_fringe(conf, PB[Quad+HOMmax], ps);
      if (!conf.Cart_Bend) {
	if (Pirho != 0e0)
	  EdgeFocus(conf, Pirho, PTx1, Pgap, ps);
      } else {
	p_rot(conf, PTx1, ps); bend_fringe(conf, Pirho, ps);
      }

      if (Pthick == thick) {
	if (!conf.Cart_Bend) {
	  // Polar coordinates.
	  h_ref = Pirho; dL = PL/PN;
	} else {
	  // Cartesian coordinates.
	  h_ref = 0e0;
	  if (Pirho == 0e0)
	    dL = PL/PN;
	  else
	    dL = 2e0/Pirho*sin(PL*Pirho/2e0)/PN;
	}

	dL1 = c_1*dL; dL2 = c_2*dL; dkL1 = d_1*dL; dkL2 = d_2*dL;

	for (seg = 1; seg <= PN; seg++) {
	  if (conf.emittance && !conf.Cavity_on) {
	    // Needs A^-1.
	    curly_dH_x += is_tps<tps>::get_curly_H(ps);
	    dI[4] += is_tps<tps>::get_dI_eta(ps);
	  }

	  Drift(conf, dL1, ps);
	  thin_kick(conf, Porder, PB, dkL1, Pirho, h_ref, ps);
	  Drift(conf, dL2, ps);
	  thin_kick(conf, Porder, PB, dkL2, Pirho, h_ref, ps);

	  if (conf.emittance && !conf.Cavity_on) {
	    // Needs A^-1.
	    curly_dH_x += 4e0*is_tps<tps>::get_curly_H(ps);
	    dI[4] += 4e0*is_tps<tps>::get_dI_eta(ps);
	  }

	  Drift(conf, dL2, ps);
	  thin_kick(conf, Porder, PB, dkL1, Pirho, h_ref, ps);
	  Drift(conf, dL1, ps);

	  if (conf.emittance && !conf.Cavity_on) {
	    // Needs A^-1.
	    curly_dH_x += is_tps<tps>::get_curly_H(ps);
	    dI[4] += is_tps<tps>::get_dI_eta(ps);
	  }
	}

	if (conf.emittance && !conf.Cavity_on) {
	  // Needs A^-1.
	  curly_dH_x /= 6e0*PN;
	  dI[1] += PL*is_tps<tps>::get_dI_eta(ps)*Pirho;
	  dI[2] += PL*sqr(Pirho);
	  dI[3] += PL*fabs(cube(Pirho));
	  dI[4] *=
	    PL*Pirho*(sqr(Pirho)+2e0*PBpar[Quad+HOMmax])
	    /(6e0*PN);
	  dI[5] += PL*fabs(cube(Pirho))*curly_dH_x;
	}
      } else
	thin_kick(conf, Porder, PB, 1e0, 0e0, 0e0, ps);

      // Fringe fields.
      if (!conf.Cart_Bend) {
	if (Pirho != 0e0)
	  EdgeFocus(conf, Pirho, PTx2, Pgap, ps);
      } else {
	bend_fringe(conf, -Pirho, ps); p_rot(conf, PTx2, ps);
      }
      if (conf.quad_fringe && (PB[Quad+HOMmax] != 0e0))
	quad_fringe(conf, -PB[Quad+HOMmax], ps);
    }
    break;

  default:
    std::cerr <<  "Mpole_Pass: Method not supported " << Name
	      <<  " method " << Pmethod <<  std::endl;
    throw ts::NotImplemented();
    break;
  }

  // now handled by LocalCoordinateElement
  // LtoG(ps, dS, dT, Pc0, Pc1, Ps1);
}
