#include <thor_scsi/core/multipoles.h>
#include <thor_scsi/core/machine.h>
#include <thor_scsi/elements/field_kick.h>
#include <thor_scsi/elements/element_helpers.h>
#include <thor_scsi/elements/utils.h>
#include <tps/tps_type.h>
#include <tps/math_pass.h>

#include <sstream>

namespace ts = thor_scsi;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

/* ==========================================================================
   Support functions
   required by other elements too ?

   If so: should go to element_helpers.h
   ========================================================================== */

namespace thor_scsi::elements {
	template<typename T>
	void edge_focus(const tsc::ConfigType &conf, const double irho,
			const double phi /* in degrees */, const double gap, ss_vect<T> &ps);

	template<typename T>
	void p_rot(const tsc::ConfigType &conf, double phi, ss_vect<T> &ps);

	template<typename T>
	void bend_fringe(const tsc::ConfigType &conf, const double hb, ss_vect<T> &ps);

	template<typename T>
	void quad_fringe(const tsc::ConfigType &conf, const double b2, ss_vect<T> &ps);

}

/**
 *
 * \verbatim embed:rst:leading-asterisk
 * .. Warning::
 *
 *     phi is in degrees
 * \endverbatim
 *
 */
template<typename T>
void tse::edge_focus(const tsc::ConfigType &conf, const double irho,
		     const double phi /* in degrees */, const double gap, ss_vect<T> &ps)
{
  ps[px_] += irho*tan(degtorad(phi))*ps[x_];
  if (!conf.dip_edge_fudge) {
    // Remark: Leads to a diverging Taylor map (see SSC-141).
    // ps[py_] -=
    //   irho*tan(degtorad(phi)-get_psi(irho, phi, gap))
    //   *ps[y_]/(1e0+ps[delta_]);
    // Leading order correction.
    ps[py_] -=
      irho*tan(degtorad(phi)-get_psi(irho, phi, gap))*ps[y_]*(1e0-ps[delta_]);
  } else
    ps[py_] -= irho*tan(degtorad(phi)-get_psi(irho, phi, gap))*ps[y_];
}

/*
 *
 * @f[\phi@f] ... dipole bend angle
 *
 * \verbatim embed:rst:leading-asterisk
 * .. Warning::
 *
 *     phi is in degrees
 *
 * .. Todo:
 *     how does it differ from PRotTransform? Should it be implemented there?
 * \endverbatim
 */
template<typename T>
void tse::p_rot(const tsc::ConfigType &conf, double phi, ss_vect<T> &ps)
{
  T  c, s, t, pz, val; // p
  ss_vect<T> ps1;

  c = cos(degtorad(phi));
  s = sin(degtorad(phi));
  t = tan(degtorad(phi));
  pz = get_p_s(conf, ps);

  if (!conf.H_exact && !conf.Cart_Bend) {
     ps[px_] = s*pz + c*ps[px_];
  } else {
    // ps1 = ps; p = c*pz - s*ps1[px_];
    // px[x_]   = ps1[x_]*pz/p; px[px_] = s*pz + c*ps1[px_];
    // px[y_]  += ps1[x_]*ps1[py_]*s/p;
    // px[ct_] += (1e0+ps1[delta_])*ps1[x_]*s/p;

    ps1 = ps; val = 1e0 - ps1[px_]*t/pz;
    ps[x_]  = ps1[x_]/(c*val);
    ps[px_] = ps1[px_]*c + s*pz;
    ps[y_]  = ps1[y_] + t*ps1[x_]*ps1[py_]/(pz*val);
    ps[ct_] = ps1[ct_] + ps1[x_]*(1e0+ps1[delta_])*t/(pz*val);
  }
}

/*
 *
 *
 * \verbatim embed:rst:leading-asterisk
 *
 * .. Todo:
 *
 *   Revisit speed of light check!
 * Throw exception?
 *
 * \endverbatim
 */
template<typename T>
void tse::bend_fringe(const tsc::ConfigType &conf, const double hb, ss_vect<T> &ps)
{
  T          coeff, u, pz, pz2, pz3;
  ss_vect<T> ps1;

  coeff = -hb/2e0;
  ps1 = ps;
  pz = get_p_s(conf, ps);
  pz2 = sqr(pz); pz3 = pz*pz2;
  u = 1e0 + 4e0*coeff*ps1[px_]*ps1[y_]*ps1[py_]/pz3;
  if (u >= 0e0) {
    ps[y_]  = 2e0*ps1[y_]/(1e0+sqrt(u));
    ps[x_]  = ps1[x_] - coeff*sqr(ps[y_])*(pz2+sqr(ps1[px_]))/pz3;
    ps[py_] = ps1[py_] + 2e0*coeff*ps1[px_]*ps[y_]/pz;
    ps[ct_] = ps1[ct_] - coeff*ps1[px_]*sqr(ps[y_])*(1e0+ps1[delta_])/pz3;
  } else {
	  std::cerr << __FILE__ << "::" <<  __FUNCTION__ << "@" <<  __LINE__
		    <<" : (bend_fringe): *** Speed of light exceeded!" << std::endl;
    std::stringstream stream;
    stream << ": (bend_fringe): *** Speed of light exceeded!" << " u = " << u;

    ps[x_] = NAN; ps[px_] = NAN; ps[y_] = NAN; ps[py_] = NAN;
    ps[delta_] = NAN; ps[ct_] = NAN;
    throw ts::PhysicsViolation(stream.str());
  }
}


template<typename T>
void tse::quad_fringe(const tsc::ConfigType &conf, const double b2, ss_vect<T> &ps)
{
  T u, p_s;

  u = b2/(12e0*(1e0+ps[delta_])); p_s = u/(1e0+ps[delta_]);
  ps[py_] /= 1e0 - 3e0*u*sqr(ps[y_]); ps[y_] -= u*tse::cube(ps[y_]);
  if (conf.Cavity_on) ps[ct_] -= p_s*tse::cube(ps[y_])*ps[py_];
  ps[px_] /= 1e0 + 3e0*u*sqr(ps[x_]);
  if (conf.Cavity_on) ps[ct_] += p_s*tse::cube(ps[x_])*ps[px_];
  ps[x_] += u*tse::cube(ps[x_]); u = u*3e0; p_s = p_s*3e0;
  ps[y_] = exp(-u*sqr(ps[x_]))*ps[y_]; ps[py_] = exp(u*sqr(ps[x_]))*ps[py_];
  ps[px_] += 2e0*u*ps[x_]*ps[y_]*ps[py_];
  if (conf.Cavity_on) ps[ct_] -= p_s*sqr(ps[x_])*ps[y_]*ps[py_];
  ps[x_] = exp(u*sqr(ps[y_]))*ps[x_]; ps[px_] = exp(-u*sqr(ps[y_]))*ps[px_];
  ps[py_] -= 2e0*u*ps[y_]*ps[x_]*ps[px_];
  if (conf.Cavity_on) ps[ct_] += p_s*sqr(ps[y_])*ps[x_]*ps[px_];
}



/* ==========================================================================
   End support functions
   ========================================================================== */



void tse::FieldKick::FieldKickForthOrder::computeIntegrationSteps(void)
{

	if(!parent){
		return;
	}
	const double Pirho = parent->getCurvature(), length = parent->getLength();
	double dL;
	auto n_steps = this->getNumberOfIntegrationSteps();

	if(parent->assumingCurvedTrajectory()){
		// along the arc
		dL = 2e0/ Pirho * sin(length * Pirho/2e0) / n_steps;
	}else{
		// along the straight line
		dL = length / n_steps;
	}
	this->splitIntegrationStep(dL, &this->dL1, &this->dL2, &this->dkL1, &this->dkL2);
}

void tse::FieldKick::FieldKickForthOrder::
splitIntegrationStep(const double dL, double *dL1, double *dL2, double *dkL1, double *dkL2)
  const
{
	const int n_steps = this->getNumberOfIntegrationSteps();

	*dL1 = this->c_1*dL;
	*dL2 = this->c_2*dL;
	*dkL1 = this->d_1*dL;
	*dkL2 = this->d_2*dL;
}


/**
 *
 * @ brief: 4th order symplectic integrator in the body
 *
 *
 * explanation:
 * Symplectic integrator
 * 2nd order
 *
 *  L/2 -> bnl -> L/2
 *
 * Here 4th order
 *
 * 4 th order
 *
 * d_1 * L -> c_1 * bnl -> d_2 * L -> c_2 * bnl -> d_1 * L
 *
 * d_1 + d_2 + d_1 = L
 *
 * d_2 negative drift
 * c_1 + c_2 = 1
 *
 *
 */
template<typename T>
inline void tse::FieldKick::FieldKickForthOrder::
_localPass(tsc::ConfigType &conf, ss_vect<T> &ps)
{


	double dL = 0.0, h_ref = 0.0;
	auto PN = this->integration_steps;
	const double length = parent->getLength(), Pirho = parent->getCurvature();

	THOR_SCSI_LOG(DEBUG) << "\n  name = " << parent->name << " N = " << PN << "\n";

	if (!conf.Cart_Bend) {
		// Polar coordinates.
		h_ref = Pirho; dL = length/PN;
	} else {
		// Cartesian coordinates.
		h_ref = 0e0;
		//if (Pirho == 0e0){
		if(!parent->assumingCurvedTrajectory()){
			// along the straight line
			dL = length/PN;
		}else{
			// along the arc
			dL = 2e0/Pirho*sin(length*Pirho/2e0)/PN;
		}
	}

	/*
	 * Calculating the individual pieces
	 */
	double dL1, dL2, dkL1, dkL2;
	this->splitIntegrationStep(dL, &dL1, &dL2, &dkL1, &dkL2);

	// std::cout <<  "local pass " <<  std::endl;
	// std::cout.flush();

	auto intp_shared_ptr = this->getFieldInterpolator();
	// std::cout << "intp_shared_ptr  " << intp_shared_ptr << " count "
	// 	  << intp_shared_ptr.use_count() << std::endl;

	auto& t_intp = *(intp_shared_ptr.get());
	auto* parent = this->parent;

	if(!parent){
		throw std::logic_error("parent was nullptr");
	}
	// computeRadiationIntegralsStart
#ifdef SYNCHROTRON_INTEGRALS
	parent->_synchrotronIntegralsInit(conf, ps);
#endif /* SYNCHROTRON_INTEGRALS */

	/* 4th order integration steps  */
	for (int seg = 1; seg <= PN; seg++) {
		const int rad_step = (seg - 1) * 4;
		// computeRadiationIntegralsStep
		THOR_SCSI_LOG(DEBUG) << "\n  seg = " << seg << "\n";


#ifdef SYNCHROTRON_INTEGRALS
		parent->_synchrotronIntegralsStep(conf, ps, rad_step);
#endif /* SYNCHROTRON_INTEGRALS */

		drift_pass(conf, dL1, ps);
		// call to radiation before thin kick
		parent->thinKickAndRadiate(conf, t_intp, dkL1, Pirho, h_ref, ps);
		drift_pass(conf, dL2, ps);
		// call to radiation before thin kick
		parent->thinKickAndRadiate(conf, t_intp, dkL2, Pirho, h_ref, ps);

		// why this step only here
		// computeRadiationIntegralsStep
#ifdef SYNCHROTRON_INTEGRALS
		parent->_synchrotronIntegralsStep(conf, ps, rad_step + 1);
#endif /* SYNCHROTRON_INTEGRALS */

		drift_pass(conf, dL2, ps);
		// call to radiation before thin kick
		parent->thinKickAndRadiate(conf, t_intp, dkL1, Pirho, h_ref, ps);
		drift_pass(conf, dL1, ps);

		// computeRadiationIntegralsStep
#ifdef SYNCHROTRON_INTEGRALS
		parent->_synchrotronIntegralsStep(conf, ps, rad_step + 2);
#endif /* SYNCHROTRON_INTEGRALS */
	}
#ifdef SYNCHROTRON_INTEGRALS
	parent->_synchrotronIntegralsFinish(conf, ps);
#endif /* SYNCHROTRON_INTEGRALS */
}


tse::FieldKick::FieldKick(const Config &config) : tse::FieldKickAPI(config)
{
	// Field interpolation type
	this->intp = nullptr;
	this->asThick(false);
	this->setIntegrationMethod(config.get<double>("Method", 4));
	this->setNumberOfIntegrationSteps(config.get<double>("N"));
	this->setLength(config.get<double>("L"));

	this->setBendingAngle(config.get<double>("T", 0.0));
	this->setEntranceAngle(config.get<double>("T1", 0.0));
	this->setExitAngle(config.get<double>("T2", 0.0));
	this->integ4O.setParent(this);

}

tse::FieldKick::FieldKick(tse::FieldKick&& O)
	: tse::FieldKickAPI(std::move(O))

{
	this->setBendingAngle(O.getBendingAngle());
	this->setEntranceAngle(O.getEntranceAngle());
	this->setExitAngle(O.getExitAngle());
	this->setCurvature(O.getCurvature());
	this->asThick(O.isThick());

	this->setNumberOfIntegrationSteps(O.getNumberOfIntegrationSteps());
	this->setIntegrationMethod(O.getIntegrationMethod());

	this->Pgap = O.Pgap;
	this->integ4O.setParent(this);
	this->rad_del = std::move(O.rad_del);
	return;

	std::cerr << __FILE__ << "::" << __FUNCTION__ << "@" << __LINE__
		  << " integration method " << this->getIntegrationMethod() << std::endl;

}

void tse::FieldKick::show(std::ostream& strm, const int level) const
{
	tse::LocalGalileanPRot::show(strm, level);
	if(level >= 1){
		/* at least intercept it with a blank */
		strm << " lens type: " <<  (this->isThick() ? "thick" : "thin") << " ";
		if(this->isThick()){
			strm << " n steps " << this->getNumberOfIntegrationSteps() << " " ;
		}
		strm << " angles: bending " << this->getBendingAngle()
		     << " entrance " << this->getEntranceAngle()
		     << " exit " << this->getExitAngle();
		strm << " curvature " << this->getCurvature()
		     << " has curvature " << std::boolalpha
		     << this->assumingCurvedTrajectory();
		if(!this->intp){
			strm << " NO interpolater set!";
		} else {
			strm << " ";
			this->intp->show(strm, level);
		}

		auto rad_del = this->getRadiationDelegate();
		if(!rad_del){
			strm << "radiation delegate=None";
		} else {
			strm << " radiation delegate=";
			rad_del->show(strm, level);
		}
		strm<< ", ";
	}
}


template<typename T>
void tse::FieldKick::_quadFringe(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps)
{

	if (!conf.quad_fringe){ return ; }

	std::cerr << __FILE__ << "::" << __FUNCTION__ << "@" << __LINE__
		  << "Code not yet tested " << std::endl;
	throw thor_scsi::NotImplemented("Quadrupole fringe implemented but code not yet" \
					" tested");

	auto muls = std::dynamic_pointer_cast<tsc::TwoDimensionalMultipoles>(this->getFieldInterpolator());
	if(!muls){
		std::cerr << __FILE__ << "::" << __FUNCTION__ << "@" << __LINE__
			  << "Quadfringe can currently only be calculated for "
			  << " PlanarMultipole interpolator " << std::endl;
		throw thor_scsi::NotImplemented("Quadfringe only works with multipoles as field interpolators");
	}
#warning "Not calling gradient but doing direct look up"
	// T Gy=0e0, Gx=0e0;
	// should be checked how to implement
	// muls->gradient(ps[x_], ps[y_], &Gx, &Gy);
	const tsc::cdbl C2 = muls->getMultipole(2);
	tse::quad_fringe(conf, C2.real(), ps);

	/*
	   if (conf.quad_fringe && (PB[Quad+HOMmax] != 0e0)){
	         tse::quad_fringe(conf, PB[Quad+HOMmax], ps);
	   }
	*/
}

template<typename T>
inline void tse::FieldKick::
thinKickAndRadiate(const thor_scsi::core::ConfigType &conf,
		   const thor_scsi::core::Field2DInterpolation& intp, const double L,
		   const double h_bend, const double h_ref, ss_vect<T> &ps)
{

	// const auto x = ps[x_];
	// const auto y = ps[y_];
 	T BxoBrho, ByoBrho;

	//intp.field(ps[x_], ps[y_], &BxoBrho, &ByoBrho);
	/*
	  if(!this->intp){
		throw std::logic_error("Interpolation object not set!");
	}
	*/

	THOR_SCSI_LOG(DEBUG) << "\n  ps = " << ps << "\n";

	intp.field(ps[x_], ps[y_], &BxoBrho, &ByoBrho);

	auto rad = this->getRadiationDelegate();
	if(rad){
		THOR_SCSI_LOG(DEBUG) <<  "\nComputing radiation:\n";
		std::array<T, 3> B = {BxoBrho, ByoBrho + h_bend, 0};
		rad->radiate(conf, ps, L, h_ref, B);
	}
	tse::thin_kick(conf, BxoBrho, ByoBrho, L, h_bend, h_ref, ps);

	THOR_SCSI_LOG(DEBUG) << "\n  ps = " << ps << "\n";
}

/**
 * @brief: thin kick: element length 0, integral kick effect
 *
 * typically used for implementing correctors
 *
 */
template<typename T>
inline void tse::FieldKick::_localPassThin(const tsc::ConfigType &conf, ss_vect<T> &ps)
{
	// length has to be zero here
	/// todo: add a check at least for debug purposes
	const bool debug = false;
	const double length = 1.0;
	if(debug){
		std::cerr << "calling thin kick " << std::endl;
	}
	// call to radiation before thin kick
	this->thinKickAndRadiate(conf, *this->intp, length, 0e0, 0e0, ps);
}


template<typename T>
inline void tse::FieldKick::_localPassBody(tsc::ConfigType &conf, ss_vect<T> &ps)
{
	/* define integration step */
	if (!this->isThick()) {
		tse::FieldKick::_localPassThin(conf, ps);
		return;
	}
	this->integ4O._localPass(conf, ps);
}

/*
 *
 *
 *
 * \verbatim embed:rst:leading-asterisk
 *
 *

 * .. Todo:
 *
 *     Review if mpole should be simplified....
 *
 *    Three function calls
 *       edge
 *       body
 *       edge
 *
 *   Revisit switch irho == 0e0
 *
 *   Split up functionality
 *   Revisit radiation calculations
 *
 *
 * \endverbatim
 */
template<typename T>
void tse::FieldKick::_localPass(tsc::ConfigType &conf, ss_vect<T> &ps)
{


	switch (Pmethod) {
	case Meth_Fourth: break;
		/*
		  Pmethod should be test when set ... thus can not pass
		  down here
		*/
	default:
		std::cerr <<  __FILE__ << "::" << __FUNCTION__<< "@" << __LINE__
			  << " element '" << this->name << "'"
			  << " : Method # "<< this->Pmethod << "  not supported "
			  << std::endl;
		throw ts::SanityCheckError("Methods checked but still ending here?");
		break;

	}

#ifdef SYNCHROTRON_INTEGRALS
	this->_synchrotronIntegralsInit(conf, ps);
#endif /* SYNCHROTRON_INTEGRALS */

	// Set start

	// Matrix method
	if (conf.mat_meth) {
		std::cerr << __FILE__ << "::" << __FUNCTION__ << "@" << __LINE__
			  << ": matrix method not implemented" << std::endl;
		throw thor_scsi::NotImplemented("Matrix method not implemented");
		/*
		if(Porder <= Quad){
			ps = mat_pass(M_elem, ps);
		}
		*/
		return;
		// if (conf.emittance && !conf.Cavity_on) &&
		// 	(PL != 0e0) && (Pirho != 0e0)) get_dI_eta_5(this);

	}

	// symplectic integration below
	// Fringe fields.
	this->_quadFringe(conf, ps);

	if (!conf.Cart_Bend) {
		if (this->assumingCurvedTrajectory()){
			tse::edge_focus(conf, Pirho, PTx1, Pgap, ps);
		}
	} else {
		// here in Carthesian coordinates

		/* horizontal focusing: purely geometric effect */
		tse::p_rot(conf, PTx1, ps);
		/* vertical focusing: leading order effect */
		tse::bend_fringe(conf, Pirho, ps);
	}

	// Body calculation delegated
	tse::FieldKick::_localPassBody(conf, ps);

	// Fringe fields.
	if (!conf.Cart_Bend) {
		if (this->assumingCurvedTrajectory()){
			tse::edge_focus(conf, Pirho, PTx2, Pgap, ps);
		}
	} else {
		tse::bend_fringe(conf, -Pirho, ps); p_rot(conf, PTx2, ps);
	}
	this->_quadFringe(conf, ps);
}


template void tse::FieldKick::_localPass(tsc::ConfigType &conf, ss_vect<double> &ps);
template void tse::FieldKick::_localPass(tsc::ConfigType &conf, ss_vect<tps> &ps);
