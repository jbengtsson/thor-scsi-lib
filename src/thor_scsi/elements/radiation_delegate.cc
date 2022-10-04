#include <thor_scsi/core/machine.h>
#include <thor_scsi/elements/radiation_delegate.h>
#include <thor_scsi/elements/element_helpers.h>
#include <thor_scsi/elements/utils.h>

namespace ts = thor_scsi;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

std::string tse::RadiationDelegateInterface::repr(void) const
{
	std::stringstream strm;
	strm << *this;
	return strm.str();
}
std::string tse::RadiationDelegateKickInterface::repr(void) const
{
	std::stringstream strm;
	strm << *this;
	return strm.str();
}

template <typename T>
inline void tse::RadiationDelegate::computeAndStoreCurlyH(const ss_vect<T> &ps)
{
	this->curly_dH_x = is_tps<T>::get_curly_H(ps);
}
template <typename T>
inline void tse::RadiationDelegate::_view(const tsc::ElemType& elem, const ss_vect<T> &ps, const enum tsc::ObservedState state, const int cnt)
{
	switch(state){
	case tsc::ObservedState::start:
		this->reset();
		this->delegator_name = elem.name;
		this->delegator_index = elem.index;
		return;
		break;
	case tsc::ObservedState::end:
		this->computeAndStoreCurlyH(ps);
		return;
		break;
	default:
		return;
	}
}

//template
//void tse::RadiationDelegate::_view(const tsc::ElemType& elem, const ss_vect<double> &ps, const enum tsc::ObservedState state, const int cnt);
//template
//void tse::RadiationDelegate::_view(const tsc::ElemType& elem, const ss_vect<tps> &ps, const enum tsc::ObservedState state, const int cnt);
void tse::RadiationDelegate::view(const tsc::ElemType& elem, const ss_vect<double> &ps, const enum tsc::ObservedState state, const int cnt){
	_view(elem, ps, state, cnt);
}
void tse::RadiationDelegate::view(const tsc::ElemType& elem, const ss_vect<tps> &ps, const enum tsc::ObservedState state, const int cnt){
	_view(elem, ps, state, cnt);
}

template<typename T>
inline void tse::RadiationDelegateKick::synchrotronIntegralsFinish(const FieldKickAPI &kick, const ss_vect<T> &ps)
{

	// Why only when cavities are not on ?
	// Needs A^-1.
	const T x = ps[x_], y = ps[y_];
	double Gx = NAN, Gy = NAN;
	kick.getFieldInterpolator()->gradient(x, y, &Gx, &Gy);
	const double PL = kick.getLength();
	const double Pirho = kick.getCurvature();
	const auto PN = kick.getNumberOfIntegrationSteps();
	this->curly_dH_x /= 6e0*PN;
	this->dI[1] += PL*is_tps<tps>::get_dI_eta(ps)*Pirho;
	this->dI[2] += PL*sqr(Pirho);
	this->dI[3] += PL*fabs(tse::cube(Pirho));
	this->dI[4] *=
		PL*Pirho*(sqr(Pirho)+2e0*Gy)
		/(6e0*PN);
	this->dI[5] += PL*fabs(tse::cube(Pirho))*curly_dH_x;

}

template<typename T>
inline void tse::RadiationDelegateKick::synchrotronIntegralsStep(const ss_vect<T> &ps)
{
	// Needs A^-1.
	this->curly_dH_x += is_tps<tps>::get_curly_H(ps);
	this->dI[4] += is_tps<tps>::get_dI_eta(ps);
}


inline void tse::RadiationDelegateKick::diffusion(const tps &B2_perp,  const tps &ds, const tps &p_s0,  const ss_vect<tps> &A)
{

	int          j;
	double       B_66;
	ss_vect<tps> A_inv;

	if (B2_perp > 0e0){
		// Fix move function to RadiationDelegateKick
		auto q_fluct = this->q_fluct;

		B_66 = (q_fluct*pow(B2_perp.cst(), 1.5)*pow(p_s0, 4)*ds).cst();
		A_inv = Inv(A);
		// D_11 = D_22 = curly_H_x,y * B_66 / 2,
		// curly_H_x,y = eta_Fl^2 + etap_Fl^2
		for (j = 0; j < 3; ++j){
			this->D_rad[j] +=
				(sqr(A_inv[j*2][delta_])+sqr(A_inv[j*2+1][delta_]))*B_66/2e0;
		}
	}
}


template <typename T>
inline void tse::RadiationDelegateKick::_view(const FieldKickAPI& kick, const ss_vect<T> &ps, const enum tsc::ObservedState state, const int cnt)
{
	switch(state){
	case tsc::ObservedState::start:
		this->reset();
		this->delegator_name = kick.name;
		this->delegator_index = kick.index;
		return;
		break;
	case tsc::ObservedState::event:
		this->synchrotronIntegralsStep(ps);
		return;
		break;
	case tsc::ObservedState::end:
		this->synchrotronIntegralsFinish(kick, ps);
		return;
		break;
	default:
		return;
	}
}
void tse::RadiationDelegate::show(std::ostream& strm, int level) const{
	strm << "RadiationDelegate for "
	     << this->delegator_name << "["<< this->delegator_index << "]"
	     <<" curly_dH_x " << this->curly_dH_x;
}

void tse::RadiationDelegateKick::show(std::ostream& strm, int level) const
{
	strm << "RadiationDelegateKick for "
	     << this->delegator_name << "["<< this->delegator_index << "]"
	     << ":, energy " << this->getEnergy()
	     << " curly_dH_x " << this->curly_dH_x;
	strm << " synchotron integrals = [";
	int cnt = 0;
	for(auto dIe : this->dI){
		if(cnt){
			strm << ", ";
		}
		strm << dIe;
		++cnt;
	}
	strm << "]";
	strm << " diffusion = [";
	cnt = 0;
	for(auto c : this->D_rad){
		if(cnt){
			strm << ", ";
		}
		strm << c;
		++cnt;
	}
	strm << "]";
}

/**
 * @brief Computing |B^2_perp| perpendicular to the arc of circle.
 */
template<typename T>
void get_B2(const double h_ref, const std::array<T,3> B, const ss_vect<T> &xp,
	    T &B2_perp, T &B2_par)
{
  // compute B_perp^2 and B_par^2
  T xn, e[3];

  xn = 1e0/sqrt(sqr(1e0+xp[x_]*h_ref)+sqr(xp[px_])+sqr(xp[py_]));
  e[X_] = xp[px_]*xn; e[Y_] = xp[py_]*xn; e[Z_] = (1e0+xp[x_]*h_ref)*xn;

  // left-handed coordinate system
  B2_perp =
    sqr(B[Y_]*e[Z_]-B[Z_]*e[Y_]) + sqr(B[X_]*e[Y_]-B[Y_]*e[X_])
    + sqr(B[Z_]*e[X_]-B[X_]*e[Z_]);

  THOR_SCSI_LOG(DEBUG)
    << "\nB2_perp = " << B2_perp;
//  B2_par = sqr(B[X_]*e[X_]+B[Y_]*e[Y_]+B[Z_]*e[Z_]);
}

void tse::RadiationDelegateKick::setEnergy(const double val)
{
        // energy in eV
	this->energy = val;

	// the equation below is for GeV
	auto energy_GeV = this->energy / 1e9;
	this->q_fluct = C_q*C_gamma/(M_PI*sqr(1e-9*m_e))*pow(energy_GeV, 5e0);

}

template<typename T>
static inline double ps_coor_ref_value(T coor);

static inline double ps_coor_ref_value(double coor)
{
	return coor;
}
static inline double ps_coor_ref_value(tps coor)
{
	return coor.cst();
}
template<typename T>
static bool check_ps_finite(ss_vect<T> ps, const double max_val = 1e3)
{
	bool check_ps_finite = true;
	for(int i=0; i < nv_tps; ++i){
		double ref_val = ps_coor_ref_value(ps[i]);
		check_ps_finite &= std::isfinite(ref_val);
		check_ps_finite &= (std::abs(ref_val) < max_val);
	}
	return check_ps_finite;
}


template<typename T>
void tse::RadiationDelegateKick::
radiate(const thor_scsi::core::ConfigType &conf, ss_vect<T> &ps, const double L,
	const double h_ref, const std::array<T, 3> B)
{
	// M. Sands "The Physics of Electron Storage Rings" SLAC-121, p. 98.
	// ddelta/d(ds) = -C_gamma*E_0^3*(1+delta)^2*(B_perp/(Brho))^2/(2*pi)
	T  p_s0, p_s1, ds, B2_perp = 0e0, B2_par = 0e0;

	const bool radiation = conf.radiation;
	const bool compute_diffusion = conf.emittance;
	if(!radiation){
	    return;
	}
#if 0
	THOR_SCSI_LOG(INFO)
	  << "\nRadiate ->:\n" << this->delegator_name << "\n" << "  ps = "
	  << ps;
#else
	THOR_SCSI_LOG(INFO) << "\nRadiate ->:\n" << "\n" << "  ps = " << ps;
#endif

	if(!check_ps_finite(ps)){
		std::stringstream strm;
		strm <<  "ps unbound "; ps.show(strm, 10, false);
		std::cerr << strm.str() << std::endl;
		THOR_SCSI_LOG(ERROR) <<  "Check radiaton" << strm.str()
				     << " \n";
		throw ts::PhysicsViolation(strm.str());
	}

	// longitudinal component
	p_s0 = get_p_s(conf, ps);
	ss_vect<T> cs = ps; //.clone()
	// Large ring: x' and y' unchanged.
	cs[px_] /= p_s0;
	cs[py_] /= p_s0;

	if(!check_ps_finite(cs)){
		std::stringstream strm;
		strm << "ps unbound "; ps.show(strm, 10, false);
		std::cerr << strm.str() << std::endl;
		THOR_SCSI_LOG(ERROR)
		  <<  "Check radiaton" << strm.str() << " \n";

		throw ts::PhysicsViolation(strm.str());
	}

	// H = -p_s => ds = H*L.
	ds = (1e0+cs[x_]*h_ref+(sqr(cs[px_])+sqr(cs[py_]))/2e0)*L;
	// Compute perpendicular reference curve for comoving frame.
	get_B2(h_ref, B, cs, B2_perp, B2_par);

	//THOR_SCSI_LOG(INFO)
	// Test on legacy code
	// Should be 1 at the end ...
	const double energy_scale = 1e0;
	const double cl_rad =
	  C_gamma * cube(this->energy * energy_scale) / (2e0 * M_PI);

	if (radiation) {
		ps[delta_] -= cl_rad*sqr(p_s0)*B2_perp*ds;
		p_s1 = get_p_s(conf, ps);
		ps[px_] = cs[px_]*p_s1;
		ps[py_] = cs[py_]*p_s1;
	}
	if(!check_ps_finite(cs)){
		std::stringstream strm;
		strm << "ps unbound "; ps.show(strm, 10, false);
		std::cerr << strm.str() << std::endl;
		THOR_SCSI_LOG(ERROR) <<  "Check radiaton" << strm.str()
				     << " \n";

		throw ts::PhysicsViolation(strm.str());
	}

	if (compute_diffusion){
		this->diffusion(B2_perp, ds, p_s0, cs);
	}

	if(!check_ps_finite(cs)){
		std::stringstream strm;
		strm << "ps unbound "; ps.show(strm, 10, false);
		std::cerr << strm.str() << std::endl;
		THOR_SCSI_LOG(ERROR) <<  "Check radiaton" << strm.str()
				     << " \n";

		throw ts::PhysicsViolation(strm.str());
	}
#if 0
	THOR_SCSI_LOG(INFO)
	  << "\nRadiate ->:\n" << this->delegator_name << "\n" << "  ps = "
	  << ps;
#else
	THOR_SCSI_LOG(INFO) << "\nRadiate ->:\n" << "\n" << "  ps = " << ps;
#endif
}



//template void tse::RadiationDelegate::_view(const FieldKickAPI& kick, const ss_vect<double> &ps, const enum tsc::ObservedState state, const int cnt);
//template void tse::RadiationDelegate::_view(const FieldKickAPI& kick, const ss_vect<tps> &ps, const enum tsc::ObservedState state, const int cnt);
void tse::RadiationDelegateKick::
view(const FieldKickAPI& kick, const ss_vect<double> &ps,
     const enum ObservedState state, const int cnt)
{
	std::cout << "Rad Del.view(ss_vect<double>) for element "
		  << kick.name << "at index" << kick.index << std::endl;
	_view(kick, ps, state, cnt);
}

void tse::RadiationDelegateKick::
view(const FieldKickAPI& kick, const ss_vect<tps> &ps,
     const enum ObservedState state, const int cnt)
{
	std::cout<< "Rad Del.view(ss_vect<tps>) for element " << kick.name
		 << "at index" << kick.index << std::endl;
	_view(kick, ps, state, cnt);
}
template void tse::RadiationDelegateKick::
radiate(const thor_scsi::core::ConfigType &conf, ss_vect<double> &ps,
	const double L, const double h_ref, const std::array<double, 3> B);
template void tse::RadiationDelegateKick::
radiate(const thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps,
	const double L, const double h_ref, const std::array<tps, 3> B);
