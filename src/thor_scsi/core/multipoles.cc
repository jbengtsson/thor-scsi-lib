#include <thor_scsi/core/multipoles.h>
#include <thor_scsi/core/math_comb.h>
#include <tps/utils.h>
#include <algorithm>

namespace tsc = thor_scsi::core;

tsc::TwoDimensionalMultipoles::TwoDimensionalMultipoles(const unsigned int h_max)
{
	if(h_max<=1){
		throw std::logic_error("max multipole must be at least 1");
	}
	this->m_max_multipole = h_max;
	this->coeffs.resize(h_max);
	const cdbl_intern zero(0.0, 0.0);
	for(auto h : this->coeffs){
		h = zero;
	}
}

#if 0
// Why do I need a copy constructor ??
tsc::TwoDimensionalMultipoles::TwoDimensionalMultipoles(const tsc::TwoDimensionalMultipoles& o): coeffs(o.coeffs){
	this->m_max_multipole = o.m_max_multipole;
}
#endif

tsc::TwoDimensionalMultipoles::TwoDimensionalMultipoles(const tsc::TwoDimensionalMultipoles&& o) : coeffs(std::move(o.coeffs))
{
	this->m_max_multipole = o.m_max_multipole;
}

tsc::TwoDimensionalMultipoles& tsc::TwoDimensionalMultipoles::operator= (const tsc::TwoDimensionalMultipoles &o)
{
	this->coeffs = o.coeffs;
	this->m_max_multipole = o.m_max_multipole;
	return *this;
}


tsc::TwoDimensionalMultipoles tsc::TwoDimensionalMultipoles::clone(void) const
{
	return tsc::TwoDimensionalMultipoles(std::vector<tsc::cdbl_intern>(this->coeffs));
}

tsc::TwoDimensionalMultipoles::TwoDimensionalMultipoles(std::vector<cdbl_intern> const coeffs)
{
	if(coeffs.size()<=1){
		throw std::logic_error("max multipole must be at least 1");
	}
	this->coeffs = coeffs;
	this->m_max_multipole = this->coeffs.size();
}


tsc::TwoDimensionalMultipoles& tsc::TwoDimensionalMultipoles::operator += (const tsc::TwoDimensionalMultipoles &other)
{

	if(this->getCoeffs().size() != other.getCoeffs().size()){
		throw std::runtime_error("Length of coeffcients does not match");
	}
	auto& c = this->getCoeffs();
	auto& oc = other.getCoeffs();

	std::transform(oc.begin(), oc.end(), c.begin(), c.begin(), std::plus<>());

	return *this;
}

tsc::TwoDimensionalMultipoles tsc::TwoDimensionalMultipoles::operator+ (const tsc::TwoDimensionalMultipoles &other) const
{

	unsigned int n = std::max(this->getCoeffs().size(), other.getCoeffs().size());
	TwoDimensionalMultipoles nh(n);
	nh += *this;
	nh += other;
	return nh;

}

void tsc::TwoDimensionalMultipoles::show(std::ostream& strm, int level) const
{
	strm << "interpolation='2D multipoles'";
	if (level > 2){
		strm << "(num=" << this->coeffs.size();
	}
	if (level  > 3){
		strm << ", muls={";
		int i = 1;
		for(auto c : this->coeffs){
			if(i > 1){
				strm << ", ";
			}
			strm << i << ":[" << c.real() << ", " << c.imag() << "]";
			++i;
		}
		strm << "}";
	}
	strm<<")";
}

void tsc::TwoDimensionalMultipoles::applyTranslation(const tsc::cdbl dzs)
{

	for (unsigned int i = 0; i < this->coeffs.size(); ++i) {
		cdbl_intern dzi(dzs.real(), dzs.imag());
		for (unsigned j = i+1; j < this->coeffs.size(); ++j) {
			this->coeffs[i] +=  double(binom(j, i)) * (this->coeffs[j] * dzi);
		}
	}
}

#if 0
std::vector<double> toTracyRepresentation(const std::vector<cdbl> &coeffs)
{
	/*
	 * Todo:
	 *    Check if this should not always be to size HOMmax
	 */
	std::vector<double> tracy_coeffs(coeffs.size() * 2);

	// #define HOMmax   21     // [a_n, b_n] <=> [-HOMmax..HOMmax].
	int i = 0;
	for(auto c : coeffs){
		/*
		 * the skew harmonics
		 * Todo:
		 *     check if tracy requires a negative flag here
		 */
		tracy_coeffs[i] = c.imag();
		++ i;
	}

	// followed by the normal ones: don't reset i they just follow afterwards
	/*
	 * Todo:
	 *    Check if 0 is normal or skew or empty
	 */
	for(auto c : coeffs){
		tracy_coeffs[i] = c.real();
	}
	return tracy_coeffs;

}
#endif

/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
