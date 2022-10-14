#include <thor_scsi/core/multipoles.h>
#include <thor_scsi/core/math_comb.h>
#include <tps/utils.h>
#include <algorithm>
#include <iostream>

namespace tsc = thor_scsi::core;

tsc::TwoDimensionalMultipoles::TwoDimensionalMultipoles(const unsigned int h_max)
	: m_max_multipole(h_max)
	, coeffs(h_max, (0.0, 0.0) )
{
	if(h_max<=1){
		throw std::logic_error("max multipole must be at least 1");
	}
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

tsc::TwoDimensionalMultipoles::TwoDimensionalMultipoles(const tsc::TwoDimensionalMultipoles&& o)
	:  m_max_multipole(o.m_max_multipole)
	 , coeffs(std::move(o.coeffs))
{
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
	: m_max_multipole( coeffs.size() )
	, coeffs(coeffs)
{
	if(coeffs.size()<=1){
		throw std::logic_error("max multipole must be at least 1");
	}
}


tsc::TwoDimensionalMultipoles& tsc::TwoDimensionalMultipoles::right_multiply(const std::vector<double> &scale, const bool begnin)
{
	auto& c = this->coeffs;
	size_t n_coeff = c.size(), n_scale = scale.size();

	if(begnin){
		if(n_coeff > n_scale){
			// Check for exception matching python length error
			std::stringstream strm;
			strm << "Received " << n_scale << " elements for scaling."
			     << " being begnin I would ignore excess elements. but I got only "
			     << n_coeff     << " .";
			throw std::runtime_error(strm.str());
		}
	} else {
		if(n_coeff != n_scale){
			// Check for exception matching python length error
			std::stringstream strm;
			strm << "Received " << n_scale << " elements for scaling."
			     << " this does not match the number of coefficients " << n_coeff
			     << " .";
			throw std::runtime_error(strm.str());
		}
	}
	std::transform(c.begin(), c.end(), scale.begin(), c.begin(), std::multiplies<>());
	return *this;

}

tsc::TwoDimensionalMultipoles& tsc::TwoDimensionalMultipoles::right_add(const TwoDimensionalMultipoles &other, const bool begnin)
{

	auto& c = this->coeffs;
	size_t n_coeff = c.size(), n_add = other.size();

	if(begnin){
		if(n_coeff >  n_add){
			// Check for exception matching python length error
			std::stringstream strm;
			strm << "Received " << n_add << " elements for adding."
			     << " being begnin I would ignore excess elements. But I got only " << n_coeff
			     << " .";
			throw std::runtime_error(strm.str());
		}
	} else {
		if(n_coeff != n_add){
			// Check for exception matching python length error
			std::stringstream strm;
			strm << "Received " << n_add << " elements for adding."
			     << " this does not match the number of coefficients " << n_coeff
			     << ".";
			throw std::runtime_error(strm.str());
		}
	}
	auto& oc = other.getCoeffs();

	std::transform(c.begin(), c.end(), oc.begin(), c.begin(), std::plus<>());

	return *this;

}

tsc::TwoDimensionalMultipoles& tsc::TwoDimensionalMultipoles::operator *= (const std::vector<double> &scale)
{
	// should be false ... leave it that way to get further with prototyping
	return this->right_multiply(scale, true);
}

tsc::TwoDimensionalMultipoles  tsc::TwoDimensionalMultipoles::operator * (const std::vector<double> &scale) const
{
	TwoDimensionalMultipoles nh = this->clone();
	nh *= scale;
	return nh;
}

tsc::TwoDimensionalMultipoles& tsc::TwoDimensionalMultipoles::operator += (const tsc::TwoDimensionalMultipoles &other)

{
	// should be false ... leave it that way to get further with prototyping
	return this->right_add(other, true);

}


tsc::TwoDimensionalMultipoles tsc::TwoDimensionalMultipoles::operator+ (const tsc::TwoDimensionalMultipoles &other) const
{

	// unsigned int n = std::max(this->getCoeffs().size(), other.getCoeffs().size());
	std::cerr << "Adding multipoles. this size " << this->getCoeffs().size()
		  << " others size" << other.getCoeffs().size()
		  << std::endl;
	TwoDimensionalMultipoles nh = this->clone();
	nh += other;
	return nh;

}


tsc::TwoDimensionalMultipoles& tsc::TwoDimensionalMultipoles::operator += (const double other)
{

	auto& c = this->getCoeffs();
	for (auto& val : c){
		val += other;
	}
	return *this;
}

tsc::TwoDimensionalMultipoles tsc::TwoDimensionalMultipoles::operator+ (const double other) const
{

	TwoDimensionalMultipoles nh = this->clone();
	std::cerr << "Adding const "<< other << "to these multipoles. this size " << this->getCoeffs().size()
		  << std::endl;
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

tsc::BegninTwoDimensionalMultipolesDelegator tsc::BegninTwoDimensionalMultipolesDelegator::clone(void) const
{
	tsc::TwoDimensionalMultipoles nh = this->delegate->clone();
        std::shared_ptr<tsc::TwoDimensionalMultipoles> ptr(nh);
	tsc::BegninTwoDimensionalMultipolesDelegator nd(ptr);
	return nd;
}

tsc::BegninTwoDimensionalMultipolesDelegator& tsc::BegninTwoDimensionalMultipolesDelegator::operator += (const tsc::TwoDimensionalMultipoles &other)
{
	this->delegate->right_add(other, true);
	return *this;

}

tsc::BegninTwoDimensionalMultipolesDelegator& tsc::BegninTwoDimensionalMultipolesDelegator::operator *= (const std::vector<double> &scale)
{
	this->delegate->right_multiply(scale, true);
	return *this;
}



tsc::BegninTwoDimensionalMultipolesDelegator tsc::BegninTwoDimensionalMultipolesDelegator::operator * (const std::vector<double> &scale) const
{
	auto nh = this->clone();
	nh *= scale;
	return nh;
}


tsc::BegninTwoDimensionalMultipolesDelegator tsc::BegninTwoDimensionalMultipolesDelegator::operator + (const tsc::TwoDimensionalMultipoles &other) const
{
	auto nh = this->clone();
	nh += other;
	return nh;
}

#endif

/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
