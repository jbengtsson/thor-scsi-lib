#include <thor_scsi/core/multipoles.h>
#include <thor_scsi/core/math_comb.h>
#include <tps/utils.h>
#include <algorithm>
#include <iostream>

namespace tsc = thor_scsi::core;

#if 0
// Why do I need a copy constructor ??
tsc::TwoDimensionalMultipolesKnobbed::TwoDimensionalMultipolesKnobbed(const tsc::TwoDimensionalMultipolesKnobbed& o): coeffs(o.coeffs){
	this->m_max_multipole = o.m_max_multipole;
}
#endif

/*
template
tsc::TwoDimensionalMultipolesKnobbed::TwoDimensionalMultipolesKnobbed<double>(const tsc::TwoDimensionalMultipolesKnobbed<double>&& o);
*/

/*
tsc::TwoDimensionalMultipolesKnobbed& tsc::TwoDimensionalMultipolesKnobbed::operator= (const tsc::TwoDimensionalMultipolesKnobbed &o)
{
	this->coeffs = o.coeffs;
	this->m_max_multipole = o.m_max_multipole;
	return *this;
}
*/
/*
tsc::TwoDimensionalMultipolesKnobbed tsc::TwoDimensionalMultipolesKnobbed::clone(void) const
{
	return tsc::TwoDimensionalMultipolesKnobbed(std::vector<tsc::cdbl_intern>(this->coeffs));
}
*/

/*
tsc::TwoDimensionalMultipolesKnobbed::TwoDimensionalMultipolesKnobbed(std::vector<cdbl_intern> const coeffs)
	: m_max_multipole( coeffs.size() )
	, coeffs(coeffs)
{
	if(coeffs.size()<=1){
		throw std::logic_error("max multipole must be at least 1");
	}
}
*/


template<typename T>
void tsc::right_multiply_helper(const std::vector<T> &scale, const bool begnin, std::vector<T>* coeffs)
{
	auto& c = *coeffs;
	size_t n_coeff = c.size(), n_scale = scale.size();

	if (begnin) {
                if (n_coeff > n_scale) {
			// Check for exception matching python length error
			std::stringstream strm;
			strm << "Received " << n_scale << " elements for scaling."
			     << " being begnin I would ignore excess elements. but I need at least "
			     << n_coeff << " .";
			throw std::runtime_error(strm.str());
                }
	} else {
                if (n_coeff != n_scale) {
			// Check for exception matching python length error
			std::stringstream strm;
			strm << "Received " << n_scale << " elements for scaling."
			     << " this does not match the number of coefficients " << n_coeff
			     << " .";
			throw std::runtime_error(strm.str());
                }
	}
	std::transform(c.begin(), c.end(), scale.begin(), c.begin(), std::multiplies<>());
}

template<typename T>
void tsc::right_add_helper(const std::vector<T>& other, const bool begnin, std::vector<T>* coeffs)
{

	auto &c = *coeffs;
	size_t n_coeff = c.size(), n_add = other.size();

	if (begnin) {
		if (n_coeff > n_add) {
                        // Check for exception matching python length error
                        std::stringstream strm;
                        strm << "Received " << n_add << " elements for adding."
                             << " being begnin I would ignore excess elements. But I got only " << n_coeff
                             << " .";
                        throw std::runtime_error(strm.str());
		}
	} else {
		if (n_coeff != n_add) {
                        // Check for exception matching python length error
                        std::stringstream strm;
                        strm << "Received " << n_add << " elements for adding."
                             << " this does not match the number of coefficients " << n_coeff
                             << ".";
                        throw std::runtime_error(strm.str());
		}
	}

	std::transform(c.begin(), c.end(), other.begin(), c.begin(), std::plus<>());
}

template<typename T>
void tsc::multipoles_show(std::ostream& strm, const std::vector<T>& coeffs, const int level)
{
	bool debug_print_made = false;
	strm << "interpolation='2D multipoles'";
	if (level > 2){
                strm << "(num=" << coeffs.size();
	}
	if (level  > 3){
                strm << ", muls={";
                int i = 1;
                for(auto c : coeffs){
			if(i > 1){
				strm << ", ";
			}
			/* deactivated for variant development */
			if(!debug_print_made) {
				debug_print_made = true;
			}
			strm << i << ":[" << c.real() << ", " << c.imag() << "]";
			++i;
                }
                strm << "}";
	}
	strm<<")";
}

template
void tsc::right_add_helper<std::complex<double>>(const std::vector<std::complex<double>>& other, const bool begnin,
						 std::vector<std::complex<double>>* coeffs);
template
void tsc::right_multiply_helper<std::complex<double>>(const std::vector<std::complex<double>>& other, const bool begnin,
						      std::vector<std::complex<double>>* coeffs);
template
void tsc::right_add_helper<gtpsa::CTpsaOrComplex>(const std::vector<gtpsa::CTpsaOrComplex>& other, const bool begnin,
						 std::vector<gtpsa::CTpsaOrComplex>* coeffs);
template
void tsc::right_multiply_helper<gtpsa::CTpsaOrComplex>(const std::vector<gtpsa::CTpsaOrComplex>& other, const bool begnin,
						 std::vector<gtpsa::CTpsaOrComplex>* coeffs);

template
void tsc::multipoles_show<std::complex<double>>(std::ostream& strm, const std::vector<std::complex<double>>& coeffs, const int level);
template
void tsc::multipoles_show<gtpsa::CTpsaOrComplex>(std::ostream& strm, const std::vector<gtpsa::CTpsaOrComplex>& coeffs, const int level);

//tsc::TwoDimensionalMultipolesKnobbed& tsc::TwoDimensionalMultipolesKnobbed::

//tsc::TwoDimensionalMultipolesKnobbed& tsc::TwoDimensionalMultipolesKnobbed::operator *= (const std::vector<double> &scale)

//tsc::TwoDimensionalMultipolesKnobbed  tsc::TwoDimensionalMultipolesKnobbed::operator * (const std::vector<double> &scale) const



//tsc::TwoDimensionalMultipolesKnobbed& tsc::TwoDimensionalMultipolesKnobbed::operator += (const tsc::TwoDimensionalMultipolesKnobbed &other)

//tsc::TwoDimensionalMultipolesKnobbed tsc::TwoDimensionalMultipolesKnobbed::operator+ (const tsc::TwoDimensionalMultipolesKnobbed &other) const


//tsc::TwoDimensionalMultipolesKnobbed& tsc::TwoDimensionalMultipolesKnobbed::operator += (const double other)

//tsc::TwoDimensionalMultipolesKnobbed tsc::TwoDimensionalMultipolesKnobbed::operator+ (const double other) const

//void tsc::TwoDimensionalMultipolesKnobbed::show(std::ostream& strm, int level) const

//void tsc::TwoDimensionalMultipolesKnobbed::applyTranslation(const tsc::cdbl dzs)




#if 0

tsc::BegninTwoDimensionalMultipolesDelegator tsc::BegninTwoDimensionalMultipolesDelegator::clone(void) const
{
	tsc::TwoDimensionalMultipolesKnobbed nh = this->delegate->clone();
        std::shared_ptr<tsc::TwoDimensionalMultipolesKnobbed> ptr(nh);
	tsc::BegninTwoDimensionalMultipolesDelegator nd(ptr);
	return nd;
}

tsc::BegninTwoDimensionalMultipolesDelegator& tsc::BegninTwoDimensionalMultipolesDelegator::operator += (const tsc::TwoDimensionalMultipolesKnobbed &other)
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


tsc::BegninTwoDimensionalMultipolesDelegator tsc::BegninTwoDimensionalMultipolesDelegator::operator + (const tsc::TwoDimensionalMultipolesKnobbed &other) const
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
