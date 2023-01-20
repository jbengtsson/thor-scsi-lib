#include <thor_scsi/core/multipoles.h>
#include <thor_scsi/core/math_comb.h>
#include <tps/utils.h>
#include <algorithm>
#include <iostream>

namespace tsc = thor_scsi::core;

namespace thor_scsi::core{

}

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

//tsc::TwoDimensionalMultipolesKnobbed& tsc::TwoDimensionalMultipolesKnobbed::right_multiply(const std::vector<double> &scale, const bool begnin)

//tsc::TwoDimensionalMultipolesKnobbed& tsc::TwoDimensionalMultipolesKnobbed::right_add(const TwoDimensionalMultipolesKnobbed &other, const bool begnin)

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
