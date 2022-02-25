#include <thor_scsi/core/harmonics.h>
#include <thor_scsi/core/math_comb.h>

namespace tsc = thor_scsi::core;


tsc::PlanarHarmonics tsc::PlanarHarmonics::operator + (tsc::PlanarHarmonics &other)
{
	unsigned int n = std::max(this->getCoeffs().size(), other.getCoeffs().size());
	unsigned int i;

	PlanarHarmonics nh(n);
	auto ncoeffs = nh.getCoeffs();

	i = 0;
	for(auto c : this->getCoeffs()){
		ncoeffs[i] = c;
		++i;
	}

	i = 0;
	for(auto c : other.getCoeffs()){
		ncoeffs[i] = c;
		++i;
	}

	return nh;
}

#include <gsl/gsl_complex.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_vector_complex_double.h>
#include <stdlib.h>
#include <iostream>
tsc::cdbl tsc::PlanarHarmonics::field_taylor(const tsc::cdbl z)
{
	tsc::cdbl_intern zp(z.real(), z.imag()), ztmp = zp,
		field(0.0, 0.0);
	
	for(size_t i = 0; i< this->coeffs.size(); ++i){
		if(i == 0){
			field += this->coeffs[i];
		} else {
			field += this->coeffs[i] * zp;
			zp *= ztmp;
		}
	}
	tsc::cdbl result(double(field.real()),double(field.imag()));
	return result;
}

tsc::cdbl tsc::PlanarHarmonics::field_gsl(const tsc::cdbl z)
{

	gsl_complex gz, *mem=NULL;
	GSL_SET_COMPLEX(&gz, z.real(), z.imag());

	auto n = this->coeffs.size();
	gsl_complex *vec = NULL;
	
	vec = (gsl_complex *) calloc(n, sizeof(gsl_vector)); 
	if(!vec){
		throw std::logic_error("No vec allocated!");
	}
	for(size_t i= 0; i<n; ++i){
		auto tmp = this->coeffs[i];
		GSL_SET_COMPLEX(&vec[i], double(tmp.real()), double(tmp.imag()));
	}

	gsl_complex gr = gsl_complex_poly_complex_eval((const gsl_complex*) vec, int(n), gz);
	free(vec);
	return tsc::cdbl(GSL_REAL(gr), GSL_IMAG(gr));
	
}

/**
 *
 *
 */
void tsc::PlanarHarmonics::applyTranslation(const tsc::cdbl dzs)
{

	for (unsigned int i = 0; i < this->coeffs.size(); ++i) {
		cdbl_intern dzi(dzs.real(), dzs.imag());
		for (unsigned j = i+1; j < this->coeffs.size(); ++j) {
			this->coeffs[i] +=  double(binom(j, i)) * this->coeffs[j] * dzi;
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
