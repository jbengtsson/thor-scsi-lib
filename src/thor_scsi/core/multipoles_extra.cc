#include <gsl/gsl_complex.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_vector_complex_double.h>
#include <thor_scsi/core/multipoles.h>

namespace tsc = thor_scsi::core;

tsc::cdbl tsc::PlanarMultipoles::field_gsl(const tsc::cdbl z)
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

tsc::cdbl tsc::PlanarMultipoles::field_taylor(const tsc::cdbl z)
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

