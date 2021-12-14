#include "harmonics.h"
#include "math_comb.h"
using namespace thor_scsi;



PlanarHarmonics PlanarHarmonics::operator + (PlanarHarmonics &other)
{
	unsigned int n = std::max(this->getCoeffs()->size(), other.getCoeffs()->size());
	unsigned int i;

	PlanarHarmonics nh(n);
	auto *ncoeffs = nh.getCoeffs();

	i = 0;
	for(auto c : *(this->getCoeffs())){
		ncoeffs->at(i) = c;
		++i;
	}

	i = 0;
	for(auto c : *(other.getCoeffs())){
		ncoeffs->at(i) = c;
		++i;
	}

	return nh;
}

/**
 *  
 *
 */
void PlanarHarmonics::applyTranslation(const cdbl dzs)
{

	auto coeffs = this->getCoeffs();

	for (unsigned int i = 0; i < coeffs->size(); ++i) {
		cdbl dzi = dzs;
		for (unsigned j = i+1; j < coeffs->size(); ++j) {
			coeffs->at(i) +=  double(binom(j, i)) * coeffs->at(j) * dzi;
			dzi *= dzs;
		}
	}
}


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


/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
