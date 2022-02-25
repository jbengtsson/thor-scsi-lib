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
