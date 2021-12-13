#ifndef THOR_SCSI_HARMONICS_H
#define THOR_SCSI_HARMONICS_H
#include <vector>
#include <complex>

namespace thor_scsi {
#include <cassert>
	const int max_harmonic = 21;

	typedef std::complex<double> cdbl;
	// pure virtual class of 2D Harmonics
	// using European convention
	class HarmonicsBase{
	};


	class Harmonics : public HarmonicsBase{
	public:

		Harmonics(const unsigned int h_max=max_harmonic){
			this->m_max_harmonic = h_max;
			this->coeffs.resize(h_max);
		};
		inline void checkHarmonicIndex(const unsigned int n){
			if(n == 0){
				// European convention
				assert(0);
			}else if (n > this->m_max_harmonic){
				assert(0);
			}
		}
		inline cdbl getHarmonic(const unsigned int n){
			this->checkHarmonicIndex(n);
			unsigned  use_n = n - 1;
			assert(use_n >=0);
			assert(use_n <this->m_max_harmonic);
			return this->coeffs[use_n];
		}

		inline void setHarmonic(const cdbl c, const unsigned int n){
			this->checkHarmonicIndex(n);
			unsigned use_n = n - 1;
			assert(use_n >=0);
			assert(use_n <this->m_max_harmonic);
			this->coeffs[use_n] = c;
		}

		// In place operation.
		inline void applyRollAngle(const double angle){
			cdbl I (0, 1);
			double phase = 0.0;
			for(unsigned int i = 0; i < this->coeffs.size(); ++i){
				phase += angle;
				cdbl scale = exp(I * phase);
				this->coeffs[i] *= scale;
			}
		}

		// In place operation. best provided as separate function
		// I consider it too long for inplace operation
		// dz = dx + I * dy
		void applyTranslation(const cdbl dz);

		inline Harmonics& operator *= (double scale){
			for(unsigned int i = 0; i<=this->coeffs.size(); ++i){
				this->coeffs.at(i) *= scale;
			}
			return *this;
		}

		inline Harmonics operator * (double scale){
			Harmonics nh = *this;
			return nh*scale;
		}
		Harmonics operator +(Harmonics &other);

		// Missing: add plus and minus etc ....
		inline unsigned int getMaxHarmonic(void){
			return this->m_max_harmonic;
		}

		//protected:
		inline std::vector<cdbl> *getCoeffs(void){
			return &this->coeffs;
		}

	private:
		unsigned int m_max_harmonic;
		std::vector<cdbl> coeffs;


	};

	// calculate tracy representation from classical multipoles
	std::vector<double> toTracyRepresentation(const std::vector<cdbl> &coeffs);


	// why shoud one need that ...
	// void toHarmonicRepresentation(std::vector<cdbl> vec);

	// consider if caching the coefficients internally is a good idea
	// If all code constantly access them, it would be a good idea
	// std::vector<<double> coeffs_tracy_representation;

}
#endif /* THOR_SCSI_HARMONICS_H */


/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
