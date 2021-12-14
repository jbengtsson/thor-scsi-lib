#ifndef THOR_SCSI_HARMONICS_H
#define THOR_SCSI_HARMONICS_H
#include <vector>
#include <complex>
#include <cassert>

namespace thor_scsi {
	const int max_harmonic = 21;

	typedef std::complex<double> cdbl;
	/// pure virtual class of 2D Harmonics
	/// using European convention
	class HarmonicsBase{
	};
        /**
	 *  Representation of planar 2D harmonics
	 *
	 *  @f[
	 *     \mathbf{B(z)} = B_y + i B_x =
	 *	\sum\limits_{n=1}^N \mathbf{C}_n
	 *	\left(\frac{\mathbf{z}}{R_\mathrm{ref}}\right)^{(n-1)}
	 *  @f]
	 *
	 *  with \f$ \mathbf{z} = x + i y \f$ and \f$ \mathbf{C}_n = B_n + i A_n \f$
	 *  and \f$R_\mathrm{ref}\f$ the reference radius.
	 *  \f$N\f$ corresponds to the maximum harmonic
	 *  Please note: the class adheres to the European convention
	 */
	class PlanarHarmonics : public HarmonicsBase{
	public:
		/**
		   Just allocates an set of zero harmonics
		 */
		inline PlanarHarmonics(const unsigned int h_max=max_harmonic){
			this->m_max_harmonic = h_max;
			this->coeffs.resize(h_max);
			const cdbl zero(0.0, 0.0);
			for(auto h : this->coeffs){
				h = zero;
			}
		};
		/** Check if harmonic index is within range of representation

		    \verbatim embed:rst

		    .. todo::
		        Raise an exception std::range_error

		    \endverbatim
		 */
		inline void checkHarmonicIndex(const unsigned int n){
			if(n == 0){
				// European convention
				assert(0);
			}else if (n > this->m_max_harmonic){
				assert(0);
			}
		}
		/** get n'th harmonic

		    \verbatim embed:rst
		    .. todo::
		        Review if the harmonic index check can be avoided as
			std::vector implements the required check
		    \endverbatim
		 */
		inline cdbl getHarmonic(const unsigned int n){
			this->checkHarmonicIndex(n);
			unsigned  use_n = n - 1;
			assert(use_n >=0);
			assert(use_n <this->m_max_harmonic);
			return this->coeffs[use_n];
		}

		/** set n'th harmonic
		 */
		inline void setHarmonic(const cdbl c, const unsigned int n){
			this->checkHarmonicIndex(n);
			unsigned use_n = n - 1;
			assert(use_n >=0);
			assert(use_n <this->m_max_harmonic);
			this->coeffs[use_n] = c;
		}

		/**
		   applys a roll angle to the harmonics

		   @f[
		   \mathbf{C´}_n =  \mathbf{C}_n \exp{(I n \alpha)}
		   @f]
		 */
		inline void applyRollAngle(const double alpha){
			cdbl I (0, 1);
			double phase = 0.0;
			for(unsigned int i = 0; i < this->coeffs.size(); ++i){
				phase += alpha;
				cdbl scale = exp(I * phase);
				this->coeffs[i] *= scale;
			}
		}

		/**
		 * Apply translation
		 * \f$ \mathbf{\Delta z}_\mathrm{s} = \mathbf{\Delta z}/R_\mathrm{rref} \f$
		 * with \f$ \mathbf{\Delta z} = \Delta x + i \Delta y \f$
		 *
		 * Implementation
		 *  @f[
		 *      \mathbf{C'}_n =  \sum_{j=n}^N \  {j \choose i} \mathbf{C}_j \ \mathbf{\Delta z}_\mathrm{s}^{j-n+1}
		 *  @f]
		 *
		 * \verbatim embed:rst:leading-asterisk		 *
		 * .. warning::
		 *        Code not yet checked
		 *
		 * \endverbatim		 *
		 */
		void applyTranslation(const cdbl dzs);

		/**
		 *  scale all harmonics by a factor scale in place
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * .. todo::
		 *        Review if complex scale makes sense
		 *
		 * \endverbatim
		 */
		inline PlanarHarmonics& operator *= (double scale){
			for(unsigned int i = 0; i<=this->coeffs.size(); ++i){
				this->coeffs.at(i) *= scale;
			}
			return *this;
		}

		/**
		 * scale harmonics and return a new object.
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * .. todo::
		 *        Implementation correct ?
		 * \endverbatim
		 */
		inline PlanarHarmonics operator * (double scale){
			PlanarHarmonics nh = *this;
			return nh*scale;
		}
		/**
		 * add
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * Add an other set of harmonics to this set.
		 *
		 * .. todo::
		 *      *  Implementation correct ?
		 *      *  Review behaviour if representations of different size
		 * \endverbatim
		 */
		PlanarHarmonics operator +(PlanarHarmonics &other);

		/**
		 * The maximum index represented.
		 */
		inline unsigned int getMaxHarmonicIndex(void){
			return this->m_max_harmonic;
		}

		/**
		 * access to the coefficents.
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * .. note::
		 *     access to the same memory. If you change the coefficients
		 *     outside you also change them here.
		 *
		 * \endverbatim
		 */
		inline std::vector<cdbl> *getCoeffs(void){
			return &this->coeffs;
		}

	private:
		unsigned int m_max_harmonic;
		std::vector<cdbl> coeffs;


	};

	/// convert representation to tracy representation
	std::vector<double> toTracyRepresentation(const std::vector<cdbl> &coeffs);


	// why should one need that ...
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
