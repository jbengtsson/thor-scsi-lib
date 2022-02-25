#ifndef _THOR_SCSI_CORE_HARMONICS_H_
#define _THOR_SCSI_CORE_HARMONICS_H_ 1
#include <vector>
#include <complex>
#include <cassert>
#include <stdexcept>

// #define THOR_SCSI_USE_F128 1
#ifdef THOR_SCSI_USE_F128
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/complex128.hpp>
#endif // THOR_SCSI_USE_F128
 
namespace thor_scsi::core {

	const int max_harmonic = 21;

	typedef std::complex<double> cdbl;
#ifdef THOR_SCSI_USE_F128
	typedef boost::multiprecision::complex128 cdbl_intern;
#else // THOR_SCSI_USE_F128	
	typedef std::complex<double> cdbl_intern;
#endif // THOR_SCSI_USE_F128
	
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
	 *
	 *  Todo: rename to multipoles
	 */
	class PlanarHarmonics : public HarmonicsBase{
	public:
		/**
		   Just allocates an set of zero harmonics
		 */
		inline PlanarHarmonics(const unsigned int h_max=max_harmonic){
			if(h_max<=1){
				throw std::logic_error("max harmonic must be at least 1");
			}
			this->m_max_harmonic = h_max;
			this->coeffs.resize(h_max);
			const cdbl_intern zero(0.0, 0.0);
			for(auto h : this->coeffs){
				h = zero;
			}
		};
		
		PlanarHarmonics(PlanarHarmonics&& o):
			coeffs(std::move(o.coeffs)){this->m_max_harmonic = o.m_max_harmonic;};

		PlanarHarmonics(PlanarHarmonics& o):
			coeffs(o.coeffs){this->m_max_harmonic = o.m_max_harmonic;};

		inline PlanarHarmonics clone(void) const {
			return PlanarHarmonics(std::vector<cdbl_intern>(this->coeffs));
		}

	private:
		/**
		 * Todo: memory handling!
		 */
		inline PlanarHarmonics(std::vector<cdbl_intern> const coeffs) {
			if(coeffs.size()<=1){
				throw std::logic_error("max harmonic must be at least 1");
			}
			this->coeffs = coeffs;
			this->m_max_harmonic = this->coeffs.size();
		}
		
	public:
		/**
		 * @brief compute the field at position z
		 *
		 * Uses Horner equation
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 *   .. todo::
		 *       check accuracy of complex calculation
		 *
		 *  Todo : cross check with GSL ...
		 * \endverbatim		 
		 */
		cdbl field_gsl(const cdbl z);
		cdbl field_taylor(const cdbl z);
		inline cdbl field3(const cdbl z){
			return field_gsl(z);
		}
		inline cdbl field2(const cdbl z){
			return field_taylor(z);
		}
		inline cdbl field(const cdbl z){
			int n = this->coeffs.size() -1;
			cdbl_intern t_field = this->coeffs[n];
			cdbl_intern z_tmp(z.real(), z.imag());
			for(int i=n - 2; i >= 0; --i){
				cdbl_intern tmp = this->coeffs[i];
				t_field = t_field * z_tmp + tmp;
			}
			cdbl result(double(t_field.real()), double(t_field.imag()));
			return result;
		}
		/*
		inline cdbl field(const double x, const double y){
			const cdbl z(x, y);
			return field(z);
		}

		inline cdbl field(const double x){
			return field(x, 0);
		}
		*/
		/** Check if harmonic index is within range of representation

		    \verbatim embed:rst

		    .. todo::
		        Raise an exception std::range_error

		    \endverbatim
		 */
		inline void checkHarmonicIndex(const unsigned int n){
			if(n <= 0){
				// European convention
				throw std::length_error("harmonics index <= 0");
				assert(0);
			}else if (n > this->m_max_harmonic){
				throw std::length_error("harmonics index >  max harmonic");
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
			assert(use_n > 0);
			assert(use_n <this->m_max_harmonic);
			cdbl_intern c = this->coeffs[use_n];
			cdbl res(double(c.real()), double(c.imag()));
			return res;
		}

		/** set n'th harmonic
		 */
		inline void setHarmonic(const unsigned int n, const cdbl c){
			this->checkHarmonicIndex(n);
			unsigned use_n = n - 1;
			assert(use_n > 0);
			assert(use_n <this->m_max_harmonic);
			this->coeffs[use_n] = cdbl_intern(c.real(), c.imag());
		}

		/**
		 * applys a roll angle to the coordinate system z
		 *
		 *  @f[
		 *  \mathbf{C´}_n =  \mathbf{C}_n \exp{(I n \alpha)}
		 *  @f]
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 *   .. todo::
		 *       check accuracy of complex calculation
		 *       check if phase should not be -alpha ....
		 * \endverbatim
		 *
		 */
		inline void applyRollAngle(const double alpha){
			cdbl_intern I (0, 1);
			for(size_t i = 0; i < this->coeffs.size(); ++i){
				const double phase = alpha * (i + 1);
				cdbl_intern scale = exp(I * phase);
				this->coeffs[i] *= scale;
			}
		}

		/**
		 * @brief translates coordiante system
		 *
		 * Apply following translation to the coordinate system of the harmonics
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
		 * .. todo::
		 *       check accuracy of complex calculation
		 *
		 * \endverbatim		 *
		 */
		void applyTranslation(const cdbl dzs);
		inline void applyTranslation(const double dx, const double dy){
			const cdbl tmp(dx, dy);
			applyTranslation(tmp);
		}
		inline void applyTranslation(const double dx){
			applyTranslation(dx, 0.0);
		}

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
		inline PlanarHarmonics& operator *= (const double scale){
			for(unsigned int i = 0; i< this->coeffs.size(); ++i){
				this->coeffs[i] *= scale;
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
		inline PlanarHarmonics operator * (const double scale) const{
			PlanarHarmonics n = this->clone();
			n *= scale;
			return n;
		}
		inline PlanarHarmonics operator / (const double scale) const{
			return *this * (1.0/scale);
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
		inline std::vector<cdbl_intern>& getCoeffs(void){
			return this->coeffs;
		}

	private:
		unsigned int m_max_harmonic;
		std::vector<cdbl_intern> coeffs;


	};

#if 0
	/// convert representation to tracy representation
	std::vector<double> toTracyRepresentation(const std::vector<cdbl> &coeffs);


	// why should one need that ...
	// void toHarmonicRepresentation(std::vector<cdbl> vec);

	// consider if caching the coefficients internally is a good idea
	// If all code constantly access them, it would be a good idea
	// std::vector<<double> coeffs_tracy_representation;
#endif
}
#endif /* _THOR_SCSI_HARMONICS_H_ */
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
