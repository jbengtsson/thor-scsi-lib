#ifndef _THOR_SCSI_CORE_MULTIPOLES_H_
#define _THOR_SCSI_CORE_MULTIPOLES_H_ 1
#include <vector>
#include <complex>
#include <cassert>
#include <stdexcept>
#include <ostream>

#include <thor_scsi/core/field_interpolation.h>
#include <thor_scsi/core/exceptions.h>

// #define THOR_SCSI_USE_F128 1
#ifdef THOR_SCSI_USE_F128
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/complex128.hpp>
#endif // THOR_SCSI_USE_F128

namespace thor_scsi::core {

	/**
	 * maximum multipole to use
	 *
	 * currently test_multipoles pass for max_multipole of 6
	 * otherwise the rather small tolerances have to be made more wide
	 */
	const int max_multipole = 6;

	typedef std::complex<double> cdbl;

#ifdef THOR_SCSI_USE_F128
	typedef boost::multiprecision::complex128 cdbl_intern;
#else // THOR_SCSI_USE_F128
	typedef std::complex<double> cdbl_intern;
#endif // THOR_SCSI_USE_F128



        /**
	 *  Representation of planar 2D harmonics / multipoles
	 *
	 *  @f[
	 *     \mathbf{B(z)} = B_y + i B_x =
	 *	\sum\limits_{n=1}^N \mathbf{C}_n
	 *	\left(\frac{\mathbf{z}}{R_\mathrm{ref}}\right)^{(n-1)}
	 *  @f]
	 *
	 *  with \f$ \mathbf{z} = x + i y \f$ and \f$ \mathbf{C}_n = B_n + i A_n \f$
	 *  and \f$R_\mathrm{ref}\f$ the reference radius.
	 *  \f$N\f$ corresponds to the maximum multipole
	 *
	 * \verbatim embed:rst:leading-asterisk
	 *
	 *  Please note: for the methods
	 *
	 *      - setMultipole
	 *      - getMultipoles
	 *
	 *  the class adheres to the European convention (i.e. dipole = 1, quadrupole = 2, ..).
	 *  All other methods do not require direct indexing.
	 *
	 *  The coefficients are internally represented by a standard vector, Thus
	 *  its index needs to be reduced by one.
 	 *
	 * .. Todo::
	 *       minimise number of coefficients that require that to be evaluated
	 *       rename to  TwoDimensionalMultipoles
	 * \endverbatim
	 */

	class PlanarMultipoles : public Field2DInterpolation {
	public:
		/**
		   Just allocates an set of zero multipoles
		 */
		inline PlanarMultipoles(const unsigned int h_max=max_multipole){
			if(h_max<=1){
				throw std::logic_error("max multipole must be at least 1");
			}
			this->m_max_multipole = h_max;
			this->coeffs.resize(h_max);
			const cdbl_intern zero(0.0, 0.0);
			for(auto h : this->coeffs){
				h = zero;
			}
		};
		virtual inline ~PlanarMultipoles(void){};
		PlanarMultipoles(PlanarMultipoles&& o):
			coeffs(std::move(o.coeffs)){this->m_max_multipole = o.m_max_multipole;};

		PlanarMultipoles(const PlanarMultipoles& o):
			coeffs(o.coeffs){this->m_max_multipole = o.m_max_multipole;};

		inline PlanarMultipoles clone(void) const {
			return PlanarMultipoles(std::vector<cdbl_intern>(this->coeffs));
		}

	private:
		/**
		 * \verbatim embed:rst:leading-asterisk
		 * Todo:
		 *     memory handling of the coefficient array
		 *
		 *  Args:
		 *     coeffs: a vector of complex coefficients
		 * \endverbatim
		 */
		inline PlanarMultipoles(std::vector<cdbl_intern> const coeffs) {
			if(coeffs.size()<=1){
				throw std::logic_error("max multipole must be at least 1");
			}
			this->coeffs = coeffs;
			this->m_max_multipole = this->coeffs.size();
		}

	public:
		/**
		 * @brief compute the field at position z
		 * Uses Horner equation
		 *
		 * \verbatim embed:rst:leading-asterisk
		 * Args:
		 *      z: position (x + I y) at which to implement the position
		 *
		 * Returns:
		 *      complex field By + I * Bx
		 * \endverbatim
		 */

#if 0
		inline cdbl field(const cdbl z) const {
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

		virtual inline void field(const double x, const double y, double *Bx, double * By) const override final{
			const cdbl z(x, y);
			const cdbl r = field(z);
			*By = r.real();
			*Bx = r.imag();
		}

#endif
		inline cdbl field(const cdbl z) const {
			double Bx=NAN, By=NAN;
			this->field(z.real(), z.imag(), &Bx, &By);
			cdbl tmp(By, Bx);
			return tmp;
		}

		template<typename T>
		inline void _field(const T x, const T y, T *Bx, T * By)  const {
			int n = this->coeffs.size() -1;
			T ByoBrho = this->coeffs[n].real(),  BxoBrho = this->coeffs[n].imag();
			for(int i=n - 2; i >= 0; --i){
				cdbl_intern tmp = this->coeffs[i];
				T ByoBrho1 = x * ByoBrho - y * BxoBrho + tmp.real();
				BxoBrho    = y * ByoBrho + x * BxoBrho + tmp.imag();
				ByoBrho = ByoBrho1;
			}
			*Bx = (BxoBrho);
			*By = ByoBrho;
			/*
			  MB[HOMmax-Order]
				ByoBrho = MB[Order+HOMmax]; BxoBrho = MB[HOMmax-Order];
    for (j = Order-1; j >= 1; j--) {
      ByoBrho1 = x0[x_]*ByoBrho - x0[y_]*BxoBrho + MB[j+HOMmax];
      BxoBrho  = x0[y_]*ByoBrho + x0[x_]*BxoBrho + MB[HOMmax-j];
      ByoBrho  = ByoBrho1;
    }
			*/

		}
		virtual inline void field(const double x, const double y, double *Bx, double * By)  const override final {
			_field(x, y, Bx, By);
		}
		virtual inline void field(const tps x, const tps y, tps *Bx, tps *By) const override final{
			// "Need to understand how to interpolate field with tps"
			_field(x, y, Bx, By);
		}
		virtual inline void gradient(const tps x, const tps y, tps *Gx, tps *Gy) const override final{
			// "Need to understand how to interpolate gradient with tps"
			throw thor_scsi::NotImplemented("Multipoles: gradient in tps not implemented");
		}
		virtual inline void gradient(const tps x, const tps y, double *Gx, double *Gy) const override final{
			// "Need to understand how to interpolate gradient with tps"
		        // throw thor_scsi::NotImplemented("Computing gradient in doubles for coordinates in tps x and tps y not implemented");
			const cdbl pos(0, 0);
			const cdbl tmp = this->gradient(pos);
			*Gy = tmp.real();
			*Gx = tmp.imag();
		}
		/**
		 *
		 * Todo: Gradient here
		 */
		inline cdbl gradient(const cdbl unused) const {
			return this->getMultipole(2);
		}
 		virtual inline void gradient(const double x, const double y, double *Gx, double *Gy) const override final{
			cdbl tmp = this->gradient(cdbl(x,y));
			*Gy = tmp.real();
			*Gx = tmp.imag();
		}
		/** @brief Check if multipole index is within range of representation

		    \verbatim embed:rst

		    .. todo::
		        Raise an exception std::range_error

		    \endverbatim
		 */
		inline void checkMultipoleIndex(const unsigned int n) const {
			if(n <= 0){
				// European convention
				throw std::length_error("multipoles index <= 0");
				assert(0);
			}else if (n > this->m_max_multipole){
				throw std::length_error("multipoles index >  max multipole");
				assert(0);
			}
		}

		/** @brief get n'th multipole

		    \verbatim embed:rst
		    .. todo::
		        Review if the multipole index check can be avoided as
			std::vector implements the required check
		    \endverbatim
		 */
		inline cdbl getMultipole(const unsigned int n) const {
			this->checkMultipoleIndex(n);
			unsigned  use_n = n - 1;
			// array indices matche the american notation
			// assert(use_n >= 0);
			assert(use_n <this->m_max_multipole);
			cdbl_intern c = this->coeffs[use_n];
			cdbl res(double(c.real()), double(c.imag()));
			return res;
		}

		/**
		 * @brief set n'th multipole
		 */
		inline void setMultipole(const unsigned int n, const cdbl c){
			this->checkMultipoleIndex(n);
			unsigned use_n = n - 1;
			// assert(use_n >= 0);
			assert(use_n <this->m_max_multipole);
			this->coeffs[use_n] = cdbl_intern(c.real(), c.imag());
		}

		/**
		 * @brief: applys a roll angle to the coordinate system z
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
		 * Apply following translation to the coordinate system of the multipoles
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
		 *       The  accuracy of complex calculation was checked (see test). It seems that
		 *       loss of accuracy is significant. It is suggested to transform the point to the
		 *       coordinate system of the multipoles instead.
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
		 *  scale all multipoles by a factor scale in place
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * .. todo::
		 *        Review if complex scale makes sense
		 *
		 * \endverbatim
		 */
		inline PlanarMultipoles& operator *= (const double scale){
			for(unsigned int i = 0; i< this->coeffs.size(); ++i){
				this->coeffs[i] *= scale;
			}
			return *this;
		}

		/**
		 * scale multipoles and return a new object.
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * .. todo::
		 *        Implementation correct ?
		 * \endverbatim
		 */
		inline PlanarMultipoles operator * (const double scale) const {
			PlanarMultipoles n = this->clone();
			n *= scale;
			return n;
		}
		inline PlanarMultipoles operator / (const double scale) const {
			return *this * (1.0/scale);
		}
		/**
		 * add
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * Add an other set of multipoles to this set.
		 *
		 * .. todo::
		 *      *  Implementation correct ?
		 *      *  Review behaviour if representations of different size
		 * \endverbatim
		 */
		PlanarMultipoles operator +(PlanarMultipoles &other) const;

		/**
		 * The maximum index represented.
		 */
		inline unsigned int getMultipoleMaxIndex(void){
			return this->m_max_multipole;
		}

		/**
		 * @brief access to the coefficents.
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
		/**
		 * @brief: access to the coefficents (const version)
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * .. todo::
		 *    is that here sufficient to make the user clear that these
		 *    coefficients should be considered constant
		 *
		 * \endverbatim
		 */
		inline const std::vector<cdbl_intern> getCoeffs(void) const {
			return this->coeffs;
		}


		/**
		 * @brief: access to coeffcs for python
		 *
		 * pybind11 also seems to have limits for disambuigity
		 */
		inline const std::vector<cdbl_intern> getCoeffsConst(void) const {
			return this->getCoeffs();
		}

virtual void show(std::ostream& strm, int level) const override final;

	private:
		unsigned int m_max_multipole;
		std::vector<cdbl_intern> coeffs;


	};

#if 0
	/// convert representation to tracy representation
	std::vector<double> toTracyRepresentation(const std::vector<cdbl> &coeffs);


	// why should one need that ...
	// void toMultipoleRepresentation(std::vector<cdbl> vec);

	// consider if caching the coefficients internally is a good idea
	// If all code constantly access them, it would be a good idea
	// std::vector<<double> coeffs_tracy_representation;
#endif
}
#endif /* _THOR_SCSI_MULTIPOLES_H_ */
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
