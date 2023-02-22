#ifndef _THOR_SCSI_CORE_MULTIPOLES_H_
#define _THOR_SCSI_CORE_MULTIPOLES_H_ 1

#include <algorithm>
#include <complex>
#include <cassert>
#include <memory>
#include <ostream>
#include <vector>
#include <stdexcept>
#include <thor_scsi/core/field_interpolation.h>
#include <thor_scsi/core/exceptions.h>
#include <thor_scsi/core/multipole_types.h>
// for binom
#include <tps/utils.h>
#include <gtpsa/utils.hpp>
#include <gtpsa/utils_tps.hpp>

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

	//typedef std::complex<double> cdbl;


    template<typename Tcc, typename Tcp> // a complex type
    inline auto honer_complex(const Tcc& coeffs, const Tcp& z) {
        size_t n = coeffs.size() -1;
        // how to handle maximum order in case of ctpsa?
        auto rB = gtpsa::clone(coeffs[n]) * (1.0 + z * 0e0);
        for(int i=n - 1; i >= 0; --i) {
            rB = z * rB + coeffs[i];
        }
        return rB;
    }

    /**
     *
     * solely here for support of tps ... to be removed with tps
     */
    inline void honer_complex_as_real(const std::vector<std::complex<double>>& coeffs,
                                      const tps& x, const tps& y, tps *Bx, tps * By) {
        size_t n = coeffs.size() -1;
        // how to handle maximum order in case of ctpsa?
        tps rBy = coeffs[n].real();
        tps rBx = coeffs[n].imag();
        for(int i=n - 2; i >= 0; --i) {
            auto tmp = coeffs[i];
            auto trBy = x * rBy - y * rBx + tmp.real();
            rBx       = y * rBy + x * rBx + tmp.imag();
            rBy  = std::move(trBy);
        }
        *Bx = rBx;
        *By = rBy;
    }

    // these helper functions required to transmit the type used for the cofficients down to the functions
    // otherwise this code would require to be part of the class definition
    template<typename T>
    void right_multiply_helper(const std::vector<T> &scale, const bool begnin, std::vector<T>* coeffs);
    template<typename T>
    void right_add_helper(const std::vector<T> &other, const bool begnin, std::vector<T>* coeffs);

    template<typename T>
    void multipoles_show(std::ostream& strm, const std::vector<T>& coeffs, const int level);

	/**
	 *  @brief: Representation of planar 2D harmonics / multipoles
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
	 *
	 *       * minimise number of coefficients that require that to be evaluated
	 *       * separate operator (+, -) interface to a begnin and a strict one
	 *         (similar to pandas data frames .loc and .iloc)
	 * \endverbatim
	 */

    template<class C, typename = typename C::complex_type, typename = typename C::double_type>
	class TwoDimensionalMultipolesKnobbed : public Field2DInterpolationKnobbed<C> {
    protected:
        using double_type = typename C::double_type;
        using complex_type = typename C::complex_type;
        using complex_intern_type = typename C::complex_intern_type;

    public:
        /**
	     *  @brief Just allocates an set of multipoles of default value
         *
         *  default value is required for coefficients of tpsa type. these need to contain a
         *  shared ptr to a a common description object.
		 */
        TwoDimensionalMultipolesKnobbed(const complex_type& default_value, const unsigned int h_max=max_multipole)
                : m_max_multipole(h_max)
                , coeffs(h_max, default_value )
        {
            if(h_max<=1){
                throw std::logic_error("max multipole must be at least 1");
            }
            //this->coeffs.resize(h_max);
            for(auto h : this->coeffs){
                h = default_value;
            }
        }

        TwoDimensionalMultipolesKnobbed(std::vector<complex_type> const coeffs)
                : m_max_multipole( coeffs.size() )
                , coeffs(coeffs)
        {
            if(coeffs.size()<=1){
                throw std::logic_error("max multipole must be at least 1");
            }
        }

        virtual inline ~TwoDimensionalMultipolesKnobbed(void){};
		// Why do I need a copy constructor ?? Did I miss an assignment operator
		// TwoDimensionalMultipolesKnobbed(const TwoDimensionalMultipolesKnobbed& o);
        TwoDimensionalMultipolesKnobbed(const TwoDimensionalMultipolesKnobbed& o)
	    :  m_max_multipole(std::move(o.m_max_multipole))
	    , coeffs(o.coeffs)
	    {}

        TwoDimensionalMultipolesKnobbed(const TwoDimensionalMultipolesKnobbed&& o)
	    :  m_max_multipole(std::move(o.m_max_multipole))
                , coeffs(std::move(o.coeffs))
        {}

	// Required e.g. for engineering tolerance studies
	TwoDimensionalMultipolesKnobbed& operator= (const  TwoDimensionalMultipolesKnobbed &o)
        {
            this->coeffs = o.coeffs;
            this->m_max_multipole = o.m_max_multipole;
            return *this;
        }


        TwoDimensionalMultipolesKnobbed clone(void) const{
            return TwoDimensionalMultipolesKnobbed(std::vector<complex_intern_type>(this->coeffs));
        }

		inline size_t size(void) const {
			return this->coeffs.size();
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
        template<typename Tc>
        inline auto _cfield(const Tc z) const {
            return honer_complex(this->coeffs, z);
        }


#if 1
        inline auto field(const std::complex<double> z){
            return this->_cfield(z);
        }

        inline void _field(const double x, const double y, double *Bx, double *By) const {
            std::complex<double> z(x, y);
            auto tmp = this->_cfield<std::complex<double>>(z);
            // in case a power series is returned
            // review _cfield and honer implementation?
            std::complex<double> cst = gtpsa::cst(tmp);
            *By = cst.real();
            *Bx = cst.imag();
        }

        inline void _field(const gtpsa::tpsa& x, const gtpsa::tpsa& y, gtpsa::tpsa *Bx, gtpsa::tpsa *By) const {
            gtpsa::ctpsa z(x, y);
            auto tmp = this->_cfield(z);
            tmp.real(By);
            tmp.imag(Bx);
        }

        inline void _field(const tps& x, const tps& y, tps *Bx, tps *By) const {
            throw std::runtime_error("_field with tps arguments not implemented for all class template types");
        }

#else
        template<typename T>
		inline void _field(const T& x, const T& y, T *Bx, T * By)  const {
             T rBy  = gtpsa::same_as_instance(x);
			T rBx  = gtpsa::same_as_instance(y);
			rBy = this->coeffs[n].real();
			rBx = this->coeffs[n].imag();
			for(int i=n - 2; i >= 0; --i){
			        complex_intern_type tmp = this->coeffs[i];
				auto trBy = x * rBy - y * rBx + tmp.real();
				rBx       = y * rBy + x * rBx + tmp.imag();
				rBy  = std::move(trBy);
			}
			*Bx = std::move(rBx);
			*By = std::move(rBy);
		}
		template<typename T>
		inline void _field(const T& x, const T& y, T *Bx, T * By)  const {
			int n = this->coeffs.size() -1;
			// T rBy   = std::move( gtpsa::same_as_instance(x) );
			// T rBx   = std::move( gtpsa::same_as_instance(y) );
			T trBy  = std::move( gtpsa::same_as_instance(y) );
			T term1 = std::move( gtpsa::same_as_instance(x) );
			T term2 = std::move( gtpsa::same_as_instance(x) );

			rBy = this->coeffs[n].real();
			rBx = this->coeffs[n].imag();
			for(int i=n - 2; i >= 0; --i){
				cdbl_intern tmp = this->coeffs[i];
				trBy = x * rBy - y * rBx + tmp.real();
				rBx  = y * rBy + x * rBx + tmp.imag();
				rBy  = std::move(trBy);
			}
			*Bx = std::move(rBx);
			*By = std::move(rBy);
		}
#endif
		virtual inline void field(const double&      x, const double&      y, double      *Bx, double      *By) const override      { _field(x, y, Bx, By); }
	    // "Need to understand how to interpolate field with tps"
		virtual inline void field(const tps&         x, const tps&        y, tps         *Bx, tps         *By) const override       { _field(x, y, Bx, By); }
		virtual inline void field(const gtpsa::tpsa& x, const gtpsa::tpsa& y, gtpsa::tpsa *Bx, gtpsa::tpsa *By) const override      { _field(x, y, Bx, By); }

		virtual inline void gradient(const tps& x, const tps&    y, tps    *Gx, tps     *Gy) const override final{
			// "Need to understand how to interpolate gradient with tps"
			throw thor_scsi::NotImplemented("Multipoles: gradient in tps not implemented");
		}
		virtual inline void gradient(const gtpsa::tpsa& x, const gtpsa::tpsa&    y, gtpsa::tpsa    *Gx, gtpsa::tpsa     *Gy) const override final{
			// "Need to understand how to interpolate gradient with tps"
			throw thor_scsi::NotImplemented("Multipoles: gradient in tps, gtpsa::tpsa not implemented");
		}
		virtual inline void gradient(const tps& x, const tps&    y, double *Gx, double *Gy) const override final{
			// "Need to understand how to interpolate gradient with tps"
		        // throw thor_scsi::NotImplemented("Computing gradient in doubles for coordinates in tps x and tps y not implemented");
			const auto tmp = this->gradient({0e0, 0e0});
			*Gy = tmp.real();
			*Gx = tmp.imag();
		}
		virtual inline void gradient(const gtpsa::tpsa& x, const gtpsa::tpsa&    y, double *Gx, double *Gy) const override final{
			// "Need to understand how to interpolate gradient with tps"
		        // throw thor_scsi::NotImplemented("Computing gradient in doubles for coordinates in tps x and tps y not implemented");
			auto tmp = this->gradient({0e0, 0e0});
			*Gy = tmp.real();
			*Gx = tmp.imag();
		}
		/**
		 *
		 * Todo: Gradient here
		 */
		 /*
        inline complex_type gradient(const complex_type unused) const {
            return this->getMultipole(2);
        }
        */
        /*
         inline complex_type gradient(const std::complex<double>& unused) const {
             return this->getMultipole(2);
         }
         */
        inline std::complex<double> gradient(const std::complex<double>& unused) const {
            auto tmp = this->getMultipole(2);
            return gtpsa::cst(tmp);
        }
 		virtual inline void gradient(const double& x, const double& y, double *Gx, double *Gy) const override final {
			auto tmp = this->gradient(std::complex<double>(x,y));
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
		inline complex_type getMultipole(const unsigned int n) const {
			this->checkMultipoleIndex(n);
			unsigned  use_n = n - 1;
			// array indices matche the american notation
			// assert(use_n >= 0);
			assert(use_n <this->m_max_multipole);
			complex_type c = this->coeffs[use_n];
            return c;
			//complex_type res(c.real(), c.imag());
			//return res;
		}

		/**
		 * @brief set n'th multipole
		 */
		inline void setMultipole(const unsigned int n, const complex_type c){
			this->checkMultipoleIndex(n);
			unsigned use_n = n - 1;
			// assert(use_n >= 0);
			assert(use_n <this->m_max_multipole);
            this->coeffs[use_n] = c;
			//this->coeffs[use_n] = complex_type(c.real(), c.imag());
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
			complex_intern_type I (0e0, 1e0);
			for(size_t i = 0; i < this->coeffs.size(); ++i){
				const double phase = alpha * (i + 1);
				complex_intern_type scale = exp(I * phase);
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
        void applyTranslation(const complex_type dzs) {
            for (unsigned int i = 0; i < this->coeffs.size(); ++i) {
                complex_intern_type dzi(dzs.real(), dzs.imag());
                for (unsigned j = i + 1; j < this->coeffs.size(); ++j) {
                    this->coeffs[i] += double(binom(j, i)) * (this->coeffs[j] * dzi);
                }
            }
        }

        inline void applyTranslation(const double dx, const double dy) {
            const complex_type tmp(dx, dy);
            applyTranslation(tmp);
        }

        inline void applyTranslation(const double dx) {
            applyTranslation(dx, 0.0);
        }

	private:
	TwoDimensionalMultipolesKnobbed& right_multiply (const std::vector<complex_type>& scale, const bool begnin) {
		right_multiply_helper<complex_intern_type>(scale, begnin, &this->coeffs);
		return *this;
	}
	TwoDimensionalMultipolesKnobbed& right_add (const TwoDimensionalMultipolesKnobbed &other, const bool begnin){
		right_add_helper<complex_intern_type>(other.coeffs, begnin, &this->coeffs);
                return *this;
            }


        public:

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
	    TwoDimensionalMultipolesKnobbed& operator *= (const complex_type scale){
		    auto& c = this->coeffs;
		    auto f = [scale](const complex_intern_type c){
			    return c * scale;
		    };
		    std::transform(c.begin(), c.end(), c.begin(), f);
			return *this;
	    }

	    /*
	     * @brief: Scaling each individual by a vector
	     *
	     * @args: bengin: ignore if surplas number of scaling vector elements are provided
	     */
	    TwoDimensionalMultipolesKnobbed& operator *= (const std::vector<complex_type>& scale){
		    // should be false ... leave it that way to get further with prototyping
		    return this->right_multiply(scale, true);
	    }
	    TwoDimensionalMultipolesKnobbed &operator += (const TwoDimensionalMultipolesKnobbed &other){
		    // should be false ... leave it that way to get further with prototyping
		    return this->right_add(other, true);

            }
	    TwoDimensionalMultipolesKnobbed &operator += (const std::vector<complex_type> &other){
		    bool benign = true;
		    right_add_helper<complex_intern_type>(other, benign, &this->coeffs);
            return *this;
            }
	    TwoDimensionalMultipolesKnobbed& operator += (const double other) {
		    auto& c = this->getCoeffs();
		    for (auto& val : c){
			    val += other;
		    }
		    return *this;
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
	    TwoDimensionalMultipolesKnobbed operator+(const TwoDimensionalMultipolesKnobbed &other) const {
		    /*
		      std::cerr << "Adding multipoles. this size " << this->getCoeffs().size()
                      << " others size" << other.getCoeffs().size()
                      << std::endl;
                      */
		    TwoDimensionalMultipolesKnobbed nh = this->clone();
		    nh += other;
		    return nh;
	    }

	    TwoDimensionalMultipolesKnobbed  operator + (const std::vector<complex_type>& offset) const{
		    TwoDimensionalMultipolesKnobbed nh = this->clone();
		    nh += offset;
		    return nh;
	    }


	    TwoDimensionalMultipolesKnobbed operator + (const double other) const {
		    TwoDimensionalMultipolesKnobbed nh = this->clone();
		    /* std::cerr << "Adding const "<< other << "to these multipoles. this size " << this->getCoeffs().size()
		       << std::endl; */
		    nh += other;
		    return nh;
	    }

	    TwoDimensionalMultipolesKnobbed  operator * (const std::vector<complex_type>& scale) const{
		    TwoDimensionalMultipolesKnobbed nh = this->clone();
		    nh *= scale;
		    return nh;
	    }

	    TwoDimensionalMultipolesKnobbed  operator * (const complex_type& scale) const{
		    TwoDimensionalMultipolesKnobbed nh = this->clone();
		    nh *= scale;
		    return nh;
	    }

	    /**
	     * The maximum index represented.
	     */
	    inline unsigned int getMultipoleMaxIndex(void) const {
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
		inline std::vector<complex_intern_type>& getCoeffs(void) {
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
		inline const std::vector<complex_intern_type>& getCoeffs(void) const {
			return this->coeffs;
		}
		/**
		 * @brief: access to coeffcs for python
		 *
		 * pybind11 also seems to have limits for disambuigity
		 */
		inline const std::vector<complex_intern_type>& getCoeffsConst(void) const {
			return this->getCoeffs();
		}

	        virtual void show(std::ostream& strm, int level) const override final {
			multipoles_show(strm, this->coeffs, level);
		}

	protected:
		unsigned int m_max_multipole;
		std::vector<complex_intern_type> coeffs;


	};

  #if 0
        /**
	 * @brief: Representation of planar 2D harmonics / multipoles, begnin for excess elements
	 *
	 * Functionally the same implementation as the base class, but will ignore excess elements
	 *
	 */
	class BegninTwoDimensionalMultipolesDelegator {
	public:
		BegninTwoDimensionalMultipolesDelegator(std::shared_ptr<TwoDimensionalMultipolesKnobbed> obj);
		BegninTwoDimensionalMultipolesDelegator clone(void) const;

		// Shall one provide a short cut to the other operators too?
		BegninTwoDimensionalMultipolesDelegator& operator *= (const std::vector<double>& scale);
		BegninTwoDimensionalMultipolesDelegator  operator * (const std::vector<double>& scale) const;

		BegninTwoDimensionalMultipolesDelegator operator +(const TwoDimensionalMultipolesKnobbed &other) const;
		BegninTwoDimensionalMultipolesDelegator& operator +=(const TwoDimensionalMultipolesKnobbed &other);

	private:
		std::shared_ptr<TwoDimensionalMultipolesKnobbed> delegate;
	};
#endif

    class TwoDimensionalMultipoles : public TwoDimensionalMultipolesKnobbed<StandardDoubleType>
    {
        using base = TwoDimensionalMultipolesKnobbed<StandardDoubleType>;
    public:
        TwoDimensionalMultipoles(const base&& o)
                : base(std::move(o))
        {}

        TwoDimensionalMultipoles(const base& o)
                : base(o)
        {}

	TwoDimensionalMultipoles(const complex_type & default_value, const unsigned int h_max=max_multipole)
            : base(default_value, h_max)
            {}

        inline auto field(const std::complex<double> z){ return base::field(z); }

        virtual inline void field(const double& x, const double& y, double *Bx, double *By) const override final {
            base::field(x, y, Bx, By);
        }
        virtual inline void field(const gtpsa::tpsa& x, const gtpsa::tpsa& y, gtpsa::tpsa *Bx, gtpsa::tpsa *By) const override final {
            base::field(x, y, Bx, By);
        }
        virtual inline void field(const tps& x, const tps& y, tps *Bx, tps *By) const override final {
            honer_complex_as_real(this->coeffs, x, y, Bx, By);
        }

        inline TwoDimensionalMultipoles& operator += (const TwoDimensionalMultipoles           & o)   { base::operator+= (o); return *this;  }
        inline TwoDimensionalMultipoles& operator += (const std::complex<double>               & o)   { base::operator+= (o); return *this;  }
	inline TwoDimensionalMultipoles& operator += (const std::vector<std::complex<double>>  & o)   { base::operator+= (o); return *this;  }

        // TwoDimensionalMultipoles& operator *= (const TwoDimensionalMultipoles           & o)   { base::operator*= (o); return *this;  }
        inline TwoDimensionalMultipoles& operator *= (const std::complex<double>               & o)   { base::operator*= (o); return *this;  }
	inline TwoDimensionalMultipoles& operator *= (const std::vector<std::complex<double>>  & o)   { base::operator*= (o); return *this;  }

        inline TwoDimensionalMultipoles operator +   (const TwoDimensionalMultipoles           & o) const  { return TwoDimensionalMultipoles( base::operator+ (o) );  }
        inline TwoDimensionalMultipoles operator +   (const std::complex<double>               & o) const  { return TwoDimensionalMultipoles( base::operator+ (o) );  }
	inline TwoDimensionalMultipoles operator +   (const std::vector<std::complex<double>>  & o) const  { return TwoDimensionalMultipoles( base::operator+ (o) );  }

        // TwoDimensionalMultipoles operator *   (const TwoDimensionalMultipoles           & o) const  { return TwoDimensionalMultipoles( base::operator+ (o) );  }
        inline TwoDimensionalMultipoles operator *   (const std::complex<double>               & o) const  { return TwoDimensionalMultipoles( base::operator* (o) );  }
	inline TwoDimensionalMultipoles operator *   (const std::vector<std::complex<double>>  & o) const  { return TwoDimensionalMultipoles( base::operator* (o) );  }

    };

    inline TwoDimensionalMultipoles operator + (const std::complex<double>& a, const TwoDimensionalMultipoles& b) {return b + a; }
    inline TwoDimensionalMultipoles operator * (const std::complex<double>& a, const TwoDimensionalMultipoles& b) {return b * a; }
    /*
    inline TwoDimensionalMultipoles& operator+ (const TwoDimensionalMultipoles& a, const double& b) {
        TwoDimensionalMultipoles()
    }
*/
    /*
    template
    TwoDimensionalMultipoles<StandardDoubleType>::;
    }
     */

    class TwoDimensionalMultipolesTpsa : public TwoDimensionalMultipolesKnobbed<TpsaVariantType>
    {
        using base = TwoDimensionalMultipolesKnobbed<TpsaVariantType>;
    public:
        TwoDimensionalMultipolesTpsa(const base&& o)
                : base(std::move(o))
        {}

        TwoDimensionalMultipolesTpsa(const base& o)
                : base(o)
        {}

	TwoDimensionalMultipolesTpsa(const complex_type & default_value, const unsigned int h_max=max_multipole)
            : base(default_value, h_max)
            {}

	inline TwoDimensionalMultipolesTpsa& operator += (const TwoDimensionalMultipolesTpsa       & o)   { base::operator+= (o); return *this;  }
        inline TwoDimensionalMultipolesTpsa& operator += (const complex_type                       & o)   { base::operator+= (o); return *this;  }
	inline TwoDimensionalMultipolesTpsa& operator += (const std::vector<complex_type>          & o)   { base::operator+= (o); return *this;  }

        // TwoDimensionalMultipolesTpsa& operator *= (const TwoDimensionalMultipolesTpsa           & o)   { base::operator*= (o); return *this;  }
        inline TwoDimensionalMultipolesTpsa& operator *= (const complex_type                       & o)   { base::operator*= (o); return *this;  }
	inline TwoDimensionalMultipolesTpsa& operator *= (const std::vector<complex_type>          & o)   { base::operator*= (o); return *this;  }

        inline TwoDimensionalMultipolesTpsa operator +   (const TwoDimensionalMultipolesTpsa       & o) const  { return TwoDimensionalMultipolesTpsa( base::operator+ (o) );  }
        inline TwoDimensionalMultipolesTpsa operator +   (const complex_type                       & o) const  { return TwoDimensionalMultipolesTpsa( base::operator+ (o) );  }
	inline TwoDimensionalMultipolesTpsa operator +   (const std::vector<complex_type>          & o) const  { return TwoDimensionalMultipolesTpsa( base::operator+ (o) );  }

        // TwoDimensionalMultipolesTpsa operator *   (const TwoDimensionalMultipolesTpsa           & o) const  { return TwoDimensionalMultipolesTpsa( base::operator+ (o) );  }
        inline TwoDimensionalMultipolesTpsa operator *   (const complex_type                       & o) const  { return TwoDimensionalMultipolesTpsa( base::operator* (o) );  }
	inline TwoDimensionalMultipolesTpsa operator *   (const std::vector<complex_type>          & o) const  { return TwoDimensionalMultipolesTpsa( base::operator* (o) );  }

    };

    inline TwoDimensionalMultipolesTpsa operator + (const gtpsa::CTpsaOrComplex& a, const TwoDimensionalMultipolesTpsa& b) {return b + a; }
    inline TwoDimensionalMultipolesTpsa operator * (const gtpsa::CTpsaOrComplex& a, const TwoDimensionalMultipolesTpsa& b) {return b * a; }

}
#endif /* _THOR_SCSI_MULTIPOLES_H_ */
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
