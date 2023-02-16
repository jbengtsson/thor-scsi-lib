#ifndef _THOR_SCSI_CORE_FIELD_INTERPOLATION_H_
#define _THOR_SCSI_CORE_FIELD_INTERPOLATION_H_ 1
#include <ostream>
#include <gtpsa/tpsa.hpp>
#include <tps/tps_type.h>
#include <thor_scsi/core/multipole_types.h>

namespace thor_scsi::core {
  	/**
	 * Pure virtual class of multipoles
	 *
	 * External code just requires to calculate multipoles
	 *
	 * Todo: move to API part?
	 */
    template<class C, typename = typename C::complex_type, typename = typename C::double_type>
    class Field2DInterpolationKnobbed{
    protected:
        using double_type = typename C::double_type;
        using complex_type = typename C::complex_type;
        using complex_intern_type = typename C::complex_intern_type;

    public:
		virtual ~Field2DInterpolationKnobbed(){};
		/**
		 * @brief interpolate field at position x and y
		 *
		 * \verbatim embed:rst:leading-asterisk
		 *
		 * Args:
		 *     x: first coordinate (typically horizontal position)
		 *     y: second coordinate  (typically vertical position)
		 *     Bx: pointer to variable where field component `x` will be stored
		 *     By: pointer to variable where field component `y` will be stored
		 *
		 * Warning:
		 *     Note: * thor_scsi* is using a left hand coordinate system.
		 *     This function here is expecting all coordinates in a left
		 *     handed coordiante system
		 *
		 * \endverbatim
		 */
		virtual inline void field(const double&      x, const double&      y, double      *Bx, double      *By) const = 0;
		virtual inline void field(const tps&         x, const tps&         y, tps         *Bx, tps         *By) const = 0;
		virtual inline void field(const gtpsa::tpsa& x, const gtpsa::tpsa& y, gtpsa::tpsa *Bx, gtpsa::tpsa *By) const = 0;

		/**
		 * @brief interpolate the gradient at the current position
		 *
		 * \verbatim embed:rst:leading-asterisk
		 * \endverbatim
		 */
		virtual void gradient(const double& x, const double& y, double *Gx, double *Gy) const = 0;
		/**
		 * @todo review interface: with tps it could be that Gx is computed but never used
		 * @todo should Gx and Gy be of type tps or type double
		 */
		virtual void gradient(const tps&         x, const tps&         y, tps         *Gx, tps         *Gy) const = 0;
		virtual void gradient(const gtpsa::tpsa& x, const gtpsa::tpsa& y, gtpsa::tpsa *Gx, gtpsa::tpsa *Gy) const = 0;
		virtual void gradient(const tps&         x, const tps&         y, double      *Gx, double      *Gy) const = 0;
		virtual void gradient(const gtpsa::tpsa& x, const gtpsa::tpsa& y, double      *Gx, double      *Gy) const = 0;
		virtual void show(std::ostream&, int level) const = 0;

		std::string prettyClassname(void) const;
		// facilitate python representation
		std::string repr(void) const;
		std::string pstr(void) const;



	};

    typedef Field2DInterpolationKnobbed<StandardDoubleType> Field2DInterpolation;
    typedef Field2DInterpolationKnobbed<TpsaVariantType>    Field2DInterpolationDependence;

    template<class C>
    inline
    std::ostream& operator<<(std::ostream& strm, const Field2DInterpolationKnobbed<C> & intp) {
        intp.show(strm, 0);
        return strm;
    }

}
#endif /* _THOR_SCSI_FIELD_INTERPOLATION_H_ */
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
