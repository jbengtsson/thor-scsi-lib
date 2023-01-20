//
// Created by mfp on 08.01.23.
//

#ifndef THOR_SCSI_CORE_MULTIPOLE_TYPES_H
#define THOR_SCSI_CORE_MULTIPOLE_TYPES_H
#include <complex>
#include <gtpsa/tpsa.hpp>
#include <gtpsa/ctpsa.hpp>
#include <gtpsa/tpsa_double_variant.hpp>
#include <gtpsa/ctpsa_complex_variant.hpp>

namespace thor_scsi::core {
    /**
     * @brief define which type of multipoles to use
     *
     * @tparam C type matching complex type
     * @tparam D type matching standard double type
     *
     * @todo rename double to float or real?
     */
    template<typename C, typename D>
    struct MultipoleType
    {
        using complex_type =  C;
        using double_type = D;
        // currently not supporting a sepearte intern complex type
        using complex_intern_type = C;
    };

#ifdef THOR_SCSI_USE_F128
    #error "currently not supporting f128 type"
	typedef boost::multiprecision::complex128 cdbl_intern;
#else // THOR_SCSI_USE_F128
    /* standard multipoles ... complex based on doubles */
    typedef MultipoleType<std::complex<double>, double> StandardDoubleType;
    //typedef std::complex<double> cdbl_intern;
#endif // THOR_SCSI_USE_F128

#if 0
    /*
     * multipoles: using ctpsa or tpsa objects
     * currently unused
     */
    typedef MultipoleType<gtpsa::ctpsa, gtpsa::tpsa> TpsaType;
#endif
    /* multipoles: using variants ... */
    typedef MultipoleType<gtpsa::CTpsaOrComplex, gtpsa::TpsaOrDouble> TpsaVariantType;
}
#endif //THOR_SCSI_CORE_MULTIPOLE_TYPES_H
