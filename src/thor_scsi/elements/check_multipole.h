#ifndef _THOR_SCSI_ELEMENTS_CHECK_MULTIPOLE_H_
#define _THOR_SCSI_ELEMENTS_CHECK_MULTIPOLE_H_
#include <memory>
#include <thor_scsi/core/multipoles.h>

#include <thor_scsi/core/multipole_types.h>

typedef typename thor_scsi::core::StandardDoubleType::complex_type cdbl;
namespace tsc = thor_scsi::core;

template<class C>
struct CheckMultipoles{
    using complex_type = typename C::complex_type;
    //void only_major_multipole_set(std::shared_ptr<C> muls, const complex_type ref, const size_t n_major);
    void only_major_multipole_set(std::shared_ptr<tsc::TwoDimensionalMultipolesKnobbed<C>> muls, const typename C::complex_type ref, const size_t n_major)
    {
        const auto& coeffs =  muls->getCoeffs();
        for(size_t i = 0; i<coeffs.size(); ++i){
            const auto& c = coeffs[i];
            /*
              thor_scsi counts multipoles following European Convention
              while C/C++ indexing follows American convention
            */
            if(i == (n_major -1)){
                if(std::abs(c.real()) > 1e-3){
                    BOOST_CHECK_CLOSE(c.real(), ref.real(), 1e-12);
                } else {
                    BOOST_CHECK_SMALL(c.real(), 1e-12);
                }

                if(std::abs(c.imag()) > 1e-3){
                    BOOST_CHECK_CLOSE(c.imag(), ref.imag(), 1e-12);
                } else {
                    BOOST_CHECK_SMALL(c.imag(), 1e-12);
                }
            } else {
                BOOST_CHECK_SMALL(c.real(), 1e-12);
                BOOST_CHECK_SMALL(c.imag(), 1e-12);
            }
        }
    }

};

typedef CheckMultipoles<tsc::StandardDoubleType> CheckMultipolesStd;

#endif
