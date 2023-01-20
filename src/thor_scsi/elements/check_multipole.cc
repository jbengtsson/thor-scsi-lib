#include <boost/test/unit_test.hpp>
#include "check_multipole.h"

namespace tsc = thor_scsi::core;
#if 0
template<class C>
void CheckMultipoles<C>::only_major_multipole_set(std::shared_ptr<C> muls, const typename C::complex_type ref, const size_t n_major)
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
#endif

