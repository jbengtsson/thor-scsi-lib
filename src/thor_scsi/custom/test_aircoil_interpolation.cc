#define BOOST_TEST_MODULE drift
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <thor_scsi/custom/aircoil_interpolation.h>

namespace tsc = thor_scsi::core;
namespace tsu = thor_scsi::custom;

BOOST_AUTO_TEST_CASE(test02_drift_zero_length)
{
    tsu::aircoil_filament_t f1 = {20e-3, 10e-3, 700e0};
    tsu::AirCoilMagneticField am({f1});

}
