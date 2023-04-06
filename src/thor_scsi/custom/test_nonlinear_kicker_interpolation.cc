#define BOOST_TEST_MODULE nonlinear_kicker
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <thor_scsi/custom/nonlinear_kicker_interpolation.h>

namespace tsc = thor_scsi::core;
namespace tsu = thor_scsi::custom;

BOOST_AUTO_TEST_CASE(test10_nlk_zero_pos)
{

    // BESSY II non linear kicker
    const double s = 14e-2, scale = 1 - s, current = 700, m_cur = current * s, cur = current * scale * scale;

    tsu::aircoil_filament_t outer  = {17e-3 * scale, 15e-3 * scale,   cur};
    tsu::aircoil_filament_t inner  = { 8e-3 * scale,  7e-3 * scale,   cur};
    tsu::aircoil_filament_t mirror = { 8e-3 * scale,  5e-3 * scale, m_cur};
    tsu:: NonlinearKicker nlk({inner, outer, mirror});

    const double x = 0e0, y = 0e0;
    double Bx=1e200, By=1e200;
    nlk.field(x, y, &Bx, &By);

    BOOST_CHECK_SMALL(Bx, 1e-12);
    BOOST_CHECK_SMALL(By, 1e-12);

}
