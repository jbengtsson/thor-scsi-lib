#define BOOST_TEST_MODULE nonlinear_kicker
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <thor_scsi/custom/nonlinear_kicker_interpolation.h>
#include <iostream>


namespace tsc = thor_scsi::core;
namespace tsu = thor_scsi::custom;

BOOST_AUTO_TEST_CASE(test10_nlk_zero_pos)
{

    // BESSY II non linear kicker
    const double s = 14e-2, scale = 1 - s, current = 700, m_cur = current * s, cur = current * scale * scale;
    const double wire_pos = 8e-3;
    tsu::aircoil_filament_t outer  = {17e-3 * scale, 15e-3 * scale,   cur};
    tsu::aircoil_filament_t inner  = { wire_pos * scale,  7e-3 * scale,   cur};
    tsu::aircoil_filament_t mirror = { wire_pos * scale,  5e-3 * scale, m_cur};
    tsu:: NonLinearKickerInterpolation nlk({inner, outer, mirror});

    // check that no field returned at scale of 0
    BOOST_CHECK_SMALL(nlk.getScale(), 0e0);
    {
        const double x = 0e0, y = 0e0;
        double Bx=1e200, By=1e200;
        nlk.field(x, y, &Bx, &By);

        BOOST_CHECK_SMALL(Bx, 1e-12);
        BOOST_CHECK_SMALL(By, 1e-12);
    }

    // check that no field returned at scale of 0 if right below the "wire"
    {
        const double x = 8e-3, y = 0e0;
        double Bx=1e200, By=1e200;

        nlk.field(x, y, &Bx, &By);

        BOOST_CHECK_SMALL(Bx, 1e-12);
        BOOST_CHECK_SMALL(By, 1e-12);
    }

    // check that field returned at scale of 1 if right below the "wire"
    const double ref_field_y = 0.0232;
    nlk.setScale(1e0);
    {
        const double x = wire_pos, y = 0e0;
        double Bx=1e200, By=1e200;

        nlk.field(x, y, &Bx, &By);

        BOOST_CHECK_SMALL(Bx,              1e-12);
        BOOST_CHECK_CLOSE(By, ref_field_y, 0.1);
    }


    auto desc = std::make_shared<gtpsa::desc>(6,2);
    {
        gtpsa::tpsa x(desc, 1), y(desc, 1), Bx(desc, 1), By(desc, 1);
        x.setVariable(0e0, 0 + 1);
        y.setVariable(0e0, 2 + 1);

        // zero for zero field at zero position if scale at zero
        nlk.setScale(0e0);

        nlk.field(x, y, &Bx, &By);


        BOOST_CHECK_SMALL(Bx.cst(), 1e-12);
        BOOST_CHECK_SMALL(By.cst(), 1e-12);

        // zero for zero field at zero position if scale at one
        nlk.setScale(1e0);

        nlk.field(x, y, &Bx, &By);

        BOOST_CHECK_SMALL(Bx.cst(), 1e-12);
        BOOST_CHECK_SMALL(By.cst(), 1e-12);

        // check if scaling works for full field and for the tpsa series
        double scale = 1.0;
        nlk.setScale(scale);
        x.set(0e0, wire_pos);
        BOOST_CHECK_CLOSE(x.cst(), wire_pos, 1e-12);

        nlk.field(x, y, &Bx, &By);

        std::cout << "scale " << scale
                  << "\nx\n" << x
                  << "y\n" << y
                  << "Bx\n" << Bx
                  << "By\n" << By
                  << std::endl;

        BOOST_CHECK_SMALL(Bx.cst(), 1e-12);
        BOOST_CHECK_CLOSE(By.cst(), ref_field_y * scale, 0.1);

        // does it work at 10 % of  the field
        scale = 0.1;
        nlk.setScale(scale);

        nlk.field(x, y, &Bx, &By);
        std::cout << "scale " << scale
                  << "\nx\n" << x
                  << "y\n" << y
                  << "Bx\n" << Bx
                  << "By\n" << By
                  << std::endl;

        BOOST_CHECK_SMALL(Bx.cst(), 1e-12);
        BOOST_CHECK_CLOSE(By.cst(), ref_field_y * scale, 0.1);

        // does it work at 0 the field
        scale = 0e0;
        nlk.setScale(scale);

        nlk.field(x, y, &Bx, &By);
        std::cout << "scale " << scale
                  << "\nx\n" << x
                  << "y\n" << y
                  << "Bx\n" << Bx
                  << "By\n" << By
                  << std::endl;

        BOOST_CHECK_SMALL(Bx.cst(), 1e-12);
        BOOST_CHECK_CLOSE(By.cst(), ref_field_y * scale, 0.1);

    }


}
