#define BOOST_TEST_MODULE nonlinear_kicker
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <thor_scsi/custom/nonlinear_kicker.h>
#include <iostream>

namespace tsc = thor_scsi::core;
namespace tsu = thor_scsi::custom;

static void check_identity_mat(const arma::mat& mat){
    arma::mat Id(mat.n_rows, mat.n_cols, arma::fill::eye);
    arma::mat diff = mat - Id;
    arma::mat chk = arma::sum(arma::sum(arma::abs(diff), 1), 1);
    double val = chk(0, 0);
    BOOST_CHECK_SMALL(val, 1e-12);
}

static void check_vec_zero(const std::vector<double>& vec)
{
    for(const auto& e: vec){
        BOOST_CHECK_SMALL(e, 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(test10_nlk_zero_pos)
{
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("Method", 4.0);
	C.set<double>("N", 1);
	C.set<double>("L", 0e0);

	tsc::ConfigType  calc_config;
	tsu::NonLinearKickerType nlk(C);

	BOOST_CHECK_EQUAL(nlk.isThick(), false);
	BOOST_CHECK_SMALL(nlk.getLength(), 1e-12);
	BOOST_CHECK_SMALL(nlk.getFieldInterpolator()->getScale(), 1e-12);

	auto desc = std::make_shared<gtpsa::desc>(6, 2);
	gtpsa::ss_vect<gtpsa::tpsa> v1 (desc,1);
    v1.set_identity();
	nlk.propagate(calc_config, v1);

    check_identity_mat( v1.jacobian());
    check_vec_zero(v1.cst());

}


BOOST_AUTO_TEST_CASE(test20_nlk_kicker)
{
    Config C;
    C.set<std::string>("name", "test");
    C.set<double>("Method", 4.0);
    C.set<double>("N", 1);
    C.set<double>("L", 0e0);

    // Non linear kicker ... a bit easier to trace
    const double current = 2e2;
    const double wire_pos = 20e-3;
    tsu::aircoil_filament_t wire = {wire_pos, wire_pos,   current};
    std::vector<tsu::aircoil_filament_t> currents {wire};

    tsc::ConfigType  calc_config;
    tsu::NonLinearKickerType nlk(C);

    auto intp = std::make_shared<tsu::NonLinearKickerInterpolation>(currents);
    nlk.setFieldInterpolator(intp);

    BOOST_CHECK_EQUAL(nlk.isThick(), false);
    BOOST_CHECK_SMALL(nlk.getLength(), 1e-12);

    auto desc = std::make_shared<gtpsa::desc>(6, 2);
    gtpsa::ss_vect<gtpsa::tpsa> v1 (desc,1);
    v1.set_identity();
    intp->setScale(0e0);
    BOOST_CHECK_SMALL(intp->getScale(), 1e-12);

    // check if off at center
    nlk.propagate(calc_config, v1);
    std::cout << "nlk off at center " <<v1 << std::endl;
    check_vec_zero(v1.cst());
    check_identity_mat(v1.jacobian());

    // check if on at center
    v1.set_identity();
    nlk.propagate(calc_config, v1);
    intp->setScale(1e0);
    BOOST_CHECK_CLOSE(intp->getScale(), 1e0, 1e-12);
    std::cout << "nlk on: but center "<< v1 << std::endl;

    check_vec_zero(v1.cst());
    check_identity_mat(v1.jacobian());

    // check if on at center
    v1.set_identity();
    v1[0].set(0e0, wire_pos);
    BOOST_CHECK_CLOSE(v1.cst()[0], wire_pos, 1e-12);
    nlk.propagate(calc_config, v1);
    intp->setScale(1e0);
    BOOST_CHECK_CLOSE(intp->getScale(), 1e0, 1e-12);
    std::cout << "nlk on: but at wire pos "<< v1 << std::endl;

    std::vector<double> vec = v1.cst();
    BOOST_CHECK_CLOSE(vec[0], wire_pos, 1e-12);

    // 4 wires, wires in corner of square
    // as it is position in x, only the wires of quadrant 2 and 3 contribute
    // 4 and 1 cancel each other
    const double px = vec[1];
    const double wp2 = wire_pos/2, r2 = wp2*wp2 + 4 * wire_pos*wire_pos;
    const double root = (1/wp2 * 16e0/100e0);
    const double expected_px     = -1e-7 * current / 2e0 /  root;
    const double expected_dx_dy  = -1e-7 * current / 2e0 / (root*root)  * 2 * wire_pos;
    BOOST_CHECK_CLOSE(px, expected_px, 1e-12);

    vec[0] = 0e0;
    vec[1] = 0e0;
    check_vec_zero(vec);
    arma::mat jac = v1.jacobian();

#warning "Missing check for derivatives"
    //BOOST_CHECK_CLOSE(jac(1, 0), expected_dx_dy, 1e-12);
    //BOOST_CHECK_CLOSE(jac(3, 2), expected_dx_dy, 1e-12);

    jac(1, 0) = 0e0;
    jac(3, 2) = 0e0;

    check_identity_mat(jac);

}


BOOST_AUTO_TEST_CASE(test30_nlk_bessyii_kicker)
{
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("Method", 4.0);
	C.set<double>("N", 1);
	C.set<double>("L", 0e0);

	// BESSY II non linear kicker
	const double s = 14e-2, scale = 1 - s, current = 700, m_cur = current * s, cur = current * scale * scale;
    const double wire_pos = 8e-3;
	tsu::aircoil_filament_t outer  = {17e-3 * scale, 15e-3 * scale,   cur};
	tsu::aircoil_filament_t inner  = { wire_pos * scale,  7e-3 * scale,   cur};
	tsu::aircoil_filament_t mirror = { wire_pos * scale,  5e-3 * scale, m_cur};
	std::vector<tsu::aircoil_filament_t> currents {inner, outer, mirror};

	tsc::ConfigType  calc_config;
	tsu::NonLinearKickerType nlk(C);

	auto intp = std::make_shared<tsu::NonLinearKickerInterpolation>(currents);
	nlk.setFieldInterpolator(intp);

	BOOST_CHECK_EQUAL(nlk.isThick(), false);
	BOOST_CHECK_SMALL(nlk.getLength(), 1e-12);

	auto desc = std::make_shared<gtpsa::desc>(6, 2);
	gtpsa::ss_vect<gtpsa::tpsa> v1 (desc,1);
	v1.set_identity();
	intp->setScale(0e0);
	BOOST_CHECK_SMALL(intp->getScale(), 1e-12);
	nlk.propagate(calc_config, v1);
	std::cout << "nlk off at center" <<v1 << std::endl;

    v1.set_identity();
    nlk.propagate(calc_config, v1);
    intp->setScale(1e0);
    BOOST_CHECK_CLOSE(intp->getScale(), 1e0, 1e-12);
    std::cout << "nlk on: but center "<< v1 << std::endl;

    v1.set_identity();
    v1[0].set(0e0, wire_pos);
    nlk.propagate(calc_config, v1);
    intp->setScale(1e0);
    BOOST_CHECK_CLOSE(intp->getScale(), 1e0, 1e-12);
    std::cout << "nlk on: but at wire pos "<< v1 << std::endl;



}
