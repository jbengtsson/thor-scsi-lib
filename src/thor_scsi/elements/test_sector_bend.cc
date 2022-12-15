#define BOOST_TEST_MODULE sector_bend
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/core/field_interpolation.h>
#include <thor_scsi/core/multipoles.h>
#include <thor_scsi/elements/bending.h>
#include "check_multipole.h"
#include <tps/enums.h>
#include <ostream>
#include <armadillo>
#include <cmath>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;


auto a_desc = std::make_shared<gtpsa::desc>(6, 1);
auto tpsa_ref = gtpsa::tpsa(a_desc, gtpsa::init::default_);

/**
 * code originally from tracy / thor_scsi on hold
 */
static arma::mat get_sbend_mat(const double length, const double b2,
			       const double phi,
			       const double delta)
{

	arma::mat mat = arma::mat(6, 6);

	const double
		phi_rad = phi / 180 * M_PI,
		L       = length,
		h       = phi_rad / L,
		K_x     = b2 + sqr(h),
		K_y     = fabs(b2),
		psi_x   = sqrt(fabs(K_x)/(1e0+delta))*L,
		psi_y   = sqrt(K_y/(1e0+delta))*L;

	mat.eye(6, 6);

	if (K_x > 0e0) {
		mat(x_, x_)      = cos(psi_x);
		mat(x_, px_)     = sin(psi_x)/sqrt(K_x*(1e0+delta));
		mat(x_, delta_)  = (1e0-cos(psi_x))*h/K_x;
		mat(px_, x_)     = -sqrt(K_x*(1e0+delta))*sin(psi_x);
		mat(px_, px_)    = cos(psi_x);
		mat(px_, delta_) = sin(psi_x)*sqrt(1e0+delta)*h/sqrt(K_x);

		if (psi_y != 0e0) {
			mat(y_, y_)    = cosh(psi_y);
			mat(y_, py_)   = sinh(psi_y)/sqrt(K_y*(1e0+delta));
			mat(py_, y_)   = sqrt(K_y*(1e0+delta))*sinh(psi_y);
			mat(py_, py_)  = cosh(psi_y);
		} else {
			mat(y_, py_)   = L/(1e0+delta);
		}
		mat(ct_, x_)     = sin(psi_x)*sqrt(1e0+delta)*h/sqrt(K_x);
		mat(ct_, px_)    = (1e0-cos(psi_x))*h/K_x;
		mat(ct_, delta_) =
			(psi_x-sin(psi_x))*sqrt(1e0+delta)*sqr(h)/pow(K_x, 3e0/2e0);
	} else if (K_x < 0e0) {
		mat(x_, x_)      = cosh(psi_x);
		mat(x_, px_)     = sinh(psi_x)/sqrt(fabs(K_x)*(1e0+delta));
		mat(x_, delta_)  = -(1e0-cosh(psi_x))*h/fabs(K_x);
		mat(px_, x_)     = sqrt(fabs(K_x)*(1e0+delta))*sinh(psi_x);
		mat(px_, px_)    = cosh(psi_x);
		mat(px_, delta_) = sinh(psi_x)*sqrt(1e0+delta)*h/sqrt(fabs(K_x));

		if (psi_y != 0e0) {
			mat(y_, y_)    = cos(psi_y);
			mat(y_, py_)   = sin(psi_y)/sqrt(K_y*(1e0+delta));
			mat(py_, y_)   = -sqrt(K_y*(1e0+delta))*sin(psi_y);
			mat(py_, py_)  = cos(psi_y);
		} else {
			mat(y_, py_)   = L/(1e0+delta);
		}
		mat(ct_, x_)     = sinh(psi_x)*sqrt(1e0+delta)*h/sqrt(fabs(K_x));
		mat(ct_, px_)    = -(1e0-cosh(psi_x))*h/fabs(K_x);
		mat(ct_, delta_) =
			- (psi_x-sinh(psi_x))*sqrt(1e0+delta)*sqr(h)/pow(fabs(K_x), 3e0/2e0);
	} else {
		// K_x = 0.
		mat(x_, px_) = L/(1e0+delta);
		mat(y_, py_) = L/(1e0+delta);
	}
	return mat;
}

BOOST_AUTO_TEST_CASE(test00_consistency)
{
	BOOST_CHECK_CLOSE(M_PI, 22e0/7e0, 1);
	BOOST_CHECK_CLOSE(M_PI, 355e0/113e0, 1e-3);
}

BOOST_AUTO_TEST_CASE(test00_sector_bend_ostream)
{
	const double length = 1e-3, b2 = 0e0, phi = 0e0, delta = 0e0;
	arma::mat mat = get_sbend_mat(length, b2, phi, delta);

	boost::test_tools::output_test_stream output;
	mat.print(output, "Analytic result");
	BOOST_CHECK( !output.is_empty( false ) );

}

BOOST_AUTO_TEST_CASE(test01_sector_bend_ostream)
{
	const double length = 1e0, b2 = 1e0, phi = 5e0, delta = 0e0;
	arma::mat mat = get_sbend_mat(length, b2, phi, delta);

	boost::test_tools::output_test_stream output;
	mat.print(output, "Analytic result");
	BOOST_CHECK( !output.is_empty( false ) );

}

BOOST_AUTO_TEST_CASE(test10_sector_bend_ostream)
{

	tsc::ConfigType calc_config;
	Config C;
	C.set<std::string>("name", "test");
	C.set<double>("K", 0.0);
	C.set<double>("L", 0e0);
	C.set<double>("N", 1);

	auto bend = tse::BendingType(C);

	boost::test_tools::output_test_stream output;
	output << bend;
	BOOST_CHECK( !output.is_empty( false ) );
}




template<typename T>
static void
symplectic_result_check_zeros_elements(T mat)
{

	BOOST_CHECK_SMALL(mat[x_][y_],         1e-12);
	BOOST_CHECK_SMALL(mat[x_][py_],        1e-12);
	BOOST_CHECK_SMALL(mat[x_][ct_],        1e-12);

	BOOST_CHECK_SMALL(mat[px_][y_],        1e-12);
	BOOST_CHECK_SMALL(mat[px_][py_],       1e-12);
	BOOST_CHECK_SMALL(mat[px_][ct_],       1e-12);

	BOOST_CHECK_SMALL(mat[y_][x_],         1e-12);
	BOOST_CHECK_SMALL(mat[y_][px_],        1e-12);
	BOOST_CHECK_SMALL(mat[y_][delta_],     1e-12);
	BOOST_CHECK_SMALL(mat[y_][ct_],        1e-12);

	BOOST_CHECK_SMALL(mat[py_][x_],        1e-12);
	BOOST_CHECK_SMALL(mat[py_][px_],       1e-12);
	BOOST_CHECK_SMALL(mat[py_][ct_],       1e-12);
	BOOST_CHECK_SMALL(mat[py_][delta_],    1e-12);

	BOOST_CHECK_SMALL(mat[ct_][y_],        1e-12);
	BOOST_CHECK_SMALL(mat[ct_][py_],       1e-12);

 	BOOST_CHECK_SMALL(mat[delta_][x_],     1e-12);
	BOOST_CHECK_SMALL(mat[delta_][px_],    1e-12);
 	BOOST_CHECK_SMALL(mat[delta_][y_],     1e-12);
	BOOST_CHECK_SMALL(mat[delta_][py_],    1e-12);
	BOOST_CHECK_SMALL(mat[delta_][ct_],    1e-12);


}


static void
symplectic_result_check_non_zeros_elements(gtpsa::ss_vect<gtpsa::tpsa> ps_arg, arma::mat mat)
{
	//arma::inplace_trans(mat);
	const double eps = 1e-6;

	arma::mat ps = ps_arg.jacobian();
	BOOST_CHECK_CLOSE( ps.at(  x_    , x_     ), mat.at(  x_    , x_     ), eps);

	BOOST_CHECK_CLOSE( ps.at(  x_    , px_    ), mat.at(  x_    , px_    ), eps);
	BOOST_CHECK_CLOSE( ps.at(  x_    , delta_ ), mat.at(  x_    , delta_ ), eps);

	BOOST_CHECK_CLOSE( ps.at( px_    , x_     ), mat.at( px_    , x_     ), eps);
	BOOST_CHECK_CLOSE( ps.at( px_    , px_    ), mat.at( px_    , px_    ), eps);
	BOOST_CHECK_CLOSE( ps.at( px_    , delta_ ), mat.at( px_    , delta_ ), eps);

	BOOST_CHECK_CLOSE( ps.at(  y_    , y_     ), mat.at(  y_    , y_     ), eps);
	BOOST_CHECK_CLOSE( ps.at(  y_    , py_    ), mat.at(  y_    , py_    ), eps);

	BOOST_CHECK_CLOSE( ps.at( py_    , y_     ), mat.at( py_    , y_     ), eps);
	BOOST_CHECK_CLOSE( ps.at( py_    , py_    ), mat.at( py_    , py_    ), eps);

	BOOST_CHECK_CLOSE( ps.at( ct_    , x_     ), mat.at( ct_    , x_     ), eps);
	BOOST_CHECK_CLOSE( ps.at( ct_    , px_    ), mat.at( ct_    , px_    ), eps);
	BOOST_CHECK_CLOSE( ps.at( ct_    , delta_ ), mat.at( ct_    , delta_ ), eps);
	BOOST_CHECK_CLOSE( ps.at( ct_    , ct_    ), mat.at( ct_    , ct_    ), eps);

	BOOST_CHECK_CLOSE( ps.at( delta_ , delta_ ), mat.at( delta_ , delta_ ), eps);

}

static arma::mat
compute_omega_matrix(void)
{
	const int n = 3;
	int k;
	arma::mat omega = arma::mat(2*n, 2*n);

	omega.zeros();
#if 1
	omega(x_,     px_   ) =  1e0;
	omega(px_,    x_    ) = -1e0;

	omega(y_,     py_   ) =  1e0;
	omega(py_,    y_    ) = -1e0;

	// columns and rows are swapped
	omega(ct_,    delta_) = -1e0;
	omega(delta_, ct_   ) =  1e0;

#else
	for (k = 0; k < 2; k++) {
	  omega(2*k, 2*k+1) = 1e0;
	  omega(2*k+1, 2*k) = -1e0;
	}
	omega(ct_, delta_) = -1e0;
	omega(delta_, ct_) = 1e0;

#endif


	// omega.print(std::cout, "Omega");
	return omega;
}

BOOST_AUTO_TEST_CASE(test03_omega_properties)
{
	arma::mat omega = compute_omega_matrix();

	BOOST_CHECK_CLOSE(arma::det(omega), 1, 1e-12);

	arma::mat tmp = arma::abs(arma::trans(omega) + omega);
	// tmp.print(std::cout, "omega test transpose");
	double omega_trans_omega = arma::accu(tmp);
	BOOST_CHECK_SMALL(omega_trans_omega, 1e-12);

	arma::mat tmp2 = arma::abs(arma::inv(omega) + omega);
	// tmp2.print(std::cout, "omega test inverse");
	double omega_inv_omega = arma::accu(tmp2);
	BOOST_CHECK_SMALL(omega_inv_omega, 1e-12);

}

static arma::mat
compute_symplectic_test(arma::mat mat, arma::mat om)
{
	arma::mat mat_t = arma::trans(mat);
	// std::cout << "Compute symplectic test" << std::endl;
	// mat.print(std::cout, "mat");
	arma::mat test = mat_t * (om * mat);
	// test.print(std::cout, "test");
	arma::mat expected_zero = test -om;
	// expected_zero.print(std::cout, "all should be zero ");
	return test;
}

static void
check_symplectisism_submats(const arma::mat mat)
{

	const double eps = 1e-12;
	arma::mat om(2,2, arma::fill::zeros);
	om(0, 1) =  1;
	om(1, 0) = -1;

	arma::mat x_px = compute_symplectic_test(mat(arma::span(0, 1), arma::span(0, 1)), om);
	double x_px_check = arma::accu(arma::abs(x_px - om));
	BOOST_CHECK_SMALL(x_px_check, eps);

	arma::mat y_py = compute_symplectic_test(mat(arma::span(2, 3), arma::span(2, 3)), om);
	double y_py_check = arma::accu(arma::abs(y_py - om));
	BOOST_CHECK_SMALL(y_py_check, eps);

	arma::mat om_d(2,2, arma::fill::zeros);
	om_d(0, 1) =  -1;
	om_d(1, 0) =   1;

	arma::mat delta_ct = compute_symplectic_test(mat(arma::span(4, 5), arma::span(4, 5)), om_d);
	double delta_ct_check = arma::accu(arma::abs(delta_ct - om_d));
	BOOST_CHECK_SMALL(delta_ct_check, eps);
}

static void
check_symplectisism(const arma::mat mat)
{
	const double eps = 1e-12;

	check_symplectisism_submats(mat);

	arma::mat omega = compute_omega_matrix();
	arma::mat a_test = compute_symplectic_test(mat, omega);

	BOOST_CHECK_CLOSE(a_test(x_,     px_    ),  1e0, eps);
	BOOST_CHECK_CLOSE(a_test(px_,    x_     ), -1e0, eps);

	BOOST_CHECK_CLOSE(a_test(y_,     py_    ),  1e0, eps);
	BOOST_CHECK_CLOSE(a_test(py_,    y_     ), -1e0, eps);

	BOOST_CHECK_CLOSE(a_test(ct_,    delta_ ), -1e0, eps);
	BOOST_CHECK_CLOSE(a_test(delta_, ct_    ),  1e0, eps);

	arma::mat expected_zero = a_test - omega;
	double a_sum = arma::accu(arma::abs(expected_zero));
	BOOST_CHECK_SMALL(a_sum, eps);

	if (a_sum > eps) {
		omega.print(std::cout, "omega");
		a_test.print(std::cout, "expected to match omega");
		expected_zero.print(std::cout, "test mat: all elements should be zero");
	}
}

BOOST_AUTO_TEST_CASE(test20_sector_analytic_result_symplectic)
{

	tsc::ConfigType calc_config;
	Config C;
	const double length = 1e0, b2 = 1e0, phi = 5e0, delta = 0e0;

	arma::mat mat = get_sbend_mat(length, b2, phi, delta);
	std::cout.setf(std::ios::scientific);
	std::cout.precision(6);
	std::cout.width(14);
	//mat.raw_print(std::cout, "Symplectic analytic result:");

	//std::cout << "Symplectic integrator:\n"
	//	  << std::setw(7) << ps << std::endl;

	// std::cout << ps[x_][px_] << " " <<  mat.at(x_, px_) << std::endl;
	// std::cout << ps[px_][x_] << " " << mat.at(px_, x_) << std::endl;
	//symplectic_result_check_zeros_elements(ps);
	// symplectic_result_check_non_zeros_elements(ps, mat);

	arma::mat  mat_rev(mat);
	// swap delta and ct
	mat_rev.col(4) = mat.col(5);
	mat_rev.col(5) = mat.col(4);
	mat_rev.row(4) = mat.row(5);
	mat_rev.row(5) = mat.row(4);


	// arma::mat omega = compute_omega_matrix();
	// arma::mat test = compute_symplectic_test(mat, omega);
	// test.print("Symplectic: test matrix");
	// (test - omega).print("Symplectic: comparsion all zeros?");
	check_symplectisism(mat);
	 //compute_symplectic_test(mat(arma::span(0, 5), arma::span(0, 5)));
}

static void extract_ps_jac(arma::mat mat, arma::mat *ps, arma::mat *jac)
{
	*jac = mat(arma::span(0, 5), arma::span(0, 5));
	*ps = mat(arma::span(6, 6), arma::span(0, 6));


}

template<typename T>
static void to_ps_jac(gtpsa::ss_vect<T>& ssv, arma::mat *ps, arma::mat *jac)
{
    *ps = ssv.cst();
    *jac = ssv.jacobian();
}

static void to_ps_jac(gtpsa::ss_vect<tps>& ssv, arma::mat *ps, arma::mat *jac)
{
    throw std::runtime_error("Not implemented");
    // extract_ps_jac(maptomat(ssv), ps, jac);
}

BOOST_AUTO_TEST_CASE(test21_sector_tps_symplectic)
{
	tsc::ConfigType calc_config;
	Config C;
	const double length = 1e0, b2 = 1e0, phi = 5e0, delta = 0e0;

	C.set<std::string>("name", "test");
	C.set<double>("K", b2);
	C.set<double>("L", length);
	C.set<double>("T", phi);
	C.set<double>("N", 100);

	auto bend = tse::BendingType(C);

	// std::cout << C << std::endl;
	// bend.show(std::cout, 4); std::cout << std::endl;

	gtpsa::ss_vect<gtpsa::tpsa> ps(tpsa_ref);
	ps.set_identity();

	const gtpsa::ss_vect<gtpsa::tpsa> ps_ref = ps.clone();
	bend.propagate(calc_config, ps);

	arma::mat ps2, jac;

	to_ps_jac(ps, &ps2, &jac);
	check_symplectisism(jac);
#if 0
	extract_ps_jac(maptomat(ps), &ps2, &jac);

	// jac.print("checking jacobian ");
#endif

}

BOOST_AUTO_TEST_CASE(test30_bend_config)
{
	tsc::ConfigType calc_config;
	Config C;
	const double length = 1e0, b2 = 1e0, phi = 5e0, delta = 0e0, T1=1.2e0, T2=0.75e0;
	const double phi_rad = phi * M_PI / 180e0;
	const double curvature = phi_rad / length;
	C.set<std::string>("name", "test");
	C.set<double>("K",  b2);
	C.set<double>("L",  length);
	C.set<double>("T",  phi);
	C.set<double>("T1", T1);
	C.set<double>("T2", T2);
	C.set<double>("N", 100);

	auto bend_ref = tse::BendingType(C);
	// std::cout << "bend " << bend_ref << std::endl;
	// std::cout << "bend "; bend_ref.show(std::cout, 4); std::cout << std::endl;

	const double eps = 1e-12;
	{
		tse::BendingType& bend = bend_ref;

		BOOST_CHECK_CLOSE(bend.getLength(),        length,    eps);
		BOOST_CHECK_CLOSE(bend.getBendingAngle(),  phi,       eps);
		BOOST_CHECK_CLOSE(bend.getEntranceAngle(), T1,        eps);
		BOOST_CHECK_CLOSE(bend.getExitAngle(),     T2,        eps);
		BOOST_CHECK_CLOSE(bend.getCurvature(),     curvature, eps);

		BOOST_CHECK(bend.assumingCurvedTrajectory() == true);


		auto mul = bend.getMultipoles()->getMultipole(2);
		BOOST_CHECK_CLOSE(mul.real(), b2, eps);
		BOOST_CHECK_SMALL(mul.imag(),     eps);
	}

	{
		tse::BendingType  bend = std::move(bend_ref);

		BOOST_CHECK_CLOSE(bend.getLength(),        length,    eps);
		BOOST_CHECK_CLOSE(bend.getBendingAngle(),  phi,       eps);
		BOOST_CHECK_CLOSE(bend.getEntranceAngle(), T1,        eps);
		BOOST_CHECK_CLOSE(bend.getExitAngle(),     T2,        eps);
		BOOST_CHECK_CLOSE(bend.getCurvature(),     curvature, eps);

		BOOST_CHECK(bend.assumingCurvedTrajectory() == true);


		auto mul = bend.getMultipoles()->getMultipole(2);
		BOOST_CHECK_CLOSE(mul.real(), b2, eps);
		BOOST_CHECK_SMALL(mul.imag(),     eps);
	}

}
