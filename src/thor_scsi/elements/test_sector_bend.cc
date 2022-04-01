#define BOOST_TEST_MODULE sector_bend
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <cmath>
#include <thor_scsi/core/field_interpolation.h>
#include <thor_scsi/core/multipoles.h>
#include <thor_scsi/elements/bending.h>
#include "check_multipole.h"
#include <tps/ss_vect.h>
#include <tps/enums.h>
#include <ostream>
#include <armadillo>
#include <cmath>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

/**
 * code originally from tracy / thor_scsi on hold
 */
arma::mat get_sbend_mat(const double length, const double b2,
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
	auto mat = get_sbend_mat(length, b2, phi, delta);
	std::cout << "Analytic result:\n" << mat << std::endl;
}

BOOST_AUTO_TEST_CASE(test01_sector_bend_ostream)
{
	const double length = 1e0, b2 = 1e0, phi = 5e0, delta = 0e0;
	auto mat = get_sbend_mat(length, b2, phi, delta);
	std::cout << "Analytic result:\n" << mat << std::endl;
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
symplectic_result_check_non_zeros_elements(ss_vect<tps> ps, arma::mat mat)
{
	//arma::inplace_trans(mat);
	const double eps = 1e-6;
	BOOST_CHECK_CLOSE((ps[x_][x_]),         (mat.at(x_, x_)),         eps);

	BOOST_CHECK_CLOSE((ps[x_][px_]),        (mat.at(x_, px_)),        eps);
	BOOST_CHECK_CLOSE((ps[x_][delta_]),     (mat.at(x_, delta_)),     eps);

	BOOST_CHECK_CLOSE((ps[px_][x_]),        (mat.at(px_, x_)),        eps);
	BOOST_CHECK_CLOSE((ps[px_][px_]),       (mat.at(px_, px_)),       eps);
	BOOST_CHECK_CLOSE((ps[px_][delta_]),    (mat.at(px_, delta_)),    eps);

	BOOST_CHECK_CLOSE((ps[y_][y_]),         (mat.at(y_, y_)),         eps);
	BOOST_CHECK_CLOSE((ps[y_][py_]),        (mat.at(y_, py_)),        eps);

	BOOST_CHECK_CLOSE((ps[py_][y_]),        (mat.at(py_, y_)),        eps);
	BOOST_CHECK_CLOSE((ps[py_][py_]),       (mat.at(py_, py_)),       eps);

	BOOST_CHECK_CLOSE((ps[ct_][x_]),        (mat.at(ct_, x_)),        eps);
	BOOST_CHECK_CLOSE((ps[ct_][px_]),       (mat.at(ct_, px_)),       eps);
	BOOST_CHECK_CLOSE((ps[ct_][delta_]),    (mat.at(ct_, delta_)),    eps);
	BOOST_CHECK_CLOSE((ps[ct_][ct_]),       (mat.at(ct_, ct_)),       eps);

	BOOST_CHECK_CLOSE((ps[delta_][delta_]), (mat.at(delta_, delta_)), eps);


}


static arma::mat
compute_omega_matrix(void)
{
	const int n = 3;
	int k;
	arma::mat omega = arma::mat(2*n, 2*n);

	omega.zeros();
	for (k = 0; k < 2; k++) {
	  omega(2*k, 2*k+1) = 1e0;
	  omega(2*k+1, 2*k) = -1e0;
	}
	omega(ct_, delta_) = 1e0;
	omega(delta_, ct_) = -1e0;

	omega.print(std::cout, "omega");

	BOOST_CHECK_CLOSE(arma::det(omega), 1, 1e-12);

	auto tmp = arma::abs(arma::trans(omega) + omega);
	double test = arma::accu(tmp);
	BOOST_CHECK_SMALL(test, 1e-12);

	arma::mat tmp2 = arma::abs(arma::inv(omega) + omega);
	double test2 = arma::sum(arma::sum(tmp2));
	BOOST_CHECK_SMALL(test2, 1e-12);

	return omega;
}

static const auto omega = compute_omega_matrix();

static void
compute_symplectic_test(arma::mat mat)
{
	arma::mat mat_t = arma::trans(mat);
	mat.print(std::cout, "compute_symplectic_test:\nmat");
	arma::mat test = mat_t * (omega * mat);
	test.print(std::cout, "test");
}
static void
check_symplectisism(arma::mat mat)
{

}


BOOST_AUTO_TEST_CASE(test20_sector_bend_compare)
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

	std::cout << C << std::endl;
	bend.show(std::cout, 4); std::cout << std::endl;

	auto mat = get_sbend_mat(length, b2, phi, delta);
	std::cout.setf(std::ios::scientific);
	std::cout.precision(6);
	std::cout.width(14);
	mat.raw_print(std::cout, "Analytic result:");

	ss_vect<tps> ps_orig;
	ps_orig.identity();

	ss_vect<tps> ps = ps_orig;

	bend.pass(calc_config, ps);
	std::cout << "Symplectic integrator:\n"
		  << std::setw(7) << ps << std::endl;

	// std::cout << ps[x_][px_] << " " <<  mat.at(x_, px_) << std::endl;
	// std::cout << ps[px_][x_] << " " << mat.at(px_, x_) << std::endl;
	//symplectic_result_check_zeros_elements(ps);
	// symplectic_result_check_non_zeros_elements(ps, mat);

	compute_symplectic_test(mat(arma::span(0, 5), arma::span(0, 5)));
}
