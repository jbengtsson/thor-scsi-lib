#define BOOST_TEST_MODULE quadrupole
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <tps/utils.h>

BOOST_AUTO_TEST_CASE(test10_factorial)
{
	BOOST_TEST(fact(0L) ==   1);
	BOOST_TEST(fact(1L) ==   1);
	BOOST_TEST(fact(2L) ==   2);
	BOOST_TEST(fact(3L) ==   6);
	BOOST_TEST(fact(4L) ==  24);

        // reference number computed by scipy.special.factorial(15, exact=True)
	BOOST_TEST(fact(15L) == 1307674368000L);

        // reference number computed by scipy.special.factorial(15, exact=True)
	BOOST_TEST(fact(20L) == 2432902008176640000L);
}


BOOST_AUTO_TEST_CASE(test10_binom_pascal_triangle)
{
	BOOST_TEST(binom(0, 1) ==   0);
	BOOST_TEST(binom(0, 2) ==   0);
	BOOST_TEST(binom(0, 3) ==   0);

	// the top
	BOOST_TEST(binom(1, 1) ==   1);

	// quadratic expansion
	BOOST_TEST(binom(2, 0) ==   1);
	BOOST_TEST(binom(2, 1) ==   2);
	BOOST_TEST(binom(2, 2) ==   1);

	// cubic expansion
	BOOST_TEST(binom(3, 0) ==   1);
	BOOST_TEST(binom(3, 1) ==   3);
	BOOST_TEST(binom(3, 2) ==   3);
	BOOST_TEST(binom(3, 3) ==   1);

	// 4th order  expansion
	BOOST_TEST(binom(4, 0) ==   1);
	BOOST_TEST(binom(4, 1) ==   4);
	BOOST_TEST(binom(4, 2) ==   6);
	BOOST_TEST(binom(4, 3) ==   4);
	BOOST_TEST(binom(4, 4) ==   1);

	// 5th order  expansion
	BOOST_TEST(binom(5, 0) ==   1);
	BOOST_TEST(binom(5, 1) ==   5);
	BOOST_TEST(binom(5, 2) ==  10);
	BOOST_TEST(binom(5, 3) ==  10);
	BOOST_TEST(binom(5, 4) ==   5);
	BOOST_TEST(binom(5, 5) ==   1);

	// 6th order  expansion
	BOOST_TEST(binom(6, 0) ==   1);
	BOOST_TEST(binom(6, 1) ==   6);
	BOOST_TEST(binom(6, 2) ==  15);
	BOOST_TEST(binom(6, 3) ==  20);
	BOOST_TEST(binom(6, 4) ==  15);
	BOOST_TEST(binom(6, 5) ==   6);
	BOOST_TEST(binom(6, 6) ==   1);
}

BOOST_AUTO_TEST_CASE(test10_binom_pascal_triangle_nok)
{
	BOOST_TEST(nok_hist(0, 1) ==   0);
	BOOST_TEST(nok_hist(0, 2) ==   0);
	BOOST_TEST(nok_hist(0, 3) ==   0);

	// the top
	BOOST_TEST(nok_hist(1, 1) ==   1);

	// quadratic expansion
	BOOST_TEST(nok_hist(2, 0) ==   1);
	BOOST_TEST(nok_hist(2, 1) ==   2);
	BOOST_TEST(nok_hist(2, 2) ==   1);

	// cubic expansion
	BOOST_TEST(nok_hist(3, 0) ==   1);
	BOOST_TEST(nok_hist(3, 1) ==   3);
	BOOST_TEST(nok_hist(3, 2) ==   3);
	BOOST_TEST(nok_hist(3, 3) ==   1);

	// 4th order  expansion
	BOOST_TEST(nok_hist(4, 0) ==   1);
	BOOST_TEST(nok_hist(4, 1) ==   4);
	BOOST_TEST(nok_hist(4, 2) ==   6);
	BOOST_TEST(nok_hist(4, 3) ==   4);
	BOOST_TEST(nok_hist(4, 4) ==   1);

	// 5th order  expansion
	BOOST_TEST(nok_hist(5, 0) ==   1);
	BOOST_TEST(nok_hist(5, 1) ==   5);
	BOOST_TEST(nok_hist(5, 2) ==  10);
	BOOST_TEST(nok_hist(5, 3) ==  10);
	BOOST_TEST(nok_hist(5, 4) ==   5);
	BOOST_TEST(nok_hist(5, 5) ==   1);

	// 6th order  expansion
	BOOST_TEST(nok_hist(6, 0) ==   1);
	BOOST_TEST(nok_hist(6, 1) ==   6);
	BOOST_TEST(nok_hist(6, 2) ==  15);
	BOOST_TEST(nok_hist(6, 3) ==  20);
	BOOST_TEST(nok_hist(6, 4) ==  15);
	BOOST_TEST(nok_hist(6, 5) ==   6);
	BOOST_TEST(nok_hist(6, 6) ==   1);
}
