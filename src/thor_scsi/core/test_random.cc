#define BOOST_TEST_MODULE random
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#include <random>
#include <ostream>

BOOST_AUTO_TEST_CASE(test01_std_random)
{
	// check the range of the generator
	// std::default_random_engine generator;
	std::mt19937_64 generator;
	std::uniform_real_distribution<double> distribution(0, 1.0);

	// the default random generator on my computer is rather small
	unsigned long checkval = (1L<<31);

	std::cerr << "Generator range: "  << generator.min()
		  << ".." << generator.max()
		  << " check value " << checkval
		  << std::endl;

	BOOST_CHECK(generator.min() >= 0);
	BOOST_CHECK(generator.max() > checkval);

	{
		double number = distribution(generator);
		std::cerr << "\t first (default) sample: "  << number << std::endl;

		BOOST_CHECK(number >= 0);
		BOOST_CHECK(number <= 1);
	}
}


BOOST_AUTO_TEST_CASE(test02_real_random_seed)
{

	std::mt19937_64 generator;
	std::uniform_real_distribution<double> distribution(0, 1.0);
	std::random_device rdev;

	unsigned long long seed1 = rdev();

	std::cerr << "From random device: seed used "  << unsigned(seed1) << std::endl;

	generator.seed(seed1);
	{
		double number = distribution(generator);
		std::cerr << "\t first sample  "  << number << std::endl;

		BOOST_CHECK(number >= 0);
		BOOST_CHECK(number <= 1);
	}
}


BOOST_AUTO_TEST_CASE(test03_normal_distribution)
{

	double mean = 1.0,  stddev = 2.0;

	std::mt19937_64 generator;
	std::normal_distribution<double> distribution(mean, stddev);


	BOOST_CHECK_CLOSE(distribution.mean(), mean, 1e-12);
	BOOST_CHECK_CLOSE(distribution.stddev(), stddev, 1e-12);

 	{
		double number = distribution(generator);
		std::cerr << "Normal distribution N(" << mean << "," << stddev
			  << "): sampled r "  << number << std::endl;

	}
}

BOOST_AUTO_TEST_CASE(test10_normal_realistic_alignment)
{
	double mean = .4,  stddev = 30e-6;

	std::mt19937_64 generator;
	std::normal_distribution<double> distribution(mean, stddev);

	BOOST_CHECK_CLOSE(distribution.mean(), mean, 1e-12);
	BOOST_CHECK_CLOSE(distribution.stddev(), stddev, 1e-12);

	{
		double number = distribution(generator);
		std::cerr << "Normal distribution N(" << mean << "," << stddev
			  << "): sampled r "  << number << std::endl;

	}
}
