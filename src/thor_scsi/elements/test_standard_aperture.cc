#define BOOST_TEST_MODULE aperture
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <thor_scsi/elements/standard_aperture.h>
#include <cmath>

namespace tse = thor_scsi::elements;

BOOST_AUTO_TEST_CASE(test10_circular_aperture_center)
{
	auto ap = tse::CircularAperture(1);

	BOOST_CHECK_CLOSE(ap.isWithin( 0, 0), 1, 1e-12);
	BOOST_CHECK_SMALL(ap.isWithin(-1, 0),    1e-12);

}

BOOST_AUTO_TEST_CASE(test11_circular_aperture_off_x)
{
	auto ap = tse::CircularAperture(1, 1, 0);

	BOOST_CHECK_CLOSE(ap.isWithin(1,  0), 1, 1e-12);
	BOOST_CHECK_SMALL(ap.isWithin(0,  0),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin(2,  0),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin(1,  1),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin(1, -1),  1e-12);

}


BOOST_AUTO_TEST_CASE(test20_rectangular_aperture)
{
	auto ap = tse::RectangularAperture(2, 2);
	// center
	BOOST_CHECK_CLOSE(ap.isWithin( 0,  0), 1, 1e-12);

	// mid border
	BOOST_CHECK_SMALL(ap.isWithin( 1,  0),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin(-1,  0),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin( 0,  1),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin( 0, -1),  1e-12);

	// corners
	BOOST_CHECK_SMALL(ap.isWithin( 1,  1),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin(-1, -1),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin(-1,  1),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin( 1, -1),  1e-12);

	// inside
	BOOST_CHECK_CLOSE(ap.isWithin( 0.9,  0),  .1, 1e-12);
	BOOST_CHECK_CLOSE(ap.isWithin( 0.9, .5),  .1, 1e-12);
	// outside
	BOOST_CHECK_CLOSE(ap.isWithin( 1.1,  0), -.1, 1e-12);
	BOOST_CHECK_CLOSE(ap.isWithin( 1.1, .5), -.1, 1e-12);
}

BOOST_AUTO_TEST_CASE(test21_rectangular_aperture)
{
	auto ap = tse::RectangularAperture(2, 2, 1, 1);
	// center
	BOOST_CHECK_CLOSE(ap.isWithin( 1,  1), 1, 1e-12);

	// mid border
	BOOST_CHECK_SMALL(ap.isWithin( 2,  1),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin( 0,  1),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin( 1,  2),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin( 1,  0),  1e-12);

	// corners
	BOOST_CHECK_SMALL(ap.isWithin( 0,  2),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin( 2,  0),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin( 2,  2),  1e-12);
	BOOST_CHECK_SMALL(ap.isWithin( 0, -0),  1e-12);

	// inside
	BOOST_CHECK_CLOSE(ap.isWithin( 1.9,  1),  .1, 1e-12);
	BOOST_CHECK_CLOSE(ap.isWithin( 1.9, .5),  .1, 1e-12);
	// outside
	BOOST_CHECK_CLOSE(ap.isWithin( 2.1,  1), -.1, 1e-12);
	BOOST_CHECK_CLOSE(ap.isWithin( 2.1, .5), -.1, 1e-12);
}
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
