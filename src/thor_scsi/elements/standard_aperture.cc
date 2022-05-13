#include <thor_scsi/elements/standard_aperture.h>
#include <algorithm>

namespace tse = thor_scsi::elements;

tse::CircularAperture::CircularAperture(const double radius, const double x, const double y)
{
	this->m_x = x;
	this->m_y = y;
	this->m_radius = radius;
}

double tse::CircularAperture::isWithin(const double x, const double y) const
{

	const double dx = x - this->m_x, dy = y - this->m_y;
	const double radius = sqrt(dx*dx + dy*dy);
	return this->m_radius - radius;
}

void tse::CircularAperture::show(std::ostream& strm, int level) const
{
	TwoDimensionalAperture::show(strm, level);
	strm << " centre (" << this->m_x << ", " << this->m_y << "), radius " << this->m_radius;
}

tse::RectangularAperture::RectangularAperture(const double width, const double height, const double x, const double y)
{
	this->m_x = x;
	this->m_y = y;
	this->m_width2 = width / 2e0;
	this->m_height2 = height / 2e0;
}

double tse::RectangularAperture::isWithin(const double x, const double y) const
{
	double x_off, y_off;
	this->distances(x, y, &x_off, &y_off);
	/*
	 * what to return ...
	 *
	 *
	 * If negative return this value ... out in this dimension
	 * If both negative return the one more negative
	 *  why ... this value has to be optimised first
	 * otherwise ... again the smaller one
	 */
	return std::min(x_off, y_off);
}

void tse::RectangularAperture::show(std::ostream& strm, int level) const
{
	TwoDimensionalAperture::show(strm, level);
	strm << " centre (" << this->m_x << ", " << this->m_y << "),"
	     << " width " << this->m_width2 * 2
	     << " height " << this->m_height2 * 2;
}

/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
