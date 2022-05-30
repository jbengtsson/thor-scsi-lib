#ifndef _THOR_SCSI_CORE_STANDARD_APERTURE_H_
#define _THOR_SCSI_CORE_STANDARD_APERTURE_H_ 1
#include <cmath>
#include <ostream>

#include <thor_scsi/core/aperture.h>

namespace thor_scsi::elements {
	/**
	 * @brief: a geometric boundary condition
	 */
	class CircularAperture : public thor_scsi::core::TwoDimensionalAperture{
	public:
		CircularAperture(const double radius, const double x=0, const double y=0);

		/*
		 *
		 * Returns distance to
		 */
		virtual double isWithin(const double x, const double y) const override final;
		virtual void show(std::ostream&, int level) const override;
		virtual const char * type_name(void) const override {
			return "CircularAperture";
		};

	private:
		double m_x = 0e0, m_y = 0e0, m_radius = 0e0;
	};


	class RectangularAperture : public thor_scsi::core::TwoDimensionalAperture {
	public:
		/*
		 * @brief: Rectangle center, width and height
		 *
		 * Could review if an angle should be added for a rotated rectangle
		 */
		RectangularAperture(const double width, const double height, const double x=0, const double y=0);

		/*
		 * @brief Returns distance of point to boundary
		 */
		inline void distances(const double x, const double y, double *x_off, double *y_off) const {
			double dx = x - this->m_x, dy = y - this->m_y;
			dx = std::abs(dx);
			dy = std::abs(dy);

			// how much off the aperture
			*x_off = this->m_width2 - dx;
			*y_off = this->m_height2 - dy;
		}

		virtual double isWithin(const double x, const double y) const override final;
		virtual void show(std::ostream&, int level) const override;
		virtual const char * type_name(void) const override {
			return "RectangularAperture";
		};

	private:
		double m_x = 0e0, m_y = 0e0, m_width2 = 0e0, m_height2 = 0e0;
	};

}
#endif /* _THOR_SCSI_STANDARD_APERTURE_H_ */
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
