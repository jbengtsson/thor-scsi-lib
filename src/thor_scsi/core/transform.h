#ifndef _THOR_SCSI_CORE_TRANSFORM_H_
#define _THOR_SCSI_CORE_TRANSFORM_H_ 1

#include <array>
#include <cmath>
#include <ostream>


/*
 * Should I separate the coefficients from the implementation ...
 */
namespace thor_scsi::core {
	/**
	 * Candidate to be replaced by quaterions
	 *
	 * Todo:
	 *    - Memory management. smart pointers?
	 *
	 * Note:
	 *      It does not compute the tranform it self, at is just here to handle
	 *      phase space (state) vectors. Implemented in an derived class
	 *
	 */
       class Galilean2DTransform {
       public:
		///< Euclidian Group: dx, dy
		inline Galilean2DTransform(void){
			setdS(0.0, 0.0);
			setRoll(0.0);
		}
		inline Galilean2DTransform(Galilean2DTransform&& O) :
			m_dS(std::move(O.m_dS)),
			m_dT(std::move(O.m_dT))
			{}

		~Galilean2DTransform(){};

		inline void setdS(const double dx, const double dy)  {
			m_dS[0] = dx;
			m_dS[1] = dy;
		}

		///< Euclidian Group: Roll angle
		inline void setRoll(const double roll)  {
			m_dT[0] = cos(roll);
			m_dT[1] = sin(roll);
		}

		///< Euclidian Group: Roll angle
		inline double getRoll(void) const {
			return atan2(m_dT[1], m_dT[0]);
		}
		inline double getDx(void) const {
			return m_dS[0];
		}
		inline double getDy(void) const {
			return m_dS[1];
		}
		inline void setDx(const double x){
			m_dS[0] = x;
		}
		inline void setDy(const double y){
			m_dS[1] = y;
		}

		///< Euclidian Group: dx, dy
		inline const std::array<double, 2>& getdS(void)  {
			return m_dS;
		}

		inline const std::array<double, 2>& getdT(void)  {
			return m_dT;
		}

		inline void show(std::ostream& strm, int level) const {
			strm << "dx = " << this->getDx() << ", "
			     << "dy = " << this->getDy() << ", "
			     << "roll = " << this->getRoll();
		}

	protected:
		std::array<double, 2>
		m_dS{0e0, 0e0},              ///< Transverse displacement.
				m_dT{0e0, 0e0};              ///< part of rotation matrix = (cos(dT), sin(dT)).
	};

	inline
	std::ostream& operator<<(std::ostream& strm, const Galilean2DTransform& tf)
	{
		tf.show(strm, 0);
		return strm;
	}

	/**
	 * This needs a proper explanation
	 *
	 * What kind of transform
	 *
	 * naming of the parameters?
	 *
	 * Reference to CAD Tool document ?
	 */
	class PRotTransform {
	public:
		inline PRotTransform(void){
			this->c0 = 0e0;
			this->c1 = 0e0;
			this->s1 = 0e0;
		}

		inline PRotTransform(PRotTransform&& O){
			this->c0 = O.c0;
			this->c1 = O.c1;
			this->s1 = O.s1;
		}

		inline void setC0(const double val){this->c0 = val;}
		inline void setC1(const double val){this->c1 = val;}
		inline void setS1(const double val){this->s1 = val;}
		inline double getC0(void) const {return this->c0;}
		inline double getC1(void) const {return this->c1;}
		inline double getS1(void) const {return this->s1;}

		inline void show(std::ostream& strm, int level) const {
			strm << "c0 = " << this->getC0() << ", "
			     << "c1 = " << this->getC1() << ", "
			     << "s1 = " << this->getS1();
		}

	private:
		double c0, c1, s1;
	};

	inline
	std::ostream& operator<<(std::ostream& strm, const PRotTransform& tf)
	{
		tf.show(strm, 0);
		return strm;
	}


}


#endif // _THOR_SCSI_CORE_TRANSFORM_H_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
