#ifndef _THOR_SCSI_CORE_TRANSFORM_H_
#define _THOR_SCSI_CORE_TRANSFORM_H_ 1

#include <array>
#include <cmath>

namespace thor_scsi{
	namespace core {
                 /**
		  * Candidate to be replaced by quaterions
		  */

		class Euclidian2DTranform {
		public:
			///< Eucledian Group: dx, dy
			inline Euclidian2DTranform(void){
				setdS(0.0, 0.0);
				setRoll(0.0);
			}
			~Euclidian2DTranform(){};

			inline void setdS(const double dx, const double dy)  {
				m_dS[0] = dx;
				m_dS[1] = dy;
			}
			///< Eucledian Group: Roll angle
			inline void setRoll(const double roll)  {
				m_dT[0] = cos(roll);
				m_dT[1] = sin(roll);
			}

			///< Eucledian Group: Roll angle
			inline const double getRoll(void)  {
				return atan2(m_dT[1], m_dT[0]);
			}
			inline const double getDx(void){
				return m_dS[0];
			}
			inline const double getDy(void){
				return m_dS[1];
			}
			inline void setDx(const double x){
				m_dS[0] = x;
			}
			inline void setDy(const double y){
				m_dS[1] = y;
			}
			///< Eucledian Group: dx, dy
			inline const std::array<double, 2>& getdS(void)  {
				return m_dS;
			}

			inline const std::array<double, 2>& getdT(void)  {
				return m_dT;
			}

			Euclidian2DTranform(Euclidian2DTranform&& O) :
				m_dS(std::move(O.m_dS)),
				m_dT(std::move(O.m_dT))
				{}

		private:
			std::array<double, 2>
			        m_dS{0e0, 0e0},              ///< Transverse displacement.
				m_dT{0e0, 0e0};              ///< dT = (cos(dT), sin(dT)).
		};
	}
}


#endif // _THOR_SCSI_CORE_TRANSFORM_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
