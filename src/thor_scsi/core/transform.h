#ifndef _THOR_SCSI_CORE_TRANSFORM_H_
#define _THOR_SCSI_CORE_TRANSFORM_H_ 1

#include <array>
#include <cmath>
#include <ostream>
#include <thor_scsi/core/multipole_types.h>

/* required for parameter study */
using gtpsa::sin;
using gtpsa::cos;

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

    template<class C, typename = typename C::double_type>
       class Galilean2DTransformKnobbed {
       protected:
           using double_type = typename C::double_type;
//#warning "Galilean transform not yet knobbed"
           std::array<
	       double_type
	       // double
                   , 2>
                   m_dS{0e0, 0e0},              ///< Transverse displacement.
           m_dT{0e0, 0e0};              ///< part of rotation matrix = (cos(dT), sin(dT)).

       public:
	       ///< Euclidian Group: dx, dy
	       inline Galilean2DTransformKnobbed(void){
		       setdS(0.0, 0.0);
			setRoll(0.0);
	       }
	       inline Galilean2DTransformKnobbed(Galilean2DTransformKnobbed&& O) :
		       m_dS(std::move(O.m_dS)),
		       m_dT(std::move(O.m_dT))
		       {}

	       virtual ~Galilean2DTransformKnobbed(){};
	       Galilean2DTransformKnobbed& operator= (const Galilean2DTransformKnobbed& O){
		       this->m_dS[0] = O.m_dS[0];
		       this->m_dS[1] = O.m_dS[1];
		       this->m_dT[0] = O.m_dT[0];
		       this->m_dT[1] = O.m_dT[1];
		       return *this;
	       }

		inline void setdS(const double_type dx, const double_type dy)  {
			m_dS[0] = dx;
			m_dS[1] = dy;
		}

		///< Euclidian Group: Roll angle
		inline void setRoll(const double_type roll)  {
			m_dT[0] = cos(roll);
			m_dT[1] = sin(roll);
		}

		///< Euclidian Group: Roll angle
		inline double_type getRoll(void) const {
			return atan2(m_dT[1], m_dT[0]);
		}
		inline double_type getDx(void) const {
			return m_dS[0];
		}
		inline double_type getDy(void) const {
			return m_dS[1];
		}
		inline void setDx(const double_type x){
			m_dS[0] = x;
		}
		inline void setDy(const double_type y){
			m_dS[1] = y;
		}

		///< Euclidian Group: dx, dy
		inline const std::array<double_type, 2>& getdS(void)  {
			return m_dS;
		}

		inline const std::array<double_type, 2>& getdT(void)  {
			return m_dT;
		}

		inline void show(std::ostream& strm, int level) const {
			strm << "dx = " << this->getDx() << ", "
			     << "dy = " << this->getDy() << ", "
			     << "roll = " << this->getRoll();
		}

       };

    template<class C>
	inline
	std::ostream& operator<<(std::ostream& strm, const Galilean2DTransformKnobbed<C>& tf)
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
    template<class C, typename = typename C::double_type>
	class PRotTransformKnobbed {
    protected:
        using double_type = typename C::double_type;
	public:
	    inline PRotTransformKnobbed(void)
		: c0(0e0)
		, c1(0e0)
		, s1(0e0)
		{}
	    inline PRotTransformKnobbed(PRotTransformKnobbed&& O)
		: c0(O.c0)
		, c1(O.c1)
		, s1(O.s1)
		{}

	    virtual ~PRotTransformKnobbed(void) {}

	    inline PRotTransformKnobbed& operator= (const PRotTransformKnobbed& O)
		{
			this->c0 = O.c0;
			this->c1 = O.c1;
			this->s1 = O.s1;
			return *this;
		}

		inline void setC0(const double_type val){this->c0 = val;}
		inline void setC1(const double_type val){this->c1 = val;}
		inline void setS1(const double_type val){this->s1 = val;}
		inline double_type getC0(void) const {return this->c0;}
		inline double_type getC1(void) const {return this->c1;}
		inline double_type getS1(void) const {return this->s1;}

		inline void show(std::ostream& strm, int level) const {
			strm << "c0 = " << this->getC0() << ", "
			     << "c1 = " << this->getC1() << ", "
			     << "s1 = " << this->getS1();
		}

// #warning "prot not yet parametrizable"
	protected:
        //double
	double_type
        c0, c1, s1;
	};

    template<class C>
	inline
	std::ostream& operator<<(std::ostream& strm, const PRotTransformKnobbed<C>& tf)
	{
		tf.show(strm, 0);
		return strm;
	}


    typedef Galilean2DTransformKnobbed<StandardDoubleType> Galilean2DTransform;
    typedef PRotTransformKnobbed<StandardDoubleType> PRotTransform;
}


#endif // _THOR_SCSI_CORE_TRANSFORM_H_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
