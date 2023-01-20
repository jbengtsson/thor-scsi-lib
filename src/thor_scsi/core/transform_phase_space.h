#ifndef _THOR_SCSI_CORE_TRANSFORM_PHASE_SPACE_H_
#define _THOR_SCSI_CORE_TRANSFORM_PHASE_SPACE_H_ 1

#include <tps/enums.h>
#include <gtpsa/utils.hpp>
#include <thor_scsi/core/transform.h>

namespace thor_scsi::core {
    template <class C>
	class PhaseSpaceGalilean2DTransformKnobbed : public Galilean2DTransformKnobbed<C> {
	public:
		/*
		 * Todo
		 *  Review if can be replaced with Eigen even is tpsa is used
		 *
		 *  Review if it can be split up for position and direction
		 */
		inline PhaseSpaceGalilean2DTransformKnobbed() : Galilean2DTransformKnobbed<C>() {}
		inline virtual ~PhaseSpaceGalilean2DTransformKnobbed() {};
		PhaseSpaceGalilean2DTransformKnobbed(PhaseSpaceGalilean2DTransformKnobbed&& O) = default;
		PhaseSpaceGalilean2DTransformKnobbed& operator= (const PhaseSpaceGalilean2DTransformKnobbed &O) {
			(Galilean2DTransformKnobbed<C>&)(*this) = (const Galilean2DTransformKnobbed<C>&) (O);
			return *this;
		}
		template<typename T>
		inline void forward(T & ps){
			forwardTranslation(ps);
			forwardRotation(ps);
		}

		template<typename T>
		inline void backward(T & ps){
			backwardRotation(ps);
			backwardTranslation(ps);
		}

		/*
		  template<typename T>
		  void GtoL(ss_vect<T> &ps, std::vector<double> &S, std::vector<double> &R,
		  const double c0, const double c1, const double s1)
		  {
		  ss_vect<T> ps1;

		  // Simplified rotated p_rot: R^-1(theta_des) prot(phi/2) R(theta_des).
		  ps[px_] += c1; ps[py_] += s1;
		  // Eucluclidian transformation:
		  //   first Translate,
		  ps[x_] -= S[X_]; ps[y_] -= S[Y_];
		  //   then Rotate.
		  ps1 = ps;
		  ps[x_]  =  R[X_]*ps1[x_]  + R[Y_]*ps1[y_];
		  ps[px_] =  R[X_]*ps1[px_] + R[Y_]*ps1[py_];
		  ps[y_]  = -R[Y_]*ps1[x_]  + R[X_]*ps1[y_];
		  ps[py_] = -R[Y_]*ps1[px_] + R[X_]*ps1[py_] ;
		  // Simplified p_rot.
		  ps[px_] -= c0;
		  // Phase space vector is now in magnet's local coordinates.
		  }
		*/
		template<typename T>
		inline void forwardTranslation(T & ps){
			ps[x_] -= this->m_dS[X_];
			ps[y_] -= this->m_dS[Y_];
		}

		template<typename T>
		/**
		 * todo: external storage ....
		 */
		inline void forwardRotation(T& ps){
		        T ps1 = clone(ps);
			auto& R =  this->m_dT;
			ps[x_]  =  R[X_] * ps1[x_]  + R[Y_] * ps1[y_];
			ps[px_] =  R[X_] * ps1[px_] + R[Y_] * ps1[py_];
			ps[y_]  = -R[Y_] * ps1[x_]  + R[X_] * ps1[y_];
			ps[py_] = -R[Y_] * ps1[px_] + R[X_] * ps1[py_] ;
		}

		/*
                template<typename T>
                void LtoG(ss_vect<T> &ps, std::vector<double> &S, std::vector<double> &R,
                	  double c0, double c1, double s1)
                {
                  ss_vect<T> ps1;

                  // Reverse of GtoL, with inverted Euclidian.
                  // Simplified p_rot.
                  ps[px_] -= c0;
                  // Inverted Eucluclidian transformation:
                  //   first Rotate.
                  ps1 = ps;
                  ps[x_]  = R[X_]*ps1[x_]  - R[Y_]*ps1[y_];
                  ps[px_] = R[X_]*ps1[px_] - R[Y_]*ps1[py_];
                  ps[y_]  = R[Y_]*ps1[x_]  + R[X_]*ps1[y_];
                  ps[py_] = R[Y_]*ps1[px_] + R[X_]*ps1[py_];
                  //   then Translate.
                  ps[x_] += S[X_]; ps[y_] += S[Y_];
                  // Rotated p_rot.
                  ps[px_] += c1; ps[py_] += s1;
                }
		*/
		template<typename T>
		inline void backwardRotation(T& ps){
		        T ps1 = ps.clone();
			auto& R = this->m_dT;
			ps[x_]  = R[X_] * ps1[x_]  - R[Y_] * ps1[y_];
			ps[px_] = R[X_] * ps1[px_] - R[Y_] * ps1[py_];
			ps[y_]  = R[Y_] * ps1[x_]  + R[X_] * ps1[y_];
			ps[py_] = R[Y_] * ps1[px_] + R[X_] * ps1[py_];
		}
		template<typename T>
		inline void backwardTranslation(T & ps){
			ps[x_] += this->m_dS[X_];
			ps[y_] += this->m_dS[Y_];
		}

	};

	/**
	 * @brief PRot Transformation of Phase Space : Mixin
	 *
	 * As it contains steps to be executed before and after the
	 * Galileian Transformation it does not seem to make sense to support it as
	 * standalone application
	 */
    template <class C>
	class PhaseSpacePRotTransformMixinKnobbed: public PRotTransformKnobbed<C> {
	public:
		inline PhaseSpacePRotTransformMixinKnobbed(void) : PRotTransformKnobbed<C>() {}
		inline virtual ~PhaseSpacePRotTransformMixinKnobbed(void) {};

		PhaseSpacePRotTransformMixinKnobbed(PhaseSpacePRotTransformMixinKnobbed&& O)  = default;

		PhaseSpacePRotTransformMixinKnobbed& operator=(const PhaseSpacePRotTransformMixinKnobbed& O){
			static_cast<PRotTransform&>(*this) = static_cast<const PRotTransform&>(O);
			return *this;
		}
        // explicit template for evaluating double vector with tpsa argument
		template<typename T>
        inline void forwardStep1(T & ps){
            // Simplified rotated p_rot: R^-1(theta_des) prot(phi/2) R(theta_des).
            ps[px_] += this->c1;
            ps[py_] += this->s1;
        }
		template<typename T>
        inline void forwardStep2(T & ps){
            // Simplified p_rot.
            ps[px_] -= this->c0;
        }
		template<typename T>
		inline void backwardStep1(T & ps){
			// Reverse of GtoL, with inverted Euclidian.
			// Simplified p_rot.
			ps[px_] -= this->c0;
		}
		template<typename T>
		inline void backwardStep2(T & ps){
			// Rotated p_rot.
			ps[px_] += this->c1;
			ps[py_] += this->s1;
		}
	};



    template<class C>
	class PhaseSpaceGalileanPRot2DTransformKnobbed : public PhaseSpaceGalilean2DTransformKnobbed<C>, PhaseSpacePRotTransformMixinKnobbed<C> {
	public:
		inline PhaseSpaceGalileanPRot2DTransformKnobbed(void) :
                PhaseSpaceGalilean2DTransformKnobbed<C>(),
			PhaseSpacePRotTransformMixinKnobbed<C>() {}
		inline virtual ~PhaseSpaceGalileanPRot2DTransformKnobbed() {};

		inline PhaseSpaceGalileanPRot2DTransformKnobbed(PhaseSpaceGalileanPRot2DTransformKnobbed&& o) :
			PhaseSpaceGalilean2DTransformKnobbed<C>(std::move(o)),
			PhaseSpacePRotTransformMixinKnobbed<C>(std::move(o))
			{}

		PhaseSpaceGalileanPRot2DTransformKnobbed& operator= (const PhaseSpaceGalileanPRot2DTransformKnobbed& o){
			static_cast<PhaseSpaceGalilean2DTransformKnobbed<C>&>(*this) = static_cast<const PhaseSpaceGalilean2DTransformKnobbed<C>&>(o);
			static_cast<PhaseSpacePRotTransformMixinKnobbed<C>&>(*this) = static_cast<const PhaseSpacePRotTransformMixinKnobbed<C>&>(o);
			return *this;
		}

		template<typename T>
		inline void forward(T & ps){
			PhaseSpacePRotTransformMixinKnobbed<C>::forwardStep1(ps);
			PhaseSpaceGalilean2DTransformKnobbed<C>::forward(ps);
			PhaseSpacePRotTransformMixinKnobbed<C>::forwardStep2(ps);
		}
		template<typename T>
		inline void backward(T & ps){
			PhaseSpacePRotTransformMixinKnobbed<C>::backwardStep1(ps);
			PhaseSpaceGalilean2DTransformKnobbed<C>::backward(ps);
			PhaseSpacePRotTransformMixinKnobbed<C>::backwardStep2(ps);
		}
	};

    typedef Galilean2DTransformKnobbed<StandardDoubleType> Galilean2DTransform;
    typedef PhaseSpaceGalilean2DTransformKnobbed<StandardDoubleType> PhaseSpaceGalilean2DTransform;

    // preparations: ss_vect uses double but parameter is tpsa
    // currently a hack...
    // to be removed when class templates are used for the elements
    typedef PhaseSpacePRotTransformMixinKnobbed<TpsaVariantType> PhaseSpacePRotTransformTpsa;
    template<> template<>
    inline void PhaseSpacePRotTransformTpsa::forwardStep1(gtpsa::ss_vect<double> & ps){
        // Simplified rotated p_rot: R^-1(theta_des) prot(phi/2) R(theta_des).
        ps[px_] += gtpsa::cst(this->c1);
        ps[py_] += gtpsa::cst(this->s1);
    }

    template<>template<>
    inline void PhaseSpacePRotTransformTpsa::forwardStep2(gtpsa::ss_vect<double> & ps){
        // Simplified p_rot.
        ps[px_] -= gtpsa::cst(this->c0);
    }


}
#endif //_THOR_SCSI_CORE_TRANSFORM_PHASE_SPACE_H_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
