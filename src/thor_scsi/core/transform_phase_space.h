#ifndef _THOR_SCSI_CORE_TRANSFORM_PHASE_SPACE_H_
#define _THOR_SCSI_CORE_TRANSFORM_PHASE_SPACE_H_ 1

#include <tps/enums.h>
#include <thor_scsi/core/transform.h>

namespace thor_scsi::core {
	class PhaseSpaceGalilean2DTransform : public Galilean2DTransform {
	public:
		/*
		 * Todo
		 *  Review if can be replaced with Eigen even is tpsa is used
		 *
		 *  Review if it can be split up for position and direction
		 */
		inline PhaseSpaceGalilean2DTransform() : Galilean2DTransform() {}
		inline virtual ~PhaseSpaceGalilean2DTransform() {};
		PhaseSpaceGalilean2DTransform(PhaseSpaceGalilean2DTransform&& O) = default;
		PhaseSpaceGalilean2DTransform& operator= (const PhaseSpaceGalilean2DTransform &O) {
			(Galilean2DTransform&)(*this) = (const Galilean2DTransform&) (O);
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
			T ps1 = ps;
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
			T ps1 = ps;
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
	class PhaseSpacePRotTransformMixin: public PRotTransform {
	public:
		inline PhaseSpacePRotTransformMixin(void) : PRotTransform() {}
		inline virtual ~PhaseSpacePRotTransformMixin(void) {};

		PhaseSpacePRotTransformMixin(PhaseSpacePRotTransformMixin&& O)  = default;

		PhaseSpacePRotTransformMixin& operator=(const PhaseSpacePRotTransformMixin& O){
			static_cast<PRotTransform&>(*this) = static_cast<const PRotTransform&>(O);
			return *this;
		}
		template<typename T>
		inline void forwardStep1(T & ps){
			// Simplified rotated p_rot: R^-1(theta_des) prot(phi/2) R(theta_des).
			ps[px_] += this->c1; ps[py_] += this->s1;
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

	class PhaseSpaceGalileanPRot2DTransform : PhaseSpaceGalilean2DTransform, PhaseSpacePRotTransformMixin {
	public:
		inline PhaseSpaceGalileanPRot2DTransform(void) :
			PhaseSpaceGalilean2DTransform(),
			PhaseSpacePRotTransformMixin() {}
		inline virtual ~PhaseSpaceGalileanPRot2DTransform() {};

		inline PhaseSpaceGalileanPRot2DTransform(PhaseSpaceGalileanPRot2DTransform&& o) :
			PhaseSpaceGalilean2DTransform(std::move(o)),
			PhaseSpacePRotTransformMixin(std::move(o))
			{}

		PhaseSpaceGalileanPRot2DTransform& operator= (const PhaseSpaceGalileanPRot2DTransform& o){
			static_cast<PhaseSpaceGalilean2DTransform&>(*this) = static_cast<const PhaseSpaceGalilean2DTransform&>(o);
			static_cast<PhaseSpacePRotTransformMixin&>(*this) = static_cast<const PhaseSpacePRotTransformMixin&>(o);
			return *this;
		}

		template<typename T>
		inline void forward(T & ps){
			PhaseSpacePRotTransformMixin::forwardStep1(ps);
			PhaseSpaceGalilean2DTransform::forward(ps);
			PhaseSpacePRotTransformMixin::forwardStep2(ps);
		}
		template<typename T>
		inline void backward(T & ps){
			PhaseSpacePRotTransformMixin::backwardStep1(ps);
			PhaseSpaceGalilean2DTransform::backward(ps);
			PhaseSpacePRotTransformMixin::backwardStep2(ps);
		}
	};
}
#endif //_THOR_SCSI_CORE_TRANSFORM_PHASE_SPACE_H_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
