#ifndef _THOR_SCSI_CORE_TRANSFORM_PHASE_SPACE_H_
#define _THOR_SCSI_CORE_TRANSFORM_PHASE_SPACE_H_ 1

#include <tps/enums.h>
#include <gtpsa/utils.hpp>
#include <thor_scsi/core/transform.h>
#include <gtpsa/ss_vect.h>

namespace thor_scsi::core {
    /**
     *
     * Deliberatly not using gtpsa::cst as other conversions could be required
     * @todo where to put these conversions ?
     */
    template <typename Ti, typename To>
    inline void to_base_type(const Ti* input, To* output) {
        if((const void *)(input) == (void *)(output)){
            return;
        }
        *output =  *input;
    }
    template <>
    inline void to_base_type(const gtpsa::tpsa* input, double *output) {  *output = input->cst(); }
    template <>
    inline void to_base_type(const gtpsa::TpsaOrDouble* input, double *output) {  *output = input->cst(); }
    template <>
    inline void to_base_type(const gtpsa::TpsaOrDouble* input, gtpsa::tpsa *output) {  *output = input->toTpsaType(*output); }
    //template <>
    //inline void to_base_type(const gtpsa::GTpsaOrBase<gtpsa::TpsaVariantDoubleTypes>* input, gtpsa::tpsa *output) {  *output = input->asTpsaType(); }


    template <class C>
	class PhaseSpaceGalilean2DTransformKnobbed : public Galilean2DTransformKnobbed<C> {
	using double_type = typename C::double_type;
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
		inline void forward(gtpsa::ss_vect<T> & ps){
			forwardTranslation(ps);
			forwardRotation(ps);
		}

		template<typename T>
		inline void backward(gtpsa::ss_vect<T> & ps){
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
		inline void forwardTranslation(gtpsa::ss_vect<T> & ps){
		    auto dx = gtpsa::same_as_instance(ps[x_]);
		    auto dy = gtpsa::same_as_instance(ps[y_]);
		    to_base_type<double_type, T>(&this->m_dS[X_], &dx);
		    to_base_type<double_type, T>(&this->m_dS[Y_], &dy);
		    ps[x_] -= dx;
		    ps[y_] -= dy;
		}

		template<typename T>
		/**
		 * todo: external storage ....
		 */
		inline void forwardRotation(gtpsa::ss_vect<T>& ps){
			const auto& rx =  this->m_dT[X_];
			const auto& ry =  this->m_dT[Y_];

			gtpsa::ss_vect<T> ps1 = clone(ps);
			auto x  =   rx * ps1[x_ ] + ry * ps1[y_ ];
			auto px =   rx * ps1[px_] + ry * ps1[py_];
			auto y  =  -ry * ps1[x_ ] + rx * ps1[y_ ];
			auto py =  -ry * ps1[px_] + rx * ps1[py_];
			to_base_type(&x , &ps[x_] );
			to_base_type(&px, &ps[px_]);
			to_base_type(&y , &ps[y_] );
			to_base_type(&py, &ps[py_]);
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
		inline void backwardRotation(gtpsa::ss_vect<T>& ps){
		        gtpsa::ss_vect<T> ps1 = ps.clone();
			auto& R = this->m_dT;
			auto x  = R[X_] * ps1[x_]  - R[Y_] * ps1[y_];
			auto px = R[X_] * ps1[px_] - R[Y_] * ps1[py_];
			auto y  = R[Y_] * ps1[x_]  + R[X_] * ps1[y_];
			auto py = R[Y_] * ps1[px_] + R[X_] * ps1[py_];

			to_base_type(&x , &ps[x_] );
			to_base_type(&px, &ps[px_]);
			to_base_type(&y , &ps[y_] );
			to_base_type(&py, &ps[py_]);

		}
		template<typename T>
		inline void backwardTranslation(gtpsa::ss_vect<T> & ps){
            auto dx = gtpsa::same_as_instance(ps[x_]);
            auto dy = gtpsa::same_as_instance(ps[y_]);
            to_base_type(&this->m_dS[X_], &dx);
            to_base_type(&this->m_dS[Y_], &dy);
			ps[x_] += dx;
			ps[y_] += dy;
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
        inline void forwardStep1(gtpsa::ss_vect<T> & ps){
            // Simplified rotated p_rot: R^-1(theta_des) prot(phi/2) R(theta_des).
            T dpx = gtpsa::same_as_instance(ps[px_]);
            T dpy = gtpsa::same_as_instance(ps[py_]);
            to_base_type(&this->c1, &dpx);
            to_base_type(&this->s1, &dpy);
            ps[px_] += dpx;
            ps[py_] += dpy;
        }
	template<typename T>
        inline void forwardStep2(gtpsa::ss_vect<T> & ps){
            // Simplified p_rot.
        T dpr = gtpsa::same_as_instance(ps[px_]);
        to_base_type(&this->c0, &dpr);
            ps[px_] -= dpr;
        }
		template<typename T>
		inline void backwardStep1(gtpsa::ss_vect<T> & ps){
			// Reverse of GtoL, with inverted Euclidian.
			// Simplified p_rot.
#warning "Investigate sign of backward step 1 versus forward step2! "
            T dpr = gtpsa::same_as_instance(ps[px_]);
            to_base_type(&this->c0, &dpr);
			ps[px_] -= dpr;
		}
		template<typename T>
		inline void backwardStep2(gtpsa::ss_vect<T> & ps){
#warning "Investigate sign of forward step 1 versus backward step2 ! "
			// Rotated p_rot.
            T dpx = gtpsa::same_as_instance(ps[px_]);
            T dpy = gtpsa::same_as_instance(ps[py_]);
            to_base_type(&this->c1, &dpx);
            to_base_type(&this->s1, &dpy);
			ps[px_] += dpx;
			ps[py_] += dpy;
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
		inline void forward(gtpsa::ss_vect<T> & ps){
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
