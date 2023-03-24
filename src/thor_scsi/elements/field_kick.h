#ifndef _THOR_SCSI_ELEMENTS_FIELD_KICK_H_
#define _THOR_SCSI_ELEMENTS_FIELD_KICK_H_

#include <thor_scsi/elements/element_local_coordinates.h>
#include <thor_scsi/elements/elements_enums.h>
// #include <thor_scsi/elements/enums.h>
#include <thor_scsi/elements/utils.h>
// move to API
#include <thor_scsi/core/exceptions.h>
#include <thor_scsi/elements/field_kick_api.h>
#include <thor_scsi/elements/radiation_delegate.h>
#include <cassert>

/* Calculate multipole kick. The kick is given by

                   e L      L delta      L x
	theta  = - --- B  + -------  -  -----  ,
             x     p    y     rho           2
                    0                    rho
                 e L
	theta  = --- B
             y   p    x
                  0
    where
                           ====
                           \                       n-1
	(B + iB  ) = B rho  >   (ia  + b ) (x + iy)
	  y    x           /       n    n
	                   ====
    where
			e      1
			-- = -----
			p    B rho
			 0
*/

/*

  Todo: ensure that irho is computed for a dipole sector bend in polar coordiantes
  (see t2lat.cc
      M->irho = (M->L != 0.0)? t*M_PI/180.0/M->L : t*M_PI/180.0;)

 */
namespace thor_scsi::elements {
    // forward declaration
    template<class>
    class FieldKickKnobbed;


    /**
 * Symplectic integrator
 * 2nd order
 *
 *  @f[ \frac{L}{2} \to bnl \to  \frac{L}{2} ]@
 *
 * Here 4th order: calculating coefficients
 *
 * 4th order:
 *
 *  @f[ d_1 L \to c_1  bnL \to d_1 L \to c_2 bnl \to d_1 L ]@
 *
 * with
 *
 * @f[ d_1 + d_2 + d_1 = L @f]
 *
 *  @f[c_1 + c_2 = 1 @f]
 *
 * d_2:  negative thus creates a negative drift
 *
 * @todo review interaction to configuration to see if integration intervals can be precomputed
 */
    template<class C>
    class FieldKickDelegate {
    public:
        FieldKickDelegate(FieldKickKnobbed<C> *parent){
            this->parent = parent;
        }
        FieldKickDelegate(void){
            this->parent = nullptr;
        }
        virtual ~FieldKickDelegate() {};

        inline void setNumberOfIntegrationSteps(const int n){
            this->integration_steps = n;
            this->computeIntegrationSteps();
        }

        inline int getNumberOfIntegrationSteps(void) const {
            return this->integration_steps;
        }

        inline double getLength(void) const {
            return parent->getLength();
        }

        void setParent(FieldKickKnobbed<C> * p){
            this->parent = p;
            this->computeIntegrationSteps();
        }

    protected:
        inline auto getFieldInterpolator(void) const {
            // required to cast parent to const ?
            const FieldKickKnobbed<C> * ptr = this->parent;
            if(!ptr){
                throw std::runtime_error("Unknown parent!");
            }
            return ptr->getFieldInterpolator();
        }

        double PL = 0.0;
        int integration_steps = 1;
        FieldKickKnobbed<C> *parent = nullptr;

    private:
        // should be not defined as pointers are available
        // FieldKickDelegate(const FieldKickDelegate& o) = delete;
        FieldKickDelegate& operator= (const FieldKickDelegate& o) = delete;

        virtual void computeIntegrationSteps(void) = 0;

    };

    template<class C>
    class FieldKickForthOrder : public FieldKickDelegate<C> {
    public:
        /**
         * @brief: calculate lengthes for drifts and kicks
         *
         * 4th order integration method is given by
         */
        void splitIntegrationStep(const double dL, double *dL1, double *dL2,
                                  double *dkL1, double *dkL2) const;

        //
        // as it is a templated function ... not defined virtual ...
        template<typename T>
        void _localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<T> &ps);

        inline std::unique_ptr<std::vector<double>> getDriftLength(void) const {
            auto res = std::make_unique<std::vector<double>>(2);
            res->at(0) = this->dL1;
            res->at(1) = this->dL2;
            return res;
        }

        inline std::unique_ptr<std::vector<double>> getKickLength(void) const {
            auto res = std::make_unique<std::vector<double>>(2);
            res->at(0) = this->dL1;
            res->at(1) = this->dL2;
            return res;
        }

        // not really useful yet ... need to push configuration into
        // objects
        // fix for Cartesian bend treatment
        // don't use this internal parametes yet
        void computeIntegrationSteps(void) override final;

    private:

        /*
         * c_1 = 1/(2*(2-2^(1/3))),    c_2 = (1-2^(1/3))/(2*(2-2^(1/3)))
         *  d_1 = 1/(2-2^(1/3)),        d_2 = -2^(1/3)/(2-2^(1/3))
         */

        const double c_1 = 1e0/(2e0*(2e0-thirdroot(2e0)));
        const double c_2 = 0.5e0 - c_1;
        const double d_1 = 2e0*c_1;
        const double d_2 = 1e0 - 2e0*d_1;

        // Consider one set for
        // an other set for Cartesian Bends
        double dL1 = 0e0, dL2 = 0e0, dkL1 = 0e0, dkL2 = 0e0;

    };

    /**
	 * Calculate multipole kick. The kick is given by
	 * @f[
	 *     \vartheta_x =  \frac{e}{p_\mathrm{0}} L B_y + \frac{L}{\rho}
	 *                  \left(\delta - \frac{x}{\rho} \right)
	 * @f]
	 * @f[
	 \vartheta_y =  \frac{e}{p_\mathrm{0}} L B_x
	 * @f]
	 * with
	 * @f[
	 * \frac{e}{p_\mathrm{0}} = \frac{1}{B\, \rho}
	 * @f]
	 *
	 * \verbatim embed:rst:leading-asterisk
	 *   The magnetic field expansion is represented by a
	 *   :any:`Field`Field2DInterpolation object:
	 *
	 * .. Warning:
	 *        Derived classes are responsible to compute h_bend and h_ref
	 *        from given fields and choice of coordinate system (polar or Cartesian)
	 *
	 * This code can treat thin lenses, i.e. componets that have a length = 0e0.
	 * Cartesian local coordinate system is used for thin lenses. For these h_bend = 0e0, and
	 * href = 0e0.
	 *
	 * .. Todo:
	 *    fix treating Porder
	 * \endverbatim
	 */
	const int HOMmax = 21;
	// typedef std::vector<double> MpoleArray;


	/**
	 *
	 * SI units are used internally
	 * apart from energy [GeV]
	 *
	 *
	 * PhasespaceChange
	 * FieldPropagate ?
	 *
	 * Assumption: The magnetic field is constant within this element.
	 *             Hence e.g. local curvature and gradient are constant
	 */
    template<class C>
	class FieldKickKnobbed : public FieldKickAPIKnobbed<C> {
    protected:
        using complex_type = typename C::complex_type;
        using double_type  = typename C::double_type;
	public:
		// using thor_scsi::core::ObservedState;
		/*
		 * @todo configuration for method: use int
		 *
		 * @todo field level
		 */
		// friend FieldKick* Mpole_Alloc(void);
		//	ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
		// void print(const std::string &);
		/// std::string repr(void);
		///< additional information of mpole for mpole
		///< more than just solely type name
		///< used by meth:`repr`
		// std::string repr_add(void);

		FieldKickKnobbed(const Config &config);
		virtual ~FieldKickKnobbed(){}
		FieldKickKnobbed(FieldKickKnobbed&& O);

		const char* type_name(void) const override { return "field_kick"; };
		virtual void show(std::ostream& strm, const int level) const override;
		void inline setIntegrationMethod(const int n){
			this->validateIntegrationMethod(n);
			this->Pmethod = n;
		}
		int  getIntegrationMethod(void) const {
			return this->Pmethod;
		}

		/*
		 * @brief: number of times the integration is repeated
		 */
		inline void setNumberOfIntegrationSteps(const int n){
			this->integ4O.setNumberOfIntegrationSteps(n);
		}

		inline int getNumberOfIntegrationSteps(void) const {
			return this->integ4O.getNumberOfIntegrationSteps();
		}

		/**
		 * @brief Use this element as thick one
		 *
		 * If not set true it will be a thin one
		 *
		 * see :meth:`isThick` for further details
		 */
		void inline asThick(const bool flag){
			this->Pthick = flag;
		}

		/**
		 * @brief true if a thick elment
		 *
		 *
		 * if true interpolated field value does not need to be
		 * multiplied with the length thus the field kick element is
		 * approximated as a very thin lens (length 0 0)
		 *
		 *
		 * \verbatim embed:rst:leading-asterisk
		 * Returns:
		 *       true if field interpolation values are used as integral
		 *       values, false if field interpolation values and length
		 *       are used for rectangular model
		 *
		 * .. Todo::
		 *
		 *     remove asIntegral as it is code duplication
		 *
		 * \endverbatim
		 *
		 *
		 */
		bool inline isThick(void) const {
			return this->Pthick;
		}

		/**
		 * @brief set length. if length == 0 interpolation will be
		 *        considered an integral one
		 *
		 * Todo: review if interface should be kept that manner
		 */
		virtual void inline setLength(const double& length) override final {
			// delegate ...
            FieldKickAPIKnobbed<C>::setLength(length);
			if(0e0 == gtpsa::cst(length)){
				this->asThick(false);
			} else {
				this->asThick(true);
			}
		}


		/*
		 * @brief: element takes geometric effects on the trajectory into account
		 *
		 * Checking that
		 * @f[ \frac{1}{\rho} = irho == 0 @f]
		 */
		inline bool assumingCurvedTrajectory(void) const {
			return !(this->Pirho == 0e0);
		}

		/**
		 * @Todo: clarify if angle should be treated in radian or degree
		 *
		 * @Todo: clarify relation to entrance and exit angle
		 */
		inline void setBendingAngle(const double angle){
			this->Pbending_angle = angle;
		}

		inline double getBendingAngle(void) const {
			return this->Pbending_angle;
		}

		inline void setEntranceAngle(const double angle) {
			this->PTx1 = angle;
		}

		inline double getEntranceAngle(void) const{
			return this->PTx1;
		}

		/**
		 * @todo: add test and check
		 */
		inline void setExitAngle(const double angle){
			this->PTx2 = angle;
		}

		inline double getExitAngle(void) const {
			return this->PTx2;
		}

		virtual void localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<double>      &ps) override final { _localPropagate(conf, ps);}
		virtual void localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<gtpsa::tpsa> &ps) override final { _localPropagate(conf, ps);}
	    /*
	        virtual void localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<tps>         &ps) override final {
			// _localPropagate(conf, ps);
			throw std::runtime_error("field kick local propagate not implemented for ss_vect<tps> ");
		}
	    */

        // Required as radiation is now handled by a delegate
        template<typename T>
        inline void thinKickAndRadiate(const thor_scsi::core::ConfigType &conf,
                                       const thor_scsi::core::Field2DInterpolationKnobbed<C>& intp,
                                       const double L, const double h_bend, const double h_ref,
                                       gtpsa::ss_vect<T> &ps);

	  private:
		template<typename T>
			void _localPropagate(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<T> &ps);

		template<typename T>
		        void _localPropagateThin(const thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<T> &ps);

		template<typename T>
		        void _localPropagateBody(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<T> &ps);


		void inline validateIntegrationMethod(const int n) const {
			switch(n){
			case Meth_Fourth:
				return;
			default:
				std::stringstream strm;
				strm << "Only implemented integration method 4 but found " << n;
				throw thor_scsi::NotImplemented(strm.str());
			}
		}

	  public:
		double
		Pbending_angle = 0e0,                     ///<  Todo: Already defined or combination of PTx1 and PTx2?
			PTx1 = 0e0,                      ///<  Bend angle [deg]:  hor. entrance angle
			PTx2 = 0e0,                      ///<  Bend angle [deg]: hor. exit angle.
			Pgap = 0e0;                      ///< Total magnet gap [m].

		/*
		 * see :any:`isThick` or :any:`asThick` for a description
		 */
    public:
        /**
         *
         * @todo should a check be made that forth order is requested ?
         */
        inline const FieldKickDelegate<C>& getFieldKickDelegator(void) const {
            return this->integ4O;
        }

        /*
         * radiation delegate part
         */
    public:
        typedef thor_scsi::elements::RadiationDelegateKickKnobbed<thor_scsi::elements::FieldKickAPIKnobbed<C>> rad_del_t;

    protected:
        std::shared_ptr<rad_del_t> rad_del;

    public:
		inline void setRadiationDelegate(std::shared_ptr<rad_del_t> p){
			this->rad_del = p;
		}
		inline auto getRadiationDelegate(void) const {
			return this->rad_del;
		}

	private:
		inline bool computeSynchrotronIntegrals(const thor_scsi::core::ConfigType &conf){
			return conf.emittance && !conf.Cavity_on;
		}
		inline auto _getRadiationDelegate(void) {
			return this->rad_del.get();
		}

    public:
		template<typename T>
		inline void _synchrotronIntegralsInit(const thor_scsi::core::ConfigType &conf,  gtpsa::ss_vect<T> &ps){
			if(this->computeSynchrotronIntegrals(conf)){
				auto ob = this->_getRadiationDelegate();
				if(ob){
					ob->view(*this, ps, thor_scsi::core::ObservedState::start, 0);
				}
			}
		}

		// first step of synchrotron integrals increment / local contribution
		// rename it to
		template<typename T>
		inline void _synchrotronIntegralsFinish(const thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<T> &ps){
			if(this->computeSynchrotronIntegrals(conf)){
				auto obj = this->_getRadiationDelegate();
				if(obj){
					obj->view(*this, ps, thor_scsi::core::ObservedState::end, 0);
				}
			}
		}

		// calculate the effect of radiation
		template<typename T>
		inline void _synchrotronIntegralsStep(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<T> &ps, const int step) {
			if(this->computeSynchrotronIntegrals(conf)){
				auto obj = this->_getRadiationDelegate();
				if(obj){
					obj->view(*this, ps, thor_scsi::core::ObservedState::event, step);
				}
			}
		}

		// calculate quadfringe if quadrupole and required
		template<typename T>
		void _quadFringe(thor_scsi::core::ConfigType &conf, gtpsa::ss_vect<T> &ps);

		FieldKickForthOrder<C> integ4O;
		int  Pmethod;                 ///< Integration Method.
		bool Pthick;                  ///< Thick or thin element

	};

	typedef FieldKickKnobbed<thor_scsi::core::StandardDoubleType> FieldKick;
	typedef FieldKickKnobbed<thor_scsi::core::TpsaVariantType> FieldKickTpsa;

} // Name space

#endif // _THOR_SCSI_ELEMENTS_FIELD_KICK_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
