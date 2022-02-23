#ifndef _THOR_SCSI_CORE_ELEMENTS_H_
#define _THOR_SCSI_CORE_ELEMENTS_H_ 1

#include <thor_scsi/core/config.h>
#include <thor_scsi/core/defines.h>
#include <thor_scsi/core/elements_basis.h>




namespace thor_scsi {
	namespace elements {

#if 0

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
		/**
		 * Calcualte multipole kick. The kick is given by
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
		 * The magnetic field expansion used here is
		 *  @f[
		 *     \mathbf{B(z)} = B_y + i B_x =
		 *	\sum\limits_{n=1}^N \mathbf{C}_n
		 *	\left(\frac{\mathbf{z}}{R_\mathrm{ref}}\right)^{(n-1)}
		 *  @f]
		 *
		 *  with \f$ \mathbf{z} = x + i y \f$ and \f$ \mathbf{C}_n = B_n + i A_n \f$
		 *  and \f$R_\mathrm{ref}\f$ the reference radius.
		 *  \f$N\f$ corresponds to the maximum harmonic
		 *
		 * Furthermore note that \f$R_\mathrm{ref}\f$ = 1
		 *
		 */
		typedef std::vector<double> MpoleArray;
		class MpoleType : public ElemType {
		public:
			int
			Pmethod,                   ///< Integration Method.
				PN,                        ///< Number of integration steps.
				Porder,                    ///< The highest order in PB.
				n_design;                  ///< multipole order (design).
			double
			PdTpar,                    ///< Roll angle [deg]: design
				PdTsys,                    //<                    systematic
				PdTrms,                    ///<                    rms
				PdTrnd,                    ///<                    random number.
				PTx1,                      ///< Bend angle [deg]:  hor. entrance angle
				PTx2,                      ///<  Bend angle [deg]: hor. exit angle.
				Pgap,                      ///< Total magnet gap [m].
				Pirho,                     ///< 1/rho [1/m].
				Pc0,                       ///< Corrections for roll error of bend.
				Pc1,
				Ps1;
			std::vector<double>
			PdSsys{0e0, 0e0},                  ///< Displacement errors [m]: systematic
				PdSrms{0e0, 0e0},          ///<                          rms
				PdSrnd{0e0, 0e0};          ///<                          random number.
			/* todo: review to implement is a complex arrays */
			MpoleArray
			PBpar,                             ///< Multipole strengths: design
				PBsys,                     ///< Multipole strengths: systematic
				PBrms,                     ///< Multipole strengths: random variation \f$ <\sigma> \f$
				PBrnd,                     ///< random number used for computing PBrms.
				PB;                        ///< Multipole strengths: total used by integrator
			pthicktype
			Pthick;
			std::vector< std::vector<double> >
			M_elem                     ///< Transport matrix & orbit.
			{{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}};


			friend MpoleType* Mpole_Alloc(void);
			ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
			void print(const std::string &);
			std::string repr(void);
			///< additional information of mpole for mpole
			///< more than just solely type name
			///< used by meth:`repr`
			std::string repr_add(void);

			/**
			 *  sets up the  calculation of the Euclidian group: translation
			 */
			void SetdS(void);
			void SetdT(void);
			/**
			 * triggers calculation of member PB from members PBrms, PBrnd, PBsys
			 */
			void SetPB(const int n);
			double GetdT(void);
			double GetPB(const int n);

			template<typename T>
			void Mpole_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps)
				{ Mpole_Pass(conf, ps); };
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps)
				{ Mpole_Pass(conf, ps); };

			inline double GetKpar(int order){
				return PBpar[order+HOMmax];
			}

			inline void SetKpar(int order, double val){
				PBpar[order+HOMmax] = val;
			}

		};
		/**
		 * RF cavity
		 */
		class CavityType : public ElemType {
		public:
			bool
			entry_focus,               ///< Edge focusing at entry.
				exit_focus;                ///< Edge focusing at exit.
			int
			PN,                        ///< Number of integration steps.
				Ph;                        ///< Harmonic number.
			double
			Pvolt,                     ///< Vrf [V].
				Pfreq,                     ///< Vrf [Hz].
				phi;                       ///< RF phase.

			friend CavityType* Cavity_Alloc(void);
			ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
			void print(const std::string &);
			std::string repr(void);

			void SetdS(void) {};
			void SetdT(void) {};
			void SetPB(const int n) {};
			double GetdT(void) { return 0e0; };
			double GetPB(const int n) { return 0e0; };

			template<typename T>
			void Cavity_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps)
				{ Cavity_Pass(conf, ps); };
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps)
				{ Cavity_Pass(conf, ps); };

			inline void setVoltage(double val){this->Pvolt = val;}
			inline double getVoltage(void){return this->Pvolt;}

			inline void setFrequency(double val){this->Pfreq = val;}
			inline double getFrequency(void){return this->Pfreq;}

		};

		/**
		 * A simple beam position monitor. can only sit on the ideal position
		 *
		 */
		class MarkerType : public ElemType {
		public:
			friend MarkerType* Marker_Alloc(void);
			ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
			void print(const std::string &);
			std::string repr(void);

			void SetdS(void) {};
			void SetdT(void) {};
			void SetPB(const int n) {};
			double GetdT(void) { return 0e0; };
			double GetPB(const int n) { return 0e0; };

			template<typename T>
			void Marker_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps)
				{ Marker_Pass(conf, ps); };
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps)
				{ Marker_Pass(conf, ps); };
		};

		/**
		 *
		 * Todo:
		 *    Check if arrays should be translated to std:vector or arma::vec
		 */
		class WigglerType : public ElemType {
		public:
			int
			Pmethod,                   ///< Integration Method.
				PN,                        ///< number of integration steps.
				n_harm,                    ///< No of harmonics.
				harm[n_harm_max],          ///< Harmonic number.
				Porder;                    ///< The highest order in PB.
			double
			PdTpar,                    ///< Roll angle [deg]: design
				PdTsys,                    ///<                   systematic
				PdTrms,                    ///<                   RMS
				PdTrnd,                    ///<                   random number.
				Lambda,                    ///< lambda.
				BoBrhoV[n_harm_max],       ///< B/Brho: ver.
				BoBrhoH[n_harm_max],       ///<         hor.
				kxV[n_harm_max],           ///< kx.
				kxH[n_harm_max],           ///< kx.
				phi[n_harm_max];           ///< phi.
			std::vector<double>
			PdSsys{0e0, 0e0},          ///< Displacement error [m]: systematic
				PdSrms{0e0, 0e0},          ///<                         RMS
				PdSrnd{0e0, 0e0};          ///<                         random number.
			MpoleArray
			PBW;

			friend WigglerType* Wiggler_Alloc(void);
			ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
			void print(const std::string &);
			std::string repr(void);

			void SetdS(void);
			void SetdT(void);
			void SetPB(const int n) {};
			double GetdT(void);
			double GetPB(const int n);

			template<typename T>
			void Wiggler_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps)
				{ Wiggler_Pass(conf, ps); };
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps)
				{ Wiggler_Pass(conf, ps); };
		};

		class InsertionType : public ElemType {
		public:
			char
			fname1[100],               ///< Filename for insertion description: 1st order.
				fname2[100];               ///< Filename for insertion description: 2nd order.
			bool
			linear,                    ///< if true linear interpolation else spline.
				firstorder,                ///< true if first order kick map loaded.
				secondorder,               ///< true if second order kick map loaded.
				long_comp;                 ///< flag for longitudinal comp.
			int
			Pmethod,                   ///< Integration Method.
				PN,                        ///< number of integration steps.
				nx,                        ///< Horizontal point number.
				nz,                        ///< Vertical point number.
				Porder;                    ///< The highest order in PB.
			double
			scaling,                   ///< static scaling factor as in BETA ESRF.
				phi,                       ///< Bend angle.
				tabx[IDXMAX],              ///< spacing in H-plane.
				tabz[IDZMAX],              ///< spacing in V-plane.
				thetax[IDZMAX][IDXMAX],
				thetax1[IDZMAX][IDXMAX],   ///< 1 for first order.
				thetaz[IDZMAX][IDXMAX],
				thetaz1[IDZMAX][IDXMAX],
				B2[IDZMAX][IDXMAX],        ///< B^2_perp.
				**tx,
				**tz,
				**f2x,
				**f2z,
				**tx1,
				**tz1,
				**f2x1,
				**f2z1,                    ///< a voir.
				*tab1,
				*tab2,                     ///< tab of x and z meshes from Radia code.
                                                           ///< Roll angle.
				PdTpar,                    ///< design [deg].
				PdTsys,                    ///< systematic [deg].
				PdTrms,                    ///< rms [deg].
				PdTrnd;                    ///< random number.
///< Strength
//  double Plperiod;           // Length Period [m].
//  int Pnperiod;              // Number of periods.
//  double PBoBrho;            // B/Brho.
//  double PKx;                // kx.
//  mpolArray PBW;
			std::vector<double>
			PdSsys{0e0, 0e0},          ///< Displacement error [m]: systematic
				PdSrms{0e0, 0e0},          ///<                         rms
				PdSrnd{0e0, 0e0};          ///<                         random number.

			friend InsertionType* Insertion_Alloc(void);
			ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
			void print(const std::string &);
			std::string repr(void);

			void SetdS(void) {};
			void SetdT(void) {};
			void SetPB(const int n) {};
			double GetdT(void) { return 0e0; };
			double GetPB(const int n) { return 0e0; };

			template<typename T>
			void Insertion_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps)
				{ Insertion_Pass(conf, ps); };
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps)
				{ Insertion_Pass(conf, ps); };
		};

		class FieldMapType : public ElemType {
		public:
			int
			n_step,                    ///< number of integration steps.
				n[3],                      ///< no of steps.
				cut;                       ///< cut in z direction.
			double
			scl,
				phi,
				x0,
				Lr,
				Lm,
				Ld,
				L1,
				dx[3],
				*x[3],                     ///< [dx, dy, dz], [x, y, z].
				***BoBrho[3],
				***BoBrho2[3],             ///< [B_x, B_y, B_z].
				***AoBrho[2],
				***AoBrho2[2];             ///< [Ax(x, y, z), Ay(x, y, z)], spline info.

			friend FieldMapType* FieldMap_Alloc(void);
			ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
			void print(const std::string &);
			std::string repr(void);

			void SetdS(void) {};
			void SetdT(void) {};
			void SetPB(const int n) {};
			double GetdT(void) { return 0e0; };
			double GetPB(const int n) { return 0e0; };

			template<typename T>
			void FieldMap_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps)
				{ FieldMap_Pass(conf, ps); };
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps)
				{ FieldMap_Pass(conf, ps); };
		};

		class SpreaderType : public ElemType {
		public:
			double
			E_max[Spreader_max];       ///< energy levels in increasing order.
			CellType
			*Cell_ptrs[Spreader_max];

			friend SpreaderType* Spreader_Alloc(void);
			ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
			void print(const std::string &);
			std::string repr(void);

			void SetdS(void) {};
			void SetdT(void) {};
			void SetPB(const int n) {};
			double GetdT(void) { return 0e0; };
			double GetPB(const int n) { return 0e0; };

			template<typename T>
			void Spreader_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps)
				{ Spreader_Pass(conf, ps); };
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps)
				{ Spreader_Pass(conf, ps); };
		};

		class RecombinerType : public ElemType {
		public:
			double
			E_min,
				E_max;

			friend RecombinerType* Recombiner_Alloc(void);
			ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
			void print(const std::string &);
			std::string repr(void);

			void SetdS(void) {};
			void SetdT(void) {};
			void SetPB(const int n) {};
			double GetdT(void) { return 0e0; };
			double GetPB(const int n) { return 0e0; };

			template<typename T>
			void Recombiner_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps)
				{ Recombiner_Pass(conf, ps); };
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps)
				{ Recombiner_Pass(conf, ps); };
		};

		class SolenoidType : public ElemType {
		public:
			int
			N;                         ///< Number of integration steps.
			double
			BoBrho,                    ///< normalized field strength.
			// Roll angle.
				dTpar,                     ///< design [deg].
				dTsys,                     ///< systematic [deg].
				dTrms,                     ///< rms [deg].
				dTrnd;                     ///< random number.
			std::vector<double>
			PdSsys{0e0, 0e0},          ///< Displacement errors [m]: systematic
				PdSrms{0e0, 0e0},          ///<                          rms
				PdSrnd{0e0, 0e0};          ///<                          random number.

			friend SolenoidType* Solenoid_Alloc(void);
			ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
			void print(const std::string &);
			std::string repr(void);

			void SetdS(void) {};
			void SetdT(void) {};
			void SetPB(const int n) {};
			double GetdT(void) { return 0e0; };
			double GetPB(const int n) { return 0e0; };

			template<typename T>
			void Solenoid_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps)
				{ Solenoid_Pass(conf, ps); };
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps)
				{ Solenoid_Pass(conf, ps); };
		};

		class MapType : public ElemType {
		public:
			double
			eta_x,
				etap_x;
			std::vector<double>
			dnu{0e0, 0e0},
				alpha{0e0, 0e0},
				beta{0e0, 0e0};
			ss_vect<tps>
			M;

			friend MapType* Map_Alloc(void);
			ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse);
			void print(const std::string &);
			std::string repr(void);

			void SetdS(void) {};
			void SetdT(void) {};
			void SetPB(const int n) {};
			double GetdT(void) { return 0e0; };
			double GetPB(const int n) { return 0e0; };

			template<typename T>
			void Map_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps)
				{ Map_Pass(conf, ps); };
			void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps)
				{ Map_Pass(conf, ps); };
		};
		DriftType* Drift_Alloc(void);           ///< Todo: Some factory implementation required?
		MarkerType* Marker_Alloc(void);		///< Todo: Some factory implementation required?
		MpoleType* Mpole_Alloc(void);		///< Todo: Some factory implementation required?
		CavityType* Cavity_Alloc(void);		///< Todo: Some factory implementation required?
		WigglerType* Wiggler_Alloc(void);	///< Todo: Some factory implementation required?
		InsertionType* Insertion_Alloc(void);	///< Todo: Some factory implementation required?
		FieldMapType* FieldMap_Alloc(void);	///< Todo: Some factory implementation required?
		SpreaderType* Spreader_Alloc(void);	///< Todo: Some factory implementation required?
		RecombinerType* Recombiner_Alloc(void);	///< Todo: Some factory implementation required?
		SolenoidType* Solenoid_Alloc(void);	///< Todo: Some factory implementation required?
		MapType* Map_Alloc(void);               ///< Todo: Some factory implementation required?

#endif

		std::vector< std::vector<double> > get_transp_mat(ElemType *elem, const double delta);

	}
}
#endif /* _THOR_SCSI_CORE_ELEMENTS_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
