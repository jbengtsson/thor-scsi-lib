#ifndef  _THOR_SCSI_CORE_CONFIG_H_
#define _THOR_SCSI_CORE_CONFIG_H_ 1

#include <vector>

namespace thor_scsi {
	namespace core {
		class ConfigType {
		public:
			bool
			trace,
				reverse_elem,
				stable,
				ErrFlag,
				Cavity_on,                    ///< if true, cavity turned on
				radiation,                    ///< if true, radiation turned on
				emittance,
				quad_fringe,                  ///< quadrupole hard-edge fringe field.
				H_exact,                      ///< "Small Ring" Hamiltonian.
				Cart_Bend,
				dip_edge_fudge,               ///< Dipole Edge fudge.
				pathlength,                   ///< Absolute Path Length.
				Aperture_on,
				EPU,
				mat_meth,                     ///< Matrix method.
				IBS,                          ///< Intrabeam Scattering.
				tuneflag,
				chromflag,
				codflag,
				mapflag,
				passflag,
				overflag,
				chambre;

			long int
			Cell_nLoc,                    ///< Number of Elements.
				Elem_nFam,                    ///< Number of Families.
				CODimax;                      /** closed Orbit Finder: max number of
								 iterations. */

			int
			bpm,                          ///< BPM Number.
				hcorr,                        ///< Corrector: Horizontal number,
				vcorr,                        ///< Corrector: Vertical number.
				qt,                           //< Corrector: Vertical corrector number. Todo: compare to vcorr
				gs,                           //< Girder: start marker,
				ge,                           //< Girder:  end marker.
				RingType,                     //< 1 if a ring (0 if transfer line).
				lossplane;                    /** lost in: horizontal    1
								 vertical      2
								 longitudinal  3 */
			double
			dPcommon,                     //< dp for numerical differentiation.
				dPparticle,                   //< Energy deviation.
				delta_RF,                     //< RF Acceptance.
				Omega,                        //< Synchrotron Frequency.
				U0,                           //< Energy Loss per turn [keV].
				Alphac,                       //< Linear Momentum Compaction.
				Energy,                       //< Beam Energy.
				dE,                           //< Energy Loss.
				CODeps,                       //< Closed Orbit precision.
				Qb,                           //< Bunch Charge.
				alpha_z,                      //< Long. alpha and beta.
				beta_z,
				beta0,                        //< Relativistic factors.
				gamma0;
			std::vector<double>
			TotalTune{0e0, 0e0, 0e0},     //< Transverse tunes.
				Chrom{0e0, 0e0},              //< Linear chromaticities.
				alpha_rad{0e0, 0e0, 0e0},     //< Damping coeffs.
				D_rad{0e0, 0e0, 0e0},         //< Diffusion coeffs (Floquet space).
				J{0e0, 0e0, 0e0},             //< Partition numbers.
				tau{0e0, 0e0, 0e0},           //< Damping times.
				D_IBS{0e0, 0e0, 0e0},         //< Diffusion matrix (Floquet space).
				eps{0e0, 0e0, 0e0},           //< Eigenemittances.
				epsp{0e0, 0e0, 0e0},          //< Trans. & Long. projected emittances.
				CODvect                       //< Closed orbit.
				{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
				wr                            //< Eigenvalues Re & Im part.
				{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
				wi
				{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
			std::vector< std::vector<double> >
			OneTurnMat
			{{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
			 {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}},
				Ascr
				{{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}},
				Ascrinv
				{{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}},
				Vr                            //< Eigenvectors: Real part,
				{{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}},

				Vi                            //<               Imaginary part.
				{{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
					{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}};
		};
	}
}
#endif /* _THOR_SCSI_CORE_CONFIG_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
