#ifndef  _THOR_SCSI_CORE_CONFIG_H_
#define _THOR_SCSI_CORE_CONFIG_H_ 1

#include <vector>
#include <armadillo>

namespace thor_scsi {
	namespace core {
		class ConfigType {
		public:
			bool
			trace = false, ///< consider to remove (handled by elements now ...)
				reverse_elem = false,
				stable = false,
				ErrFlag  = false,
				Cavity_on  = false,           ///< if true, cavity turned on
				radiation  = false,           ///< if true, radiation turned on
				emittance  = false,
				quad_fringe  = false,        ///< quadrupole hard-edge fringe field.
				H_exact  = false,            ///< "Small Ring" Hamiltonian.
				Cart_Bend  = false,
				dip_edge_fudge  = true,               ///< Dipole Edge fudge.
				pathlength  = false,                   ///< Absolute Path Length.
				Aperture_on = false,                  ///< Aperture limitation used ?
				EPU  = false,
				mat_meth  = false,                     ///< Matrix method.
				IBS  = false,                          ///< Intrabeam Scattering.
				tuneflag  = false,
				chromflag  = false,
				codflag  = false,
				mapflag  = false,
				passflag  = false,
				overflag = false,
				chambre  = false;

			long int
			Cell_nLoc = -1,                    ///< Number of Elements. I thnk not required any more
				Elem_nFam = -1,            ///< Number of Families.
				CODimax =  -1;             /** closed Orbit Finder: max number of
							       iterations. */

			int
			bpm = -1,                          ///< BPM Number.
				hcorr = -1,                        ///< Corrector: Horizontal number,
				vcorr = -1,                        ///< Corrector: Vertical number.
				qt = -1,                           //< Corrector: Vertical corrector number. Todo: compare to vcorr
				gs = -1,                           //< Girder: start marker,
				ge = -1,                           //< Girder:  end marker.
				RingType = 1,                     //< 1 if a ring (0 if transfer line).
				lossplane = 0;                    /** lost in: horizontal    1
								 vertical      2
								 longitudinal  3 */
			double
			dPcommon = 0e0,                     //< dp for numerical differentiation.
				dPparticle  = 0e0,                   //< Energy deviation.
				delta_RF = 0e0,                     //< RF Acceptance.
				Omega = 0e0,                        //< Synchrotron Frequency.
				U0 = 0e0,                           //< Energy Loss per turn [keV].
				Alphac = 0e0,                       //< Linear Momentum Compaction.
				Energy = 1.7e9,                     //< Beam Energy in eV.
				dE = 0e0,                           //< Energy Loss.
				CODeps = 1e-6,                       //< Closed Orbit precision.
				Qb= 0e0,                           //< Bunch Charge.
				alpha_z = 0e0,                      //< Long. alpha and beta.
				beta_z = 0e0,
				beta0 = 0e0,                        //< Relativistic factors.
				gamma0 = 0e0;
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
