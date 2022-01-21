#ifndef _THOR_SCSI_COMPAT_GLOBALS_H_
#define _THOR_SCSI_COMPAT_GLOBALS_H_ 1
/**

  Objects defined as global objects to be accessible for former tracy code
 */
#include <cassert>
#include <thor_scsi/core/config.h>
#include <thor_scsi/core/lattice.h>
#include <thor_scsi/process/lsoc.h>
#include <thor_scsi/process/t2ring.h>
#include <thor_scsi/obsolete/nsls-ii_lib.h>
#include <thor_scsi/tweak/errors.h>
#include <thor_scsi/tweak/param.h>
#include <thor_scsi/exceptions.h>

namespace thor_scsi {
	namespace compat {
		/**

		   Global lattice definition. required for implementation of compatabiltiy layer
		 */
		thor_scsi::core::LatticeType lat = thor_scsi::core::LatticeType(); ///< global definition of lattice
		thor_scsi::core::ConfigType  &globval = lat.conf;                  ///< global definition of config, was named globval before
		std::vector<thor_scsi::elements::ElemType*> &Cell = lat.elems;     ///< global definition of elements

		///< Cell_Pass(0, globval.Cell_nLoc, M, lastpos) used to be called like that before
		template<typename T>
		inline void Cell_Pass(const long i0, const long i1, ss_vect<T> &ps, long &lastpos){
			lat.Cell_Pass(i0, i1, ps, lastpos);
		}

		inline bool getcod(double dp, long &lastpos){
			return lat.getcod(dp, lastpos);
		}

		inline long int ElemIndex(const std::string &name){
			return lat.ElemIndex(name);
		}

		inline long int Elem_GetPos(const int Fnum, const int Knum){
			return lat.Elem_GetPos(Fnum, Knum);
		}
		inline int GetnKid(const int Fnum){
			return lat.GetnKid(Fnum);
		}

		inline void get_eps_x(double &eps_x, double &sigma_delta, double &U_0,
				      std::vector<double> &J, std::vector<double> &tau,
				      std::vector<double> &I, const bool prt){
			return lat.get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, prt);
		}

		inline void Ring_GetTwiss(bool chroma, double dp){
			return lat.Ring_GetTwiss(chroma, dp);
		}

		inline double get_code_compat(const thor_scsi::elements::ElemType *Cell){
			return get_code(globval, *Cell);
		}

		/**
		 * Todo: reimplement
		 *
		 * Implemented in the nsls_lib:
		 *
		 * Functionality should be in lattice or family class:
		 */
		inline void no_sxt(void){
			throw thor_scsi::NotImplemented();
			// thor_scsi::no_sxt(lat);
		}

		inline void printglob(void){
			lat.print("");
		}

		inline void misalign_rms_type(const int type, const double dx_rms, const double dy_rms,
					      const double dr_rms, const bool new_rnd){
			thor_scsi::misalign_rms_type(lat, type, dx_rms, dy_rms, dr_rms, new_rnd);
		}

		/**
		 * Todo: reimplement
		 *
		 * Implemented in the nsls_lib?
		 *
		 * Functionality should be in lattice or family class:
		 */
		inline bool orb_corr(const int n_orbit){
			throw thor_scsi::NotImplemented();
			// return thor_scsi::orb_corr(lat, n_orbit);
		}

		/**
		 * Todo:
		 *      should one call Lat_Init?
		 */
		inline void Read_Lattice(const std::string &filename, bool verbose=true){
			lat.Lat_Read(filename, verbose);
			lat.Lat_Init();
		}

		inline void rdmfile(const std::string &filename){
			lat.rdmfile(filename);
		}
		inline void prtmfile(const std::string &filename){
			lat.prtmfile(filename);
		}

		inline void set_map(const int Fnum, const double dnu[]){
			thor_scsi::set_map(lat, Fnum, dnu);
		}
		inline void GDiag(int n, double C, arma::mat &A, arma::mat &Ainv, arma::mat &R,
			   arma::mat &M, double &Omega, double &alphac){
			lat.GDiag(n, C, A, Ainv, R, M, Omega, alphac);
		}
		inline void gcmat(const int bpm, const int corr, const int plane){
			thor_scsi::gcmat(lat, bpm, corr, plane);
		}
		inline void get_alphac2(void){
			thor_scsi::get_alphac2(lat);
		}
		inline void GetEmittance(const int Fnum, const bool prt){
			lat.GetEmittance(Fnum, prt);
		}
		inline double get_dynap(const double delta, const int n_aper,
					 const int n_track, const bool cod ){
			return thor_scsi::get_dynap(lat, delta, n_aper, n_track, cod);
		}

		/*
		inline void ttwiss(const std::vector<double> &alpha, const std::vector<double> &beta,
			    const std::vector<double> &eta, const std::vector<double> &etap,
			    const double dp){
			lat.ttwiss(alpha, beta, eta, etap, dp);
		}
		*/
		inline void ttwiss(const Vector2 &alpha, const Vector2 &beta,
			    const Vector2 &eta, const Vector2 &etap,
			    const double dp){
			const std::vector<double>
				av(alpha.begin(), alpha.end()),
				bv(beta.begin(), beta.end()),
				ev(eta.begin(), eta.end()),
				epv(etap.begin(), etap.end());
			lat.ttwiss(av, bv, ev, epv, dp);
		}


		inline void prt_chrom_lat(void){
			lat.prt_chrom_lat();
		}
		inline void prt_lat(std::string fname, const bool all){
			lat.prt_lat2();
		}
		inline void prt_lat(std::string fname, const bool all, const int n){
			lat.prt_lat4();
		}


	}; /* compat */
}; /* thor_scsi */
#endif /* _THOR_SCSI_COMPAT_GLOBALS_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
