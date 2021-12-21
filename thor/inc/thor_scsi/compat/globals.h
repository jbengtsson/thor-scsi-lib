#ifndef _THOR_SCSI_COMPAT_GLOBALS_H_
#define _THOR_SCSI_COMPAT_GLOBALS_H_ 1
/**

  Objects defined as global objects to be accessible for former tracy code
 */
#include <thor_scsi/core/config.h>
#include <thor_scsi/core/lattice.h>
#include <thor_scsi/process/t2ring.h>

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


	}; /* compat */
}; /* thor_scsi */
#endif /* _THOR_SCSI_COMPAT_GLOBALS_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
