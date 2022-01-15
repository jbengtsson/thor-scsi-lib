#ifndef _THOR_SCSI_COMPAT_TYPEDEFS_H_
#define _THOR_SCSI_COMPAT_TYPEDEFS_H_ 1
#include <armadillo>
#include <array>
#include <vector>
#include <tps/ss_vect.h>
#include <tps/tpsa_lin.h>
#include <thor_scsi/compat/constants.h>
#include <thor_scsi/exceptions.h>

namespace thor_scsi {
	namespace compat {
		typedef ss_vect<double>  psVector;
		// typedef psVector         Matrix[ss_dim];
		typedef arma::mat        Matrix;
		typedef std::array<double, 2> Vector2;

		/// Warning: requires to be implemented
		inline void _prtmat(const int n, const Matrix &A){
			throw thor_scsi::NotImplemented();
		}

		inline void prtmat(const int n, const std::vector< std::vector<double> > &stlmat){
			_prtmat(n, stlmattomat(stlmat));
		}
		inline void TpMat(const int n, Matrix &A){
			throw thor_scsi::NotImplemented();
		}

		inline double DetMat(const int n, const Matrix &A_){
			throw thor_scsi::NotImplemented();
		}

		inline double DetMat(const int n, const std::vector< std::vector<double> > &stlmat){
			return DetMat(n, stlmattomat(stlmat));
		}
	};
};
#endif /* _THOR_SCSI_COMPAT_TYPEDEFS_H_  */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
