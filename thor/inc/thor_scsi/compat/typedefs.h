#ifndef _THOR_SCSI_COMPAT_TYPEDEFS_H_
#define _THOR_SCSI_COMPAT_TYPEDEFS_H_ 1
#include <armadillo>
#include <array>
#include <vector>
#include <tps/ss_vect.h>
#include <tps/tpsa_lin.h>
#include <thor_scsi/compat/constants.h>

namespace thor_scsi {
	namespace compat {
		typedef ss_vect<double>  psVector;
		// typedef psVector         Matrix[ss_dim];
		typedef arma::mat        Matrix;
		typedef std::array<double, 2> Vector2;

		/// Warning: requires to be implemented
		void _prtmat(const int n, const Matrix &A);

		inline void prtmat(const int n, const std::vector< std::vector<double> > &stlmat){
			_prtmat(n, stlmattomat(stlmat));
		}

		void TpMat(const int n, Matrix &A);
		double DetMat(const int n, const Matrix &A_);
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
