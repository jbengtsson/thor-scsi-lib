#ifndef _THOR_SCSI_MATH_INTERPOLATION_H_
#define _THOR_SCSI_MATH_INTERPOLATION_H_ 1

#include <thor_scsi/core/elements_basis.h>

namespace thor_scsi {
	namespace math {
		template<typename T>
		void LinearInterpolation2(T &X, T &Z, T &TX, T &TZ, T &B2,
					  thor_scsi::elements::ElemType *Cell,
					  bool &out, int order);

		template<typename T>
		void SplineInterpolation2(T &X, T &Z, T &thetax, T &thetaz,
					  thor_scsi::elements::ElemType *Cell, bool &out);

		void splie2(double x1a[], double x2a[], double **ya,
			    int m, int n, double **y2a);

	}
}

#endif /* _THOR_SCSI_MATH_INTERPOLATION_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
