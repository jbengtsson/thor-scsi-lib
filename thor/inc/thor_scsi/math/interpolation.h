#ifndef _THOR_SCSI_MATH_INTERPOLATION_H_
#define _THOR_SCSI_MATH_INTERPOLATION_H_ 1

#include <thor_scsi/core/cells.h>

namespace thor_scsi {
  namespace math {
    template<typename T>
      void LinearInterpolation2(T &X, T &Z, T &TX, T &TZ, T &B2,
				thor_scsi::elements::CellType &Cell,
				bool &out, int order);
  }


}


#endif /* _THOR_SCSI_MATH_INTERPOLATION_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
