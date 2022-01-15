#ifndef _THOR_SCSI_COMPAT_FAMILIES_H_
#define _THOR_SCSI_COMPAT_FAMILIES_H_ 1
#include <thor_scsi/exceptions.h>

#warning "Functions in these headers require to be implemented"

namespace thor_scsi {
	namespace compat {
		/**
		 * Todo:
		 *     Review definition of families!
		 *     Appropriate functionality should be implemented in
		 *     thor_scsi::elements::ElemFamType
		 */
		inline void set_bn_design_fam(const int Fnum,
					      const int n, const double bn, const double an){
			throw thor_scsi::NotImplemented();
		}

		inline void set_dbnL_design_fam(const int Fnum,
						const int n, const double dbnL, const double danL){
			throw thor_scsi::NotImplemented();
		}

		inline void get_bn_design_elem(const int Fnum, const int Knum,
					       const int n, double &bn, double &an){
			throw thor_scsi::NotImplemented();
		}

		inline void set_dbnL_design_elem(const int Fnum, const int Knum,
						 const int n, const double dbnL, const double danL){
			throw thor_scsi::NotImplemented();
		}
	} /* namespace compat */
} /* namespace thor_scsi */


#endif /* _THOR_SCSI_COMPAT_FAMILIES_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
