#ifndef _THOR_SCSI_CORE_ELEMENTS_ENUMS_H_
#define _THOR_SCSI_CORE_ELEMENTS_ENUMS_H_ 1

namespace thor_scsi {
	namespace elements {
		//typedef std::vector<double> MpoleArray;


		/**
		 * Todo:
		 *    Check if obsoleted by class inheritance
		 */
		enum PlaneKind
		{
			Horizontal = 1,
			Vertical   = 2
		};

		enum IntMethKind
		{
			Meth_Linear = 0,
			Meth_First  = 1,
			Meth_Second = 2,
			Meth_Fourth = 4
		};
	}
}
#endif /* _THOR_SCSI_CORE_ELEMENTS_ENUMS_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
