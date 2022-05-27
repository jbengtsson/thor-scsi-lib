#ifndef _THOR_SCSI_CORE_ELEMENTS_ENUMS_H_
#define _THOR_SCSI_CORE_ELEMENTS_ENUMS_H_ 1

namespace thor_scsi {
	namespace elements {
		//typedef std::vector<double> MpoleArray;


		/**
		 * Todo:
		 *    Check if obsoleted by class inheritance
		 */
#if 0
		enum PartsKind
		{
			drift      = 0,
			Wigl       = 1,
			Mpole      = 2,
			Cavity     = 3,
			marker     = 4,
			undef      = 5,
			Insertion  = 6,
			FieldMap   = 7,
			Spreader   = 8,
			Recombiner = 9,
			Solenoid   = 10,
			Map        = 11 };


		enum MpoleKind
		{
			All   = 0,
			Dip   = 1,
			Quad  = 2,
			Sext  = 3,
			Oct   = 4,
			Dec   = 5,
			Dodec = 6
		};


#endif
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

#if 0
		const int n_harm_max   = 10,
			Spreader_max = 10;

		const double
			max_ampl = 10.0; // [m]
#endif

	}
}
#endif /* _THOR_SCSI_CORE_ELEMENTS_ENUMS_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
