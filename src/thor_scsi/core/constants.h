#ifndef _THOR_SCSI_CORE_CONSTANTS_H_
#define  _THOR_SCSI_CORE_CONSTANTS_H_ 1
namespace thor_scsi::core {
	enum phase_space_index {
		x_ = 0,
		px_ = 1,
		y_ = 2,
		py_ = 3,
		delta_ = 4,
		ct_ = 5
	};

	enum spatial_index {
		X_ = 0,
		Y_ = 1,
		Z_ = 2
	};
} // namespace thor_scsi::core

#endif /* _THOR_SCSI_CORE_CONSTANTS_H_ */
