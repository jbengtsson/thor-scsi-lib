#ifndef _THOR_SCSI_MATH_COMB_H_
#define _THOR_SCSI_MATH_COMB_H_
namespace thor_scsi::core {
	/**
	 * calculate binomial coefficient
	 *
	 * Todo: check if in standard library or in GSL
	 */
	unsigned long binom (unsigned long const n, unsigned long const k);
}
#endif /* THOR_SCSI_MATH_COMB_H */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
