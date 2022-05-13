#ifndef _THOR_SCSI_STD_MACHINE_ACCELERATOR_
#define _THOR_SCSI_STD_MACHINE_ACCELERATOR_

#include <thor_scsi/core/machine.h>
#include <tps/tps_type.h>
#include <tps/ss_vect.h>

namespace thor_scsi {
	typedef ss_vect<tps> ss_vect_tps;
	typedef ss_vect<double> ss_vect_dbl;

	class Accelerator : public thor_scsi::core::Machine {
	public:
		Accelerator(const Config &conf);
		/** @brief pass the given state through the machine
		 *
		 * @param conf Configuration of calculation
		 * @param ps state of calculation
		 * @param start The index of the first Element the state will pass through
		 * @param max The maximum number of elements through which the state will be passed
		 * @returns last element passed (check config type for lost plane)
		 * @throws std::exception sub-classes for various errors.
		 *         If an exception is thrown then the state of S is undefined.
		 *
		 * @todo proper interface design!
		 */
		template <typename T>
		int _propagate(thor_scsi::core::ConfigType& conf, ss_vect<T>& ps, size_t start, int max);

		int propagate(thor_scsi::core::ConfigType&, ss_vect_tps &ps,
			       size_t start=0,
			       int max_elements=std::numeric_limits<int>::max());
		int propagate(thor_scsi::core::ConfigType&, ss_vect_dbl &ps,
			       size_t start=0,
			       int max_elements=std::numeric_limits<int>::max());
	};
  }
#endif // _THOR_SCSI_STD_MACHINE_ACCELERATOR_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
