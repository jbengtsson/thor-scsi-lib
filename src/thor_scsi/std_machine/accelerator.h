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
		/** @brief pass the given state through the machine
		 *
		 * @param conf Configuration of calculation
		 * @param ps state of calculation
		 * @param start The index of the first Element the state will pass through
		 * @param max The maximum number of elements through which the state will be passed
		 * @throws std::exception sub-classes for various errors.
		 *         If an exception is thrown then the state of S is undefined.
		 */
		Accelerator(const Config &conf);
		template <typename T>
		void _propagate(thor_scsi::core::ConfigType& conf, ss_vect<T>& ps, size_t start, int max);

		void propagate(thor_scsi::core::ConfigType&, ss_vect_tps &ps,
			       size_t start=0,
			       int max=INT_MAX);
		void propagate(thor_scsi::core::ConfigType&, ss_vect_dbl &ps,
			       size_t start=0,
			       int max=INT_MAX);
	};
  }
#endif // _THOR_SCSI_STD_MACHINE_ACCELERATOR_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
