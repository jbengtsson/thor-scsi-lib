#ifndef _THOR_SCSI_STD_MACHINE_ACCELERATOR_
#define _THOR_SCSI_STD_MACHINE_ACCELERATOR_

#include <thor_scsi/core/machine.h>
#include <tps/tps_type.h>
// #include <tps/ss_vect.h>

namespace thor_scsi {

	typedef gtpsa::ss_vect<tps>         ss_vect_tps;
	typedef gtpsa::ss_vect<gtpsa::tpsa> ss_vect_tpsa;
	typedef gtpsa::ss_vect<double>      ss_vect_dbl;

	/**
	 * @brief propagation result of accelerator
	 *
	 * Motivated by how one should report lost elements ...
	 * ConfigType is not the place to store it ...
	 */
	class PropagationResult {
		int last_element;
		int loss_plane;
		bool success;
	};

	template<class C>
	class AcceleratorKnobbable : public thor_scsi::core::Machine {
	public:

		/**
		 * @brief initalise accelerator with a configuration file
		 *
		 * @param conf configuration for the accelerator (aka parsed lattise file)
		 * @param add_marker_at_start add a marker at the start of the lattice
		 *
		 * @warning consider if the marker is not better added manually to the lattice
		 *          file
		 */
		AcceleratorKnobbable(const Config &conf, bool add_marker_at_start=false);

		/**
		 * @brief inititsialse accelerator with a list of elements
		 *
		 * @param elements list of elements
		 * @param add_marker_at_start add a marker at the start of the lattice
		 *
		 * @warning consider if the marker is not better added manually to the lattice
		 *          file
		 * @param conf
		 */
		/*
		*/
		AcceleratorKnobbable(std::vector<std::shared_ptr<thor_scsi::core::CellVoid>> elements,
				     bool add_marker_at_start=false);

		AcceleratorKnobbable(const std::vector<std::shared_ptr<thor_scsi::core::ElemTypeKnobbed>> elements,
				     bool add_marker_at_start=false);
		/** @brief pass the given state through the machine
		 *
		 * @param conf Configuration of calculation
		 * @param ps state of calculation
		 * @param start The index of the first Element the state will pass through
		 * @param max The maximum number of elements through which the state will be passed
		 * @param tracy compatible indexing: start to refer to first element with 1 instead of zero
		 * @param add_marker_at_start add a marker at the start of the lattice
		 *
		 * @returns last element passed (check config type for lost plane)
		 *
		 * @throws std::exception sub-classes for various errors.
		 *         If an exception is thrown then the state of S is undefined.
		 *
		 * @warning  tracy compatible indexing will be removed soon as it is not consistent with global indexing
		 *           consider if the marker is not better added manually to the lattice
		 *           file
		 *
		 * @todo proper interface design!
		 */
		template <typename T>
		int _propagate(thor_scsi::core::ConfigType& conf, gtpsa::ss_vect<T>& ps, size_t start, int max, size_t n_turns, bool tracy_compatible_indexing = false);

	    /*
		int propagate(thor_scsi::core::ConfigType&, ss_vect_tps  &ps,
			      size_t start=0,
			      int max_elements=std::numeric_limits<int>::max(), size_t n_turns=1, bool tracy_compatible_indexing = false);
	    */
		int propagate(thor_scsi::core::ConfigType&, ss_vect_tpsa &ps,
			       size_t start=0,
			      int max_elements=std::numeric_limits<int>::max(), size_t n_turns=1, bool tracy_compatible_indexing = false);
		int propagate(thor_scsi::core::ConfigType&, ss_vect_dbl  &ps,
			       size_t start=0,
			      int max_elements=std::numeric_limits<int>::max(), size_t n_turns=1, bool tracy_compatible_indexing = false);
	private:
		/**
		 * @brief add a marker at the beginning of the lattice if the lattice does not start with one
		 *
		 * @warning rather than using this feature consider adding it manually to the lattice file.
		 *          this here is plagued to add more and more in evolved processing lines
		 */
		void addMarkerAtStartIfRequired(void);
		void addMarkerAtStart(void);
		template <typename T>
		int _propagate(thor_scsi::core::ConfigType& conf, ss_vect<T>& ps, size_t start, int max, size_t n_turns, bool tracy_compatible_indexing);
	};

    typedef class AcceleratorKnobbable<thor_scsi::core::StandardDoubleType> Accelerator;
    typedef class AcceleratorKnobbable<thor_scsi::core::TpsaVariantType> AcceleratorTpsa;
  }
#endif // _THOR_SCSI_STD_MACHINE_ACCELERATOR_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
