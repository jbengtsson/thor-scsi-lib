#ifndef _THOR_SCSI_CORE_ELEMENTS_BASIS_H_
#define _THOR_SCSI_CORE_ELEMENTS_BASIS_H_ 1

/**
   Definitions common for all elements

 */
#include <vector>
#include <string>
#include <tps/ss_vect.h>
#include <tps/tps_type.h>
#include <thor_scsi/core/cells.h>
#include <thor_scsi/core/elements_enums.h>
#include <thor_scsi/core/config.h>

namespace thor_scsi {
	namespace elements {
		//< Element virtual base class.
		class ElemType : public CellType {
		public:
			std::string
			Name;                      ///< Element name.
			bool
			Reverse;                   ///< reverse elements: rearange the elements in reveresed order
			double
			PL;                        ///< Length[m].
			PartsKind
			Pkind;                     ///<  Enumeration for magnet types.


			ElemType(const Config & config) : CellType(config)
				{};
			std::string repr_elem(void);     ///< auxilliary function providing a string of common information
			                                 ///< required for the different elements
			virtual std::string repr(void) = 0;  ///< representation similar to prt_elem but a bit more pythonic
			                                 ///< used by python interface to generate the information for
			                                 ///< :meth:`__repr__`

			void prt_elem(const std::string &);

			/**
			 * Todo:
			 *    Check if that is still missing an overloaded method?
			 *
			 * If understood coorectly one should review if a element factory is requireed.
			 */
			virtual ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse)
				{ return NULL; };

			/**
			 * Todo: implement taking stream or as ostream operator ....
			 */
			virtual void print(const std::string &) {};

			virtual void SetdS(void) {}; ///< Eucledian Group: dx, dy
			virtual void SetdT(void) {}; ///< Eucledian Group: Roll angle
			virtual void SetPB(const int n) {}; ///< Multipoles (total numbers)
			virtual double GetdT(void) { return 0e0; };
			virtual double GetPB(const int n) { return 0e0; };

			// C++ templates not supported for virtual functions.

			/**
			 * Propagater  step for phase space.
			 *
			 * Args:
			 *    ps : phase space
			 */
			virtual void pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps) {};
			virtual void pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps) {};

			template<typename T>
			void pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);

			///< Candidate to be moved towards aperture
			template<typename T>
			  bool CheckAmpl(thor_scsi::core::ConfigType &conf, const ss_vect<T> &x);

		};


		///< Index for lattice families & elements.
		class ElemFamType {
		public:
			ElemType 		*ElemF;
			int
			nKid,                      ///< No of kids.
				NoDBN;
			std::vector<int> KidList;   ///< Todo: position number in lattice ??
			std::vector<std::string>    DBNlist; ///< For control system. Todo: but what ?
		};
	}
}
#endif /*  _THOR_SCSI_CORE_ELEMENTS_BASIS_H_  */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
