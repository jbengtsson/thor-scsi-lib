#ifndef _THOR_SCSI_CORE_ELEMENTS_BASIS_H_
#define _THOR_SCSI_CORE_ELEMENTS_BASIS_H_ 1

#include <vector>
#include <string>
#include <tps/ss_vect.h>
#include <tps/tps_type.h>
#include <thor_scsi/core/cells.h>
#include <thor_scsi/core/elements_enums.h>
#include <thor_scsi/core/config.h>

namespace thor_scsi {
	namespace elements {
		// Element virtual base class.
		class ElemType : public CellType {
		public:
			std::string
			Name;                      // Element name.
			bool
			Reverse;                   // Reverse element.
			double
			PL;                        // Length[m].
			PartsKind
			Pkind;                     // Enumeration for magnet types.


			// representation similar to prt_elem but a bit more pythonic
			std::string repr_elem(void);
			virtual std::string repr() = 0;

			void prt_elem(const std::string &);

			/*
			 * Todo:
			 *    Check if that is still missing an overloaded method?
			 */
			virtual ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse)
				{ return NULL; };
			virtual void print(const std::string &) {};

			virtual void SetdS(void) {};
			virtual void SetdT(void) {};
			virtual void SetPB(const int n) {};
			virtual double GetdT(void) { return 0e0; };
			virtual double GetPB(const int n) { return 0e0; };

			// C++ templates not supported for virtual functions.
			virtual void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps) {};
			virtual void Elem_Pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps) {};

			template<typename T>
			  bool CheckAmpl(thor_scsi::core::ConfigType &conf, const ss_vect<T> &x);
			template<typename T>
			void Cell_Pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
		};


		// Index for lattice families & elements.
		class ElemFamType {
		public:
			ElemType 		*ElemF;
			int
			nKid,                      // No of kids.
				NoDBN;
			std::vector<int> KidList;
			std::vector<std::string>    DBNlist; // For control system.
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
