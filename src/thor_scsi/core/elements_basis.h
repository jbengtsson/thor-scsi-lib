#ifndef _THOR_SCSI_CORE_ELEMENTS_BASIS_H_
#define _THOR_SCSI_CORE_ELEMENTS_BASIS_H_ 1

/**
   Definitions common for all elements

 */
#include <vector>
#include <string>
#include <tps/ss_vect.h>
#include <tps/tps_type.h>
// #include <thor_scsi/core/cells.h>
#include <thor_scsi/core/internals.h>
#include <thor_scsi/core/cell_void.h>
// #include <thor_scsi/core/elements_enums.h>
#include <thor_scsi/core/config.h>
#include <thor_scsi/core/aperture.h>


namespace thor_scsi {
	namespace core {
		//< Element virtual base class.
		class ElemType : public CellVoid {
		public:
			bool
			Reverse = true;                   ///< reverse elements: rearange the elements in reveresed order
#if 0
			// should not be required any more
			PartsKind
			Pkind;                     ///<  Enumeration for magnet types.
#endif
			/**
			 * @brief basic element type
			 *
			 * "length" is option, and is 0.0 if omitted.
			 *
			 */
			inline ElemType(const Config & config) : CellVoid(config) {
				const double l = config.get<double>("L", 0.0);
				this->setLength(l);
			};

			ElemType(ElemType&& o) : CellVoid(std::move(o)), PL(std::move(o.PL)) {
				// std::cerr << "ElemType move ctor " << this->name <<  std::endl;
			}

			virtual inline double getLength(void) const final { return this->PL;};
			/**
			 * @brief Set the length of the element [m]
			 *
			 * field interpolation treats length 0 special.
			 *
			 */
			virtual inline void setLength(const double length) {
				this->PL = length;
			}

			/**
			 * Todo: implement taking stream or as ostream operator ....
			 */
			virtual void show(std::ostream& strm, int level) const override;
			// C++ templates not supported for virtual functions.
			/**
			 * @brief Propagator step for phase space.
			 *
			 * Args:
			 *    ps : phase space
			 *
			 * Todo:
			 *      make config constant (after config has been reworked)
			 *      rename to propagate!
			 */
			virtual void pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps) = 0;
			virtual void pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps) = 0;
			/*
			 * the non linear tps part ... to be made
			 */
			// virtual void pass(thor_scsi::core::ConfigType &conf, ss_vect<tps_nlin> &ps) = 0;

			// why is that required here?
			// template<typename T>
			// void pass(thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);

			/**
			 * @brief: check if the particles amplitude is within aperture
			 *
			 * @returns: not lost (true if within aperture)
			 *
			 * Currently  it only checks the 2D amplitude at the current spot
			 * given that an aperture is set
			 *
			 * returns always true if no aperture registered
			 */
			template<typename T>
			inline bool checkAmplitude(const ss_vect<T> &ps){
				if(!this->m_aperture){
					return true;
				}
				double x = is_double<T>::cst(ps[x_]);
				double y = is_double<T>::cst(ps[y_]);

				double d = this->m_aperture->isWithin(x, y);

				if (d < 0e0){
					return false;
				}
				return true;

			}

#if 0
			virtual void show(std::ostream&, int level) const;
#if 0
			std::string repr_elem(void);     ///< auxilliary function providing a string of common information
			                                 ///< required for the different elements
			virtual std::string repr(void) = 0;  ///< representation similar to prt_elem but a bit more pythonic
			                                 ///< used by python interface to generate the information for
			                                 ///< :meth:`__repr__`
#endif
			/**
			 * Todo:
			 *    Check if that is still missing an overloaded method?
			 *
			 * If understood coorectly one should review if a element factory is requireed.
			 */
			virtual ElemType* Elem_Init(const thor_scsi::core::ConfigType &conf, const bool reverse)
				{ return NULL; };
			void prt_elem(const std::string &);
			virtual void print(const std::string &) {};

			virtual void SetdS(void) {}; ///< Eucledian Group: dx, dy
			virtual void SetdT(void) {}; ///< Eucledian Group: Roll angle
			virtual void SetPB(const int n) {}; ///< Multipoles (total numbers)
			virtual double GetdT(void) { return 0e0; };
			virtual double GetPB(const int n) { return 0e0; };
#endif

			inline auto getAperture(void) const {
				return std::const_pointer_cast<thor_scsi::core::TwoDimensionalAperture>(this->m_aperture);
			}

			void setAperture(std::shared_ptr<thor_scsi::core::TwoDimensionalAperture> ap){
				this->m_aperture = ap;
			}

		protected:
			double PL = 0.0;                        ///< Length[m].

		private:
			// currently only implementing 2D apertures
			std::shared_ptr<thor_scsi::core::TwoDimensionalAperture> m_aperture;

		};

#if 0
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
#endif
	}
}
#endif /*  _THOR_SCSI_CORE_ELEMENTS_BASIS_H_  */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
