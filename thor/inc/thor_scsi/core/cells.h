#ifndef _THOR_SCSI_CORE_CELL_H_
#define _THOR_SCSI_CORE_CELL_H_
#include <thor_scsi/exceptions.h>
#include <flame/core/base.h>
#include <thor_scsi/core/transform.h>

#include <vector>
namespace thor_scsi {
	namespace elements {

		/**  LEGO lattice structure.
		 *
		 * Now derived from flame elementvoid
		 *
		 * Todo:
		 *      Clear up interface
		 *      double implementation ???
		 *
		 *      * propagate / pass
		 */
		class CellType : public ElementVoidBase {
		public:
			inline CellType(const Config& config): ElementVoidBase(config){
			}
	    int
	      Fnum,                      ///< Element Family #.
	      Knum;                      ///< Element Kid #.
	    double
	      S,                         ///< position in the ring. S coordinate of the desgin orbit.
	      curly_dH_x;                ///< Contribution to curly_H_x.
	    thor_scsi::core::Euclidian2DTranform transform;
	    std::vector<double>
	      dI                         ///< Contribution to I[1..5].
		{0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
	      Eta{0e0, 0e0},             ///< Eta & eta' (redundant).
	      Etap{0e0, 0e0},
	      Alpha{0e0, 0e0},           ///< Twiss parameters (redundant).
	      Beta{0e0, 0e0},
	      Nu{0e0, 0e0},
	     BeamPos                    ///< Last position of the beam this cell.
             {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};

		  /* todo: review if these should be  arma arrays */
	     std::vector< std::vector<double> >
	       maxampl                    ///< Hor & ver physical aperture.
		 {{0e0, 0e0},
		     {0e0, 0e0}},
	       A                          ///< Floquet space to phase space transformation.
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}},
		 sigma                      ///< sigma matrix (redundant).
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0}};
		 CellType
		   *next_ptr;               ///< pointer to next cell (for tracking).

		 bool trace = false;

		 inline void setTraceFlag(bool flag){
		   this->trace = flag;
		 }
		  /**
		   * Euclidian Group: element 0, actually used value
		   * \f$ rms \cdot random number \f$
		   * horizontal coordinate
		   */
		 inline const double getDx(void){
			 return this->transform.getDx();
		 }
		  /**
		   * Euclidian Group: element 0, actually used value
		   * \f$ rms \cdot random number \f$
		   * vertical coordinate
		   */
		 inline const double getDy(void){
			 return this->transform.getDy();
		 }
		  /**
		   * Warning: value should be typically applied to rms and random value
		   *          be aware that it is overriden
		   */
		 inline void setDx(double value){
			 this->transform.setDx(value);
		 }
		  /**
		   * Warning: value should be typically applied to rms and random value
		   *          be aware that it is overriden
		   */
		 inline void setDy(double value){
			 this->transform.setDy(value);
		 }
		 inline double getRoll(void){
			 return this->transform.getRoll();
		 }
		 inline void setRoll(double roll){
			 this->transform.setRoll(roll);
		 }
	  };
	}
}

#endif /*  _THOR_SCSI_CORE_CELL_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
