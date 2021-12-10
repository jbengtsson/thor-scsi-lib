#ifndef _THOR_SCSI_CORE_CELL_H_
#define _THOR_SCSI_CORE_CELL_H_
#include <thor_scsi/exceptions.h>

#include <vector>
namespace thor_scsi {
	namespace elements {

	  // LEGO lattice structure.
	  class CellType {
	  public:
	    int
	      Fnum,                      // Element Family #.
	      Knum;                      // Element Kid #.
	    double
	      S,                         // Position in the ring.
	      curly_dH_x;                // Contribution to curly_H_x.
	    std::vector<double>
	      dI                         // Contribution to I[1..5].
		{0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
	      dS{0e0, 0e0},              // Transverse displacement.
              dT{0e0, 0e0},              // dT = (cos(dT), sin(dT)).
	      Eta{0e0, 0e0},             // Eta & eta' (redundant).
	      Etap{0e0, 0e0},
	      Alpha{0e0, 0e0},           // Twiss parameters (redundant).
	      Beta{0e0, 0e0},
	      Nu{0e0, 0e0},
	     BeamPos                    // Last position of the beam this cell.
             {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};

	     std::vector< std::vector<double> >
	       maxampl                    // Hor & ver physical aperture.
		 {{0e0, 0e0},
		     {0e0, 0e0}},
	       A                          // Floquet space to phase space transformation.
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0}},
		 sigma                      // sigma matrix (redundant).
      {{0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
       {0e0, 0e0, 0e0, 0e0, 0e0, 0e0}};
		 CellType
		   *next_ptr;               // pointer to next cell (for tracking).

		 bool trace = false;

		 inline void setTraceFlag(bool flag){
		   this->trace = flag;
		 }
		 inline const double getDx(void){
		   return this->dS[0];
		 }
		 inline const double getDy(void){
		   return this->dS[1];
		 }
		 inline void setDx(double value){
		   this->dS[0] = value;
		 }
		 inline void setDy(double value){
		   this->dS[1] = value;
		 }
		 inline double getRoll(void){
		   // how needs that to be implemented
		   throw thor_scsi::NotImplemented();
		   return -1e100;
		 }
		 inline void setRoll(double roll){
		   // how needs that to be implemented
		   throw thor_scsi::NotImplemented();
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
