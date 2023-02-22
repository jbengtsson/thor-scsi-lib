#ifndef _THOR_SCSI_ELEMENTS_BENDING_H_
#define _THOR_SCSI_ELEMENTS_BENDING_H_

#include <thor_scsi/elements/classical_magnet.h>

namespace thor_scsi::elements {
	template<class C>
	class BendingTypeWithKnob : public ClassicalMagnetWithKnob<C> {
		/**
		 *1
		 * Can derive from classical magnet -> mpole -> field kick as this is currently deriving from
		 * LocalGalileanProt
		 *
		 */

		/*
		     <name>: Bending,
		     L      = <length>, ( [m] )
		     T      = <bending angle>, ( [degree] )
		     T1     = <entrance angle>, ( [degree] )
		     T2     = <exit angle>, ( [degree] )
		     gap    = <total magnet gap>, ( [m] )
		     K      = <K-value>, ( [m^-2] )
                     ( K > 0 : focusing in horizontal plane )
                     ( K < 0 : defocusing in horizontal plane )
		     N      = <# of kicks>,
		     method = <method>, ( 2 or 4. The method to divide Q into slices.)
		     ( The detail of <method> will be discussed later.)
		     Default value is 2.
		     Roll   = <roll angle>, ( [deg], design roll angle )
		     HOM    = (i, <Bi>, <Ai>, ( multipole compoments (power series) )
		     j, <Bj>, <Aj>, ( Systematic error Only )
		     ............   ( Random errors are assigned )
                     n, <Bn>, <An>); ( in a Program File using procedures )

		     Example

		     B: bending, L=0.70, T=10.0, T1:=5.0, T2:=5.0, K=-1.0, N=8, Method=2;
		 */
	public:
		inline BendingTypeWithKnob(const Config &config) : ClassicalMagnetWithKnob<C>(config){
			const double gradient = config.get<double>("K");
			this->getMultipoles()->setMultipole(2, gradient);

			/* set pirho */
			const double phi = degtorad(this->getBendingAngle());
			const double length = this->PL;
			if(length == 0.0){
				/* nothing should be required any more */
				this->setCurvature(0.0);
				return;
			}
			this->setCurvature(phi / length);
		}

		inline int getMainMultipoleNumber(void) const override final {
			return 1;
		};
		inline bool isSkew(void) const override final {
			return false;
		};
		const char* type_name(void) const override final { return "Bending"; };
	};

	typedef BendingTypeWithKnob<thor_scsi::core::StandardDoubleType> BendingType;

} // Name space

#endif // _THOR_SCSI_ELEMENTS_BENDING_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
