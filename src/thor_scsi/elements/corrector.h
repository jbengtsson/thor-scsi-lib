#ifndef _THOR_SCSI_ELEMENTS_CORRECTOR_H_
#define _THOR_SCSI_ELEMENTS_CORRECTOR_H_

#include <thor_scsi/elements/classical_magnet.h>

namespace thor_scsi::elements {
	/**
	 *
	 * @todo: consider defining it as piggy pack power converter
	 */

	template<class C>
	class CorrectorTypeWithKnob : public ClassicalMagnetWithKnob<C> {
	public:
		inline CorrectorTypeWithKnob(const Config &config)
			: ClassicalMagnetWithKnob<C>(config)
			{ /* field kick insists that it is specified ... */ }
		const char* type_name(void) const override { return "Corrector"; };
	};

	template<class C>
	class SteererTypeWithKnob: public CorrectorTypeWithKnob<C> {
	public:
		inline SteererTypeWithKnob (const Config &config)
			: CorrectorTypeWithKnob<C>(config)
			{}
		inline int getMainMultipoleNumber(void) const override final {
			return 1;
		};
	};

	template<class C>
	class HorizontalSteererTypeWithKnob : public SteererTypeWithKnob<C> {
	  public:
		inline HorizontalSteererTypeWithKnob (const Config &config)
			: SteererTypeWithKnob<C>(config)
			{}
		const char* type_name(void) const override final { return "HorizontalSteerer"; };
		inline bool isSkew(void) const override final {	return false; };

	};

	template<class C>
	class VerticalSteererTypeWithKnob : public SteererTypeWithKnob<C> {
	  public:
		inline VerticalSteererTypeWithKnob (const Config &config)
			: SteererTypeWithKnob<C>(config)
			{}
		const char* type_name(void) const override final { return "VerticalSteerer"; };
		inline bool isSkew(void) const override final { return true; };
	};

	typedef HorizontalSteererTypeWithKnob<thor_scsi::core::StandardDoubleType> HorizontalSteererType;
	typedef VerticalSteererTypeWithKnob<thor_scsi::core::StandardDoubleType>   VerticalSteererType;

} // Name space

#endif // _THOR_SCSI_ELEMENTS_CORRECTOR_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
