#ifndef _THOR_SCSI_ELEMENTS_CORRECTOR_H_
#define _THOR_SCSI_ELEMENTS_CORRECTOR_H_

#include <thor_scsi/elements/classical_magnet.h>

namespace thor_scsi::elements {
	/**
	 *
	 * @todo: consider defining it as piggy pack power converter
	 */

	class CorrectorType : public ClassicalMagnet {
	public:
		inline CorrectorType(const Config &config) : ClassicalMagnet(config)  {
			// field kick insists that it is specified ...
		}
		const char* type_name(void) const override { return "Corrector"; };
	};

	class SteererType: public CorrectorType {
	public:
		inline SteererType (const Config &config) : CorrectorType(config){
		}
		inline int getMainMultipoleNumber(void) const override final {
			return 1;
		};
	};

	class HorizontalSteererType : public SteererType {
	  public:
		inline HorizontalSteererType (const Config &config) : SteererType(config){}
		const char* type_name(void) const override final { return "HorizontalSteerer"; };
		inline bool isSkew(void) const override final {	return false; };

	};

	class VerticalSteererType : public SteererType {
	  public:
		inline VerticalSteererType (const Config &config) : SteererType(config){}
		const char* type_name(void) const override final { return "VerticalSteerer"; };
		inline bool isSkew(void) const override final { return true; };
	};

} // Name space

#endif // _THOR_SCSI_ELEMENTS_CORRECTOR_H_
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
