#ifndef _THOR_SCSI_RADIATION_OBSERVER_API_H_
#define _THOR_SCSI_RADIATION_OBSERVER_API_H_ 1
#include <thor_scsi/elements/field_kick_api.h>
#include <thor_scsi/core/elements_basis.h>
#include <gtpsa/tpsa.hpp>
#include <gtpsa/ss_vect.h>
#include <tps/tps_type.h>

namespace thor_scsi::elements {
	using thor_scsi::core::ObservedState;
    template<class C>
	class RadiationDelegateInterfaceKnobbed {
	public:
	        virtual ~RadiationDelegateInterfaceKnobbed(void){}
		virtual void view(const thor_scsi::core::ElemTypeKnobbed<C>& elem, const gtpsa::ss_vect<double>      &ps, const enum ObservedState, const int cnt) = 0;
	    // virtual void view(const thor_scsi::core::ElemTypeKnobbed<C>& elem, const gtpsa::ss_vect<tps>         &ps, const enum ObservedState, const int cnt) = 0;
		virtual void view(const thor_scsi::core::ElemTypeKnobbed<C>& elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum ObservedState, const int cnt) = 0;
		virtual void show(std::ostream& strm, int level) const{
			strm << "RadiationDelegateInterface";
		}
		// support for python
		std::string repr(void) const;

	};
    typedef RadiationDelegateInterfaceKnobbed<thor_scsi::core::StandardDoubleType> RadiationDelegateInterface;

    template<class C>
	class RadiationDelegateKickInterfaceKnobbed {
	public:
	        virtual ~RadiationDelegateKickInterfaceKnobbed(void){}
		virtual void view(const FieldKickAPIKnobbed<C>& kick, const gtpsa::ss_vect<double>      &ps, const enum ObservedState, const int cnt) = 0;
		virtual void view(const FieldKickAPIKnobbed<C>& kick, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum ObservedState, const int cnt) = 0;
	    // virtual void view(const FieldKickAPIKnobbed<C>& kick, const gtpsa::ss_vect<tps>         &ps, const enum ObservedState, const int cnt) = 0;
		virtual void show(std::ostream& strm, int level) const{
			strm << "RadiationDelegateKickInterface";
		}
		// support for python
		std::string repr(void) const;

	};
	/**
	 * @brief Radiation effect due to local field
	 *
	 * @f[
	 *     \frac{\mathrm{d}\delta}{\mathrm{d}(ds)} =
	 *           -C_{\gamma} \, E_0^3 \, (1 + \delta)^2 \, \left( \frac{B_{perp}}{B \rho}\right)^2
	 *                \frac{1}{2\pi}
	 * @f]
	 *
	 * @todo: depends on energy of ring .... currently taken from config ...
	 *
	 * M. Sands "The hysics of Electron Storage Rings" SLAC-121, p. 98.
	 */
	//inline
    template<typename C>
	std::ostream& operator<<(std::ostream& strm, const RadiationDelegateKickInterfaceKnobbed<C>& rd)
	{
		rd.show(strm, 0);
		return strm;
	}

    template<class C>
	std::ostream& operator<<(std::ostream& strm, const RadiationDelegateInterfaceKnobbed<C>& rd)
	{
		rd.show(strm, 0);
		return strm;
	}

} // namespace thor_scsi::elements

#endif /* _THOR_SCSI_RADIATION_OBSERVER_API_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
