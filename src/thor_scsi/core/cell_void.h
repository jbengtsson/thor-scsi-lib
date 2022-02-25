#ifndef _THOR_SCSI_CORE_CELL_VOID_H_
#define _THOR_SCSI_CORE_CELL_VOID_H_

/**
 * Based heavly on the code of FLAME
 * Or better said:: just a little modified from flame
 */

#include <flame/core/config.h>
#include <flame/core/util.h>
#include <tps/ss_vect.h>
#include <tps/tps_type.h>


namespace thor_scsi::core {
	struct CellVoid;
        struct Machine;
        /**
	 * @brief Allow inspection of intermediate State
	 *
	 * Use with ElementVoid::set_observer() to associate an observer with an Element.
	 * During Machine::propagate() a call will be made to Observer::view()
	 * with the element's output State.
	 */
	struct Observer : public boost::noncopyable
	{
		virtual ~Observer() {}
		//! Called from within Machine::propagate()

		template <typename State>
		void view(const CellVoid* elem, const State * state);

		// view for state vector of double
		virtual void view(const CellVoid* elem, const ss_vect<double> &ps) = 0;
		// view for state vector of truncated power series
		virtual void view(const CellVoid* elem, const ss_vect<tps> &ps) = 0;
	};



       /**
	* @brief Base class for all simulated elements ... but here agnostic from advance
	*
	* Sub-classes of ElementVoid must be registered with Machine::registerElement
	* before they will be found by Machine::Machine().c
	*/
	struct CellVoid : public boost::noncopyable
	{
		/**
		 * @brief Construct this element using the provided Config.
		 *
		 * Base class ctor makes use of Config parameters "name"
		 * "name" is required.
		 *
		 * Sub-classes are allowed to require certain parameters to be provided.
		 *
		 * @throws KeyError       If a required parameter is missing
		 * @throws std::bad_variant_access If a parameter exists, but has the wrong value type
		 */
		CellVoid(const Config& conf);
		virtual ~CellVoid();

		/** Sub-classes must provide an approprate short description string.
		 *  Must match the type name passed to Machine::registerElement().
		 */
		virtual const char* type_name(void) const =0;

		//! Propogate the given State through this Element
		//! Still dicussing which name to use
		// virtual void advance(StateBase& s) =0;

		//! The Config used to construct this element.
		inline const Config& conf() const {return p_conf;}

		const std::string name; //!< Name of this element (unique in its Machine)
		size_t index; //!< Index of this element (unique in its Machine)


		//! The current observer, or NULL
		Observer *observer() const { return p_observe; }
		/** Add Observer which will inspect the output State of this Element.
		 *  Observer instance musy outlive the Element.
		 * @param o A new Observer or NULL, will replace any existing pointer.
		 */
		void set_observer(Observer *o) { p_observe = o; }

		//! Print information about the element.
		//! level is a hint as to the verbosity expected by the caller.
		virtual void show(std::ostream&, int level) const;

		//! Used by Machine::reconfigure() to avoid re-alloc (and iterator invalidation)
		//! Assumes other has the same type.
		//! Sub-classes must call base class assign()
		//! Come c++11 this can be replaced with a move ctor.
		// virtual void assign(const ElementVoid* other ) =0;

	private:
		Observer *p_observe;
		Config p_conf;
	        friend class Machine;
	};

	inline
	std::ostream& operator<<(std::ostream& strm, const CellVoid& cell_void)
	{
		cell_void.show(strm, 0);
		return strm;
	}
}



#endif /* _THOR_SCSI_CORE_CELL_VOID_H_ */
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
