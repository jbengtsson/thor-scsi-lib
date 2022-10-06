#ifndef _THOR_SCSI_CORE_CELL_VOID_H_
#define _THOR_SCSI_CORE_CELL_VOID_H_

/**
 * Based heavly on the code of FLAME
 * Or better said:: just a little modified from flame
 */

#include <flame/core/config.h>
#include <flame/core/util.h>
#include <gtpsa/ss_vect.h>
#include <gtpsa/tpsa.hpp>
#include <tps/tps_type.h>


namespace thor_scsi::core {
	struct CellVoid;
        struct Machine;

	/**
	 *
	 *
	 * Experimenting with the observer.
	 * influenced by blueske
	 */
	enum ObservedState
	{
		/**
		 * a standard unspecified event, e.g. one integration step
		 * typically issued within the _pass or _localPass method
		 *  of an element
		 */
		event = 0,
		/**
		 * observation at start, e.g. called before
		 * called into the elements _pass method
		 */
		start,
		/**
		 * observation at end, e.g. called before
		 * called into the elements _pass method
		 */
		end,
		/**
		 * typically should never reach an observer
		 * Don't know if this state will be used
		 */
		undefined,
		/**
		 * an failure occured during processing. can be useful
		 * sometimes to record the state
		 */
		failure
	};
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

		//! Called from within the different elements ...

		// view for state vector of double
		virtual void view(std::shared_ptr<const CellVoid> elem, const gtpsa::ss_vect<double>      &ps, const enum ObservedState, const int cnt) = 0;
		// view for state vector of truncated power series
		virtual void view(std::shared_ptr<const CellVoid> elem, const gtpsa::ss_vect<tps>         &ps, const enum ObservedState,  const int cnt) = 0;
	        virtual void view(std::shared_ptr<const CellVoid> elem, const gtpsa::ss_vect<gtpsa::tpsa> &ps, const enum ObservedState,  const int cnt) = 0;

		virtual void show(std::ostream& strm, int level) const = 0;

		///< support for python __repr__ uses show
		std::string repr(void);
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
		CellVoid(CellVoid&& o);
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


		//! The current observer as shared ptr (required for python)
		std::shared_ptr<Observer> observer() const { return p_observe; }
		/** Add Observer which will inspect the output State of this Element.
		 *  Observer instance musy outlive the Element.
		 * @param o A new Observer or NULL, will replace any existing pointer.
		 */
		void set_observer(std::shared_ptr<Observer> o) { p_observe = o; }

		//! Print information about the element.
		//! level is a hint as to the verbosity expected by the caller.
		virtual void show(std::ostream&, int level) const;

		virtual std::string prettyClassname(void) const;
		// facilite python representation
		std::string repr(void) const;
		std::string pstr(void) const;
		//! Used by Machine::reconfigure() to avoid re-alloc (and iterator invalidation)
		//! Assumes other has the same type.
		//! Sub-classes must call base class assign()
		//! Come c++11 this can be replaced with a move ctor.
		// virtual void assign(const ElementVoid* other ) =0;

	private:
		std::shared_ptr<Observer> p_observe;
		Config p_conf;
	        friend class Machine;
	};

	inline
	std::ostream& operator<<(std::ostream& strm, const CellVoid& cell_void)
	{
		cell_void.show(strm, 0);
		return strm;
	}

	inline
	std::ostream& operator<<(std::ostream& strm, const Observer& observer)
	{
		observer.show(strm, 0);
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
