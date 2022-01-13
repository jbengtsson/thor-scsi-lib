#ifndef _THOR_SCSI_ELEMENTS_FILTERS_H_
#define _THOR_SCSI_ELEMENTS_FILTERS_H_ 1

/**
 * Todo:
 *    automatically include std::ranges if available
 */
// #include <algorithm>
// #ifdef __cpp_lib_ranges
// #error "ranges"
// #include <ranges>
// #else
// #error "No ranges"
// #endif

#include <thor_scsi/core/cpp_version.h>

#include <vector>
#include <set>
#include <unordered_set>
#include <cassert>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>

#include <thor_scsi/core/elements.h>

/*
 * Todo:
 *     Consider use of templates...
 */
namespace thor_scsi{
	namespace elements {
		/**
		 * Check if element is a drift (if the dynamic cast is successful)
		 */
		inline bool is_drift_type(thor_scsi::elements::ElemType *elem){
			return (dynamic_cast<thor_scsi::elements::DriftType *>(elem) != nullptr);
		}
                /**
		 * Cast Elem to drift type
		 *
		 * Report and abort if cast fails. Forseen to be only used on
		 * ranges::views::transform after elements have been filtered.
		 *
		 * Warning:
		 *       Not part of user API
		 */
		inline thor_scsi::elements::DriftType *
			as_drift_type( thor_scsi::elements::ElemType * elem)
		{
			DriftType *d = nullptr;
			d = dynamic_cast< thor_scsi::elements::DriftType *>(elem);
			if(!d){
				std::cerr << "Could not cast " << elem->Name <<
					" to DriftType" <<std::endl;
			}
			assert(d);
			return d;
		}
		/**
		 *  a view on the drits ... as drifts
		 *
		 * Typical use case
		 *
		 * ..::
		 *
		 *      for(auto : filter_drifts(lat.elems)){
		 *            do_something
		 *      }
		 */
		inline auto filter_drifts(std::vector< thor_scsi::elements::ElemType *> elems){
			return (elems | ranges::views::filter(is_drift_type)
				      | ranges::views::transform(as_drift_type));
		}
                /**
		 *  check if it is a Mpole type
		 */
		inline bool is_mpole_type(thor_scsi::elements::ElemType *elem){
			return dynamic_cast<thor_scsi::elements::MpoleType *>(elem) != nullptr;
		}
#if 0
                /**
		 *  a view on the Mpole Types ... as ElemType only
		 *
		 * returns a view on the Mpoles as ElementType
		 *
		 * Warning:
		 *      Not part of the API
		 */
		inline auto filter_mpole_as_elem(std::vector<ElemType *> elems){
			return (elems | ranges::views::filter(is_derived));
		}
#endif

                /**
		 * Cast Elem to mpole type
		 *
		 * Report and abort if cast fails. Forseen to be only used on
		 * ranges::views::transform after elements have been filtered.
		 *
		 * Warning:
		 *       Not part of user API
		 */
		inline thor_scsi::elements::MpoleType *
			as_mpole_type( thor_scsi::elements::ElemType * elem)
		{
			MpoleType *d = nullptr;
			d = dynamic_cast< thor_scsi::elements::MpoleType *>(elem);
			if(!d){
				std::cerr << "Could not cast " << elem->Name <<
					" to Mpole" <<std::endl;
			}
			assert(d);
			return d;
		}

               /**
		*  a view on the Mpole Types ... as MpoleTypes
		*
		* Typical use case
		*
		* ..::
		*
		*      for(auto : lat.elems){
		*            do_something
		*      }
		*/
                //inline ranges::transform_view<ranges::filter_view<ranges::ref_view<std::vector<Element*> >, bool (*)(Element*)>, DerivedElement* (*)(Element*)>
		inline auto filter_mpole_types(std::vector< thor_scsi::elements::ElemType *> elems){
			return (elems | ranges::views::filter(is_mpole_type)
				      | ranges::views::transform(as_mpole_type));
		}

		/**
		 * if a mpole is a sextupole
		 * Todo: consider make it a method of mpole
		 */
		inline bool mpole_is_sextupole(thor_scsi::elements::MpoleType * a_mpole){
			return (a_mpole->n_design == Sext);
		}
		inline bool elem_is_sextupole(thor_scsi::elements::ElemType * a_elem){
			if (!is_mpole_type(a_elem)){
				return false;
			}
			return mpole_is_sextupole(as_mpole_type(a_elem));
		}
		/**
		 * a view on the sextupoles
		 */
		inline auto filter_sextupoles(std::vector< thor_scsi::elements::ElemType *> elems){
			return (filter_mpole_types(elems) | ranges::views::filter(mpole_is_sextupole));
		}

		/*
		 * is the mpole a bend
		 *
		 * Todo:
		 *     is the check oof Pirho appropriate (floating point comparison)
		 *     should it be a method of mpole
		 */
		inline bool mpole_is_bend(thor_scsi::elements::MpoleType * a_mpole){
			return (a_mpole->Pirho != 0e0);
		}
		/**
		 * a view on (design) bends)
		 */

		inline auto filter_bends(std::vector< thor_scsi::elements::ElemType *> elems){
			return (filter_mpole_types(elems) | ranges::views::filter(mpole_is_bend));
		}


	}; /* elements */
}; /* thor_scsi */
#endif /*  _THOR_SCSI_ELEMENTS_FILTERS_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
