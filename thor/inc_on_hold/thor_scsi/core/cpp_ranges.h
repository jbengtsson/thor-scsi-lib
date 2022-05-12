#ifndef _THOR_SCSI_CPP_RANGES_VERSION_H_
#define _THOR_SCSI_CPP_RANGES_VERSION_H_ 1

/**
 * Use std::ranges if compiler has support for it
 * Otherwise use range/v3 see e.g. https://github.com/ericniebler/range-v3
 *
 */
#include <thor_scsi/core/cpp_version.h>
#include <version>

#ifdef __cpp_lib_ranges
#warning "using std::ranges"
#include <ranges>
// namespace ranges=std::ranges;
#else
#warning "using range-v3"
#error "using range-v3"
#include <vector>
#include <set>
#include <unordered_set>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#endif

#endif /*  _THOR_SCSI_CPP_RANGES_VERSION_H_ */
