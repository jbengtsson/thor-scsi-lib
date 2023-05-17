#define BOOST_TEST_MODULE standard_observer
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <thor_scsi/elements/standard_observer.h>

namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;

BOOST_AUTO_TEST_CASE(test20_standard_observer_tpsa)
{
    auto desc = std::make_shared<gtpsa::desc> (6, 2);
    gtpsa::ss_vect<gtpsa::tpsa> ps(desc, 2);
    tse::StandardObserver observer;

    ps[1].set(0, 2);
    BOOST_CHECK_CLOSE(ps.cst()[1], 2, 1e-12);

    // counting on that the element is not accessed ...
    observer.view(nullptr, ps, tsc::ObservedState::end, 1);
    ps[5].set(0, 7);
    BOOST_CHECK_CLOSE(ps.cst()[1], 2, 1e-12);
    BOOST_CHECK_CLOSE(ps.cst()[5], 7, 1e-12);

    // Observer should have allocated a new object
    auto test_ps = observer.getTruncatedPowerSeriesA();
    BOOST_CHECK_NE(&ps, test_ps.get());
    BOOST_CHECK_CLOSE(test_ps->cst()[1], 2, 1e-12);
    BOOST_CHECK_CLOSE(test_ps->cst()[5], 0, 1e-12);

}
