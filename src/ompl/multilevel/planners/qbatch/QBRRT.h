

#ifndef OMPL_MULTILEVEL_PLANNERS_BUNDLESPACE_QBRRT_
#define OMPL_MULTILEVEL_PLANNERS_BUNDLESPACE_QBRRT_
#include <ompl/multilevel/datastructures/BundleSpaceSequence.h>
#include <ompl/multilevel/planners/qrrt/QBRRTImpl.h>

namespace ompl
{
    namespace multilevel
    {

        using QBRRT = BundleSpaceSequence<QBRRTImpl>;

    }  // namespace multilevel
}  // namespace ompl

#endif
