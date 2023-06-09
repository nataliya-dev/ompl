

#ifndef OMPL_MULTILEVEL_PLANNERS_BundleSpace_QBRRTIMPL_
#define OMPL_MULTILEVEL_PLANNERS_BundleSpace_QBRRTIMPL_
#include <ompl/multilevel/datastructures/BundleSpaceGraph.h>
#include <ompl/datastructures/PDF.h>

namespace ompl
{
    namespace base
    {
        OMPL_CLASS_FORWARD(OptimizationObjective);
    }
    namespace multilevel
    {
        /** \brief Implementation of BundleSpace Rapidly-Exploring Random Trees Algorithm*/
        class QBRRTImpl : public ompl::multilevel::BundleSpaceGraph
        {
            using BaseT = BundleSpaceGraph;

        public:
            QBRRTImpl(const ompl::base::SpaceInformationPtr &si, BundleSpace *parent_);
            virtual ~QBRRTImpl() override;

            /** \brief One iteration of RRT with adjusted sampling function */
            virtual void grow() override;
        };
    }  // namespace multilevel
}  // namespace ompl

#endif
