

#include <ompl/multilevel/planners/qbatch/QBRRTImpl.h>
#include <ompl/multilevel/datastructures/graphsampler/GraphSampler.h>
#include <ompl/tools/config/SelfConfig.h>
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

ompl::multilevel::QBRRTImpl::QBRRTImpl(const base::SpaceInformationPtr &si, BundleSpace *parent_) : BaseT(si, parent_)
{
    setName("QBRRTImpl" + std::to_string(id_));
    setImportance("exponential");
    setGraphSampler("randomvertex");
    getGraphSampler()->disableSegmentBias();
}

ompl::multilevel::QBRRTImpl::~QBRRTImpl()
{
}

void ompl::multilevel::QBRRTImpl::grow()
{
    //(0) If first run, add start configuration
    if (firstRun_)
    {
        init();
        firstRun_ = false;

        findSection();
    }

    //(1) Get Random Sample
    sampleBundleGoalBias(xRandom_->state);

    //(2) Get Nearest in Tree
    const Configuration *xNearest = nearest(xRandom_);

    //(3) Connect Nearest to Random (within range)
    Configuration *xNext = extendGraphTowards_Range(xNearest, xRandom_);

    //(4) If extension was successful, check if we reached goal
    if (xNext && !hasSolution_)
    {
        if (isDynamic())
        {
            double dist;
            bool satisfied = getGoalPtr()->isSatisfied(xNext->state, &dist);
            if (dist < bestCost_.value())
            {
                bestCost_ = base::Cost(dist);
            }
            if (satisfied)
            {
                goalConfigurations_.push_back(xNext);
                hasSolution_ = true;
            }
        }
        else
        {
            bool satisfied = getGoalPtr()->isSatisfied(xNext->state);
            if (satisfied)
            {
                goalConfigurations_.push_back(xNext);
                hasSolution_ = true;
            }
        }
    }
}
