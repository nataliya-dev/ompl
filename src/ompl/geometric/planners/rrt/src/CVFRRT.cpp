#include <utility>

#include "ompl/geometric/planners/rrt/CVFRRT.h"
#include "ompl/base/goals/GoalSampleableRegion.h"

ompl::geometric::CVFRRT::CVFRRT(const base::SpaceInformationPtr &si, VectorField vf, double lambda,
                                unsigned int update_freq)
  : RRT(si), vf_(std::move(vf)), lambda_(lambda), nth_step_(update_freq)
{
    setName("CVFRRT");
    maxDistance_ = si->getStateValidityCheckingResolution();
}

ompl::geometric::CVFRRT::~CVFRRT() = default;

void ompl::geometric::CVFRRT::clear()
{
    RRT::clear();
    step_ = 0;
}

void ompl::geometric::CVFRRT::setup()
{
    RRT::setup();
    vfdim_ = si_->getStateSpace()->getValueLocations().size();
}

Eigen::VectorXd ompl::geometric::CVFRRT::getNewDirection(const base::State *qnear, const base::State *qrand)
{
    OMPL_INFORM("========== getNewDirection");

    // Set vrand to be the normalized vector from qnear to qrand
    Eigen::VectorXd vrand(vfdim_);
    Eigen::VectorXd vnear(vfdim_);
    for (unsigned int i = 0; i < vfdim_; i++)
    {
        vrand[i] = *si_->getStateSpace()->getValueAddressAtIndex(qrand, i) -
                   *si_->getStateSpace()->getValueAddressAtIndex(qnear, i);
        vnear[i] = *si_->getStateSpace()->getValueAddressAtIndex(qnear, i);
    }
    vrand /= si_->distance(qnear, qrand);

    // Get the vector at qnear, and normalize
    Eigen::VectorXd vfield = vf_(qnear);
    // const double lambdaScale = vfield.norm();
    // // In the case where there is no vector field present, vfield.norm() == 0,
    // // return the direction of the random state.
    // if (lambdaScale < std::numeric_limits<float>::epsilon())
    //     return vrand;
    // vfield /= lambdaScale;

    // get vrand and vf weights
    updateGain();
    double alpha = 1.0 - std::exp(lambda_ * step_);
    double beta = 1.0 - alpha;

    if (alpha > 0.999)
    {
        alpha = 1.0;
        beta = 0.0;
    }
    OMPL_INFORM("step: %d", step_);
    OMPL_INFORM("alpha: %f", alpha);
    OMPL_INFORM("beta: %f", beta);

    // calculate vnew
    Eigen::VectorXd vnew = vfield + vrand;

    OMPL_INFORM("VRAND");
    for (std::size_t i = 0; i < 7; i++)
    {
        std::cout << vrand[i] << ", ";
    }
    std::cout << std::endl;

    OMPL_INFORM("VNEAR");
    for (std::size_t i = 0; i < 7; i++)
    {
        std::cout << vnear[i] << ", ";
    }
    std::cout << std::endl;

    OMPL_INFORM("VFIELD");
    for (std::size_t i = 0; i < 7; i++)
    {
        std::cout << vfield[i] << ", ";
    }
    std::cout << std::endl;

    OMPL_INFORM("VNEW");
    for (std::size_t i = 0; i < 7; i++)
    {
        std::cout << vnew[i] << ", ";
    }
    std::cout << std::endl;

    return vnew;
}

void ompl::geometric::CVFRRT::updateGain()
{
    step_++;
}

ompl::geometric::CVFRRT::Motion *ompl::geometric::CVFRRT::extendTree(Motion *m, base::State *rstate,
                                                                     const Eigen::VectorXd &v)
{
    base::State *newState = si_->allocState();
    si_->copyState(newState, m->state);

    double d = si_->distance(m->state, rstate);
    if (d > maxDistance_)
        d = maxDistance_;

    const base::StateSpacePtr &space = si_->getStateSpace();
    for (unsigned int i = 0; i < vfdim_; i++)
        *space->getValueAddressAtIndex(newState, i) += d * v[i];
    if (!v.hasNaN() && si_->checkMotion(m->state, newState))
    {
        OMPL_INFORM("Valid motion.");
        auto *motion = new Motion();
        motion->state = newState;
        motion->parent = m;
        nn_->add(motion);
        return motion;
    }
    else
    {
        OMPL_INFORM("Motion invalid.");
        si_->freeState(newState);
        return nullptr;
    }
}

ompl::base::PlannerStatus ompl::geometric::CVFRRT::solve(const base::PlannerTerminationCondition &ptc)
{
    checkValidity();
    base::Goal *goal = pdef_->getGoal().get();
    auto *goal_s = dynamic_cast<base::GoalSampleableRegion *>(goal);

    if (!sampler_)
        sampler_ = si_->allocStateSampler();

    while (const base::State *st = pis_.nextStart())
    {
        auto *motion = new Motion(si_);
        si_->copyState(motion->state, st);
        nn_->add(motion);
    }

    if (nn_->size() == 0)
    {
        OMPL_ERROR("%s: There are no valid initial states!", getName().c_str());
        return base::PlannerStatus::INVALID_START;
    }

    OMPL_INFORM("%s: Starting planning with %u states already in datastructure", getName().c_str(), nn_->size());

    Motion *solution = nullptr;
    Motion *approxsol = nullptr;
    double approxdif = std::numeric_limits<double>::infinity();
    auto *rmotion = new Motion(si_);
    base::State *rstate = rmotion->state;
    base::State *xstate = si_->allocState();

    while (ptc == false)
    {
        // Sample random state (with goal biasing)
        if (goal_s && rng_.uniform01() < goalBias_ && goal_s->canSample())
            goal_s->sampleGoal(rstate);
        else
            sampler_->sampleUniform(rstate);

        // Find closest state in the tree
        Motion *nmotion = nn_->nearest(rmotion);

        // Modify direction based on vector field before extending
        Motion *motion = extendTree(nmotion, rstate, getNewDirection(nmotion->state, rstate));
        if (!motion)
            continue;

        // Check if we can connect to the goal
        double dist = 0;
        bool sat = goal->isSatisfied(motion->state, &dist);
        if (sat)
        {
            approxdif = dist;
            solution = motion;
            break;
        }
        if (dist < approxdif)
        {
            approxdif = dist;
            approxsol = motion;
        }
    }

    bool solved = false;
    bool approximate = false;
    if (solution == nullptr)
    {
        solution = approxsol;
        approximate = true;
    }

    if (solution != nullptr)
    {
        lastGoalMotion_ = solution;

        // Construct the solution path
        std::vector<Motion *> mpath;
        while (solution != nullptr)
        {
            mpath.push_back(solution);
            solution = solution->parent;
        }

        // Set the solution path
        auto path(std::make_shared<PathGeometric>(si_));
        for (int i = mpath.size() - 1; i >= 0; --i)
            path->append(mpath[i]->state);
        pdef_->addSolutionPath(path, approximate, approxdif, name_);
        solved = true;
    }

    si_->freeState(xstate);
    if (rmotion->state)
        si_->freeState(rmotion->state);
    delete rmotion;

    OMPL_INFORM("%s: Created %u states", getName().c_str(), nn_->size());

    return {solved, approximate};
}
