/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2015, Caleb Voss and Wilson Beebe
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Willow Garage nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

/* Authors: Caleb Voss, Wilson Beebe */

#include <utility>

#include "ompl/geometric/planners/rrt/ContactVFRRT.h"
#include "ompl/base/goals/GoalSampleableRegion.h"

namespace ompl
{
    namespace magic
    {
        /// Number of sampler to determine mean vector field norm in \ref gContactVFRRT
        static const unsigned int ContactVFRRT_MEAN_NORM_SAMPLES = 1000;
    }  // namespace magic
}  // namespace ompl

ompl::geometric::ContactVFRRT::ContactVFRRT(const base::SpaceInformationPtr &si, VectorField vf, double exploration,
                                double initial_lambda, unsigned int update_freq)
  : RRT(si), vf_(std::move(vf)), explorationSetting_(exploration), lambda_(initial_lambda), nth_step_(update_freq)
{
    setName("ContactVFRRT");
    maxDistance_ = si->getStateValidityCheckingResolution();
}

ompl::geometric::ContactVFRRT::~ContactVFRRT() = default;

void ompl::geometric::ContactVFRRT::clear()
{
    RRT::clear();
    efficientCount_ = 0;
    inefficientCount_ = 0;
    explorationInefficiency_ = 0.;
    step_ = 0;
}

void ompl::geometric::ContactVFRRT::setup()
{
    RRT::setup();
    vfdim_ = si_->getStateSpace()->getValueLocations().size();
}

Eigen::VectorXd ompl::geometric::ContactVFRRT::getNewDirection(const base::State *qnear, const base::State *qrand)
{
    OMPL_INFORM("========== getNewDirection");
    // Set vrand to be the normalized vector from qnear to qrand
    Eigen::VectorXd vrand(vfdim_);
    Eigen::VectorXd vnear(vfdim_);
    for (unsigned int i = 0; i < vfdim_; i++)
    {
        vrand[i] = *si_->getStateSpace()->getValueAddressAtIndex(qrand, i);
        vnear[i] = *si_->getStateSpace()->getValueAddressAtIndex(qnear, i);
    }
    // vrand /= si_->distance(qnear, qrand);

    // Get the vector at qnear, and normalize
    Eigen::VectorXd vfield = vf_(qrand);

    const double lambdaScale = vfield.norm();
    OMPL_INFORM("lambdaScale: %f", lambdaScale);

    // In the case where there is no vector field present, vfield.norm() == 0,
    // return the direction of the random state.
    // if (lambdaScale < std::numeric_limits<float>::epsilon())
    //     return vrand;
    // vfield /= lambdaScale;

    // // Sample a weight from the distribution
    // const double omega = biasedSampling(vrand, vfield, lambdaScale);
    // OMPL_INFORM("omega: %f", omega);

    // Determine updated direction
    // return computeAlphaBeta(omega, vrand, vfield);

    Eigen::VectorXd vnew(vfdim_);
    vnew = vfield + vrand;

    OMPL_INFORM("VRAND");
    for (std::size_t i = 0; i < vfdim_; i++)
    {
        std::cout << vrand[i] << ", ";
    }
    std::cout << std::endl;

    OMPL_INFORM("VNEAR");
    for (std::size_t i = 0; i < vfdim_; i++)
    {
        std::cout << vnear[i] << ", ";
    }
    std::cout << std::endl;

    OMPL_INFORM("VFIELD");
    for (std::size_t i = 0; i < vfdim_; i++)
    {
        std::cout << vfield[i] << ", ";
    }
    std::cout << std::endl;

    OMPL_INFORM("VNEW");
    for (std::size_t i = 0; i < vfdim_; i++)
    {
        std::cout << vnew[i] << ", ";
    }
    std::cout << std::endl;

    return vnew;
}

void ompl::geometric::ContactVFRRT::updateGain()
{
    if (step_ == nth_step_)
    {
        OMPL_INFORM("========== updateGain");

        OMPL_INFORM("before lambda_: %f", lambda_);
        OMPL_INFORM("explorationInefficiency_: %f", explorationInefficiency_);
        OMPL_INFORM("explorationSetting_: %f", explorationSetting_);

        lambda_ = lambda_ * (1 - explorationInefficiency_ + explorationSetting_);
        OMPL_INFORM("after lambda_: %f", lambda_);
        efficientCount_ = inefficientCount_ = 0;
        explorationInefficiency_ = 0;
        step_ = 0;
    }
    else
        step_++;
}

ompl::geometric::ContactVFRRT::Motion *ompl::geometric::ContactVFRRT::extendTree(Motion *nmotion, base::State *rstate,
                                                                     const Eigen::VectorXd &v)
{
    OMPL_INFORM("========== extendTree");
    // base::State *newState = si_->allocState();
    base::State *newState = rstate;
    base::State *xstate = si_->allocState();
    // si_->copyState(newState, nmotion->state);

    // double d = si_->distance(nmotion->state, rstate);
    // if (d > maxDistance_)
    // {
    //     d = maxDistance_;
    // }
    // OMPL_INFORM("d: %f", d);

    const base::StateSpacePtr &space = si_->getStateSpace();
    for (std::size_t i = 0; i < vfdim_; i++)
    {
        // *space->getValueAddressAtIndex(newState, i) += d * v[i];
        *space->getValueAddressAtIndex(newState, i) = v[i];
        // *space->getValueAddressAtIndex(newState, i) = *space->getValueAddressAtIndex(rstate, i);

        // std::cout << i << ": " << *space->getValueAddressAtIndex(newState, i) << ", ";
    }
    // std::cout << std::endl;

    double d = si_->distance(nmotion->state, newState);
    maxDistance_ = 0.0174533;
    OMPL_INFORM("d: %f", d);
    OMPL_INFORM("maxDistance_: %f", maxDistance_);
    OMPL_INFORM("maxDistance_ / d: %f", maxDistance_ / d);

    if (d > maxDistance_)
    {
        si_->getStateSpace()->interpolate(nmotion->state, newState, maxDistance_ / d, xstate);
        newState = xstate;
    }

    if (!v.hasNaN() && si_->checkMotion(nmotion->state, newState))
    {
        auto *motion = new Motion(si_);
        // motion->state = newState;
        si_->copyState(motion->state, newState);
        motion->parent = nmotion;
        // updateExplorationEfficiency(motion);
        nn_->add(motion);

        return motion;
    }
    else
    {
        si_->freeState(newState);
        si_->freeState(xstate);
        // inefficientCount_++;
        return nullptr;
    }
}

ompl::geometric::ContactVFRRT::Motion *ompl::geometric::ContactVFRRT::extendTreeRRT(Motion *nmotion, base::State *rstate)
{
    OMPL_INFORM("========== extendTreeRRT");
    base::State *xstate = si_->allocState();
    base::State *dstate = rstate;

    /* find state to add */
    double d = si_->distance(nmotion->state, rstate);
    OMPL_INFORM("d: %f", d);
    OMPL_INFORM("maxDistance_: %f", maxDistance_);
    OMPL_INFORM("maxDistance_ / d: %f", maxDistance_ / d);

    if (d > maxDistance_)
    {
        si_->getStateSpace()->interpolate(nmotion->state, rstate, maxDistance_ / d, xstate);
        dstate = xstate;
    }

    if (si_->checkMotion(nmotion->state, dstate))
    {
        auto *motion = new Motion(si_);
        si_->copyState(motion->state, dstate);
        motion->parent = nmotion;
        nn_->add(motion);
        nmotion = motion;
        return motion;
    }
    else
    {
        si_->freeState(xstate);
        return nullptr;
    }
}

void ompl::geometric::ContactVFRRT::updateExplorationEfficiency(Motion *m)
{
    Motion *near = nn_->nearest(m);
    if (distanceFunction(m, near) < si_->getStateValidityCheckingResolution())
        inefficientCount_++;
    else
        efficientCount_++;

    explorationInefficiency_ = inefficientCount_ / (double)(efficientCount_ + inefficientCount_);
}

ompl::base::PlannerStatus ompl::geometric::ContactVFRRT::solve(const base::PlannerTerminationCondition &ptc)
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

    while (!ptc)
    {
        // Sample random state (with goal biasing)
        if ((goal_s != nullptr) && rng_.uniform01() < goalBias_ && goal_s->canSample())
            goal_s->sampleGoal(rstate);
        else
            sampler_->sampleUniform(rstate);

        // Find closest state in the tree
        Motion *nmotion = nn_->nearest(rmotion);

        // Modify direction based on vector field before extending
        Motion *motion = extendTree(nmotion, rstate, getNewDirection(nmotion->state, rstate));

        // Motion *motion = extendTreeRRT(nmotion, rstate);

        if (!motion)
            continue;

        // Check if we can connect to the goal
        double dist = 0.0;
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

    if (rmotion->state)
        si_->freeState(rmotion->state);
    delete rmotion;

    OMPL_INFORM("%s: Created %u states", getName().c_str(), nn_->size());

    return {solved, approximate};
}
