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

#include "ompl/geometric/planners/rrt/CVFRRT.h"
#include "ompl/base/goals/GoalSampleableRegion.h"

namespace ompl
{
    namespace magic
    {
        /// Number of sampler to determine mean vector field norm in \ref gCVFRRT
        static const unsigned int CVFRRT_MEAN_NORM_SAMPLES = 1000;
    }  // namespace magic
}  // namespace ompl

ompl::geometric::CVFRRT::CVFRRT(const base::SpaceInformationPtr &si, VectorField vf, double exploration,
                                double initial_lambda, unsigned int update_freq)
  : RRT(si), vf_(std::move(vf)), explorationSetting_(exploration), lambda_(initial_lambda), nth_step_(update_freq)
{
    setName("CVFRRT");
    maxDistance_ = si->getStateValidityCheckingResolution();
}

ompl::geometric::CVFRRT::~CVFRRT() = default;

void ompl::geometric::CVFRRT::clear()
{
    RRT::clear();
    efficientCount_ = 0;
    inefficientCount_ = 0;
    explorationInefficiency_ = 0.;
    step_ = 0;
}

void ompl::geometric::CVFRRT::setup()
{
    RRT::setup();
    vfdim_ = si_->getStateSpace()->getValueLocations().size();
}

double ompl::geometric::CVFRRT::determineMeanNorm()
{
    ompl::base::State *rstate = si_->allocState();
    double sum = 0.;
    for (unsigned int i = 0; i < magic::CVFRRT_MEAN_NORM_SAMPLES; i++)
    {
        sampler_->sampleUniform(rstate);
        sum += vf_(rstate).norm();
    }
    si_->freeState(rstate);
    return sum / magic::CVFRRT_MEAN_NORM_SAMPLES;
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
    vnew = vrand + vfield;

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

double ompl::geometric::CVFRRT::biasedSampling(const Eigen::VectorXd &vrand, const Eigen::VectorXd &vfield,
                                               double lambdaScale)
{
    OMPL_INFORM("========== biasedSampling");
    double sigma = .25 * (vrand - vfield).squaredNorm();
    OMPL_INFORM("sigma: %f", sigma);
    updateGain();
    OMPL_INFORM("lambdaScale: %f", lambdaScale);
    OMPL_INFORM("meanNorm_: %f", meanNorm_);
    OMPL_INFORM("lambda_: %f", lambda_);
    double scaledLambda = lambda_ * lambdaScale / meanNorm_;
    OMPL_INFORM("scaledLambda: %f", scaledLambda);

    double phi = scaledLambda / (1. - std::exp(-2. * scaledLambda));
    OMPL_INFORM("phi: %f", phi);
    double z = -std::log(1. - sigma * scaledLambda / phi) / scaledLambda;
    OMPL_INFORM("z: %f", z);
    return std::sqrt(2. * z);
}

void ompl::geometric::CVFRRT::updateGain()
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

Eigen::VectorXd ompl::geometric::CVFRRT::computeAlphaBeta(double omega, const Eigen::VectorXd &vrand,
                                                          const Eigen::VectorXd &vfield)
{
    OMPL_INFORM("========== computeAlphaBeta");
    double w2 = omega * omega;
    OMPL_INFORM("w2: %f", w2);

    double c = vfield.dot(vrand);
    OMPL_INFORM("c: %f", c);

    double cc_1 = c * c - 1.;
    OMPL_INFORM("cc_1: %f", cc_1);

    double root = std::sqrt(cc_1 * w2 * (w2 - 4.));
    OMPL_INFORM("root: %f", root);

    double beta = -root / (2. * cc_1);
    OMPL_INFORM("beta: %f", beta);

    double sign = (beta < 0.) ? -1. : 1.;
    OMPL_INFORM("sign: %f", sign);

    beta *= sign;
    OMPL_INFORM("beta: %f", beta);

    double alpha = (sign * c * root + cc_1 * (2. - w2)) / (2. * cc_1);
    OMPL_INFORM("alpha: %f", alpha);

    OMPL_INFORM("VRAND");
    for (std::size_t i = 0; i < 7; i++)
    {
        std::cout << vrand[i] << ", ";
    }
    std::cout << std::endl;
    OMPL_INFORM("VFIELD");
    for (std::size_t i = 0; i < 7; i++)
    {
        std::cout << vfield[i] << ", ";
    }
    std::cout << std::endl;

    Eigen::VectorXd vnew = alpha * vfield + beta * vrand;
    OMPL_INFORM("VNEW");
    for (std::size_t i = 0; i < 7; i++)
    {
        std::cout << vnew[i] << ", ";
    }
    std::cout << std::endl;
    return alpha * vfield + beta * vrand;
}

ompl::geometric::CVFRRT::Motion *ompl::geometric::CVFRRT::extendTree(Motion *m, base::State *rstate,
                                                                     const Eigen::VectorXd &v)
{
    OMPL_INFORM("========== extendTree");
    base::State *newState = si_->allocState();
    si_->copyState(newState, m->state);

    double d = si_->distance(m->state, rstate);
    if (d > maxDistance_)
        d = maxDistance_;

    OMPL_INFORM("d: %f", d);

    const base::StateSpacePtr &space = si_->getStateSpace();
    for (unsigned int i = 0; i < vfdim_; i++)
        *space->getValueAddressAtIndex(newState, i) += d * v[i];
    if (!v.hasNaN() && si_->checkMotion(m->state, newState))
    {
        auto *motion = new Motion();
        motion->state = newState;
        motion->parent = m;
        updateExplorationEfficiency(motion);
        nn_->add(motion);
        return motion;
    }
    else
    {
        si_->freeState(newState);
        inefficientCount_++;
        return nullptr;
    }
}

void ompl::geometric::CVFRRT::updateExplorationEfficiency(Motion *m)
{
    Motion *near = nn_->nearest(m);
    if (distanceFunction(m, near) < si_->getStateValidityCheckingResolution())
        inefficientCount_++;
    else
        efficientCount_++;

    explorationInefficiency_ = inefficientCount_ / (double)(efficientCount_ + inefficientCount_);
}

ompl::base::PlannerStatus ompl::geometric::CVFRRT::solve(const base::PlannerTerminationCondition &ptc)
{
    checkValidity();
    base::Goal *goal = pdef_->getGoal().get();
    auto *goal_s = dynamic_cast<base::GoalSampleableRegion *>(goal);

    if (!sampler_)
        sampler_ = si_->allocStateSampler();

    meanNorm_ = determineMeanNorm();
    OMPL_INFORM("meanNorm_: %f", meanNorm_);

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
