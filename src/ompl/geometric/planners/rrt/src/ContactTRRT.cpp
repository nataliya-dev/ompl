/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2008, Willow Garage, Inc.
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

/* Author: Dave Coleman, Ryan Luna */

#include "ompl/geometric/planners/rrt/ContactTRRT.h"
#include "ompl/base/objectives/MechanicalWorkOptimizationObjective.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/tools/config/SelfConfig.h"
#include "ompl/tools/config/MagicConstants.h"
#include <limits>

ompl::geometric::ContactTRRT::ContactTRRT(const base::SpaceInformationPtr &si, VectorField vf)
  : base::Planner(si, "ContactTRRT"), vf_(std::move(vf))
{
    // Standard RRT Variables
    specs_.approximateSolutions = true;
    specs_.directed = true;

    Planner::declareParam<double>("range", this, &ContactTRRT::setRange, &ContactTRRT::getRange, "0.:1.:10000.");
    Planner::declareParam<double>("goal_bias", this, &ContactTRRT::setGoalBias, &ContactTRRT::getGoalBias, "0.:.05:1.");
    // Planner::declareParam<double>("temp_change_factor", this, &ContactTRRT::setTempChangeFactor,
    //                               &ContactTRRT::getTempChangeFactor, "0.:.1:1.");
    // Planner::declareParam<double>("init_temperature", this, &ContactTRRT::setInitTemperature,
    //                               &ContactTRRT::getInitTemperature);
    Planner::declareParam<double>("frontier_threshold", this, &ContactTRRT::setFrontierThreshold,
                                  &ContactTRRT::getFrontierThreshold);
    Planner::declareParam<double>("frontier_node_ratio", this, &ContactTRRT::setFrontierNodeRatio,
                                  &ContactTRRT::getFrontierNodeRatio);
    // Planner::declareParam<double>("cost_threshold", this, &ContactTRRT::setCostThreshold,
    //                               &ContactTRRT::getCostThreshold);
}

ompl::geometric::ContactTRRT::~ContactTRRT()
{
    freeMemory();
}

void ompl::geometric::ContactTRRT::clear()
{
    Planner::clear();
    sampler_.reset();
    freeMemory();
    if (nearestNeighbors_)
        nearestNeighbors_->clear();
    lastGoalMotion_ = nullptr;

    // Clear ContactTRRT specific variables ---------------------------------------------------------
    nonfrontierCount_ = 1;
    frontierCount_ = 1;  // init to 1 to prevent division by zero error
    if (opt_)
        bestCost_ = worstCost_ = opt_->identityCost();
}

void ompl::geometric::ContactTRRT::setup()
{
    Planner::setup();
    tools::SelfConfig selfConfig(si_, getName());

    // ContactTRRT Specific Variables
    frontierNodeRatio_ = 0.1;  // 1/10, or 1 nonfrontier for every 10 frontier
    maxDistance_ = 0.5;
    frontierThreshold_ = 0.1;

    vfdim_ = si_->getStateSpace()->getValueLocations().size();

    if (!pdef_ || !pdef_->hasOptimizationObjective())
    {
        OMPL_INFORM("%s: No optimization objective specified.  Defaulting to mechanical work minimization.",
                    getName().c_str());
        opt_ = std::make_shared<base::MechanicalWorkOptimizationObjective>(si_);
    }
    else
        opt_ = pdef_->getOptimizationObjective();

    // Set maximum distance a new node can be from its nearest neighbor
    if (maxDistance_ < std::numeric_limits<double>::epsilon())
    {
        selfConfig.configurePlannerRange(maxDistance_);
        maxDistance_ *= magic::COST_MAX_MOTION_LENGTH_AS_SPACE_EXTENT_FRACTION;
    }

    // Set the threshold that decides if a new node is a frontier node or non-frontier node
    if (frontierThreshold_ < std::numeric_limits<double>::epsilon())
    {
        frontierThreshold_ = si_->getMaximumExtent() * 0.01;
        OMPL_DEBUG("%s: Frontier threshold detected to be %lf", getName().c_str(), frontierThreshold_);
    }

    // Create the nearest neighbor function the first time setup is run
    if (!nearestNeighbors_)
        nearestNeighbors_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Motion *>(this));

    // Set the distance function
    nearestNeighbors_->setDistanceFunction([this](const Motion *a, const Motion *b) { return distanceFunction(a, b); });

    // Setup ContactTRRT specific variables ---------------------------------------------------------
    nonfrontierCount_ = 1;
    frontierCount_ = 1;  // init to 1 to prevent division by zero error
    bestCost_ = worstCost_ = opt_->identityCost();
    initDataFile();
}

void ompl::geometric::ContactTRRT::freeMemory()
{
    // Delete all motions, states and the nearest neighbors data structure
    if (nearestNeighbors_)
    {
        std::vector<Motion *> motions;
        nearestNeighbors_->list(motions);
        for (auto &motion : motions)
        {
            if (motion->state)
                si_->freeState(motion->state);
            delete motion;
        }
    }
}

ompl::base::PlannerStatus
ompl::geometric::ContactTRRT::solve(const base::PlannerTerminationCondition &plannerTerminationCondition)
{
    // Basic error checking
    checkValidity();

    // Goal information
    base::Goal *goal = pdef_->getGoal().get();
    auto *goalRegion = dynamic_cast<base::GoalSampleableRegion *>(goal);

    // Input States ---------------------------------------------------------------------------------

    // Loop through valid input states and add to tree
    while (const base::State *state = pis_.nextStart())
    {
        // Allocate memory for a new start state motion based on the "space-information"-size
        auto *motion = new Motion(si_);

        // Copy destination <= source
        si_->copyState(motion->state, state);

        // Set cost for this start state
        motion->cost = opt_->stateCost(motion->state);

        if (nearestNeighbors_->size() == 0)  // do not overwrite best/worst from previous call to solve
            worstCost_ = bestCost_ = motion->cost;

        // Add start motion to the tree
        nearestNeighbors_->add(motion);
    }

    // Check that input states exist
    if (nearestNeighbors_->size() == 0)
    {
        OMPL_ERROR("%s: There are no valid initial states!", getName().c_str());
        return base::PlannerStatus::INVALID_START;
    }

    // Create state sampler if this is ContactTRRT's first run
    if (!sampler_)
        sampler_ = si_->allocStateSampler();

    // Debug
    OMPL_INFORM("%s: Starting planning with %u states already in datastructure", getName().c_str(),
                nearestNeighbors_->size());

    // Solver variables ------------------------------------------------------------------------------------

    // the final solution
    Motion *solution = nullptr;
    // the approximate solution, returned if no final solution found
    Motion *approxSolution = nullptr;
    // track the distance from goal to closest solution yet found
    double approxDifference = std::numeric_limits<double>::infinity();

    // distance between states - the intial state and the interpolated state (may be the same)
    double randMotionDistance;

    // Create random motion and a pointer (for optimization) to its state
    auto *randMotion = new Motion(si_);
    Motion *nearMotion;

    // STATES
    // The random state
    base::State *randState = randMotion->state;
    // The new state that is generated between states *to* and *from*
    base::State *interpolatedState = si_->allocState();  // Allocates "space information"-sized memory for a state
    // The chosen state btw rand_state and interpolated_state
    base::State *newState;

    std::size_t num_goal_samples = 0;
    std::size_t num_uniform_samples = 0;

    // Begin sampling --------------------------------------------------------------------------------------
    while (plannerTerminationCondition() == false)
    {
        OMPL_INFORM("=== New state ===");
        // I.

        // Sample random state (with goal biasing probability)
        if (goalRegion && rng_.uniform01() < goalBias_ && goalRegion->canSample())
        {
            // Bias sample towards goal
            goalRegion->sampleGoal(randState);
            num_goal_samples++;
        }
        else
        {
            // Uniformly Sample
            sampler_->sampleUniform(randState);
            num_uniform_samples++;
        }

        OMPL_INFORM("num_goal_samples %u: ", num_goal_samples);
        OMPL_INFORM("num_uniform_samples %u: ", num_uniform_samples);

        // II.

        // Find closest state in the tree
        nearMotion = nearestNeighbors_->nearest(randMotion);

        // III.

        // Distance from near state q_n to a random state
        randMotionDistance = si_->distance(nearMotion->state, randState);

        // OMPL_INFORM("randMotionDistance: %f", randMotionDistance);
        // OMPL_INFORM("maxDistance_: %f", maxDistance_);

        // Check if the rand_state is too far away
        if (randMotionDistance > maxDistance_)
        {
            // OMPL_INFORM("maxDistance_ / randMotionDistance: %f", maxDistance_ / randMotionDistance);

            // Computes the state that lies at time t in [0, 1] on the segment that connects *from* state to *to* state.
            // The memory location of *state* is not required to be different from the memory of either *from* or *to*.
            si_->getStateSpace()->interpolate(nearMotion->state, randState, maxDistance_ / randMotionDistance,
                                              interpolatedState);

            // Update the distance between near and new with the interpolated_state
            randMotionDistance = si_->distance(nearMotion->state, interpolatedState);
            // OMPL_INFORM("randMotionDistance: %f", randMotionDistance);

            // Use the interpolated state as the new state
            newState = interpolatedState;
        }
        else  // Random state is close enough
            newState = randState;

        // OMPL_INFORM("randMotionDistance: %f", randMotionDistance);
        sampleNum_++;

        // IV.
        // this stage integrates collision detections in the presence of obstacles and checks for collisions
        if (!si_->checkMotion(nearMotion->state, newState))
        {
            OMPL_INFORM("Collision state");
            continue;  // try a new sample
        }

        // base::Cost childCost = opt_->stateCost(newState);
        // base::Cost childCost = opt_->motionCost(nearMotion->state, newState);
        // OMPL_INFORM("childCost.value(): %f", childCost.value());

        goal->isSatisfied(newState, &distToGoal_);

        bool is_valid = perLinkTransitionTestWedighted(nearMotion, newState);
        // bool is_valid = perBranchTransitionTest(nearMotion, randMotionDistance, childCost);
        if (!is_valid)
        {
            saveData();
            temp_ = nearMotion->temp_;
            setMaxTemp(temp_);
            continue;
        }
        // Minimum Expansion Control
        // A possible side effect may appear when the tree expansion toward unexplored regions remains slow, and the
        // new nodes contribute only to refine already explored regions.
        // if (!minExpansionControl(randMotionDistance))
        //     continue;  // give up on this one and try a new sample

        // V.

        // Create a motion
        auto *motion = new Motion(si_);
        si_->copyState(motion->state, newState);
        motion->parent = nearMotion;  // link q_new to q_near as an edge
        // motion->cost = childCost;

        motion->vthresh = nearMotion->vthresh;
        motion->vnumfail = nearMotion->vnumfail;
        motion->vtemp = nearMotion->vtemp;

        motion->nFail_ = nearMotion->nFail_;
        motion->temp_ = nearMotion->temp_;

        // Add motion to data structure
        nearestNeighbors_->add(motion);

        // OMPL_INFORM("bestCost_.value(): %f", bestCost_.value());
        // OMPL_INFORM("worstCost_.value(): %f", worstCost_.value());
        // OMPL_INFORM("motion->cost.value(): %f", motion->cost.value());

        if (opt_->isCostBetterThan(motion->cost, bestCost_))  // motion->cost is better than the existing best
            bestCost_ = motion->cost;
        if (opt_->isCostBetterThan(worstCost_, motion->cost))  // motion->cost is worse than the existing worst
            worstCost_ = motion->cost;

        OMPL_INFORM("+++ State added +++");

        // VI.

        // Check if this motion is the goal
        double distToGoal = 0.0;
        bool isSatisfied = goal->isSatisfied(motion->state, &distToGoal);
        OMPL_INFORM("distToGoal: %f", distToGoal);

        setMinDistToGoal(distToGoal);
        setMaxTemp(temp_);
        saveData();
        if (isSatisfied)
        {
            approxDifference = distToGoal;  // the tolerated error distance btw state and goal
            solution = motion;              // set the final solution
            break;
        }

        // Is this the closest solution we've found so far
        if (distToGoal < approxDifference)
        {
            approxDifference = distToGoal;
            approxSolution = motion;
        }

    }  // end of solver sampling loop

    // Finish solution processing --------------------------------------------------------------------

    bool solved = false;
    bool approximate = false;

    // Substitute an empty solution with the best approximation
    if (solution == nullptr)
    {
        solution = approxSolution;
        approximate = true;
    }

    // Generate solution path for real/approx solution
    if (solution != nullptr)
    {
        lastGoalMotion_ = solution;

        // construct the solution path
        std::vector<Motion *> mpath;
        while (solution != nullptr)
        {
            mpath.push_back(solution);
            solution = solution->parent;
        }

        // set the solution path
        auto path(std::make_shared<PathGeometric>(si_));
        for (int i = mpath.size() - 1; i >= 0; --i)
            path->append(mpath[i]->state);

        pdef_->addSolutionPath(path, approximate, approxDifference, getName());
        solved = true;
    }

    // Clean up ---------------------------------------------------------------------------------------

    si_->freeState(interpolatedState);
    if (randMotion->state)
        si_->freeState(randMotion->state);
    delete randMotion;

    OMPL_INFORM("%s: Created %u states", getName().c_str(), nearestNeighbors_->size());

    return {solved, approximate};
}

void ompl::geometric::ContactTRRT::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);

    std::vector<Motion *> motions;
    if (nearestNeighbors_)
        nearestNeighbors_->list(motions);

    if (lastGoalMotion_)
        data.addGoalVertex(base::PlannerDataVertex(lastGoalMotion_->state));

    for (auto &motion : motions)
    {
        if (motion->parent == nullptr)
            data.addStartVertex(base::PlannerDataVertex(motion->state));
        else
            data.addEdge(base::PlannerDataVertex(motion->parent->state), base::PlannerDataVertex(motion->state));
    }
}

bool ompl::geometric::ContactTRRT::perBranchTransitionTest(Motion *parentMotion, double dist,
                                                           const base::Cost &childCost)
{
    // const base::State *parentState = parentMotion->state;
    base::Cost parentCost = parentMotion->cost;
    OMPL_INFORM("childCost.value(): %f", childCost.value());
    OMPL_INFORM("parentCost.value(): %f", parentCost.value());
    // OMPL_INFORM("costThreshold_.value(): %f", costThreshold_.value());

    // Disallow any cost that is not better than the cost threshold
    // if (!opt_->isCostBetterThan(childCost, costThreshold_))
    // {
    //     OMPL_INFORM("Cost is above threshold.");
    //     return false;
    // }

    // always allow cost that is better than parent
    if (opt_->isCostBetterThan(childCost, parentCost))
    {
        OMPL_INFORM("Motion cost is better than parent cost");
        return true;
    }

    // Always accept if the cost is near or below zero
    if (childCost.value() < 1e-4)
    {
        OMPL_INFORM("Cost is near or below zero.");
        return true;
    }

    double dCost = childCost.value() - parentCost.value();
    OMPL_INFORM("dCost: %f", dCost);
    OMPL_INFORM("dist: %f", dist);
    double dCostDist = dCost / dist;
    OMPL_INFORM("dCostDist: %f", dCostDist);

    OMPL_INFORM("temp_: %f", parentMotion->temp_);
    // OMPL_INFORM("K_: %f", K_);

    double tranProb = exp(-1.0 * dCostDist / (parentMotion->temp_ * parentMotion->K_));
    OMPL_INFORM("tranProb: %f", tranProb);

    double randProb = (double)rand() / RAND_MAX;
    OMPL_INFORM("randProb: %f", randProb);
    if (tranProb > randProb)
    {
        parentMotion->temp_ /= parentMotion->tempChangeFactor_;
        parentMotion->nFail_ = 0;
        OMPL_INFORM("tranProb > randProb: temp_: %f", parentMotion->temp_);
        return true;
    }
    else
    {
        if (parentMotion->nFail_ > parentMotion->nFailMax_)
        {
            parentMotion->temp_ *= parentMotion->tempChangeFactor_;
            OMPL_INFORM("nFail > nFailMax_: temp_: %f", parentMotion->temp_);
            parentMotion->nFail_ = 0;
        }
        else
        {
            parentMotion->nFail_++;
            OMPL_INFORM("nFail_: %d", parentMotion->nFail_);
        }

        return false;
    }
}

bool ompl::geometric::ContactTRRT::perLinkTransitionTestWedighted(Motion *parentMotion, base::State *newState)
{
    base::State *nearState = parentMotion->state;
    Eigen::VectorXd vfield = vf_(nearState);
    const double lambdaScale = vfield.norm();
    OMPL_INFORM("lambdaScale: %f", lambdaScale);
    if (lambdaScale < std::numeric_limits<float>::epsilon())
    {
        return true;
    }

    Eigen::VectorXd vnew(vfdim_);
    Eigen::VectorXd vnear(vfdim_);

    for (std::size_t i = 0; i < vfdim_; i++)
    {
        vnew[i] = *si_->getStateSpace()->getValueAddressAtIndex(newState, i);
        vnear[i] = *si_->getStateSpace()->getValueAddressAtIndex(nearState, i);
    }

    Eigen::VectorXd vdiff(vnew - vnear);

    Eigen::VectorXd &vtemp = parentMotion->vtemp;
    Eigen::VectorXd &vnumfail = parentMotion->vnumfail;

    bool is_valid = true;
    double prev_cost = 0;
    std::size_t link_to_inc = 0;
    double min_temp = vtemp[0];

    for (std::size_t i = 0; i < vfdim_; i++)
    {
        Eigen::VectorXd subfield = vfield.head(i + 1);
        Eigen::VectorXd subnew = vdiff.head(i + 1);
        double cost = subfield.squaredNorm() - subfield.dot(subnew) - prev_cost;
        prev_cost += std::abs(cost);

        double &temp = vtemp[i];
        double &numfail = vnumfail[i];

        double tranProb = exp(-1.0 * std::abs(cost) / (temp * parentMotion->K_));

        OMPL_INFORM("%ld cost: %f", i, cost);
        OMPL_INFORM("%ld prev_cost: %f", i, prev_cost);
        OMPL_INFORM("%ld temp: %f", i, temp);
        OMPL_INFORM("tranProb: %f", tranProb);

        if (tranProb < 0.5)
        {
            is_valid = false;
            if (vtemp[i] < min_temp)
            {
                link_to_inc = i;
            }
        }
        else
        {
            temp /= 1.2;
            numfail = 0;
        }
    }

    if (!is_valid)
    {
        double &temp = vtemp[link_to_inc];
        double &numfail = vnumfail[link_to_inc];
        temp *= 1.9;
        numfail++;
        OMPL_INFORM("State is not valid.");
        return false;
    }
    else
    {
        return true;
    }
}

bool ompl::geometric::ContactTRRT::perLinkTransitionTest(Motion *parentMotion, base::State *newState)
{
    base::State *nearState = parentMotion->state;
    Eigen::VectorXd vfield = vf_(nearState);
    const double lambdaScale = vfield.norm();
    OMPL_INFORM("lambdaScale: %f", lambdaScale);
    if (lambdaScale < std::numeric_limits<float>::epsilon())
    {
        return true;
    }

    // vfield.normalize();

    Eigen::VectorXd vnew(vfdim_);
    Eigen::VectorXd vnear(vfdim_);

    for (std::size_t i = 0; i < vfdim_; i++)
    {
        vnew[i] = *si_->getStateSpace()->getValueAddressAtIndex(newState, i);
        vnear[i] = *si_->getStateSpace()->getValueAddressAtIndex(nearState, i);
    }

    Eigen::VectorXd vdiff(vnew - vnear);

    Eigen::VectorXd &vtemp = parentMotion->vtemp;
    Eigen::VectorXd &vnumfail = parentMotion->vnumfail;

    bool is_valid = true;
    double prev_cost = 0;

    for (std::size_t i = 0; i < vfdim_; i++)
    {
        Eigen::VectorXd subfield = vfield.head(i + 1);
        Eigen::VectorXd subnew = vdiff.head(i + 1);
        double cost = subfield.squaredNorm() - subfield.dot(subnew) - prev_cost;
        prev_cost += std::abs(cost);

        // double cost = vfield[i];

        double &temp = vtemp[i];
        double &numfail = vnumfail[i];

        // OMPL_INFORM("subs size: %ld", subfield.size());
        // OMPL_INFORM("%ld cost: %f", i, cost);
        // OMPL_INFORM("%ld prev_cost: %f", i, prev_cost);
        // OMPL_INFORM("%ld temp: %f", i, temp);
        // OMPL_INFORM("%ld numfail: %f", i, numfail);
        // OMPL_INFORM("%ld distToGoal_: %f", i, distToGoal_);

        double tranProb = exp(-1.0 * std::abs(cost) / (temp * parentMotion->K_));
        // OMPL_INFORM("tranProb: %f", tranProb);

        // double randProb = (double)rand() / RAND_MAX;
        // OMPL_INFORM("randProb: %f", randProb);

        if (tranProb > 0.5)
        {
            temp /= 1.25;
            OMPL_INFORM("tranProb > randProb: %f", temp);

            numfail = 0;
            continue;
        }
        else
        {
            temp *= 1.3;
            // if (numfail > parentMotion->nFailMax_)
            // {
            //     OMPL_INFORM("nFail > nFailMax_: temp: %f", temp);
            //     parentMotion->nFail_ = 0;
            // }
            // else
            // {
            //     numfail++;
            //     OMPL_INFORM("numfail: %d", numfail);
            // }
            is_valid = false;
            break;
        }
    }

    if (!is_valid)
    {
        OMPL_INFORM("State is not valid.");
        return false;
    }
    else
    {
        return true;
    }
}

bool ompl::geometric::ContactTRRT::minExpansionControl(double randMotionDistance)
{
    OMPL_INFORM("frontierThreshold_: %f", frontierThreshold_);

    if (randMotionDistance > frontierThreshold_)
    {
        // participates in the tree expansion
        ++frontierCount_;

        OMPL_INFORM("frontierCount_: %f", (double)frontierCount_);
        OMPL_INFORM("nonfrontierCount_: %f", (double)nonfrontierCount_);

        return true;
    }
    else
    {
        // participates in the tree refinement

        // check our ratio first before accepting it
        if ((double)nonfrontierCount_ / (double)frontierCount_ > frontierNodeRatio_)
        {
            OMPL_INFORM("frontierNodeRatio_: %f", (double)frontierNodeRatio_);
            OMPL_INFORM("frontierCount_: %f", (double)frontierCount_);
            OMPL_INFORM("nonfrontierCount_: %f", (double)nonfrontierCount_);
            // reject this node as being too much refinement
            return false;
        }

        ++nonfrontierCount_;
        OMPL_INFORM("frontierCount_: %f", (double)frontierCount_);
        OMPL_INFORM("nonfrontierCount_: %f", (double)nonfrontierCount_);

        return true;
    }
}

void ompl::geometric::ContactTRRT::setMinDistToGoal(double dist)
{
    if (dist < minDistTGoal_)
    {
        minDistTGoal_ = dist;
    }
}

void ompl::geometric::ContactTRRT::setMaxTemp(double temp)
{
    if (temp > maxTemp_)
    {
        maxTemp_ = temp;
    }
}

void ompl::geometric::ContactTRRT::saveData()
{
    std::fstream file("/home/nataliya/action_ws/src/tacbot/scripts/contactTRRT.csv", std::ios::out | std::ios::app);
    if (file.is_open())
    {
        file << sampleNum_;
        file << ",";
        file << distToGoal_;
        file << ",";
        file << minDistTGoal_;
        file << ",";
        file << worstCost_.value();
        file << ",";
        file << temp_;
        file << ",";
        file << maxTemp_;
        file << "\n";
        file.close();
    }
    else
    {
        OMPL_ERROR("Unable to open file for writing.");
    }
}

void ompl::geometric::ContactTRRT::initDataFile()
{
    std::fstream file("/home/nataliya/action_ws/src/tacbot/scripts/contactTRRT.csv", std::ios::out | std::ios::trunc);
    if (file.is_open())
    {
        file << "sampleNumber";
        file << ",";
        file << "distToGoal";
        file << ",";
        file << "minDistToGoal";
        file << ",";
        file << "worstCost";
        file << ",";
        file << "temp";
        file << ",";
        file << "maxTemp";
        file << "\n";
        file.close();
    }
    else
    {
        OMPL_ERROR("Unable to open file for writing.");
    }
}
