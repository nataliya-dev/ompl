/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2010, Rice University
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
 *   * Neither the name of the Rice University nor the names of its
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

/* Author: Ioan Sucan */

#include "ompl/base/DiscreteMotionValidator.h"
#include "ompl/util/Exception.h"
#include <queue>

#include <kdl/chainiksolverpos_nr_jl.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/frames.hpp>
#include <kdl/jntarray.hpp>
#include <kdl_parser/kdl_parser.hpp>
#include <kdl/chainiksolverpos_nr_jl.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/frames.hpp>
#include <kdl/jntarray.hpp>
#include <kdl_parser/kdl_parser.hpp>

#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
py::scoped_interpreter python{};
auto trajectoryClassifier = py::module::import("trajectory_model").attr("spilled");
auto quatToEuler = py::module::import("trajectory_model").attr("quat_to_euler");


void ompl::base::DiscreteMotionValidator::defaultSettings()
{
    stateSpace_ = si_->getStateSpace().get();
    if (stateSpace_ == nullptr)
        throw Exception("No state space for motion validator");
}

bool ompl::base::DiscreteMotionValidator::checkMotion(const State *s1, const State *s2,
                                                      std::pair<State *, double> &lastValid) const
{
    /* assume motion starts in a valid configuration so s1 is valid */

    bool result = true;
    int nd = stateSpace_->validSegmentCount(s1, s2);

    if (nd > 1)
    {
        /* temporary storage for the checked state */
        State *test = si_->allocState();

        for (int j = 1; j < nd; ++j)
        {
            stateSpace_->interpolate(s1, s2, (double)j / (double)nd, test);
            if (!si_->isValid(test))
            {
                lastValid.second = (double)(j - 1) / (double)nd;
                if (lastValid.first != nullptr)
                    stateSpace_->interpolate(s1, s2, lastValid.second, lastValid.first);
                result = false;
                break;
            }
        }
        si_->freeState(test);
    }

    if (result)
        if (!si_->isValid(s2))
        {
            lastValid.second = (double)(nd - 1) / (double)nd;
            if (lastValid.first != nullptr)
                stateSpace_->interpolate(s1, s2, lastValid.second, lastValid.first);
            result = false;
        }

    if (result)
        valid_++;
    else
        invalid_++;

    return result;
}

std::array<double, 7> fk(std::array<double, 7> q) {
  /*
   * Calculate End Effector Pose from Joint Angles using KDL
   * @param q: joint angles
   * @param urdf_path: path to urdf file
   * @return: end effector pose
   */

  KDL::Tree my_tree;
  kdl_parser::treeFromFile("/home/ava/npm/panda.urdf", my_tree);
  KDL::Chain chain;

  my_tree.getChain("panda_link0", "panda_link8", chain);

  KDL::JntArray jointpositions = KDL::JntArray(chain.getNrOfJoints());

  for (int i = 0; i < 7; i++) {
    jointpositions(i) = q[i];
  }

  KDL::Frame cartpos;

  KDL::ChainFkSolverPos_recursive fksolver = KDL::ChainFkSolverPos_recursive(chain);

  bool kinematics_status;
  kinematics_status = fksolver.JntToCart(jointpositions, cartpos);

  if (kinematics_status >= 0) {
    std::array<double, 7> ee_pose;
    for (int i = 0; i < 3; i++) {
      ee_pose[i] = cartpos.p[i];
    }
    for (int i = 0; i < 4; i++) {
      ee_pose[i + 3] = cartpos.M.data[i];
    }
    return ee_pose;
  } else {
    throw ompl::Exception("Could not calculate forward kinematics");
  }
}

bool ompl::base::DiscreteMotionValidator::checkMotion(const State *s1, const State *s2) const
{
    // a, b, c, w = {0.725851,-0.687294,-0.0277177,-0.687799}
    // euler = [ x: -116.8784625, y: 37.9156883, z: 63.044746 ]

    std::vector<double> joint_values;
    si_->getStateSpace()->copyToReals(joint_values, s2);
    std::array<double, 7> jv;
    std::copy_n(joint_values.begin(), 7, jv.begin());
    std::array<double, 7> cartesian_pos = fk(jv);

    std::array<double, 4> quaternion = {cartesian_pos[3], cartesian_pos[4], cartesian_pos[5], cartesian_pos[6]};
    py::list quat_py = py::cast(quaternion);
    auto resultobj = quatToEuler(quat_py);
    std::vector<double> euler_angles = resultobj.cast<std::vector<double>>();

    std::array<double, 3> stable_angles = {-118.682, 42.0002, -61.2353};
    // the absolute value of the difference between the current and stable angles 
    // should be less than 5 degrees
    for (int i = 0; i < 3; i++) {
        if (abs(euler_angles[i] - stable_angles[i]) > 45) {
            // std::cout<<"Invalid euler angles"<<std::endl;
            return false;
        }
    }



    /* assume motion starts in a valid configuration so s1 is valid */
    if (!si_->isValid(s2))
    {
        invalid_++;
        return false;
    }

    bool result = true;
    int nd = stateSpace_->validSegmentCount(s1, s2);

    /* initialize the queue of test positions */
    std::queue<std::pair<int, int>> pos;
    if (nd >= 2)
    {
        pos.emplace(1, nd - 1);

        /* temporary storage for the checked state */
        State *test = si_->allocState();

        /* repeatedly subdivide the path segment in the middle (and check the middle) */
        while (!pos.empty())
        {
            std::pair<int, int> x = pos.front();

            int mid = (x.first + x.second) / 2;
            stateSpace_->interpolate(s1, s2, (double)mid / (double)nd, test);

            if (!si_->isValid(test))
            {
                result = false;
                break;
            }

            pos.pop();

            if (x.first < mid)
                pos.emplace(x.first, mid - 1);
            if (x.second > mid)
                pos.emplace(mid + 1, x.second);
        }

        si_->freeState(test);
    }

    if (result)
        valid_++;
    else
        invalid_++;

    return result;
}



bool ompl::base::DiscreteMotionValidator::checkTrajectorySoFar(std::vector<base::State *> trajectory_so_far) const
{
    OMPL_INFORM("Check Trajectory So Far");
    // return true;
    if (trajectory_so_far.size() <= 1) {
        return true;
    }

    std::list<Eigen::VectorXd> points;
    for (int i=0; i<trajectory_so_far.size(); i++){
        base::State *state = trajectory_so_far[i];
        std::vector<double> joint_values;
        si_->getStateSpace()->copyToReals(joint_values, state);
        std::array<double, 7> jv;
        std::copy_n(joint_values.begin(), 7, jv.begin());
        Eigen::VectorXd new_point(7);
        new_point << jv[0], jv[1], jv[2], jv[3], jv[4] ,jv[5] ,jv[6];
        points.push_back(new_point);
    }

    Eigen::VectorXd goal_point(7);
    goal_point << 1.67791, 0.949664, -1.31687, -1.48235, 0.907626, 1.77346, 1.2502;
    points.push_back(goal_point);


    double resample_dt_ = 0.001;
    double path_tolerance = 0.09;
    Eigen::VectorXd max_velocity(7);
    Eigen::VectorXd max_acceleration(7);
    max_velocity << 1.0875, 1.0875, 1.0875, 1.0875, 1.305, 1.305, 1.305;
    max_acceleration << 1.875, 0.9375, 1.25, 1.5625, 1.875, 2.5, 2.5; // TODO clean
    ompl::trajectory_processing::Trajectory parameterized(ompl::trajectory_processing::Path(points, path_tolerance), max_velocity, max_acceleration, 0.001);
    size_t sample_count = std::ceil(parameterized.getDuration() / resample_dt_);

    std::vector<std::array<double, 7>> trajectory_in_cartesian;
    for (size_t sample = 0; sample <= sample_count; ++sample)
    {
        double t = std::min(parameterized.getDuration(), sample * resample_dt_);
        Eigen::VectorXd position = parameterized.getPosition(t);
        std::array<double, 7> jv = {position[0], position[1], position[2], position[3], position[4], position[5], position[6]};
        std::array<double, 7> cartesian_pos = fk(jv);
        trajectory_in_cartesian.push_back(cartesian_pos);
    }


    py::list py_traj_cartesian = py::cast(trajectory_in_cartesian);
    auto resultobj = trajectoryClassifier(py_traj_cartesian);
    double result = resultobj.cast<double>();
    std::cout<<"Probability of spilling is: "<<result<<std::endl;
    return result <= 0.5;
    // return true;
}
