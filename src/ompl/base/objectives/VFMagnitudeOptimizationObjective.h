
#ifndef OMPL_BASE_OBJECTIVES_VF_MAGNITUDE_OPTIMIZATION_OBJECTIVE_
#define OMPL_BASE_OBJECTIVES_VF_MAGNITUDE_OPTIMIZATION_OBJECTIVE_

#include <utility>

#include "ompl/base/OptimizationObjective.h"
#include "ompl/geometric/planners/rrt/VFRRT.h"

namespace ompl
{
    namespace base
    {
        /**
         * Optimization objective that computes the cost of a robot state based on the added magnitude of the vectors
         * acting on it.
         */
        class VFMagnitudeOptimizationObjective : public ompl::base::OptimizationObjective
        {
        public:
            /** Constructor. */
            VFMagnitudeOptimizationObjective(const ompl::base::SpaceInformationPtr &si,
                                             geometric::VFRRT::VectorField vf)
              : ompl::base::OptimizationObjective(si), vf_(std::move(vf))
            {
                description_ = "Vector Field Magnitude";
            }

            /** Assume we can always do better. */
            bool isSatisfied(ompl::base::Cost) const override
            {
                return false;
            }

            /** \brief Returns a cost with a value of 0. */
            Cost stateCost(const State *) const override
            {
                return Cost(0.);
            }

            bool isSymmetric() const override
            {
                return false;
            }

            double sum(const Eigen::VectorXd &vec) const
            {
                double total = 0.0;
                for (std::size_t i = 0; i < vec.size(); i++)
                {
                    total += vec[i];
                }
                return total;
            }

            /** Compute upstream criterion between two states. */
            ompl::base::Cost motionCost(const State *s1, const State *s2) const override
            {
                const base::StateSpacePtr &space = si_->getStateSpace();
                std::size_t vfdim = space->getValueLocations().size();
                std::size_t numSegments = space->validSegmentCount(s1, s2);
                std::vector<ompl::base::State *> interp;
                si_->getMotionStates(s1, s2, interp, numSegments - 1, true, true);
                double cost = 0;
                for (std::size_t i = 0; i < interp.size() - 1; i++)
                {
                    Eigen::VectorXd f = vf_(interp[i]);
                    cost += sum(f);
                    si_->freeState(interp[i]);
                }
                // std::cout << "cost: " << cost << std::endl;
                si_->freeState(interp[interp.size() - 1]);
                return ompl::base::Cost(cost);
            }

        protected:
            /** VectorField associated with the space. */
            geometric::VFRRT::VectorField vf_;
        };
    }  // namespace base
}  // namespace ompl

#endif