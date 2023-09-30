
#ifndef OMPL_BASE_OBJECTIVES_MINIMIZE_CONTACT_OBJECTIVE_
#define OMPL_BASE_OBJECTIVES_MINIMIZE_CONTACT_OBJECTIVE_

#include <utility>

#include "ompl/base/OptimizationObjective.h"
#include "ompl/base/samplers/informed/PathLengthDirectInfSampler.h"

namespace ompl
{
    namespace base
    {

        using ContactCalc = std::function<double(const base::State *)>;

        /**
         * Optimization objective that computes the cost of a robot state based on the added magnitude of the vectors
         * acting on it.
         */
        class MinimizeContactObjective : public ompl::base::OptimizationObjective
        {
        public:
            /** Constructor. */
            MinimizeContactObjective(const ompl::base::SpaceInformationPtr &si, ContactCalc contactCalc)
              : ompl::base::OptimizationObjective(si), contactCalc_(std::move(contactCalc))
            {
                description_ = "Minimize Contact Objective";
            }

            /** Assume we can always do better. */
            bool isSatisfied(ompl::base::Cost) const override
            {
                return false;
            }

            /** \brief Returns a cost with a value of 0. */
            Cost stateCost(const State *state) const override
            {
                return Cost(contactCalc_(state));
            }

            bool isSymmetric() const override
            {
                return false;
            }

            // ompl::base::Cost motionCostHeuristic(const State *s1, const State *s2) const
            // {
            //     return ompl::base::Cost(si_->distance(s1, s2));
            // }

            // ompl::base::InformedSamplerPtr allocInformedStateSampler(const ProblemDefinitionPtr &probDefn,
            //                                                          unsigned int maxNumberCalls) const
            // {
            //     // Make the direct path-length informed sampler and return. If OMPL was compiled with Eigen, a direct
            //     // version is available, if not a rejection-based technique can be used
            //     return std::make_shared<PathLengthDirectInfSampler>(probDefn, maxNumberCalls);
            // }

            ompl::base::Cost motionCost(const State *s1, const State *s2) const override
            {
                const base::StateSpacePtr &space = si_->getStateSpace();
                std::size_t vfdim = space->getValueLocations().size();
                std::size_t numSegments = space->validSegmentCount(s1, s2);
                std::vector<ompl::base::State *> interp;
                si_->getMotionStates(s1, s2, interp, numSegments - 1, true, true);
                double total_cost = 0;
                double dist = si_->distance(s1, s2);
                for (std::size_t i = 0; i < interp.size() - 1; i++)
                {
                    double unit_cost = contactCalc_(interp[i]);
                    total_cost += unit_cost;
                    si_->freeState(interp[i]);
                }
                si_->freeState(interp[interp.size() - 1]);

                total_cost += si_->distance(s1, s2);

                return ompl::base::Cost(total_cost);
            }

        protected:
            /** ContactCalc associated with the space. */
            ContactCalc contactCalc_;
        };
    }  // namespace base
}  // namespace ompl

#endif