#ifndef OMPL_GEOMETRIC_PLANNERS_RRT_CVFRRT_
#define OMPL_GEOMETRIC_PLANNERS_RRT_CVFRRT_

#include <limits>

#include <Eigen/Core>

#include <ompl/geometric/planners/rrt/RRT.h>

namespace ompl
{
    namespace geometric
    {

        class CVFRRT : public RRT
        {
        public:
            using VectorField = std::function<Eigen::VectorXd(const base::State *)>;

            /** Constructor. */
            CVFRRT(const base::SpaceInformationPtr &si, VectorField vf, double lambda, unsigned int update_freq);

            /** Destructor. */
            ~CVFRRT() override;

            /** Reset internal data. */
            void clear() override;

            /** Use the vector field to alter the direction of a sample. */
            Eigen::VectorXd getNewDirection(const base::State *qnear, const base::State *qrand);

            /**
             * Every nth time this function is called (where the nth step is the update
             * frequency given in the constructor) the value of lambda is updated and
             * the counts of efficient and inefficient samples added to the tree are
             * reset to 0. The measure for exploration inefficiency is also reset to 0.
             */
            void updateGain();

            /**
             * This attempts to extend the tree from the motion m to a new motion
             * in the direction specified by the vector v.
             */
            Motion *extendTree(Motion *m, base::State *rstate, const Eigen::VectorXd &v);

            /** Solve the planning problem. */
            base::PlannerStatus solve(const base::PlannerTerminationCondition &ptc) override;

            void setup() override;

        private:
            /** Vector field for the environemnt. */
            const VectorField vf_;

            /** User-specified lambda parameter. */
            double lambda_;

            /** The number of steps until lambda is updated and efficiency metrics are reset. */
            unsigned int nth_step_;

            /** Current number of steps since lambda was updated/initialized. */
            unsigned int step_{0u};

            /** Dimensionality of vector field */
            unsigned int vfdim_{0u};
        };
    }  // namespace geometric
}  // namespace ompl
#endif
