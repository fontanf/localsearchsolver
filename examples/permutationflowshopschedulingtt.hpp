#pragma once

/**
 * Permutation flow shop scheduling problem, Total tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/permutationflowshopschedulingtt.hpp
 *
 */

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/permutationflowshopschedulingtt.hpp"

namespace localsearchsolver
{

namespace permutationflowshopschedulingtt
{

using namespace orproblems::permutationflowshopschedulingtt;

class LocalScheme
{

public:

    /** Global cost: <Number of jobs, Total tardiness>; */
    using GlobalCost = std::tuple<JobPos, Time>;

    inline JobPos&       number_of_jobs(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time&        total_tardiness(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline JobPos  number_of_jobs(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time   total_tardiness(const GlobalCost& global_cost) { return std::get<0>(global_cost); }

    /*
     * Solutions.
     */

    struct Solution
    {
        std::vector<JobId> sequence;
        std::vector<Time> times;
        Time total_tardiness = 0;
    };

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_block_maximum_length = 2;
            sequencing_parameters.swap_block_maximum_length = 1;
            sequencing_parameters.reverse = false;
            sequencing_parameters.shift_reverse_block_maximum_length = 0;
            sequencing_parameters.double_bridge_number_of_perturbations = 0;
            sequencing_parameters.ruin_and_recreate_number_of_perturbations = 4;
            sequencing_parameters.ruin_and_recreate_number_of_elements_removed = 4;
        }

        sequencing2::Parameters sequencing_parameters;
    };

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters) { }

    LocalScheme(const LocalScheme& local_scheme):
        LocalScheme(local_scheme.instance_, local_scheme.parameters_) { }

    virtual ~LocalScheme() { }

    /*
     * Initial solutions.
     */

    inline Solution empty_solution() const
    {
        Solution solution;
        solution.times = std::vector<Time>(instance_.number_of_machines(), 0);
        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            -solution.sequence.size(),
            solution.total_tardiness,
        };
    }

    /*
     * Methods required by sequencing2::LocalScheme.
     */

    inline JobPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline GlobalCost bound(const Solution& solution) const
    {
        return {
            -instance_.number_of_jobs(),
            solution.total_tardiness,
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            -instance_.number_of_jobs(),
            value,
        };
    }

    inline void append(
            Solution& solution,
            JobId j) const
    {
        MachineId m = instance_.number_of_machines();
        // Update sequence.
        solution.sequence.push_back(j);
        // Update times.
        solution.times[0] = solution.times[0] + instance_.job(j).processing_times[0];
        for (MachineId i = 1; i < m; ++i) {
            if (solution.times[i - 1] > solution.times[i]) {
                solution.times[i] = solution.times[i - 1] + instance_.job(j).processing_times[i];
            } else {
                solution.times[i] = solution.times[i] + instance_.job(j).processing_times[i];
            }
        }
        // Update total tardiness.
        if (solution.times[m - 1] > instance_.job(j).due_date)
            solution.total_tardiness += (solution.times[m - 1] - instance_.job(j).due_date);
    }

private:

    /*
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

};

}

}

