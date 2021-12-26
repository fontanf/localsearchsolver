#pragma once

/**
 * Permutation flow shop scheduling problem, Total completion_time.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/permutationflowshopschedulingtct.hpp
 *
 */

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/permutationflowshopschedulingtct.hpp"

namespace localsearchsolver
{

namespace permutationflowshopschedulingtct
{

using namespace orproblems::permutationflowshopschedulingtct;

class LocalScheme
{

public:

    /** Global cost: <Number of jobs, Total completion time>; */
    using GlobalCost = std::tuple<JobPos, Time>;

    inline JobPos&             number_of_jobs(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time&        total_completion_time(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline JobPos        number_of_jobs(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time   total_completion_time(const GlobalCost& global_cost) { return std::get<0>(global_cost); }

    /*
     * Solutions.
     */

    struct Solution
    {
        std::vector<JobId> sequence;
        std::vector<Time> times;
        Time total_completion_time = 0;
    };

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_block_maximum_length = 7;
            sequencing_parameters.swap_block_maximum_length = 2;
            sequencing_parameters.reverse = true;
            sequencing_parameters.shift_reverse_block_maximum_length = 4;
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
            solution.total_completion_time,
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
            solution.total_completion_time,
        };
    }

    inline void append(
            Solution& solution,
            JobId j) const
    {
        MachineId m = instance_.number_of_machines();
        JobId n = instance_.number_of_jobs();
        JobId n_cur = solution.sequence.size();
        Time t_prec = solution.times[m - 1];
        // Update sequence.
        solution.sequence.push_back(j);
        // Update times.
        solution.times[0] = solution.times[0] + instance_.processing_time(j, 0);
        for (MachineId i = 1; i < m; ++i) {
            if (solution.times[i - 1] > solution.times[i]) {
                solution.times[i] = solution.times[i - 1] + instance_.processing_time(j, i);
            } else {
                solution.times[i] = solution.times[i] + instance_.processing_time(j, i);
            }
        }
        // Update total completion_time.
        solution.total_completion_time += (n - n_cur) * (solution.times[m - 1] - t_prec);
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

