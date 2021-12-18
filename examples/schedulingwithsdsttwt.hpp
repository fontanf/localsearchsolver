#pragma once

/**
 * Single machine scheduling problem with sequence-dependent setup times, Total
 * weighted Tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/schedulingwithsdsttwt.hpp
 *
 */


#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/schedulingwithsdsttwt.hpp"

namespace localsearchsolver
{

namespace schedulingwithsdsttwt
{

using namespace orproblems::schedulingwithsdsttwt;

class LocalScheme
{

public:

    /** Global cost: <Total weighted tardiness>; */
    using GlobalCost = std::tuple<Weight>;

    inline Weight&       total_weighted_tardiness(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight  total_weighted_tardiness(const GlobalCost& global_cost) { return std::get<0>(global_cost); }

    /*
     * Solutions.
     */

    struct Solution
    {
        std::vector<JobId> sequence;
        Time time = 0;
        Weight total_weighted_tardiness = 0;
    };

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_block_maximum_length = 13;
            sequencing_parameters.swap_block_maximum_length = 1;
            sequencing_parameters.shuffle_neighborhood_order = true;
            sequencing_parameters.number_of_perturbations = 10;
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
        solution.total_weighted_tardiness = std::numeric_limits<Weight>::max();
        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            solution.total_weighted_tardiness,
        };
    }

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline JobPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline void append(
            Solution& solution,
            JobId j) const
    {
        if (solution.sequence.size() == 0)
            solution.total_weighted_tardiness = 0;
        // Update time.
        JobId j_prev = (solution.sequence.size() > 0)?
            solution.sequence.back():
            instance_.number_of_jobs();
        solution.time += instance_.setup_time(j_prev, j);
        solution.time += instance_.job(j).processing_time;
        // Update jobs.
        solution.sequence.push_back(j);
        // Update total weighted tardiness.
        if (solution.time > instance_.job(j).due_date)
            solution.total_weighted_tardiness
                += instance_.job(j).weight
                * (solution.time - instance_.job(j).due_date);
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

