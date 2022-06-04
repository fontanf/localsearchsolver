#pragma once

/**
 * Single machine scheduling problem with sequence-dependent setup times, Total
 * weighted tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/schedulingwithsdsttwt.hpp
 *
 */


#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing.hpp"

#include "orproblems/schedulingwithsdsttwt.hpp"

namespace localsearchsolver
{

namespace schedulingwithsdsttwt
{

using namespace orproblems::schedulingwithsdsttwt;

class LocalScheme
{

public:

    using ElementId = sequencing::ElementId;
    using ElementPos = sequencing::ElementPos;

    /** Global cost: <Number of jobs, Total weighted tardiness>; */
    using GlobalCost = std::tuple<JobPos, Weight>;

    inline JobPos&                 number_of_jobs(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight&       total_weighted_tardiness(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline JobPos            number_of_jobs(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight  total_weighted_tardiness(const GlobalCost& global_cost) { return std::get<1>(global_cost); }

    /*
     * Sequence.
     */

    struct SequenceData
    {
        JobPos number_of_jobs = 0;
        JobId j_last = -1;
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
            sequencing_parameters.swap_block_maximum_length = 3;
            sequencing_parameters.shuffle_neighborhood_order = true;

            sequencing_parameters.double_bridge_number_of_perturbations = 0;
            sequencing_parameters.ruin_and_recreate_number_of_perturbations = 10;
        }

        sequencing::Parameters sequencing_parameters;
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
     * Sequence properties.
     */

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            -sequence_data.number_of_jobs,
            sequence_data.total_weighted_tardiness,
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            -instance_.number_of_jobs(),
            value,
        };
    }

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            -instance_.number_of_jobs(),
            sequence_data.total_weighted_tardiness,
        };
    }

    inline void append(
            SequenceData& sequence_data,
            ElementId j) const
    {
        // Update number_of_jobs.
        sequence_data.number_of_jobs++;
        // Update time.
        JobId j_prev = (sequence_data.j_last != -1)?
            sequence_data.j_last:
            instance_.number_of_jobs();
        sequence_data.time += instance_.setup_time(j_prev, j);
        sequence_data.time += instance_.job(j).processing_time;
        // Update total weighted tardiness.
        if (sequence_data.time > instance_.job(j).due_date)
            sequence_data.total_weighted_tardiness
                += instance_.job(j).weight
                * (sequence_data.time - instance_.job(j).due_date);
        // Update j_prev.
        sequence_data.j_last = j;
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

