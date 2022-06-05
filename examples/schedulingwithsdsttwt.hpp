/**
 * Single machine scheduling problem with sequence-dependent setup times, Total
 * weighted tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/schedulingwithsdsttwt.hpp
 *
 */

#pragma once

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

    /**
     * Global cost:
     * - Number of jobs
     * - Total weighted tardiness
     */
    using GlobalCost = std::tuple<JobPos, Weight>;

    struct SequenceData
    {
        JobPos number_of_jobs = 0;
        JobId j_last = -1;
        Time time = 0;
        Weight total_weighted_tardiness = 0;
    };

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_block_maximum_length = 13;
            sequencing_parameters.swap_block_maximum_length = 3;
            sequencing_parameters.shuffle_neighborhood_order = true;

            sequencing_parameters.double_bridge_number_of_perturbations = 0;
            sequencing_parameters.ruin_and_recreate_number_of_perturbations = 10;
            sequencing_parameters.ruin_and_recreate_number_of_elements_removed = 2;
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

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            -instance_.number_of_jobs(),
            sequence_data.total_weighted_tardiness,
        };
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId j) const
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

    const Instance& instance_;
    Parameters parameters_;

};

}

}

