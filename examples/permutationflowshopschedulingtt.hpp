/**
 * Permutation flow shop scheduling problem, Total tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/permutationflowshopschedulingtt.hpp
 *
 */

#pragma once

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing.hpp"

#include "orproblems/permutationflowshopschedulingtt.hpp"

namespace localsearchsolver
{

namespace permutationflowshopschedulingtt
{

using namespace orproblems::permutationflowshopschedulingtt;

class LocalScheme
{

public:

    /**
     * Global cost:
     * - Number of jobs
     * - Total tardiness
     */
    using GlobalCost = std::tuple<JobPos, Time>;

    inline JobPos&       number_of_jobs(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time&        total_tardiness(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline JobPos  number_of_jobs(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time   total_tardiness(const GlobalCost& global_cost) { return std::get<0>(global_cost); }

    struct SequenceData
    {
        JobPos number_of_jobs = 0;
        std::vector<Time> times;
        Time total_tardiness = 0;
    };

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

    inline SequenceData empty_sequence_data(sequencing::SequenceId) const
    {
        SequenceData sequence_data;
        sequence_data.times = std::vector<Time>(instance_.number_of_machines(), 0);
        return sequence_data;
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            -sequence_data.number_of_jobs,
            sequence_data.total_tardiness,
        };
    }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            -instance_.number_of_jobs(),
            sequence_data.total_tardiness,
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
            SequenceData& sequence_data,
            sequencing::ElementId j) const
    {
        MachineId m = instance_.number_of_machines();
        // Update times.
        sequence_data.times[0] = sequence_data.times[0] + instance_.job(j).processing_times[0];
        for (MachineId i = 1; i < m; ++i) {
            if (sequence_data.times[i - 1] > sequence_data.times[i]) {
                sequence_data.times[i] = sequence_data.times[i - 1] + instance_.job(j).processing_times[i];
            } else {
                sequence_data.times[i] = sequence_data.times[i] + instance_.job(j).processing_times[i];
            }
        }
        // Update total tardiness.
        if (sequence_data.times[m - 1] > instance_.job(j).due_date)
            sequence_data.total_tardiness += (sequence_data.times[m - 1] - instance_.job(j).due_date);
        // Update number_of_jobs.
        sequence_data.number_of_jobs++;
    }

private:

    const Instance& instance_;
    Parameters parameters_;

};

}

}

