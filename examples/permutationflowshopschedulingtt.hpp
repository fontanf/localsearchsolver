/**
 * Permutation flow shop scheduling problem, Total tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/permutationflowshopschedulingtt.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/permutationflowshopschedulingtt.hpp"

namespace localsearchsolver
{

namespace permutationflowshopschedulingtt
{

using namespace orproblems::permutationflowshopschedulingtt;

class SequencingScheme
{

public:

    /**
     * Global cost:
     * - Total tardiness
     */
    using GlobalCost = std::tuple<Time>;

    struct SequenceData
    {
        std::vector<Time> times;
        Time total_tardiness = 0;
    };

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 2;
        parameters.swap_block_maximum_length = 1;
        parameters.reverse = false;
        parameters.shift_reverse_block_maximum_length = 0;

        parameters.ruin_and_recreate_number_of_perturbations = 4;
        parameters.ruin_and_recreate_number_of_elements_removed = 4;

        return parameters;
    }

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline SequenceData empty_sequence_data(sequencing::SequenceId) const
    {
        SequenceData sequence_data;
        sequence_data.times = std::vector<Time>(instance_.number_of_machines(), 0);
        return sequence_data;
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {sequence_data.total_tardiness};
    }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {sequence_data.total_tardiness};
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {value};
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
    }

private:

    const Instance& instance_;

};

}

}

