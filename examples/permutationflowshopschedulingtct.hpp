/**
 * Permutation flow shop scheduling problem, Total completion_time.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/permutationflowshopschedulingtct.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/permutationflowshopschedulingtct.hpp"

namespace localsearchsolver
{

namespace permutationflowshopschedulingtct
{

using namespace orproblems::permutationflowshopschedulingtct;

class SequencingScheme
{

public:

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 2;
        parameters.swap_block_maximum_length = 1;
        parameters.reverse = false;
        parameters.shift_reverse_block_maximum_length = 0;

        parameters.ruin_and_recreate_number_of_perturbations = 4;
        parameters.ruin_number_of_elements_removed = 4;

        return parameters;
    }

    /**
     * Global cost:
     * - Total completion time
     */
    using GlobalCost = std::tuple<Time>;

    struct SequenceData
    {
        JobPos number_of_jobs = 0;
        std::vector<Time> times;
        Time total_completion_time = 0;
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline SequenceData empty_sequence_data(sequencing::SequenceId) const
    {
        SequenceData sequence_data;
        sequence_data.times = std::vector<Time>(instance_.number_of_machines(), 0);
        return sequence_data;
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {sequence_data.total_completion_time};
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {sequence_data.total_completion_time};
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId element_id) const
    {
        MachineId m = instance_.number_of_machines();
        JobPos n = instance_.number_of_jobs();
        Time t_prec = sequence_data.times[m - 1];
        // Update times.
        sequence_data.times[0] = sequence_data.times[0]
            + instance_.processing_time(element_id, 0);
        for (MachineId machine_id = 1; machine_id < m; ++machine_id) {
            if (sequence_data.times[machine_id - 1] > sequence_data.times[machine_id]) {
                sequence_data.times[machine_id] = sequence_data.times[machine_id - 1]
                    + instance_.processing_time(element_id, machine_id);
            } else {
                sequence_data.times[machine_id] = sequence_data.times[machine_id]
                    + instance_.processing_time(element_id, machine_id);
            }
        }
        // Update total completion_time.
        sequence_data.total_completion_time
            += (n - sequence_data.number_of_jobs)
            * (sequence_data.times[m - 1] - t_prec);
        // Update number_of_jobs.
        sequence_data.number_of_jobs++;
    }

private:

    const Instance& instance_;

};

}

}

