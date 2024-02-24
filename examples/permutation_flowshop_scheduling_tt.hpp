/**
 * Permutation flow shop scheduling problem, total tardiness
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/permutation_flowshop_scheduling_tt.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/scheduling/permutation_flowshop_scheduling_tt.hpp"

namespace localsearchsolver
{

namespace permutation_flowshop_scheduling_tt
{

using namespace orproblems::permutation_flowshop_scheduling_tt;

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
     * - Total tardiness
     */
    using GlobalCost = std::tuple<Time>;

    struct SequenceData
    {
        std::vector<Time> times;
        Time total_tardiness = 0;
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
        return {sequence_data.total_tardiness};
    }

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
            sequencing::ElementId element_id) const
    {
        MachineId m = instance_.number_of_machines();
        // Update times.
        sequence_data.times[0] = sequence_data.times[0]
            + instance_.job(element_id).processing_times[0];
        for (MachineId machine_id = 1; machine_id < m; ++machine_id) {
            if (sequence_data.times[machine_id - 1] > sequence_data.times[machine_id]) {
                sequence_data.times[machine_id] = sequence_data.times[machine_id - 1]
                    + instance_.job(element_id).processing_times[machine_id];
            } else {
                sequence_data.times[machine_id] = sequence_data.times[machine_id]
                    + instance_.job(element_id).processing_times[machine_id];
            }
        }
        // Update total tardiness.
        if (sequence_data.times[m - 1] > instance_.job(element_id).due_date) {
            sequence_data.total_tardiness += (sequence_data.times[m - 1]
                    - instance_.job(element_id).due_date);
        }
    }

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Permutation flow shop scheduling problem, total tardiness" << std::endl;
        instance_.format(os, verbosity_level);
    }

private:

    /** Instance. */
    const Instance& instance_;

};

}

}

