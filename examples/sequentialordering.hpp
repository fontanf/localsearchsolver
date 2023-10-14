/**
 * Sequential ordering problem
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/sequentialordering.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/sequentialordering.hpp"

namespace localsearchsolver
{

namespace sequentialordering
{

using namespace orproblems::sequentialordering;

class SequencingScheme
{

public:

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 8;
        parameters.swap_block_maximum_length = 2;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 3;

        parameters.ruin_and_recreate_number_of_perturbations = 10;

        return parameters;
    }

    /**
     * Global cost:
     * - Number of precedence violations
     * - Total distance
     */
    using GlobalCost = std::tuple<LocationPos, Distance>;

    struct SequenceData
    {
        sequencing::ElementId element_id_last = -1;
        Distance length = 0;
        LocationPos number_of_precedence_violations = 0;
        std::vector<uint8_t> contains;
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_locations(); }

    inline SequenceData empty_sequence_data(sequencing::SequenceId) const
    {
        SequenceData sequence_data;
        sequence_data.contains = std::vector<uint8_t>(
                instance_.number_of_locations(), 0);
        return sequence_data;
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.number_of_precedence_violations,
            sequence_data.length,
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            0,
            value,
        };
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.number_of_precedence_violations,
            sequence_data.length,
        };
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId element_id) const
    {
        // Update number_of_precedence_violations.
        for (LocationId location_id_pred: instance_.predecessors(element_id))
            if (!sequence_data.contains[location_id_pred])
                sequence_data.number_of_precedence_violations++;
        // Update time.
        if (sequence_data.element_id_last != -1) {
            Distance d = instance_.distance(sequence_data.element_id_last, element_id);
            if (d == std::numeric_limits<Distance>::max())
                d = instance_.distance(element_id, sequence_data.element_id_last);
            sequence_data.length += d;
        }
        // Update contains.
        sequence_data.contains[element_id] = 1;
        // Update element_id_last.
        sequence_data.element_id_last = element_id;
    }

private:

    const Instance& instance_;

};

}

}

