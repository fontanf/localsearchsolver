/**
 * Sequential Ordering Problem.
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
    using GlobalCost = std::tuple<VertexPos, Distance>;

    struct SequenceData
    {
        VertexId j_last = -1;
        Distance length = 0;
        VertexPos number_of_precedence_violations = 0;
        std::vector<uint8_t> contains;
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_vertices(); }

    inline SequenceData empty_sequence_data(sequencing::SequenceId) const
    {
        SequenceData sequence_data;
        sequence_data.contains = std::vector<uint8_t>(
                instance_.number_of_vertices(), 0);
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
            sequencing::ElementId j) const
    {
        // Update number_of_precedence_violations.
        for (VertexId j: instance_.predecessors(j))
            if (!sequence_data.contains[j])
                sequence_data.number_of_precedence_violations++;
        // Update time.
        if (sequence_data.j_last != -1) {
            Distance d = instance_.distance(sequence_data.j_last, j);
            if (d == std::numeric_limits<Distance>::max())
                d = instance_.distance(j, sequence_data.j_last);
            sequence_data.length += d;
        }
        // Update contains.
        sequence_data.contains[j] = 1;
        // Update j_last.
        sequence_data.j_last = j;
    }

private:

    const Instance& instance_;

};

}

}

