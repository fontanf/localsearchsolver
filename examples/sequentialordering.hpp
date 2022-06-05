/**
 * Sequential Ordering Problem.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/sequentialordering.hpp
 *
 */

#pragma once

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing.hpp"

#include "orproblems/sequentialordering.hpp"

namespace localsearchsolver
{

namespace sequentialordering
{

using namespace orproblems::sequentialordering;

class LocalScheme
{

public:

    /**
     * Global cost:
     * - Number of vertices
     * - Number of precedence violations
     * - Total distance
     */
    using GlobalCost = std::tuple<VertexPos, VertexPos, Distance>;

    struct SequenceData
    {
        VertexPos number_of_vertices = 0;
        VertexId j_last = -1;
        Distance length = 0;
        VertexPos number_of_precedence_violations = 0;
        std::vector<uint8_t> contains;
    };

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_block_maximum_length = 8;
            sequencing_parameters.swap_block_maximum_length = 2;
            sequencing_parameters.reverse = true;
            sequencing_parameters.shift_reverse_block_maximum_length = 3;
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
            -sequence_data.number_of_vertices,
            sequence_data.number_of_precedence_violations,
            sequence_data.length,
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            -instance_.number_of_vertices(),
            0,
            value,
        };
    }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_vertices(); }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            -instance_.number_of_vertices(),
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
        if (sequence_data.number_of_vertices > 0) {
            Distance d = instance_.distance(sequence_data.j_last, j);
            if (d == std::numeric_limits<Distance>::max())
                d = instance_.distance(j, sequence_data.j_last);
            sequence_data.length += d;
        }
        // Update contains.
        sequence_data.contains[j] = 1;
        // Update j_last.
        sequence_data.j_last = j;
        // Update number_of_vertices.
        sequence_data.number_of_vertices++;
    }

private:

    const Instance& instance_;
    Parameters parameters_;

};

}

}

