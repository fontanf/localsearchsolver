#pragma once

/**
 * Sequential Ordering Problem..
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/sequentialordering.hpp
 *
 */


#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/sequentialordering.hpp"

namespace localsearchsolver
{

namespace sequentialordering
{

using namespace orproblems::sequentialordering;

class LocalScheme
{

public:

    using ElementId = sequencing2::ElementId;
    using ElementPos = sequencing2::ElementPos;

    /** Global cost: <Number of vertices, Number of precedence violations, Total distance>; */
    using GlobalCost = std::tuple<ElementPos, ElementPos, Distance>;

    inline ElementPos&                    number_of_vertices(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline ElementPos&       number_of_precedence_violations(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Distance&                                  length(GlobalCost& global_cost) { return std::get<2>(global_cost); }
    inline ElementPos               number_of_vertices(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline ElementPos  number_of_precedence_violations(const GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Distance                             length(const GlobalCost& global_cost) { return std::get<2>(global_cost); }

    /*
     * Sequences.
     */

    struct Sequence
    {
        std::vector<ElementId> sequence;
        Distance length = 0;
        ElementPos number_of_precedence_violations = 0;
        std::vector<uint8_t> contains;
    };

    /*
     * Constructors and destructor.
     */

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
     * Initial sequences.
     */

    inline Sequence empty_sequence(sequencing2::SequenceId) const
    {
        Sequence sequence;
        sequence.contains = std::vector<uint8_t>(
                instance_.number_of_vertices(), 0);
        return sequence;
    }

    /*
     * Sequence properties.
     */

    inline GlobalCost global_cost(const Sequence& sequence) const
    {
        return {
            -sequence.sequence.size(),
            sequence.number_of_precedence_violations,
            sequence.length,
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

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline ElementPos number_of_elements() const { return instance_.number_of_vertices(); }

    inline GlobalCost bound(const Sequence& sequence) const
    {
        return {
            -instance_.number_of_vertices(),
            sequence.number_of_precedence_violations,
            sequence.length,
        };
    }

    inline void append(
            Sequence& sequence,
            ElementPos j) const
    {
        // Update number_of_precedence_violations.
        for (ElementId j: instance_.predecessors(j))
            if (!sequence.contains[j])
                sequence.number_of_precedence_violations++;
        // Update time.
        if (sequence.sequence.size() > 0) {
            ElementId j_prev = sequence.sequence.back();
            Distance d = instance_.distance(j_prev, j);
            if (d == std::numeric_limits<Distance>::max())
                d = instance_.distance(j, j_prev);
            sequence.length += d;
        }
        // Update contains.
        sequence.contains[j] = 1;
        // Update sequence.
        sequence.sequence.push_back(j);
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

