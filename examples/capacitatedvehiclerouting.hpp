/**
 * Capacitated vehicle routing problem.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/capacitatedvehiclerouting.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/capacitatedvehiclerouting.hpp"

namespace localsearchsolver
{

namespace capacitatedvehiclerouting
{

using namespace orproblems::capacitatedvehiclerouting;

class SequencingScheme
{

public:

    /**
     * Global cost:
     * - Overcapacity
     * - Total distance
     */
    using GlobalCost = std::tuple<Demand, Distance>;

    struct SequenceData
    {
        LocationId j_first = -1;
        LocationId j_last = -1;
        Demand demand = 0;
        Distance total_distance = 0;  // Without depot -> j_first.
    };

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 3;
        parameters.swap_block_maximum_length = 3;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 3;

        parameters.inter_shift_block_maximum_length = 3;
        parameters.inter_swap_block_maximum_length = 3;
        parameters.swap_tails = true;
        parameters.split = true;
        parameters.inter_shift_reverse_block_maximum_length = 3;
        //parameters.inter_swap_star = true;

        parameters.ruin_and_recreate_number_of_perturbations = 128;
        parameters.ruin_number_of_elements_removed = 10;
        parameters.ruin_adjacent_string_removal_weight = 1.0;

        parameters.crossover_ox_weight = 1.0;

        return parameters;
    }

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::SequencePos number_of_sequences() const
    {
        // Safety margin: 30% + 3 more vehicles than the trivial bin packing LB
        return std::ceil(1.3 * instance_.total_demand() / instance_.capacity()) + 3;
    }

    inline sequencing::ElementPos number_of_elements() const
    {
        // -1 since we don't schedule the depot.
        return instance_.number_of_locations() - 1;
    }

    inline double distance(sequencing::ElementId j1, sequencing::ElementId j2) const
    {
        return instance_.distance(j1 + 1, j2 + 1);
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        if (sequence_data.j_first == -1)
            return {0, 0};
        return {
            std::max((Demand)0, sequence_data.demand - instance_.capacity()),
            instance_.distance(0, sequence_data.j_first + 1)
                + sequence_data.total_distance
                + instance_.distance(sequence_data.j_last + 1, 0),
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
            std::max((Demand)0, sequence_data.demand - instance_.capacity()),
            sequence_data.total_distance
                + instance_.distance(0, sequence_data.j_first + 1),
        };
    }

    inline SequenceData sequence_data_init(
            sequencing::ElementId j) const
    {
        SequenceData sequence_data;
        // Uppdate j_first.
        sequence_data.j_first = j;
        // Update demand.
        sequence_data.demand = instance_.demand(j + 1);
        // Update j_last.
        sequence_data.j_last = j;
        return sequence_data;
    }

    inline bool concatenate(
            SequenceData& sequence_data,
            const SequenceData& sequence_data_2) const
    {
        // Update total_completion_time.
        sequence_data.total_distance
            += instance_.distance(
                    sequence_data.j_last + 1,
                    sequence_data_2.j_first + 1)
            + sequence_data_2.total_distance;
        // Update demand.
        sequence_data.demand += sequence_data_2.demand;
        // Update j_last.
        sequence_data.j_last = sequence_data_2.j_last;
        return true;
    }

private:

    const Instance& instance_;

};

}

}

