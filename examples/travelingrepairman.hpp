/**
 * Traveling Repairman Problem.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/travelingrepairman.hpp
 *
 */

#pragma once

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing.hpp"

#include "orproblems/travelingrepairman.hpp"

namespace localsearchsolver
{

namespace travelingrepairman
{

using namespace orproblems::travelingrepairman;

class SequencingScheme
{

public:

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 4;
        parameters.swap_block_maximum_length = 2;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 4;

        parameters.ruin_and_recreate_number_of_perturbations = 10;
        parameters.ruin_number_of_elements_removed = 10;

        parameters.crossover_ox_weight = 1;

        return parameters;
    }

    /**
     * Global cost:
     * - Total completion time
     */
    using GlobalCost = std::tuple<Time>;

    struct SequenceData
    {
        LocationPos number_of_locations = 0;
        LocationId j_first = -1;
        LocationId j_last = -1;
        Time time = 0;  // Without depot -> j_first.
        Time total_completion_time = 0;  // Without depot -> j_first.
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const
    {
        // -1 since we don't schedule the depot.
        return instance_.number_of_locations() - 1;
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.total_completion_time
                + sequence_data.number_of_locations
                * instance_.travel_time(0, sequence_data.j_first + 1),
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {value};
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.total_completion_time
                + sequence_data.number_of_locations
                * instance_.travel_time(0, sequence_data.j_first + 1),
        };
    }

    inline SequenceData sequence_data_init(LocationId j) const
    {
        SequenceData sequence_data;
        // Uppdate j_first.
        sequence_data.j_first = j;
        // Update j_last.
        sequence_data.j_last = j;
        // Update number_of_locations.
        sequence_data.number_of_locations = 1;
        return sequence_data;
    }

    bool concatenate(
            SequenceData& sequence_data,
            const SequenceData& sequence_data_2) const
    {
        Time sij = instance_.travel_time(sequence_data.j_last + 1, sequence_data_2.j_first + 1);
        // Update total_completion_time.
        sequence_data.total_completion_time
            += (sequence_data_2.total_completion_time
                    + sequence_data_2.number_of_locations
                    * (sequence_data.time + sij));
        // Update time.
        sequence_data.time += (sij + sequence_data_2.time);
        // Update j_last.
        sequence_data.j_last = sequence_data_2.j_last;
        // Update number_of_locations.
        sequence_data.number_of_locations += sequence_data_2.number_of_locations;
        return true;
    }

private:

    const Instance& instance_;

};

}

}

