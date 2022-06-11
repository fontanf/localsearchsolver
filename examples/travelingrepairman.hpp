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

    /**
     * Global cost:
     * - Number of locations
     * - Total completion time
     */
    using GlobalCost = std::tuple<LocationPos, Time>;

    struct SequenceData
    {
        LocationPos number_of_locations = 0;
        LocationId j_first = -1;
        LocationId j_last = -1;
        Time time = 0;  // Without depot -> j_first.
        Time total_completion_time = 0;  // Without depot -> j_first.
    };

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_block_maximum_length = 4;
            sequencing_parameters.swap_block_maximum_length = 2;
            sequencing_parameters.reverse = true;
            sequencing_parameters.shift_reverse_block_maximum_length = 4;

            sequencing_parameters.ruin_and_recreate_number_of_perturbations = 10;
            sequencing_parameters.ruin_and_recreate_number_of_elements_removed = 10;
        }

        sequencing::Parameters sequencing_parameters;
    };

    SequencingScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters) { }

    SequencingScheme(const SequencingScheme& sequencing_scheme):
        SequencingScheme(sequencing_scheme.instance_, sequencing_scheme.parameters_) { }

    virtual ~SequencingScheme() { }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            -sequence_data.number_of_locations,
            sequence_data.total_completion_time
                + sequence_data.number_of_locations
                * instance_.travel_time(0, sequence_data.j_first + 1),
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            -instance_.number_of_locations(),
            value,
        };
    }

    inline sequencing::ElementPos number_of_elements() const
    {
        // -1 since we don't schedule the depot.
        return instance_.number_of_locations() - 1;
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            -instance_.number_of_locations(),
            sequence_data.total_completion_time
                + sequence_data.number_of_locations
                * instance_.travel_time(0, sequence_data.j_first + 1),
        };
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId j) const
    {
        if (sequence_data.number_of_locations == 0) {
            // Uppdate j_first.
            sequence_data.j_first = j;
        } else {
            // Update time.
            sequence_data.time += instance_.travel_time(sequence_data.j_last + 1, j + 1);
        }
        // Update total_completion_time.
        sequence_data.total_completion_time += sequence_data.time;
        // Update j_last.
        sequence_data.j_last = j;
        // Update number_of_locations.
        sequence_data.number_of_locations++;
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
    Parameters parameters_;

};

}

}

