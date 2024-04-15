/**
 * Traveling salesman problem with release dates
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/traveling_salesman_with_release_dates.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/routing/traveling_salesman_with_release_dates.hpp"

namespace localsearchsolver
{
namespace traveling_salesman_with_release_dates
{

using namespace orproblems::traveling_salesman_with_release_dates;

class SequencingScheme
{

public:

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 4;
        parameters.swap_block_maximum_length = 2;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 3;

        parameters.shift_change_mode_block_maximum_length = 1;
        parameters.mode_swap = true;
        parameters.swap_with_modes = true;

        parameters.ruin_and_recreate_number_of_perturbations = 10;
        parameters.ruin_number_of_elements_removed = 10;

        return parameters;
    }

    /**
     * Global cost:
     * - Total duration
     */
    using GlobalCost = std::tuple<Time>;

    struct SequenceData
    {
        LocationId location_id_last = -1;
        Time current_trip_start = 0;
        Time current_trip_duration = 0;
        Time time_full = 0;
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const
    {
        // -1 since we don't schedule the depot.
        return instance_.number_of_locations() - 1;
    }

    inline sequencing::Mode number_of_modes(sequencing::ElementId) const
    {
        // Mode 0: Visit next location without returning to the depot.
        //         The last departure might be delayed.
        // Mode 1: Return to the depot before visiting next location.
        return 2;
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {sequence_data.time_full};
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {value};
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {sequence_data.time_full};
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId element_id,
            sequencing::Mode mode) const
    {
        // Update last_start and time.
        Time rj = instance_.release_date(element_id + 1);
        if (mode == 0) {  // No return to depot.
            if (sequence_data.current_trip_start < rj)
                sequence_data.current_trip_start = rj;
            sequence_data.current_trip_duration += instance_.travel_time(
                        sequence_data.location_id_last + 1,
                        element_id + 1);
        } else {  // Return to depot.
            sequence_data.current_trip_start = sequence_data.time_full;
            if (sequence_data.current_trip_start < rj)
                sequence_data.current_trip_start = rj;
            sequence_data.current_trip_duration = instance_.travel_time(0, element_id + 1);
        }
        sequence_data.time_full = sequence_data.current_trip_start
            + sequence_data.current_trip_duration
            + instance_.travel_time(element_id + 1, 0);
        // Update location_id_last.
        sequence_data.location_id_last = element_id;
    }

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Traveling salesman problem with release dates" << std::endl;
        instance_.format(os, verbosity_level);
    }

private:

    /** Instance. */
    const Instance& instance_;

};

}
}
