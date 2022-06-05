/**
 * Traveling salesman problem with release dates.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/travelingsalesmanwithreleasedates.hpp
 *
 */

#pragma once

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing.hpp"

#include "orproblems/travelingsalesmanwithreleasedates.hpp"

namespace localsearchsolver
{

namespace travelingsalesmanwithreleasedates
{

using namespace orproblems::travelingsalesmanwithreleasedates;

class LocalScheme
{

public:

    using ElementId = sequencing::ElementId;
    using ElementPos = sequencing::ElementPos;
    using Mode = sequencing::Mode;

    /** Global cost: <Number of locations, Total duration>; */
    using GlobalCost = std::tuple<ElementPos, Time>;

    /*
     * SequenceDatas.
     */

    struct SequenceData
    {
        LocationPos number_of_locations = 0;
        LocationId j_last = -1;
        Time current_trip_start = 0;
        Time current_trip_duration = 0;
        Time time_full = 0;
    };

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shuffle_neighborhood_order = true;

            sequencing_parameters.shift_block_maximum_length = 4;
            sequencing_parameters.swap_block_maximum_length = 2;
            sequencing_parameters.reverse = true;
            sequencing_parameters.shift_reverse_block_maximum_length = 3;

            sequencing_parameters.shift_change_mode = true;
            sequencing_parameters.mode_swap = true;
            sequencing_parameters.swap_with_modes = true;

            sequencing_parameters.double_bridge_number_of_perturbations = 0;
            sequencing_parameters.ruin_and_recreate_number_of_perturbations = 10;
            sequencing_parameters.ruin_and_recreate_number_of_elements_removed = 10;
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

    /*
     * SequenceData properties.
     */

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            -sequence_data.number_of_locations,
            sequence_data.time_full,
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            -instance_.number_of_locations(),
            value,
        };
    }

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline ElementPos number_of_elements() const
    {
        // -1 since we don't schedule the depot.
        return instance_.number_of_locations() - 1;
    }

    inline Mode number_of_modes(ElementId) const
    {
        // Mode 0: Visit next location without returning to the depot.
        //         The last departure might be delayed.
        // Mode 1: Return to the depot before visiting next location.
        return 2;
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            -instance_.number_of_locations(),
            sequence_data.time_full,
        };
    }

    inline void append(
            SequenceData& sequence_data,
            ElementId j,
            Mode mode) const
    {
        // Update last_start and time.
        Time rj = instance_.release_date(j + 1);
        if (mode == 0) {  // No return to depot.
            if (sequence_data.current_trip_start < rj)
                sequence_data.current_trip_start = rj;
            sequence_data.current_trip_duration += instance_.travel_time(sequence_data.j_last + 1, j + 1);
        } else {  // Return to depot.
            sequence_data.current_trip_start = sequence_data.time_full;
            if (sequence_data.current_trip_start < rj)
                sequence_data.current_trip_start = rj;
            sequence_data.current_trip_duration = instance_.travel_time(0, j + 1);
        }
        sequence_data.time_full = sequence_data.current_trip_start
            + sequence_data.current_trip_duration
            + instance_.travel_time(j + 1, 0);
        // Update j_last.
        sequence_data.j_last = j;
        // Update number_of_vertices.
        sequence_data.number_of_locations++;
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

