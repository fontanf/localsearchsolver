/**
 * Vehicle routing problem with time windows.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/vehicleroutingwithtimewindows.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/vehicleroutingwithtimewindows.hpp"

namespace localsearchsolver
{

namespace vehicleroutingwithtimewindows
{

using namespace orproblems::vehicleroutingwithtimewindows;

class SequencingScheme
{

public:

    /**
     * Global cost:
     * - Number of locations
     * - Reversed time
     * - Overcapacity
     * - Total travel time
     */
    using GlobalCost = std::tuple<LocationPos, Time, Demand, Time>;

    struct SequenceData
    {
        LocationPos number_of_locations = 0;
        LocationId j_first = -1;
        LocationId j_last = -1;
        Demand demand = 0;

        Time duration = 0;
        Time earliest_start = 0;
        Time latest_start = 0;
        Time reversed_time = 0;

        Time total_travel_time = 0;  // Without depot -> j_first.
    };

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 2;
        parameters.swap_block_maximum_length = 2;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 2;

        parameters.inter_shift_block_maximum_length = 1;
        parameters.inter_swap_block_maximum_length = 1;
        parameters.inter_two_opt = true;
        parameters.inter_shift_reverse_block_maximum_length = 0;
        parameters.inter_swap_star = true;

        parameters.ruin_and_recreate_number_of_perturbations = 10;
        parameters.ruin_and_recreate_number_of_elements_removed = 10;

        parameters.crossover_srex1_weight = 1;

        return parameters;
    }

    SequencingScheme(
            const Instance& instance):
        instance_(instance) { }

    SequencingScheme(const SequencingScheme& sequencing_scheme):
        SequencingScheme(sequencing_scheme.instance_) { }

    virtual ~SequencingScheme() { }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        if (sequence_data.number_of_locations == 0)
            return {0, 0, 0, 0};
        SequenceData sd;
        append(sd, -1);
        concatenate(sd, sequence_data);
        append(sd, -1);
        return {
            -sequence_data.number_of_locations,
            sd.reversed_time,
            std::max((Demand)0, sequence_data.demand - instance_.capacity()),
            sd.total_travel_time,
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            -instance_.number_of_locations(),
            0,
            0,
            value,
        };
    }

    inline sequencing::SequencePos number_of_sequences() const { return instance_.number_of_vehicles(); }

    inline sequencing::ElementPos number_of_elements() const
    {
        // -1 since we don't schedule the depot.
        return instance_.number_of_locations() - 1;
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            -instance_.number_of_locations(),
            sequence_data.reversed_time,
            std::max((Demand)0, sequence_data.demand - instance_.capacity()),
            sequence_data.total_travel_time
                + instance_.travel_time(0, sequence_data.j_first + 1),
        };
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId j) const
    {
        Time rj = instance_.location(j + 1).release_date;
        Time dj = instance_.location(j + 1).deadline;
        Time tj = instance_.location(j + 1).service_time;
        if (sequence_data.number_of_locations == 0) {
            // Uppdate j_first.
            sequence_data.j_first = j;
            // Update time windows.
            sequence_data.duration = tj;
            sequence_data.earliest_start = rj;
            sequence_data.latest_start = dj;
        } else {
            Time tij = instance_.travel_time(sequence_data.j_last + 1, j + 1);
            // Update time windows.
            Time delta = sequence_data.duration - sequence_data.reversed_time + tij;
            Time delta_wt = std::max(rj - delta - sequence_data.latest_start, (Time)0);
            Time delta_tw = std::max(sequence_data.earliest_start + delta - dj, (Time)0);
            sequence_data.duration += tj + tij + delta_wt;
            sequence_data.earliest_start = std::max(
                    rj - delta,
                    sequence_data.earliest_start) - delta_wt;
            sequence_data.latest_start = std::min(
                    dj - delta,
                    sequence_data.latest_start) + delta_tw;
            sequence_data.reversed_time += delta_tw;
            // Update total_travel_time.
            sequence_data.total_travel_time += tij;
        }
        // Update demand.
        sequence_data.demand += instance_.location(j + 1).demand;
        // Update j_last.
        sequence_data.j_last = j;
        // Update number_of_locations.
        sequence_data.number_of_locations++;
    }

    bool concatenate(
            SequenceData& sequence_data,
            const SequenceData& sequence_data_2) const
    {
        Time tij = instance_.travel_time(sequence_data.j_last + 1, sequence_data_2.j_first + 1);
        // Update time windows.
        Time delta = sequence_data.duration - sequence_data.reversed_time + tij;
        Time delta_wt = std::max(sequence_data_2.earliest_start - delta - sequence_data.latest_start, (Time)0);
        Time delta_tw = std::max(sequence_data.earliest_start + delta - sequence_data_2.latest_start, (Time)0);
        sequence_data.duration += sequence_data_2.duration + tij + delta_wt;
        sequence_data.earliest_start = std::max(
                sequence_data_2.earliest_start - delta,
                sequence_data.earliest_start) - delta_wt;
        sequence_data.latest_start = std::min(
                sequence_data_2.latest_start - delta,
                sequence_data.latest_start) + delta_tw;
        sequence_data.reversed_time += sequence_data_2.reversed_time + delta_tw;
        // Update total_completion_time.
        sequence_data.total_travel_time += tij + sequence_data_2.total_travel_time;
        // Update demand.
        sequence_data.demand += sequence_data_2.demand;
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

