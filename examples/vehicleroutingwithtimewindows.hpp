/**
 * Vehicle routing problem with time windows
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

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 2;
        parameters.swap_block_maximum_length = 2;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 2;

        parameters.inter_shift_block_maximum_length = 2;
        parameters.inter_shift_reverse_block_maximum_length = 2;
        parameters.inter_swap_block_maximum_length = 2;
        parameters.swap_tails = true;
        //parameters.inter_swap_star = true;

        parameters.ruin_and_recreate_number_of_perturbations = 128;
        parameters.ruin_number_of_elements_removed = 10;
        parameters.ruin_adjacent_string_removal_weight = 1.0;
        parameters.recreate_best_weight = 1.0;

        parameters.selective_route_exchange_crossover_1_weight = 1.0;

        return parameters;
    }

    /**
     * Global cost:
     * - Reversed time
     * - Overcapacity
     * - Total travel time
     */
    using GlobalCost = std::tuple<Time, Demand, Time>;

    struct SequenceData
    {
        sequencing::ElementId element_id_first = -1;
        sequencing::ElementId element_id_last = -1;
        Demand demand = 0;

        Time duration = 0;
        Time earliest_start = 0;
        Time latest_start = 0;
        Time reversed_time = 0;

        Time total_travel_time = 0;  // Without depot -> element_id_first.
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::SequencePos number_of_sequences() const { return instance_.number_of_vehicles(); }

    inline sequencing::ElementPos number_of_elements() const
    {
        // -1 since we don't schedule the depot.
        return instance_.number_of_locations() - 1;
    }

    inline double distance(
            sequencing::ElementId element_id_1,
            sequencing::ElementId element_id_2) const
    {
        Time r1 = instance_.location(element_id_1 + 1).release_date;
        Time d1 = instance_.location(element_id_1 + 1).deadline;
        Time s1 = instance_.location(element_id_1 + 1).service_time;
        Time r2 = instance_.location(element_id_2 + 1).release_date;
        Time d2 = instance_.location(element_id_2 + 1).deadline;
        Time t12 = instance_.travel_time(element_id_1 + 1, element_id_2 + 1);
        Time wt = r2 - (d1 + t12 + s1);
        Time tw = (r1 + s1 + t12) - d2;
        return (double)t12 + 1 * std::max((Time)0, wt) + 0.2 * std::max((Time)0, tw);
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        if (sequence_data.element_id_first == -1)
            return {0, 0, 0};
        SequenceData sd = sequence_data_init(-1);
        concatenate(sd, sequence_data);
        concatenate(sd, sequence_data_init(-1));
        return {
            sd.reversed_time,
            std::max((Demand)0, sequence_data.demand - instance_.capacity()),
            sd.total_travel_time,
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            0,
            0,
            value,
        };
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.reversed_time,
            std::max((Demand)0, sequence_data.demand - instance_.capacity()),
            sequence_data.total_travel_time
                + instance_.travel_time(0, sequence_data.element_id_first + 1),
        };
    }

    inline SequenceData sequence_data_init(
            sequencing::ElementId element_id) const
    {
        SequenceData sequence_data;
        // Uppdate element_id_first.
        sequence_data.element_id_first = element_id;
        // Update time windows.
        sequence_data.duration = instance_.location(element_id + 1).service_time;
        sequence_data.earliest_start = instance_.location(element_id + 1).release_date;
        sequence_data.latest_start = instance_.location(element_id + 1).deadline;
        // Update demand.
        sequence_data.demand += instance_.location(element_id + 1).demand;
        // Update element_id_last.
        sequence_data.element_id_last = element_id;
        return sequence_data;
    }

    inline bool concatenate(
            SequenceData& sequence_data,
            const SequenceData& sequence_data_2) const
    {
        Time tij = instance_.travel_time(
                sequence_data.element_id_last + 1,
                sequence_data_2.element_id_first + 1);
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
        // Update element_id_last.
        sequence_data.element_id_last = sequence_data_2.element_id_last;
        return true;
    }

private:

    const Instance& instance_;

};

}

}

