/**
 * Time-dependent orienteering problem
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/routing/team_orienteering.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/routing/team_orienteering.hpp"

namespace localsearchsolver
{
namespace team_orienteering
{

using namespace orproblems::team_orienteering;

class SequencingScheme
{

public:

    using Time = int64_t;

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 7;
        parameters.swap_block_maximum_length = 5;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 6;

        parameters.add_remove = true;
        parameters.replace = true;

        parameters.inter_shift_block_maximum_length = 3;
        parameters.inter_swap_block_maximum_length = 3;
        parameters.swap_tails = true;
        parameters.split = true;
        parameters.inter_shift_reverse_block_maximum_length = 3;

        parameters.force_add = true;

        return parameters;
    }

    /**
     * Global cost:
     * - Overtime
     * - Profit
     * - Total time
     */
    using GlobalCost = std::tuple<Time, Profit, Time>;

    struct SequenceData
    {
        sequencing::ElementId element_id_first = -1;
        sequencing::ElementId element_id_last = -1;
        Time duration = 0;  // Without depot -> element_id_first.
        Profit profit = 0;
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::SequencePos number_of_sequences() const { return instance_.number_of_vehicles(); }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_locations() - 2; }

    inline Time travel_time(
            LocationId location_id_1,
            LocationId location_id_2) const
    {
        return 1e4 * instance_.travel_time(location_id_1, location_id_2);
    }

    inline Time maximum_duration() const
    {
        return 1e4 * instance_.maximum_duration();
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        if (sequence_data.element_id_first == -1)
            return {0, 0, 0};
        Time total_duration = travel_time(0, sequence_data.element_id_first + 1)
                + sequence_data.duration
                + travel_time(
                        sequence_data.element_id_last + 1,
                        instance_.number_of_locations() - 1);
        return {
            std::max((Time)0, total_duration - maximum_duration()),
            -sequence_data.profit,
            total_duration,
        };
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            std::max((Time)0, sequence_data.duration - maximum_duration()),
            std::numeric_limits<Profit>::lowest(),
            sequence_data.duration,
        };
    }

    inline SequenceData sequence_data_init(
            sequencing::ElementId element_id) const
    {
        SequenceData sequence_data;
        // Uppdate element_id_first.
        sequence_data.element_id_first = element_id;
        // Update profit.
        sequence_data.profit = instance_.location(element_id + 1).profit;
        // Update element_id_last.
        sequence_data.element_id_last = element_id;
        return sequence_data;
    }

    inline bool concatenate(
            SequenceData& sequence_data,
            const SequenceData& sequence_data_2) const
    {
        // Update total_completion_time.
        sequence_data.duration
            += travel_time(
                    sequence_data.element_id_last + 1,
                    sequence_data_2.element_id_first + 1)
            + sequence_data_2.duration;
        // Update profit.
        sequence_data.profit += sequence_data_2.profit;
        // Update element_id_last.
        sequence_data.element_id_last = sequence_data_2.element_id_last;
        return true;
    }

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Team orienteering problem" << std::endl;
        instance_.format(os, verbosity_level);
    }

    void solution_write(
            const sequencing::LocalScheme<SequencingScheme>::Solution& solution,
            const std::string& certificate_path)
    {
        if (certificate_path.empty())
            return;
        std::ofstream file(certificate_path);
        if (!file.good()) {
            throw std::runtime_error(
                    "Unable to open file \"" + certificate_path + "\".");
        }

        for (VehicleId vehicle_id = 0;
                vehicle_id < instance_.number_of_vehicles();
                ++vehicle_id) {
            file << solution.sequences[vehicle_id].elements.size() << std::endl;
            for (auto se: solution.sequences[vehicle_id].elements)
                file << " " << se.element_id + 1;
            file << std::endl;
        }
    }

private:

    /** Instance. */
    const Instance& instance_;

};

}
}
