/**
 * Time-dependent orienteering problem.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/timedependentorienteering.hpp
 *
 */

#pragma once

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing.hpp"

#include "orproblems/timedependentorienteering.hpp"

namespace localsearchsolver
{

namespace timedependentorienteering
{

using namespace orproblems::timedependentorienteering;

class LocalScheme
{

public:

    /**
     * Global cost:
     * - Overtime
     * - Profit
     * - Total time
     */
    using GlobalCost = std::tuple<Time, Profit, Time>;

    struct SequenceData
    {
        LocationId j_last = -1;
        Time time_cur = 0;
        Time time_full = 0;
        Profit profit = 0;
    };

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_block_maximum_length = 7;
            sequencing_parameters.swap_block_maximum_length = 5;
            sequencing_parameters.reverse = true;
            sequencing_parameters.shift_reverse_block_maximum_length = 6;
            sequencing_parameters.add_remove = true;

            sequencing_parameters.double_bridge_number_of_perturbations = 0;
            sequencing_parameters.ruin_and_recreate_number_of_perturbations = 0;
            sequencing_parameters.force_add = true;
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

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            std::max((Time)0, sequence_data.time_full - instance_.maximum_duration()),
            -sequence_data.profit,
            sequence_data.time_full,
        };
    }

    inline LocationPos number_of_elements() const { return instance_.number_of_locations() - 2; }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            std::max((Time)0, sequence_data.time_full - instance_.maximum_duration()),
            std::numeric_limits<Profit>::lowest(),
            sequence_data.time_full,
        };
    }

    inline void append(
            SequenceData& sequence_data,
            LocationId j) const
    {
        // Update time_cur.
        sequence_data.time_cur = instance_.arrival_time(
                sequence_data.j_last + 1, j + 1, sequence_data.time_cur);
        // Update profit.
        sequence_data.profit += instance_.location(j + 1).profit;
        // Update time_full.
        LocationId jn = instance_.number_of_locations() - 1;
        sequence_data.time_full = instance_.arrival_time(j + 1, jn, sequence_data.time_cur);
        // Update j_last.
        sequence_data.j_last = j;
    }

private:

    const Instance& instance_;
    Parameters parameters_;

};

}

}

