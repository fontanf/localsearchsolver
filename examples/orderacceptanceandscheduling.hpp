/**
 * Single machine order acceptance and scheduling problem with time windows and
 * sequence_data-dependent setup times, Total weighted tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/orderacceptanceandscheduling.hpp
 *
 */

#pragma once

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing.hpp"

#include "orproblems/orderacceptanceandscheduling.hpp"

namespace localsearchsolver
{

namespace orderacceptanceandscheduling
{

using namespace orproblems::orderacceptanceandscheduling;

class LocalScheme
{

public:

    /**
     * Global cost:
     * - Reversed time
     * - Total weighted tardiness - Profit
     */
    using GlobalCost = std::tuple<Time, Weight>;

    struct SequenceData
    {
        JobId j_last = -1;
        Time time = 0;
        Time reversed_time_curr = 0;
        Time reversed_time_full = 0;
        Weight total_weighted_tardiness_curr = 0;
        Weight total_weighted_tardiness_full = 0;
        Profit profit = 0;
    };

    /*
     * Constructors and destructor.
     */

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

    inline SequenceData empty_sequence_data() const
    {
        SequenceData sequence_data;
        sequence_data.profit += instance_.job(0).profit;
        sequence_data.profit += instance_.job(instance_.number_of_jobs() - 1).profit;
        return sequence_data;
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.reversed_time_full,
            sequence_data.total_weighted_tardiness_full - sequence_data.profit,
        };
    }

    inline JobPos number_of_elements() const { return instance_.number_of_jobs() - 2; }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.reversed_time_curr,
            std::numeric_limits<Profit>::lowest(),
        };
    }

    inline void append(
            SequenceData& sequence_data,
            JobId j) const
    {
        // Update time.
        Time rj = instance_.job(j + 1).release_date;
        if (sequence_data.time < rj)
            sequence_data.time = rj;
        sequence_data.time += instance_.setup_time(sequence_data.j_last + 1, j + 1);
        sequence_data.time += instance_.job(j + 1).processing_time;
        // Update reversed_time.
        Time dj = instance_.job(j + 1).deadline;
        if (sequence_data.time > dj) {
            sequence_data.reversed_time_curr += (sequence_data.time - dj);
            sequence_data.time = dj;
        }
        // Update total weighted tardiness.
        if (sequence_data.time > instance_.job(j + 1).due_date)
            sequence_data.total_weighted_tardiness_curr
                += instance_.job(j + 1).weight
                * (sequence_data.time - instance_.job(j + 1).due_date);
        // Update profit.
        sequence_data.profit += instance_.job(j + 1).profit;
        // Update reversed_time_full and total_weighted_tardiness_full.
        sequence_data.reversed_time_full = sequence_data.reversed_time_curr;
        sequence_data.total_weighted_tardiness_full = sequence_data.total_weighted_tardiness_curr;
        Time time_full = sequence_data.time;
        JobId jn = instance_.number_of_jobs() - 1;
        Time rjn = instance_.job(jn).release_date;
        if (time_full < rjn)
            time_full = rjn;
        time_full += instance_.setup_time(j + 1, jn);
        time_full += instance_.job(jn).processing_time;
        Time djn = instance_.job(jn).deadline;
        if (time_full > djn) {
            sequence_data.reversed_time_full += (time_full - djn);
            time_full = djn;
        }
        if (time_full > instance_.job(jn).due_date)
            sequence_data.total_weighted_tardiness_full
                += instance_.job(jn).weight
                * (time_full - instance_.job(jn).due_date);
        // Update j_last.
        sequence_data.j_last = j;
    }

private:

    const Instance& instance_;
    Parameters parameters_;

};

}

}

