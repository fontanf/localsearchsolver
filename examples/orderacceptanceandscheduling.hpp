/**
 * Single machine order acceptance and scheduling problem with time windows and
 * sequence_data-dependent setup times, total weighted tardiness
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/orderacceptanceandscheduling.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/orderacceptanceandscheduling.hpp"

namespace localsearchsolver
{

namespace orderacceptanceandscheduling
{

using namespace orproblems::orderacceptanceandscheduling;

class SequencingScheme
{

public:

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 7;
        parameters.swap_block_maximum_length = 5;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 6;
        parameters.add_remove = true;
        parameters.replace = true;

        parameters.force_add = true;

        return parameters;
    }

    /**
     * Global cost:
     * - Reversed time
     * - Total weighted tardiness - Profit
     */
    using GlobalCost = std::tuple<Time, Weight>;

    struct SequenceData
    {
        sequencing::ElementId element_id_last = -1;
        Time time = 0;
        Time reversed_time_curr = 0;
        Time reversed_time_full = 0;
        Weight total_weighted_tardiness_curr = 0;
        Weight total_weighted_tardiness_full = 0;
        Profit profit = 0;
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs() - 2; }

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

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.reversed_time_curr,
            std::numeric_limits<Profit>::lowest(),
        };
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId element_id) const
    {
        // Update time.
        Time rj = instance_.job(element_id + 1).release_date;
        if (sequence_data.time < rj)
            sequence_data.time = rj;
        sequence_data.time += instance_.setup_time(
                sequence_data.element_id_last + 1,
                element_id + 1);
        sequence_data.time += instance_.job(element_id + 1).processing_time;
        // Update reversed_time.
        Time dj = instance_.job(element_id + 1).deadline;
        if (sequence_data.time > dj) {
            sequence_data.reversed_time_curr += (sequence_data.time - dj);
            sequence_data.time = dj;
        }
        // Update total weighted tardiness.
        if (sequence_data.time > instance_.job(element_id + 1).due_date)
            sequence_data.total_weighted_tardiness_curr
                += instance_.job(element_id + 1).weight
                * (sequence_data.time - instance_.job(element_id + 1).due_date);
        // Update profit.
        sequence_data.profit += instance_.job(element_id + 1).profit;
        // Update reversed_time_full and total_weighted_tardiness_full.
        sequence_data.reversed_time_full = sequence_data.reversed_time_curr;
        sequence_data.total_weighted_tardiness_full = sequence_data.total_weighted_tardiness_curr;
        Time time_full = sequence_data.time;
        JobId jn = instance_.number_of_jobs() - 1;
        Time rjn = instance_.job(jn).release_date;
        if (time_full < rjn)
            time_full = rjn;
        time_full += instance_.setup_time(element_id + 1, jn);
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
        // Update element_id_last.
        sequence_data.element_id_last = element_id;
    }

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Single machine order acceptance and scheduling problem with time windows and sequence_data-dependent setup times, total weighted tardiness" << std::endl;
        instance_.format(os, verbosity_level);
    }

private:

    /** Instdance. */
    const Instance& instance_;

};

}

}

