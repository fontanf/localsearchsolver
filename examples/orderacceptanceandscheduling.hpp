#pragma once

/**
 * Single machine order acceptance and scheduling problem with time windows and
 * sequence-dependent setup times, Total weighted tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/orderacceptanceandscheduling.hpp
 *
 */


#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/orderacceptanceandscheduling.hpp"

namespace localsearchsolver
{

namespace orderacceptanceandscheduling
{

using namespace orproblems::orderacceptanceandscheduling;

class LocalScheme
{

public:

    /** Global cost: <Reversed time, Total weighted tardiness - Profit>; */
    using GlobalCost = std::tuple<Time, Weight>;

    inline Time&        reverse_time(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight&         objective(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Time  reversed_time(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight    objective(const GlobalCost& global_cost) { return std::get<1>(global_cost); }

    /*
     * Solutions.
     */

    struct Sequence
    {
        std::vector<JobId> sequence;
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

        sequencing2::Parameters sequencing_parameters;
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
     * Initial solutions.
     */

    inline Sequence empty_sequence() const
    {
        Sequence sequence;
        sequence.profit += instance_.job(0).profit;
        sequence.profit += instance_.job(instance_.number_of_jobs() - 1).profit;
        return sequence;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Sequence& sequence) const
    {
        return {
            sequence.reversed_time_full,
            sequence.total_weighted_tardiness_full - sequence.profit,
        };
    }

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline JobPos number_of_elements() const { return instance_.number_of_jobs() - 2; }

    inline GlobalCost bound(const Sequence& sequence) const
    {
        return {
            sequence.reversed_time_curr,
            std::numeric_limits<Profit>::lowest(),
        };
    }

    inline void append(
            Sequence& sequence,
            JobId j) const
    {
        // Update time.
        Time rj = instance_.job(j + 1).release_date;
        if (sequence.time < rj)
            sequence.time = rj;
        JobId j_prev = (sequence.sequence.size() > 0)?
            sequence.sequence.back():
            -1;
        sequence.time += instance_.setup_time(j_prev + 1, j + 1);
        sequence.time += instance_.job(j + 1).processing_time;
        // Update reversed_time.
        Time dj = instance_.job(j + 1).deadline;
        if (sequence.time > dj) {
            sequence.reversed_time_curr += (sequence.time - dj);
            sequence.time = dj;
        }
        // Update jobs.
        sequence.sequence.push_back(j);
        // Update total weighted tardiness.
        if (sequence.time > instance_.job(j + 1).due_date)
            sequence.total_weighted_tardiness_curr
                += instance_.job(j + 1).weight
                * (sequence.time - instance_.job(j + 1).due_date);
        // Update profit.
        sequence.profit += instance_.job(j + 1).profit;

        sequence.reversed_time_full = sequence.reversed_time_curr;
        sequence.total_weighted_tardiness_full = sequence.total_weighted_tardiness_curr;
        Time time_full = sequence.time;
        JobId jn = instance_.number_of_jobs() - 1;
        Time rjn = instance_.job(jn).release_date;
        if (time_full < rjn)
            time_full = rjn;
        time_full += instance_.setup_time(j + 1, jn);
        time_full += instance_.job(jn).processing_time;
        Time djn = instance_.job(jn).deadline;
        if (time_full > djn) {
            sequence.reversed_time_full += (time_full - djn);
            time_full = djn;
        }
        if (time_full > instance_.job(jn).due_date)
            sequence.total_weighted_tardiness_full
                += instance_.job(jn).weight
                * (time_full - instance_.job(jn).due_date);
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

