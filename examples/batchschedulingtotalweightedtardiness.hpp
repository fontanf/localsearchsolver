/**
 * Single machine batch scheduling problem, Total weighted tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/batchschedulingtotalweightedtardiness.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/batchschedulingtotalweightedtardiness.hpp"

namespace localsearchsolver
{

namespace batchschedulingtotalweightedtardiness
{

using namespace orproblems::batchschedulingtotalweightedtardiness;

class SequencingScheme
{

public:

    /**
     * Global cost:
     * - Overcapacity
     * - Total weighted tardiness
     */
    using GlobalCost = std::tuple<Size, Weight>;

    struct SequenceData
    {
        std::vector<JobId> current_batch_jobs = {};
        Time current_batch_start = 0;
        Time current_batch_duration = 0;
        Time current_batch_end = 0;
        Size current_batch_size = 0;
        Size overcapacity = 0;
        Weight total_weighted_tardiness = 0;
    };

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 8;
        parameters.swap_block_maximum_length = 4;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 4;

        parameters.shift_change_mode = true;
        parameters.mode_swap = true;
        parameters.swap_with_modes = true;

        parameters.ruin_and_recreate_number_of_perturbations = 10;
        parameters.ruin_number_of_elements_removed = 4;

        return parameters;
    }

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.overcapacity,
            sequence_data.total_weighted_tardiness,
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            0,
            value,
        };
    }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline sequencing::Mode number_of_modes(sequencing::ElementId) const
    {
        // Mode 0: Add next job to the current batch.
        // Mode 1: Add next job in a new batch.
        return 2;
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.overcapacity,
            sequence_data.total_weighted_tardiness,
        };
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId j,
            sequencing::Mode mode) const
    {
        if (mode == 1) {  // New batch.
            sequence_data.current_batch_jobs.clear();
            sequence_data.current_batch_start = sequence_data.current_batch_end;
            sequence_data.current_batch_duration = 0;
            sequence_data.current_batch_end = sequence_data.current_batch_end;
            sequence_data.current_batch_size = 0;
        }
        // Update current_batch_start.
        if (sequence_data.current_batch_start < instance_.job(j).release_date)
            sequence_data.current_batch_start = instance_.job(j).release_date;
        // Update current_batch_duration.
        if (sequence_data.current_batch_duration < instance_.job(j).processing_time)
            sequence_data.current_batch_duration = instance_.job(j).processing_time;
        // Update total_weighted_tardiness from jobs of the current batch.
        Time new_end = sequence_data.current_batch_start + sequence_data.current_batch_duration;
        Time diff = new_end - sequence_data.current_batch_end;
        for (JobId j2: sequence_data.current_batch_jobs) {
            Time dj2 = instance_.job(j2).due_date;
            Weight wj2 = instance_.job(j2).weight;
            if (sequence_data.current_batch_end >= dj2) {
                sequence_data.total_weighted_tardiness += wj2 * diff;
            } else if (new_end <= dj2) {
            } else {
                sequence_data.total_weighted_tardiness += wj2 * (new_end - dj2);
            }
        }
        // Update current_batch_end.
        sequence_data.current_batch_end = new_end;
        // Update overcapacity.
        Size sj = instance_.job(j).size;
        Size c = instance_.capacity();
        if (sequence_data.current_batch_size >= c) {
            sequence_data.overcapacity += sj;
        } else if (sequence_data.current_batch_size + sj <= c) {
        } else {
            sequence_data.overcapacity += sequence_data.current_batch_size + sj - c;
        }
        // Update current_batch_size.
        sequence_data.current_batch_size += instance_.job(j).size;
        // Update total_weighted_tardiness.
        Time dj = instance_.job(j).due_date;
        Weight wj = instance_.job(j).weight;
        if (sequence_data.current_batch_end >= dj)
            sequence_data.total_weighted_tardiness += wj * (sequence_data.current_batch_end - dj);
        // Update current_batch_jobs.
        sequence_data.current_batch_jobs.push_back(j);
    }

private:

    const Instance& instance_;

};

}

}

