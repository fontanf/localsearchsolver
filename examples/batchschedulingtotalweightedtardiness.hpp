/**
 * Single machine batch scheduling problem, Total weighted tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/batchschedulingtotalweightedtardiness.hpp
 *
 */

#pragma once

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing.hpp"

#include "orproblems/batchschedulingtotalweightedtardiness.hpp"

namespace localsearchsolver
{

namespace batchschedulingtotalweightedtardiness
{

using namespace orproblems::batchschedulingtotalweightedtardiness;

class LocalScheme
{

public:

    using ElementId = sequencing::ElementId;
    using ElementPos = sequencing::ElementPos;
    using Mode = sequencing::Mode;

    /** Global cost: <Number of jobs, Overcapacity, Total weighted tardiness>; */
    using GlobalCost = std::tuple<JobPos, Size, Weight>;

    /*
     * SequenceDatas.
     */

    struct SequenceData
    {
        JobPos number_of_jobs = 0;
        std::vector<JobId> current_batch_jobs = {};
        Time current_batch_start = 0;
        Time current_batch_duration = 0;
        Time current_batch_end = 0;
        Size current_batch_size = 0;
        Size overcapacity = 0;
        Weight total_weighted_tardiness = 0;
    };

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shuffle_neighborhood_order = true;

            sequencing_parameters.shift_block_maximum_length = 8;
            sequencing_parameters.swap_block_maximum_length = 4;
            sequencing_parameters.reverse = true;
            sequencing_parameters.shift_reverse_block_maximum_length = 4;

            sequencing_parameters.shift_change_mode = true;
            sequencing_parameters.mode_swap = true;
            sequencing_parameters.swap_with_modes = true;

            sequencing_parameters.double_bridge_number_of_perturbations = 0;
            sequencing_parameters.ruin_and_recreate_number_of_perturbations = 10;
            sequencing_parameters.ruin_and_recreate_number_of_elements_removed = 4;
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
            -sequence_data.number_of_jobs,
            sequence_data.overcapacity,
            sequence_data.total_weighted_tardiness,
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            -instance_.number_of_jobs(),
            0,
            value,
        };
    }

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline ElementPos number_of_elements() const
    {
        return instance_.number_of_jobs();
    }

    inline Mode number_of_modes(ElementId) const
    {
        // Mode 0: Add next job to the current batch.
        // Mode 1: Add next job in a new batch.
        return 2;
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            -instance_.number_of_jobs(),
            sequence_data.overcapacity,
            sequence_data.total_weighted_tardiness,
        };
    }

    inline void append(
            SequenceData& sequence_data,
            ElementId j,
            Mode mode) const
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
        // Update number_of_vertices.
        sequence_data.number_of_jobs++;
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

