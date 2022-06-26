/**
 * Single machine scheduling problem with sequence-dependent setup times, Total
 * weighted tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/schedulingwithsdsttwt.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/schedulingwithsdsttwt.hpp"

namespace localsearchsolver
{

namespace schedulingwithsdsttwt
{

using namespace orproblems::schedulingwithsdsttwt;

class SequencingScheme
{

public:

    /**
     * Global cost:
     * - Total weighted tardiness
     */
    using GlobalCost = std::tuple<Weight>;

    struct SequenceData
    {
        JobId j_last = -1;
        Time time = 0;
        Weight total_weighted_tardiness = 0;
    };

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 13;
        parameters.swap_block_maximum_length = 3;

        parameters.ruin_and_recreate_number_of_perturbations = 10;
        parameters.ruin_and_recreate_number_of_elements_removed = 2;

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
        return {sequence_data.total_weighted_tardiness};
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {value};
    }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {sequence_data.total_weighted_tardiness};
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId j) const
    {
        // Update time.
        JobId j_prev = (sequence_data.j_last != -1)?
            sequence_data.j_last:
            instance_.number_of_jobs();
        sequence_data.time += instance_.setup_time(j_prev, j);
        sequence_data.time += instance_.job(j).processing_time;
        // Update total weighted tardiness.
        if (sequence_data.time > instance_.job(j).due_date)
            sequence_data.total_weighted_tardiness
                += instance_.job(j).weight
                * (sequence_data.time - instance_.job(j).due_date);
        // Update j_prev.
        sequence_data.j_last = j;
    }

private:

    const Instance& instance_;

};

}

}

