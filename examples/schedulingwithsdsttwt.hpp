/**
 * Single machine scheduling problem with sequence-dependent setup times, total
 * weighted tardiness
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

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 13;
        parameters.swap_block_maximum_length = 3;

        parameters.ruin_and_recreate_number_of_perturbations = 10;
        parameters.ruin_number_of_elements_removed = 2;

        return parameters;
    }

    /**
     * Global cost:
     * - Total weighted tardiness
     */
    using GlobalCost = std::tuple<Weight>;

    struct SequenceData
    {
        sequencing::ElementId element_id_last = -1;
        Time time = 0;
        Weight total_weighted_tardiness = 0;
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {sequence_data.total_weighted_tardiness};
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {value};
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {sequence_data.total_weighted_tardiness};
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId element_id) const
    {
        // Update time.
        JobId job_id_prev = (sequence_data.element_id_last != -1)?
            sequence_data.element_id_last:
            instance_.number_of_jobs();
        sequence_data.time += instance_.setup_time(job_id_prev, element_id);
        sequence_data.time += instance_.job(element_id).processing_time;
        // Update total weighted tardiness.
        if (sequence_data.time > instance_.job(element_id).due_date)
            sequence_data.total_weighted_tardiness
                += instance_.job(element_id).weight
                * (sequence_data.time - instance_.job(element_id).due_date);
        // Update job_id_prev.
        sequence_data.element_id_last = element_id;
    }

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Single machine scheduling problem with sequence-dependent setup times, total weighted tardines" << std::endl;
        instance_.format(os, verbosity_level);
    }

private:

    /** Instance. */
    const Instance& instance_;

};

}

}

