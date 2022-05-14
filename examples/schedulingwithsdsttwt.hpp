#pragma once

/**
 * Single machine scheduling problem with sequence-dependent setup times, Total
 * weighted tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/schedulingwithsdsttwt.hpp
 *
 */


#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/schedulingwithsdsttwt.hpp"

namespace localsearchsolver
{

namespace schedulingwithsdsttwt
{

using namespace orproblems::schedulingwithsdsttwt;

class LocalScheme
{

public:

    using ElementId = sequencing2::ElementId;
    using ElementPos = sequencing2::ElementPos;

    /** Global cost: <Number of jobs, Total weighted tardiness>; */
    using GlobalCost = std::tuple<ElementPos, Weight>;

    inline ElementPos&             number_of_jobs(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight&       total_weighted_tardiness(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline ElementPos        number_of_jobs(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight  total_weighted_tardiness(const GlobalCost& global_cost) { return std::get<1>(global_cost); }

    /*
     * Sequence.
     */

    struct Sequence
    {
        std::vector<ElementId> sequence;
        Time time = 0;
        Weight total_weighted_tardiness = 0;
    };

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_block_maximum_length = 13;
            sequencing_parameters.swap_block_maximum_length = 3;
            sequencing_parameters.shuffle_neighborhood_order = true;
            sequencing_parameters.double_bridge_number_of_perturbations = 0;
            sequencing_parameters.ruin_and_recreate_number_of_perturbations = 10;
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
     * Sequence properties.
     */

    inline GlobalCost global_cost(const Sequence& sequence) const
    {
        return {
            -sequence.sequence.size(),
            sequence.total_weighted_tardiness,
        };
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            -instance_.number_of_jobs(),
            value,
        };
    }

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline GlobalCost bound(const Sequence& sequence) const
    {
        return {
            -instance_.number_of_jobs(),
            sequence.total_weighted_tardiness,
        };
    }

    inline void append(
            Sequence& sequence,
            ElementId j) const
    {
        ElementPos ns = sequence.sequence.size();
        // Update time.
        ElementId j_prev = (ns)?
            sequence.sequence.back():
            instance_.number_of_jobs();
        sequence.time += instance_.setup_time(j_prev, j);
        sequence.time += instance_.job(j).processing_time;
        // Update jobs.
        sequence.sequence.push_back(j);
        // Update total weighted tardiness.
        if (sequence.time > instance_.job(j).due_date)
            sequence.total_weighted_tardiness
                += instance_.job(j).weight
                * (sequence.time - instance_.job(j).due_date);
        if (ns > 2
                && instance_.job(sequence.sequence[ns - 1]).weight > 0
                && instance_.job(sequence.sequence[ns - 2]).weight == 0)
            sequence.total_weighted_tardiness += 1000000;
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

