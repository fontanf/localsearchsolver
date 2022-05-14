#pragma once

/**
 * Permutation flow shop scheduling problem, Total completion_time.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/permutationflowshopschedulingtct.hpp
 *
 */

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/permutationflowshopschedulingtct.hpp"

namespace localsearchsolver
{

namespace permutationflowshopschedulingtct
{

using namespace orproblems::permutationflowshopschedulingtct;

class LocalScheme
{

public:

    using ElementId = sequencing2::ElementId;
    using ElementPos = sequencing2::ElementPos;

    /** Global cost: <Number of jobs, Total completion time>; */
    using GlobalCost = std::tuple<ElementPos, Time>;

    inline ElementPos&         number_of_jobs(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time&        total_completion_time(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline ElementPos    number_of_jobs(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time   total_completion_time(const GlobalCost& global_cost) { return std::get<0>(global_cost); }

    /*
     * Sequences.
     */

    struct Sequence
    {
        std::vector<ElementId> sequence;
        std::vector<Time> times;
        Time total_completion_time = 0;
    };

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_block_maximum_length = 2;
            sequencing_parameters.swap_block_maximum_length = 1;
            sequencing_parameters.reverse = false;
            sequencing_parameters.shift_reverse_block_maximum_length = 0;
            sequencing_parameters.double_bridge_number_of_perturbations = 0;
            sequencing_parameters.ruin_and_recreate_number_of_perturbations = 4;
            sequencing_parameters.ruin_and_recreate_number_of_elements_removed = 4;
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
     * Initial sequences.
     */

    inline Sequence empty_sequence(sequencing2::SequenceId) const
    {
        Sequence sequence;
        sequence.times = std::vector<Time>(instance_.number_of_machines(), 0);
        return sequence;
    }

    /*
     * Sequence properties.
     */

    inline GlobalCost global_cost(const Sequence& sequence) const
    {
        return {
            -sequence.sequence.size(),
            sequence.total_completion_time,
        };
    }

    /*
     * Methods required by sequencing2::LocalScheme.
     */

    inline ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline GlobalCost bound(const Sequence& sequence) const
    {
        return {
            -instance_.number_of_jobs(),
            sequence.total_completion_time,
        };
    }

    inline void append(
            Sequence& sequence,
            ElementId j) const
    {
        MachineId m = instance_.number_of_machines();
        ElementId n = instance_.number_of_jobs();
        ElementId n_cur = sequence.sequence.size();
        Time t_prec = sequence.times[m - 1];
        // Update sequence.
        sequence.sequence.push_back(j);
        // Update times.
        sequence.times[0] = sequence.times[0] + instance_.processing_time(j, 0);
        for (MachineId i = 1; i < m; ++i) {
            if (sequence.times[i - 1] > sequence.times[i]) {
                sequence.times[i] = sequence.times[i - 1] + instance_.processing_time(j, i);
            } else {
                sequence.times[i] = sequence.times[i] + instance_.processing_time(j, i);
            }
        }
        // Update total completion_time.
        sequence.total_completion_time += (n - n_cur) * (sequence.times[m - 1] - t_prec);
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

