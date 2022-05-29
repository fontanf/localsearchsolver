#pragma once

/**
 * Permutation flow shop scheduling problem, Total completion_time.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/distributedpfsstct.hpp
 *
 */

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/distributedpfsstct.hpp"

namespace localsearchsolver
{

namespace distributedpfsstct
{

using namespace orproblems::distributedpfsstct;

class LocalScheme
{

public:

    using ElementId = sequencing2::ElementId;
    using ElementPos = sequencing2::ElementPos;
    using SequencePos = sequencing2::SequencePos;

    /** Global cost: <Number of jobs, Total completion time>; */
    using GlobalCost = std::tuple<ElementPos, Time>;

    inline JobPos&            number_of_jobs(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time&       total_completion_time(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline JobPos       number_of_jobs(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time  total_completion_time(const GlobalCost& global_cost) { return std::get<0>(global_cost); }

    /*
     * SequenceDatas.
     */

    struct SequenceData
    {
        JobPos number_of_jobs = 0;
        JobId j_last = -1;
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
            sequencing_parameters.shift_block_maximum_length = 3;
            sequencing_parameters.swap_block_maximum_length = 2;
            sequencing_parameters.reverse = true;
            sequencing_parameters.shift_reverse_block_maximum_length = 2;

            sequencing_parameters.inter_shift_block_maximum_length = 3;
            sequencing_parameters.inter_swap_block_maximum_length = 2;
            sequencing_parameters.inter_two_opt = true;
            sequencing_parameters.inter_shift_reverse_block_maximum_length = 2;

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
     * Initial sequence_datas.
     */

    inline SequenceData empty_sequence_data(sequencing2::SequenceId) const
    {
        SequenceData sequence_data;
        sequence_data.times = std::vector<Time>(instance_.number_of_machines(), 0);
        return sequence_data;
    }

    /*
     * SequenceData properties.
     */

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            -sequence_data.number_of_jobs,
            sequence_data.total_completion_time,
        };
    }

    /*
     * Methods required by sequencing2::LocalScheme.
     */

    inline SequencePos number_of_sequences() const { return instance_.number_of_factories(); }

    inline ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            -instance_.number_of_jobs(),
            sequence_data.total_completion_time,
        };
    }

    inline void append(
            SequenceData& sequence_data,
            ElementId j) const
    {
        MachineId m = instance_.number_of_machines();
        // Update times.
        sequence_data.times[0] = sequence_data.times[0] + instance_.processing_time(j, 0);
        for (MachineId i = 1; i < m; ++i) {
            if (sequence_data.times[i - 1] > sequence_data.times[i]) {
                sequence_data.times[i] = sequence_data.times[i - 1] + instance_.processing_time(j, i);
            } else {
                sequence_data.times[i] = sequence_data.times[i] + instance_.processing_time(j, i);
            }
        }
        // Update total completion_time.
        sequence_data.total_completion_time += sequence_data.times[m - 1];
        // Update number_of_jobs.
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

