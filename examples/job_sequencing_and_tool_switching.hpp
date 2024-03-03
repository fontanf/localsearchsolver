/**
 * Permutation flow shop scheduling problem, total completion_time
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/job_sequencing_and_tool_switching.hpp
 *
 */

#pragma once

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/scheduling/job_sequencing_and_tool_switching.hpp"

namespace localsearchsolver
{
namespace job_sequencing_and_tool_switching
{

using namespace orproblems::job_sequencing_and_tool_switching;

class SequencingScheme
{

public:

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 4;
        parameters.swap_block_maximum_length = 3;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 3;

        parameters.ruin_and_recreate_number_of_perturbations = 4;
        parameters.ruin_number_of_elements_removed = 4;
        parameters.order_crossover_weight = 1;

        return parameters;
    }

    /**
     * Global cost:
     * - Number of switches
     */
    using GlobalCost = std::tuple<ToolId>;

    struct SequenceData
    {
        /** Sequence of jobs. */
        std::vector<JobId> jobs;

        /** Number of switches. */
        ToolId number_of_switches = 0;

        /** For each tool, the positions at which it is required. */
        std::vector<std::vector<JobPos>> tool_positions;

        /** Current positions in tool_positions for each tool. */
        std::vector<JobPos> tool_positions_cur;

        // Magazine.
        optimizationtools::IndexedBinaryHeap<JobPos> magazine;

        /** Next job position to try. */
        JobPos job_pos = 0;
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline SequenceData empty_sequence_data(
            sequencing::SequenceId) const
    {
        SequenceData sequence_data;
        sequence_data.tool_positions = std::vector<std::vector<JobPos>>(instance_.number_of_tools());
        sequence_data.tool_positions_cur = std::vector<JobPos>(instance_.number_of_tools(), 1);
        sequence_data.magazine = optimizationtools::IndexedBinaryHeap<JobPos>(instance_.number_of_tools());
        return sequence_data;
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {sequence_data.number_of_switches};
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {sequence_data.number_of_switches};
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId element_id) const
    {
        // Update sequence_data.jobs.
        sequence_data.jobs.push_back(element_id);

        // Update sequence_data.tool_positions.
        for (ToolId tool_id: instance_.tools(element_id)) {
            sequence_data.tool_positions[tool_id].push_back(sequence_data.jobs.size() - 1);

            if (sequence_data.magazine.contains(tool_id)
                    && sequence_data.magazine.key(tool_id) == -instance_.number_of_jobs()) {
                sequence_data.magazine.update_key(tool_id, -(sequence_data.jobs.size() - 1));
            }
        }

        // All jobs have been scheduled.
        for (;
                sequence_data.job_pos < (JobPos)sequence_data.jobs.size();
                ++sequence_data.job_pos) {

            // If the magazine contains a job for which we don't know the
            // next position at which it will be required, then stop.
            if (!sequence_data.magazine.empty()
                    && sequence_data.magazine.top().second
                    == -instance_.number_of_jobs()) {
                break;
            }

            JobId job_id = sequence_data.jobs[sequence_data.job_pos];

            // Add the tools of the current job in the magasin.
            // Set the key to ensure that they remain there during the next
            // step.
            for (ToolId tool_id: instance_.tools(job_id)) {
                if (!sequence_data.magazine.contains(tool_id))
                    sequence_data.number_of_switches++;
                sequence_data.magazine.update_key(tool_id, 0);
            }

            // Remove tools from the magazine until the magazine magazine_capacity
            // constraint is satisfied.
            while (sequence_data.magazine.size() > instance_.magazine_capacity())
                sequence_data.magazine.pop();

            // Set the right key for the tools of the current job.
            for (ToolId tool_id: instance_.tools(job_id)) {
                JobPos pos = sequence_data.tool_positions_cur[tool_id];
                // If the tool is not neaded anymore, remove it from the heap.
                if (pos == (JobPos)instance_.jobs(tool_id).size()) {
                    sequence_data.magazine.update_key(
                            tool_id,
                            -(instance_.number_of_jobs() + 1));
                    sequence_data.magazine.pop();
                } else if (pos == (JobPos)sequence_data.tool_positions[tool_id].size()) {
                    // If the tool will still be needed, but we don't know
                    // when yet, add it back with a speical key value.
                    sequence_data.magazine.update_key(
                            tool_id,
                            -instance_.number_of_jobs());
                    sequence_data.tool_positions_cur[tool_id]++;
                } else {
                    // Otherwise, the its right key value.
                    sequence_data.magazine.update_key(
                            tool_id,
                            -sequence_data.tool_positions[tool_id][pos]);
                    sequence_data.tool_positions_cur[tool_id]++;
                }
            }
        }
    }

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Job sequencing and tool switching" << std::endl;
        instance_.format(os, verbosity_level);
    }

private:

    /** Instance. */
    const Instance& instance_;

};

}
}
