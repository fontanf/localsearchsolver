/**
 * Permutation flow shop scheduling problem, total completion_time
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/job_sequencing_and_tool_switching.hpp
 *
 */

#include "read_args.hpp"

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/scheduling/job_sequencing_and_tool_switching.hpp"

using namespace localsearchsolver;
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
        parameters.shift_reverse_block_maximum_length = 4;

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

    SequencingScheme(const Instance& instance):
        instance_(instance)
    {
        // Compute the distances between jobs.
        distances_ = std::vector<std::vector<ToolId>>(
                instance.number_of_jobs(),
                std::vector<ToolId>(instance.number_of_jobs()));
        optimizationtools::IndexedSet tools_tmp(instance.number_of_tools());
        for (JobId job_id = 0;
                job_id < instance.number_of_jobs();
                ++job_id) {
            tools_tmp.clear();
            for (ToolId tool_id: instance.tools(job_id))
                tools_tmp.add(tool_id);
            for (JobId job_id_2 = 0;
                    job_id_2 < instance.number_of_jobs();
                    ++job_id_2) {
                for (ToolId tool_id: instance.tools(job_id_2)) {
                    if (!tools_tmp.contains(tool_id)) {
                        distances_[job_id][job_id_2]++;
                        distances_[job_id_2][job_id]++;
                    }
                }
            }
        }
    }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline double distance(
            sequencing::ElementId element_id_1,
            sequencing::ElementId element_id_2) const
    {
        return distances_[element_id_1][element_id_2];
    }

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

    std::vector<std::vector<ToolId>> distances_;

};

int main(int argc, char *argv[])
{
    // Create command line options.
    boost::program_options::options_description desc = setup_args();
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;;
        throw "";
    }
    try {
        boost::program_options::notify(vm);
    } catch (const boost::program_options::required_option& e) {
        std::cout << desc << std::endl;;
        throw "";
    }

    // Create instance.
    InstanceBuilder instance_builder;
    instance_builder.read(
            vm["input"].as<std::string>(),
            vm["format"].as<std::string>());
    const Instance instance = instance_builder.build();

    // Create local scheme.
    SequencingScheme sequencing_scheme(instance);
    auto sequencing_parameters = read_sequencing_args<SequencingScheme>(vm);
    sequencing::LocalScheme<SequencingScheme> local_scheme(
            sequencing_scheme,
            sequencing_parameters);

    // Run algorithm.
    std::string algorithm = vm["algorithm"].as<std::string>();
    auto output =
        (algorithm == "multi-start-local-search")?
        run_multi_start_local_search(local_scheme, vm):
        (algorithm == "iterated-local-search")?
        run_iterated_local_search(local_scheme, vm):
        (algorithm == "best-first-local-search")?
        run_best_first_local_search(local_scheme, vm):
        run_genetic_local_search(local_scheme, vm);

    // Run checker.
    if (vm["print-checker"].as<int>() > 0
            && vm["certificate"].as<std::string>() != "") {
        std::cout << std::endl
            << "Checker" << std::endl
            << "-------" << std::endl;
        instance.check(
                vm["certificate"].as<std::string>(),
                std::cout,
                vm["print-checker"].as<int>());
    }

    return 0;
}
