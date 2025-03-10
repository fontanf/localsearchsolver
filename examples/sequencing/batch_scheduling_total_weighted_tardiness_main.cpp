/**
 * Single machine batch scheduling problem, total weighted tardiness
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/batch_scheduling_total_weighted_tardiness.hpp
 *
 */

#include "read_args.hpp"

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/scheduling/batch_scheduling_total_weighted_tardiness.hpp"

using namespace localsearchsolver;
using namespace orproblems::batch_scheduling_total_weighted_tardiness;

class SequencingScheme
{

public:

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 8;
        parameters.swap_block_maximum_length = 4;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 4;

        parameters.shift_change_mode_block_maximum_length = 8;
        parameters.mode_swap = true;
        parameters.swap_with_modes = true;

        parameters.ruin_and_recreate_number_of_perturbations = 10;
        parameters.ruin_number_of_elements_removed = 4;

        return parameters;
    }

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

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline sequencing::Mode number_of_modes(sequencing::ElementId) const
    {
        // Mode 0: Add next job to the current batch.
        // Mode 1: Add next job in a new batch.
        return 2;
    }

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

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.overcapacity,
            sequence_data.total_weighted_tardiness,
        };
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId element_id,
            sequencing::Mode mode) const
    {
        const Job& job = instance_.job(element_id);
        if (mode == 1) {  // New batch.
            sequence_data.current_batch_jobs.clear();
            sequence_data.current_batch_start = sequence_data.current_batch_end;
            sequence_data.current_batch_duration = 0;
            sequence_data.current_batch_end = sequence_data.current_batch_end;
            sequence_data.current_batch_size = 0;
        }
        // Update current_batch_start.
        if (sequence_data.current_batch_start < job.release_date)
            sequence_data.current_batch_start = job.release_date;
        // Update current_batch_duration.
        if (sequence_data.current_batch_duration < job.processing_time)
            sequence_data.current_batch_duration = job.processing_time;
        // Update total_weighted_tardiness from jobs of the current batch.
        Time new_end = sequence_data.current_batch_start
            + sequence_data.current_batch_duration;
        Time diff = new_end - sequence_data.current_batch_end;
        for (JobId job_id_2: sequence_data.current_batch_jobs) {
            const Job& job_2 = instance_.job(job_id_2);
            if (sequence_data.current_batch_end >= job_2.due_date) {
                sequence_data.total_weighted_tardiness
                    += job_2.weight * diff;
            } else if (new_end <= job_2.due_date) {
            } else {
                sequence_data.total_weighted_tardiness
                    += job_2.weight * (new_end - job_2.due_date);
            }
        }
        // Update current_batch_end.
        sequence_data.current_batch_end = new_end;
        // Update overcapacity.
        Size c = instance_.capacity();
        if (sequence_data.current_batch_size >= c) {
            sequence_data.overcapacity += job.size;
        } else if (sequence_data.current_batch_size + job.size <= c) {
        } else {
            sequence_data.overcapacity += sequence_data.current_batch_size
                + job.size - c;
        }
        // Update current_batch_size.
        sequence_data.current_batch_size += job.size;
        // Update total_weighted_tardiness.
        if (sequence_data.current_batch_end >= job.due_date)
            sequence_data.total_weighted_tardiness += job.weight
                * (sequence_data.current_batch_end - job.due_date);
        // Update current_batch_jobs.
        sequence_data.current_batch_jobs.push_back(element_id);
    }

    void instance_format(
            std::ostream& os,
            int verbosity_level)
    {
        os << "Single machine batch scheduling problem, total weighted tardiness" << std::endl;
        instance_.format(os, verbosity_level);
    }

    void solution_write(
            const sequencing::LocalScheme<SequencingScheme>::Solution& solution,
            const std::string& certificate_path)
    {
        if (certificate_path.empty())
            return;
        std::ofstream file(certificate_path);
        if (!file.good()) {
            throw std::runtime_error(
                    "Unable to open file \"" + certificate_path + "\".");
        }
        std::vector<std::vector<JobId>> batches;
        for (auto se: solution.sequences[0].elements) {
            if (batches.empty() || se.mode == 1)
                batches.push_back({});
            batches.back().push_back(se.element_id);
        }
        for (const auto& batch: batches) {
            file << batch.size() << std::endl;
            for (JobId job_id: batch)
                file << job_id << " ";
            file << std::endl;
        }
    }

private:

    /** Instance. */
    const Instance& instance_;

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
