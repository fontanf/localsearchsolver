/**
 * Permutation flow shop scheduling problem, total tardiness
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/permutation_flowshop_scheduling_tt.hpp
 *
 */

#include "read_args.hpp"

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/scheduling/permutation_flowshop_scheduling_tt.hpp"

using namespace localsearchsolver;
using namespace orproblems::permutation_flowshop_scheduling_tt;

class SequencingScheme
{

public:

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 2;
        parameters.swap_block_maximum_length = 1;
        parameters.reverse = false;
        parameters.shift_reverse_block_maximum_length = 0;

        parameters.ruin_and_recreate_number_of_perturbations = 4;
        parameters.ruin_number_of_elements_removed = 4;

        return parameters;
    }

    /**
     * Global cost:
     * - Total tardiness
     */
    using GlobalCost = std::tuple<Time>;

    struct SequenceData
    {
        std::vector<Time> times;
        Time total_tardiness = 0;
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_jobs(); }

    inline SequenceData empty_sequence_data(sequencing::SequenceId) const
    {
        SequenceData sequence_data;
        sequence_data.times = std::vector<Time>(instance_.number_of_machines(), 0);
        return sequence_data;
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {sequence_data.total_tardiness};
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {sequence_data.total_tardiness};
    }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {value};
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId element_id) const
    {
        MachineId m = instance_.number_of_machines();
        // Update times.
        sequence_data.times[0] = sequence_data.times[0]
            + instance_.job(element_id).processing_times[0];
        for (MachineId machine_id = 1; machine_id < m; ++machine_id) {
            if (sequence_data.times[machine_id - 1] > sequence_data.times[machine_id]) {
                sequence_data.times[machine_id] = sequence_data.times[machine_id - 1]
                    + instance_.job(element_id).processing_times[machine_id];
            } else {
                sequence_data.times[machine_id] = sequence_data.times[machine_id]
                    + instance_.job(element_id).processing_times[machine_id];
            }
        }
        // Update total tardiness.
        if (sequence_data.times[m - 1] > instance_.job(element_id).due_date) {
            sequence_data.total_tardiness += (sequence_data.times[m - 1]
                    - instance_.job(element_id).due_date);
        }
    }

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Permutation flow shop scheduling problem, total tardiness" << std::endl;
        instance_.format(os, verbosity_level);
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
