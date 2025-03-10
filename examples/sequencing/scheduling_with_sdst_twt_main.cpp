/**
 * Single machine scheduling problem with sequence-dependent setup times, total
 * weighted tardiness
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/scheduling_with_sdst_twt.hpp
 *
 */

#include "read_args.hpp"

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/scheduling/scheduling_with_sdst_twt.hpp"

using namespace localsearchsolver;
using namespace orproblems::scheduling_with_sdst_twt;

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
