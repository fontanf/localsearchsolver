/**
 * Time-dependent orienteering problem
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/time_dependent_orienteering.hpp
 *
 */

#include "read_args.hpp"

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/routing/time_dependent_orienteering.hpp"

using namespace localsearchsolver;
using namespace orproblems::time_dependent_orienteering;

class SequencingScheme
{

public:

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 7;
        parameters.swap_block_maximum_length = 5;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 6;
        parameters.add_remove = true;
        parameters.replace = true;

        parameters.force_add = true;

        return parameters;
    }

    /**
     * Global cost:
     * - Overtime
     * - Profit
     * - Total time
     */
    using GlobalCost = std::tuple<Time, Profit, Time>;

    struct SequenceData
    {
        sequencing::ElementId element_id_last = -1;
        Time time_cur = 0;
        Time time_full = 0;
        Profit profit = 0;
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_locations() - 2; }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            std::max((Time)0, sequence_data.time_full - instance_.maximum_duration()),
            -sequence_data.profit,
            sequence_data.time_full,
        };
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return {
            std::max((Time)0, sequence_data.time_full - instance_.maximum_duration()),
            std::numeric_limits<Profit>::lowest(),
            sequence_data.time_full,
        };
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId element_id) const
    {
        // Update time_cur.
        sequence_data.time_cur = instance_.arrival_time(
                sequence_data.element_id_last + 1,
                element_id + 1,
                sequence_data.time_cur);
        // Update profit.
        sequence_data.profit += instance_.location(element_id + 1).profit;
        // Update time_full.
        sequence_data.time_full = instance_.arrival_time(
                element_id + 1,
                instance_.number_of_locations() - 1,
                sequence_data.time_cur);
        // Update element_id_last.
        sequence_data.element_id_last = element_id;
    }

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Time-dependent orienteering problem" << std::endl;
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

        for (auto se: solution.sequences[0].elements)
            file << se.element_id + 1 << std::endl;
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
