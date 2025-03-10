/**
 * Sequential ordering problem
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/sequential_ordering.hpp
 *
 */

#include "read_args.hpp"

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/scheduling/sequential_ordering.hpp"

using namespace localsearchsolver;
using namespace orproblems::sequential_ordering;

class SequencingScheme
{

public:

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 8;
        parameters.swap_block_maximum_length = 2;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 3;

        parameters.ruin_and_recreate_number_of_perturbations = 10;

        return parameters;
    }

    /**
     * Global cost:
     * - Number of precedence violations
     * - Total distance
     */
    using GlobalCost = std::tuple<LocationPos, Distance>;

    struct SequenceData
    {
        sequencing::ElementId element_id_last = -1;
        Distance length = 0;
        LocationPos number_of_precedence_violations = 0;
        std::vector<uint8_t> contains;
    };

    SequencingScheme(const Instance& instance): instance_(instance) { }

    inline sequencing::ElementPos number_of_elements() const { return instance_.number_of_locations(); }

    inline SequenceData empty_sequence_data(sequencing::SequenceId) const
    {
        SequenceData sequence_data;
        sequence_data.contains = std::vector<uint8_t>(
                instance_.number_of_locations(), 0);
        return sequence_data;
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return {
            sequence_data.number_of_precedence_violations,
            sequence_data.length,
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
            sequence_data.number_of_precedence_violations,
            sequence_data.length,
        };
    }

    inline void append(
            SequenceData& sequence_data,
            sequencing::ElementId element_id) const
    {
        // Update number_of_precedence_violations.
        for (LocationId location_id_pred: instance_.predecessors(element_id))
            if (!sequence_data.contains[location_id_pred])
                sequence_data.number_of_precedence_violations++;
        // Update time.
        if (sequence_data.element_id_last != -1) {
            Distance d = instance_.distance(sequence_data.element_id_last, element_id);
            if (d == std::numeric_limits<Distance>::max())
                d = instance_.distance(element_id, sequence_data.element_id_last);
            sequence_data.length += d;
        }
        // Update contains.
        sequence_data.contains[element_id] = 1;
        // Update element_id_last.
        sequence_data.element_id_last = element_id;
    }

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Sequential ordering problem" << std::endl;
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

        for (auto it = solution.sequences[0].elements.begin() + 1;
                it != solution.sequences[0].elements.end();
                ++it) {
            file << it->element_id << std::endl;
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
