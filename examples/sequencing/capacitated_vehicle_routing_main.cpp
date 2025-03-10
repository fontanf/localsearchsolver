/**
 * Capacitated vehicle routing problem
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/capacitated_vehicle_routing.hpp
 *
 */

#include "read_args.hpp"

#include "localsearchsolver/sequencing.hpp"

#include "orproblems/routing/capacitated_vehicle_routing.hpp"

using namespace localsearchsolver;
using namespace orproblems::capacitated_vehicle_routing;

template <typename Distances>
class SequencingScheme
{

public:

    static sequencing::Parameters sequencing_parameters()
    {
        sequencing::Parameters parameters;

        parameters.shift_block_maximum_length = 3;
        parameters.swap_block_maximum_length = 3;
        parameters.reverse = true;
        parameters.shift_reverse_block_maximum_length = 3;

        parameters.inter_shift_block_maximum_length = 3;
        parameters.inter_swap_block_maximum_length = 3;
        parameters.swap_tails = true;
        parameters.split = true;
        parameters.inter_shift_reverse_block_maximum_length = 3;
        //parameters.inter_swap_star = true;

        parameters.ruin_and_recreate_number_of_perturbations = 128;
        parameters.ruin_number_of_elements_removed = 10;
        parameters.ruin_adjacent_string_removal_weight = 1.0;
        parameters.recreate_best_weight = 1.0;

        //parameters.order_crossover_weight = 1.0;
        parameters.selective_route_exchange_crossover_1_weight = 1.0;

        return parameters;
    }

    /**
     * Global cost:
     * - Overcapacity
     * - Total distance
     */
    using GlobalCost = std::tuple<Demand, Distance>;

    struct SequenceData
    {
        sequencing::ElementId element_id_first = -1;
        sequencing::ElementId element_id_last = -1;
        Demand demand = 0;
        Distance total_distance = 0;  // Without depot -> element_id_first.
    };

    SequencingScheme(
            const Instance& instance,
            const Distances& distances):
        instance_(instance),
        distances_(distances) { }

    inline sequencing::SequencePos number_of_sequences() const
    {
        // Safety margin: 30% + 3 more vehicles than the trivial bin packing LB
        return std::ceil(1.3 * instance_.total_demand() / instance_.capacity()) + 3;
    }

    inline sequencing::ElementPos number_of_elements() const
    {
        // -1 since we don't schedule the depot.
        return instance_.number_of_locations() - 1;
    }

    inline double distance(
            sequencing::ElementId element_id_1,
            sequencing::ElementId element_id_2) const
    {
        return distances_.distance(
                element_id_1 + 1,
                element_id_2 + 1);
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        if (sequence_data.element_id_first == -1)
            return {0, 0};
        return {
            std::max((Demand)0, sequence_data.demand - instance_.capacity()),
            distances_.distance(0, sequence_data.element_id_first + 1)
                + sequence_data.total_distance
                + distances_.distance(sequence_data.element_id_last + 1, 0),
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
            std::max((Demand)0, sequence_data.demand - instance_.capacity()),
            sequence_data.total_distance
                + distances_.distance(0, sequence_data.element_id_first + 1),
        };
    }

    inline SequenceData sequence_data_init(
            sequencing::ElementId element_id) const
    {
        SequenceData sequence_data;
        // Uppdate element_id_first.
        sequence_data.element_id_first = element_id;
        // Update demand.
        sequence_data.demand = instance_.demand(element_id + 1);
        // Update element_id_last.
        sequence_data.element_id_last = element_id;
        return sequence_data;
    }

    inline bool concatenate(
            SequenceData& sequence_data,
            const SequenceData& sequence_data_2) const
    {
        // Update total_completion_time.
        sequence_data.total_distance
            += distances_.distance(
                    sequence_data.element_id_last + 1,
                    sequence_data_2.element_id_first + 1)
            + sequence_data_2.total_distance;
        // Update demand.
        sequence_data.demand += sequence_data_2.demand;
        // Update element_id_last.
        sequence_data.element_id_last = sequence_data_2.element_id_last;
        return true;
    }

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Capacitated vehicle routing problem" << std::endl;
        instance_.format(os, verbosity_level);
    }

    void solution_write(
            const typename sequencing::LocalScheme<SequencingScheme<Distances>>::Solution& solution,
            std::string certificate_path) const
    {
        if (certificate_path.empty())
            return;
        std::ofstream file(certificate_path);
        if (!file.good()) {
            throw std::runtime_error(
                    "Unable to open file \"" + certificate_path + "\".");
        }

        file << solution.sequences.size() << std::endl;
        for (const auto& sequence: solution.sequences) {
            file << sequence.elements.size() << std::endl;
            for (const auto& se: sequence.elements)
                file << se.element_id + 1 << " ";
            file << std::endl;
        }
    }

private:

    /** Instance. */
    const Instance& instance_;

    /** Distances. */
    const Distances& distances_;

};

template <typename Distances>
void run(
        const Distances& distances,
        const boost::program_options::variables_map& vm,
        const Instance& instance)
{
    // Create local scheme.
    SequencingScheme<Distances> sequencing_scheme(instance, distances);
    auto sequencing_parameters = read_sequencing_args<SequencingScheme<Distances>>(vm);
    sequencing::LocalScheme<SequencingScheme<Distances>> local_scheme(
            sequencing_scheme,
            sequencing_parameters);

    // Run algorithm.
    std::string algorithm = vm["algorithm"].as<std::string>();
    if (algorithm == "multi-start-local-search") {
        run_multi_start_local_search(local_scheme, vm);
    } else if (algorithm == "iterated-local-search") {
        run_iterated_local_search(local_scheme, vm);
    } else if (algorithm == "best-first-local-search") {
        run_best_first_local_search(local_scheme, vm);
    } else {
        run_genetic_local_search(local_scheme, vm);
    }
}

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
    instance.distances().compute_distances_explicit();

    FUNCTION_WITH_DISTANCES(
            run,
            instance.distances(),
            vm,
            instance);

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
