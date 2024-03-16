#include "examples/capacitated_vehicle_routing.hpp"
#include "localsearchsolver/best_first_local_search.hpp"

using namespace localsearchsolver;
using namespace capacitated_vehicle_routing;

template <typename Distances>
void run(
        const Distances& distances,
        const Instance& instance)
{
    // Create local scheme.
    SequencingScheme<Distances> sequencing_scheme(instance, distances);
    auto sequencing_parameters = SequencingScheme<Distances>::sequencing_parameters();;
    sequencing::LocalScheme<SequencingScheme<Distances>> local_scheme(sequencing_scheme, sequencing_parameters);

    using LocalScheme = sequencing::LocalScheme<SequencingScheme<Distances>>;
    BestFirstLocalSearchParameters<LocalScheme> bfls_parameters;
    bfls_parameters.verbosity_level = 0;
    bfls_parameters.initial_solution_ids = {3};
    bfls_parameters.new_solution_callback
        = [&sequencing_scheme](
                const Output<LocalScheme>& output)
        {
            const typename LocalScheme::Solution& solution = output.solution_pool.best();
            sequencing::SequenceId route_id = 0;
            for (sequencing::SequenceId sequence_id = 0;
                    sequence_id < sequencing_scheme.number_of_sequences();
                    ++sequence_id) {
                if (solution.sequences[sequence_id].elements.size() == 0)
                    continue;
                std::cout << "Route #" << route_id << ":";
                for (auto se: solution.sequences[sequence_id].elements)
                    std::cout << " " << se.element_id + 1;
                std::cout << std::endl;
                route_id++;
            }
            Distance obj = std::get<1>(solution.global_cost);
            std::cout << "Cost " << (double)obj * 2 << std::endl;
        };
    best_first_local_search(local_scheme, bfls_parameters);
}

int main(int argc, char *argv[])
{
    (void)argc;
    std::string instance_path = argv[1];
    InstanceBuilder instance_builder;
    instance_builder.read(
            instance_path);
    const Instance instance = instance_builder.build();
    instance.distances().compute_distances_explicit();

    FUNCTION_WITH_DISTANCES(
            run,
            instance.distances(),
            instance);

    return 0;
}
