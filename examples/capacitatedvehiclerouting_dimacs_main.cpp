#include "examples/capacitatedvehiclerouting.hpp"
#include "localsearchsolver/best_first_local_search.hpp"

using namespace localsearchsolver;
using namespace capacitatedvehiclerouting;

int main(int argc, char *argv[])
{
    (void)argc;
    std::string instance_path = argv[1];
    InstanceBuilder instance_builder;
    instance_builder.read(
            instance_path);
    Instance instance = instance_builder.build();

    // Create local scheme.
    SequencingScheme sequencing_scheme(instance);
    auto sequencing_parameters = SequencingScheme::sequencing_parameters();;
    sequencing::LocalScheme<SequencingScheme> local_scheme(sequencing_scheme, sequencing_parameters);

    using LocalScheme = sequencing::LocalScheme<SequencingScheme>;
    BestFirstLocalSearchParameters<LocalScheme> bfls_parameters;
    bfls_parameters.verbosity_level = 0;
    bfls_parameters.initial_solution_ids = {3};
    bfls_parameters.new_solution_callback
        = [&sequencing_scheme](
                const Output<LocalScheme>& output)
        {
            const LocalScheme::Solution& solution = output.solution_pool.best();
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
    const BestFirstLocalSearchOutput<LocalScheme> output = best_first_local_search(local_scheme, bfls_parameters);

    return 0;
}
