#include "examples/capacitatedvehiclerouting.hpp"
#include "localsearchsolver/best_first_local_search.hpp"

using namespace localsearchsolver;
using namespace capacitatedvehiclerouting;

int main(int argc, char *argv[])
{
    (void)argc;
    std::string instance_path = argv[1];
    Instance instance(instance_path);

    // Create local scheme.
    SequencingScheme sequencing_scheme(instance);
    auto sequencing_parameters = SequencingScheme::sequencing_parameters();;
    sequencing::LocalScheme<SequencingScheme> local_scheme(sequencing_scheme, sequencing_parameters);

    BestFirstLocalSearchOptionalParameters<sequencing::LocalScheme<SequencingScheme>> parameters_best_first_local_search;
    parameters_best_first_local_search.info.set_verbose(false);
    parameters_best_first_local_search.initial_solution_ids = {3};
    parameters_best_first_local_search.new_solution_callback
        = [&sequencing_scheme](
                const sequencing::LocalScheme<SequencingScheme>::Solution& solution)
        {
            sequencing::SequenceId route_id = 0;
            for (sequencing::SequenceId i = 0; i < sequencing_scheme.number_of_sequences(); ++i) {
                if (solution.sequences[i].elements.size() == 0)
                    continue;
                std::cout << "Route #" << route_id << ":";
                for (auto se: solution.sequences[i].elements)
                    std::cout << " " << se.j + 1;
                std::cout << std::endl;
                route_id++;
            }
            Distance obj = std::get<1>(solution.global_cost);
            std::cout << "Cost " << (double)obj * 2 << std::endl;
        };
    auto output = best_first_local_search(local_scheme, parameters_best_first_local_search);
    best_first_local_search(local_scheme, parameters_best_first_local_search);

    return 0;
}
