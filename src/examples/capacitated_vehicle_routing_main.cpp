#include "localsearchsolver/examples/capacitated_vehicle_routing.hpp"
#include "read_args.hpp"

using namespace localsearchsolver;
using namespace capacitated_vehicle_routing;

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
