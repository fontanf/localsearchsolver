#include "examples/vehicleroutingwithtimewindows.hpp"
#include "localsearchsolver/read_args.hpp"

using namespace localsearchsolver;
using namespace vehicleroutingwithtimewindows;

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
    Instance instance = instance_builder.build();

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
