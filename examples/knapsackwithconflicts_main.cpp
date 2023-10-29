#include "examples/knapsackwithconflicts.hpp"
#include "localsearchsolver/read_args.hpp"

using namespace localsearchsolver;
using namespace knapsackwithconflicts;

inline LocalScheme::Parameters read_local_scheme_args(
        const std::vector<char*> argv)
{
    LocalScheme::Parameters parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        //("seed,s", boost::program_options::value<Seed>(&parameters.seed), "")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line((Counter)argv.size(), argv.data(), desc), vm);
    try {
        boost::program_options::notify(vm);
    } catch (const boost::program_options::required_option& e) {
        std::cout << desc << std::endl;;
        throw "";
    }
    return parameters;
}

int main(int argc, char *argv[])
{
    MainArgs main_args;
    main_args.algorithm = "best-first-local-search";
    read_args(argc, argv, main_args);
    auto& os = main_args.info.os();

    // Create instance.
    Instance instance(main_args.instance_path, main_args.format);
    if (main_args.print_instance > 0) {
        os
            << "Instance" << std::endl
            << "--------" << std::endl;
        instance.print(os, main_args.print_instance);
        os << std::endl;
    }

    // Create local scheme.
    auto parameters = read_local_scheme_args(main_args.local_scheme_argv);
    LocalScheme local_scheme(instance, parameters);

    // Run algorithm.
    auto solution_pool =
        (strcmp(main_args.algorithm_argv[0], "multi-start-local-search") == 0)?
        run_multi_start_local_search(main_args, local_scheme, main_args.info):
        (strcmp(main_args.algorithm_argv[0], "iterated-local-search") == 0)?
        run_iterated_local_search(main_args, local_scheme, main_args.info):
        (strcmp(main_args.algorithm_argv[0], "best-first-local-search") == 0)?
        run_best_first_local_search(main_args, local_scheme, main_args.info):
        run_genetic_local_search(main_args, local_scheme, main_args.info);

    // Write solution.
    std::string certificate_path = main_args.info.output->certificate_path;
    local_scheme.write(solution_pool.best(), certificate_path);
    if (main_args.print_solution > 0) {
        os << std::endl
            << "Solution" << std::endl
            << "--------" << std::endl;
        local_scheme.print(
                os,
                solution_pool.best(),
                main_args.print_solution);
    }

    // Run checker.
    if (main_args.print_checker > 0
            && certificate_path != "") {
        os << std::endl
            << "Checker" << std::endl
            << "-------" << std::endl;
        instance.check(
                certificate_path,
                os,
                main_args.print_checker);
    }

    return 0;
}
