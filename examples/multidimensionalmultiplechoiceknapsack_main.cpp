#include "examples/multidimensionalmultiplechoiceknapsack.hpp"
#include "localsearchsolver/read_args.hpp"

using namespace localsearchsolver;
using namespace multidimensionalmultiplechoiceknapsack;

int main(int argc, char *argv[])
{
    MainArgs main_args;
    main_args.algorithm = "genetic_local_search";
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
    LocalScheme local_scheme(instance);

    // Run algorithm.
    auto solution_pool =
        (strcmp(main_args.algorithm_argv[0], "restarting_local_search") == 0)?
        run_restarting_local_search(main_args, local_scheme, main_args.info):
        run_genetic_local_search(main_args, local_scheme, main_args.info);

    // Write solution.
    local_scheme.write(solution_pool.best(), main_args.info.output->certificate_path);
    if (main_args.print_solution > 0) {
        os << std::endl
            << "Solution" << std::endl
            << "--------" << std::endl;
        local_scheme.print(os, solution_pool.best(), main_args.print_solution);
    }

    // Run checker.
    if (main_args.print_checker > 0
            && main_args.info.output->certificate_path != "") {
        os << std::endl
            << "Checker" << std::endl
            << "-------" << std::endl;
        instance.check(
                main_args.info.output->certificate_path,
                os,
                main_args.print_checker);
    }

    return 0;
}
