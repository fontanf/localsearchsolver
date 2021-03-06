#include "examples/distributedpfsstct.hpp"
#include "localsearchsolver/read_args.hpp"

using namespace localsearchsolver;
using namespace distributedpfsstct;

int main(int argc, char *argv[])
{
    MainArgs main_args;
    main_args.algorithm = "best_first_local_search";
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
    SequencingScheme sequencing_scheme(instance);
    auto sequencing_parameters = read_sequencing_args<SequencingScheme>(main_args.sequencing_argv);
    sequencing::LocalScheme<SequencingScheme> local_scheme(sequencing_scheme, sequencing_parameters);

    // Run algorithm.
    auto solution_pool =
        (strcmp(main_args.algorithm_argv[0], "restarting_local_search") == 0)?
        run_restarting_local_search(main_args, local_scheme, main_args.info):
        (strcmp(main_args.algorithm_argv[0], "iterated_local_search") == 0)?
        run_iterated_local_search(main_args, local_scheme, main_args.info):
        (strcmp(main_args.algorithm_argv[0], "best_first_local_search") == 0)?
        run_best_first_local_search(main_args, local_scheme, main_args.info):
        run_genetic_local_search(main_args, local_scheme, main_args.info);

    // Write solution.
    local_scheme.write(solution_pool.best(), main_args.info.output->certificate_path);
    if (main_args.print_solution) {
        std::cout << std::endl
            << "Solution" << std::endl
            << "--------" << std::endl;
        local_scheme.print(std::cout, solution_pool.best());
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
