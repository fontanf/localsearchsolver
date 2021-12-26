#include "examples/schedulingwithsdsttwt.hpp"
#include "localsearchsolver/read_args.hpp"

using namespace localsearchsolver;
using namespace schedulingwithsdsttwt;

inline LocalScheme::Parameters read_local_scheme_args(
        const std::vector<char*> argv)
{
    LocalScheme::Parameters parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("shift-block-maximum-length,", boost::program_options::value<JobPos>(&parameters.sequencing_parameters.shift_block_maximum_length), "")
        ("swap-block-maximum-length,", boost::program_options::value<JobPos>(&parameters.sequencing_parameters.swap_block_maximum_length), "")
        ("double-bridge-number-of-pertubrations", boost::program_options::value<JobPos>(&parameters.sequencing_parameters.double_bridge_number_of_perturbations), "")
        ("ruin-and-recreate-number-of-pertubrations", boost::program_options::value<JobPos>(&parameters.sequencing_parameters.ruin_and_recreate_number_of_perturbations), "")
        ("ruin-and-recreate-number-of-elements-removed", boost::program_options::value<JobPos>(&parameters.sequencing_parameters.ruin_and_recreate_number_of_elements_removed), "")
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
    main_args.algorithm = "best_first_local_search";
    read_args(argc, argv, main_args);

    // Create instance.
    Instance instance(main_args.instance_path, main_args.format);
    if (main_args.print_instance)
        std::cout << instance << std::endl;

    // Create local scheme.
    auto parameters_local_scheme_0 = read_local_scheme_args(main_args.local_scheme_argv);
    LocalScheme local_scheme_0(instance, parameters_local_scheme_0);
    sequencing2::LocalScheme<LocalScheme> local_scheme(local_scheme_0, parameters_local_scheme_0.sequencing_parameters);

    // Run algorithm.
    auto solution_pool =
        (strcmp(main_args.algorithm_argv[0], "restarting_local_search") == 0)?
        run_restarting_local_search(main_args.algorithm_args, local_scheme, main_args.info):
        (strcmp(main_args.algorithm_argv[0], "iterated_local_search") == 0)?
        run_iterated_local_search(main_args.algorithm_args, local_scheme, main_args.info):
        (strcmp(main_args.algorithm_argv[0], "best_first_local_search") == 0)?
        run_best_first_local_search(main_args.algorithm_args, local_scheme, main_args.info):
        run_genetic_local_search(main_args.algorithm_args, local_scheme, main_args.info);

    // Write solution.
    local_scheme.write(solution_pool.best(), main_args.info.output->certificate_path);
    if (main_args.print_solution)
        local_scheme.print(std::cout, solution_pool.best());

    // Run checker.
    if (main_args.info.output->certificate_path != "") {
        std::cout << std::endl;
        instance.check(main_args.info.output->certificate_path);
    }

    return 0;
}
