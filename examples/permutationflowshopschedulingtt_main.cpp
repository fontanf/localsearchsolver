#include "examples/permutationflowshopschedulingtt.hpp"
#include "localsearchsolver/read_args.hpp"

using namespace localsearchsolver;
using namespace permutationflowshopschedulingtt;

inline LocalScheme::Parameters read_local_scheme_args(
        const std::vector<char*> argv)
{
    LocalScheme::Parameters parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("shift-bloc-maximum-length,", boost::program_options::value<permutationflowshopschedulingtt::JobPos>(&parameters.sequencing_parameters.shift_bloc_maximum_length), "")
        ("swap-bloc-maximum-length,", boost::program_options::value<permutationflowshopschedulingtt::JobPos>(&parameters.sequencing_parameters.swap_bloc_maximum_length), "")
        ("crossover-ox-weight,", boost::program_options::value<double>(&parameters.sequencing_parameters.crossover_ox_weight), "")
        ("crossover-sjox-weight,", boost::program_options::value<double>(&parameters.sequencing_parameters.crossover_sjox_weight), "")
        ("crossover-sbox-weight,", boost::program_options::value<double>(&parameters.sequencing_parameters.crossover_sbox_weight), "")
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
    auto main_args = read_args(argc, argv);

    // Create instance.
    Instance instance(main_args.instance_path, main_args.format);
    if (main_args.print_instance)
        std::cout << instance << std::endl;

    // Create local scheme.
    auto parameters_local_scheme_0 = read_local_scheme_args(main_args.local_scheme_argv);
    LocalScheme local_scheme_0(instance, parameters_local_scheme_0);
    sequencing::LocalScheme<LocalScheme> local_scheme(local_scheme_0, parameters_local_scheme_0.sequencing_parameters);

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
