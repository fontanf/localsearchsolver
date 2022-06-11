#include "examples/travelingrepairman.hpp"
#include "localsearchsolver/read_args.hpp"

using namespace localsearchsolver;
using namespace travelingrepairman;

inline SequencingScheme::Parameters read_local_scheme_args(
        const std::vector<char*> argv)
{
    SequencingScheme::Parameters parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("shift-block-maximum-length,", boost::program_options::value<LocationPos>(&parameters.sequencing_parameters.shift_block_maximum_length), "")
        ("swap-block-maximum-length,", boost::program_options::value<LocationPos>(&parameters.sequencing_parameters.swap_block_maximum_length), "")
        ("reverse,", boost::program_options::value<bool>(&parameters.sequencing_parameters.reverse), "")
        ("shift-reverse-block-maximum-length,", boost::program_options::value<LocationPos>(&(parameters.sequencing_parameters.shift_reverse_block_maximum_length)), "")
        ("double-bridge-number-of-perturbations,", boost::program_options::value<LocationPos>(&parameters.sequencing_parameters.double_bridge_number_of_perturbations), "")
        ("ruin-and-recreate-number-of-perturbations,", boost::program_options::value<LocationPos>(&parameters.sequencing_parameters.ruin_and_recreate_number_of_perturbations), "")
        ("ruin-and-recreate-number-of-elements-removed,", boost::program_options::value<LocationPos>(&parameters.sequencing_parameters.ruin_and_recreate_number_of_elements_removed), "")
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
    MainArgs main_args;
    main_args.algorithm = "best_first_local_search";
    main_args.initial_solution_ids = {2};
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
    auto parameters_sequencing_scheme = read_local_scheme_args(main_args.local_scheme_argv);
    SequencingScheme sequencing_scheme(instance, parameters_sequencing_scheme);
    sequencing::LocalScheme<SequencingScheme> local_scheme(sequencing_scheme, parameters_sequencing_scheme.sequencing_parameters);

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
    std::string certificate_path = main_args.info.output->certificate_path;
    if (!certificate_path.empty()) {
        std::ofstream file(certificate_path);
        if (!file.good()) {
            throw std::runtime_error(
                    "Unable to open file \"" + certificate_path + "\".");
        }
        std::vector<std::vector<LocationId>> solution;
        for (auto se: solution_pool.best().sequences[0].elements)
            file << se.j + 1 << std::endl;
    }

    if (main_args.print_solution > 0) {
        os << std::endl
            << "Solution" << std::endl
            << "--------" << std::endl;
        local_scheme.print(os, solution_pool.best());
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
