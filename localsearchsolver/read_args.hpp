#pragma once

#include "localsearchsolver/restarting_local_search.hpp"
#include "localsearchsolver/iterated_local_search.hpp"
#include "localsearchsolver/best_first_local_search.hpp"
#include "localsearchsolver/genetic_local_search.hpp"
#include "localsearchsolver/sequencing.hpp"

#include <boost/program_options.hpp>

namespace localsearchsolver
{

template <typename LocalScheme>
inline RestartingLocalSearchOptionalParameters<LocalScheme> read_restarting_local_search_args(
        const std::vector<char*> argv)
{
    RestartingLocalSearchOptionalParameters<LocalScheme> parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("seed,s", boost::program_options::value<Seed>(&parameters.seed), "")
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

template <typename LocalScheme>
inline IteratedLocalSearchOptionalParameters<LocalScheme> read_iterated_local_search_args(
        const std::vector<char*> argv)
{
    IteratedLocalSearchOptionalParameters<LocalScheme> parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("seed,s", boost::program_options::value<Seed>(&parameters.seed), "")
        ("maximum-number-of-restarts,r", boost::program_options::value<Counter>(&parameters.maximum_number_of_restarts), "")
        ("maximum-number-of-iterations,n", boost::program_options::value<Counter>(&parameters.maximum_number_of_iterations), "")
        ("minimum-number-of-perturbations,p", boost::program_options::value<Counter>(&parameters.minimum_number_of_perturbations), "")
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

template <typename LocalScheme>
inline BestFirstLocalSearchOptionalParameters<LocalScheme> read_best_first_local_search_args(
        const std::vector<char*> argv)
{
    BestFirstLocalSearchOptionalParameters<LocalScheme> parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("seed,s", boost::program_options::value<Seed>(&parameters.seed), "")
        (",x", boost::program_options::value<Seed>(&parameters.number_of_threads_1), "")
        (",y", boost::program_options::value<Seed>(&parameters.number_of_threads_2), "")
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

template <typename LocalScheme>
inline GeneticLocalSearchOptionalParameters<LocalScheme> read_genetic_local_search_args(
        const std::vector<char*> argv)
{
    GeneticLocalSearchOptionalParameters<LocalScheme> parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("seed,s", boost::program_options::value<Seed>(&parameters.seed), "")
        (",t", boost::program_options::value<Seed>(&parameters.number_of_threads), "")
        ("population-size,", boost::program_options::value<Counter>(&parameters.maximum_size_of_the_population), "")
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

struct MainArgs
{
    std::string instance_path = "";
    std::string format = "";
    std::string algorithm = "iterated_local_search";
    std::vector<std::string> algorithm_args;
    std::vector<char*> algorithm_argv;
    std::string local_scheme = "local_scheme";
    std::vector<std::string> local_scheme_args;
    std::vector<char*> local_scheme_argv;
    std::string sequencing = "sequencing";
    std::vector<std::string> sequencing_args;
    std::vector<char*> sequencing_argv;
    optimizationtools::Info info = optimizationtools::Info();
    int print_instance = 1;
    int print_solution = 1;
    int print_checker = 1;
    bool has_goal = false;
    double goal = 0;
    std::vector<Counter> initial_solution_ids = {0};
};

MainArgs read_args(int argc, char *argv[], MainArgs& main_args)
{
    std::string output_path = "";
    std::string certificate_path = "";
    double time_limit = std::numeric_limits<double>::infinity();
    int verbosity_level = 1;

    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i", boost::program_options::value<std::string>(&main_args.instance_path)->required(), "set input path (required)")
        ("output,o", boost::program_options::value<std::string>(&output_path), "set JSON output path")
        ("certificate,c", boost::program_options::value<std::string>(&certificate_path), "set certificate path")
        ("format,f", boost::program_options::value<std::string>(&main_args.format), "set input file format (default: orlibrary)")
        ("algorithm,a", boost::program_options::value<std::string>(&main_args.algorithm), "set algorithm")
        ("local-scheme,s", boost::program_options::value<std::string>(&main_args.local_scheme), "set localscheme parameters")
        ("sequencing,", boost::program_options::value<std::string>(&main_args.sequencing), "set sequencing parameters")
        ("initial-solutions,", boost::program_options::value<std::vector<Counter>>(&main_args.initial_solution_ids)->multitoken(), "")
        ("time-limit,t", boost::program_options::value<double>(&time_limit), "Time limit in seconds\n  ex: 3600")
        ("goal,g", boost::program_options::value<double>(&main_args.goal), "set goal")
        ("only-write-at-the-end,e", "Only write output and certificate files at the end")
        ("verbosity-level,v", boost::program_options::value<int>(&verbosity_level), "set verbosity level")
        ("print-instance", boost::program_options::value<int>(&main_args.print_instance), "print instance")
        ("print-solution", boost::program_options::value<int>(&main_args.print_solution), "print solution")
        ("print-checker", boost::program_options::value<int>(&main_args.print_checker), "print checker")
        ;
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

    main_args.has_goal = (vm.count("goal"));

    main_args.algorithm_args = boost::program_options::split_unix(main_args.algorithm);
    for (std::string& s: main_args.algorithm_args)
        main_args.algorithm_argv.push_back(const_cast<char*>(s.c_str()));

    main_args.local_scheme_args = boost::program_options::split_unix(main_args.local_scheme);
    for (std::string& s: main_args.local_scheme_args)
        main_args.local_scheme_argv.push_back(const_cast<char*>(s.c_str()));

    main_args.sequencing_args = boost::program_options::split_unix(main_args.sequencing);
    for (std::string& s: main_args.sequencing_args)
        main_args.sequencing_argv.push_back(const_cast<char*>(s.c_str()));

    main_args.info = optimizationtools::Info()
        .set_verbosity_level(verbosity_level)
        .set_time_limit(time_limit)
        .set_certificate_path(certificate_path)
        .set_json_output_path(output_path)
        .set_only_write_at_the_end(vm.count("only-write-at-the-end"))
        .set_only_write_at_the_end(true)
        .set_sigint_handler()
        ;

    return main_args;
}

template <typename SequencingScheme>
inline sequencing::Parameters read_sequencing_args(
        const std::vector<char*> argv)
{
    sequencing::Parameters parameters = SequencingScheme::sequencing_parameters();
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("shift-block-maximum-length,", boost::program_options::value<sequencing::ElementPos>(&parameters.shift_block_maximum_length), "")
        ("swap-block-maximum-length,", boost::program_options::value<sequencing::ElementPos>(&parameters.swap_block_maximum_length), "")
        ("reverse,", boost::program_options::value<bool>(&parameters.reverse), "")
        ("shift-reverse-block-maximum-length,", boost::program_options::value<sequencing::ElementPos>(&(parameters.shift_reverse_block_maximum_length)), "")
        ("double-bridge-number-of-perturbations,", boost::program_options::value<sequencing::ElementPos>(&parameters.double_bridge_number_of_perturbations), "")
        ("ruin-and-recreate-number-of-perturbations,", boost::program_options::value<sequencing::ElementPos>(&parameters.ruin_and_recreate_number_of_perturbations), "")
        ("ruin-and-recreate-number-of-elements-removed,", boost::program_options::value<sequencing::ElementPos>(&parameters.ruin_number_of_elements_removed), "")
        ("crossover-ox-weight,", boost::program_options::value<double>(&parameters.crossover_ox_weight), "")
        ("crossover-sjox-weight,", boost::program_options::value<double>(&parameters.crossover_sjox_weight), "")
        ("crossover-sbox-weight,", boost::program_options::value<double>(&parameters.crossover_sbox_weight), "")
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

template <typename LocalScheme>
SolutionPool<LocalScheme> run_restarting_local_search(
        const MainArgs& main_args,
        LocalScheme& local_scheme,
        const optimizationtools::Info& info)
{
    auto parameters = read_restarting_local_search_args<LocalScheme>(main_args.algorithm_argv);
    parameters.info = info;
    parameters.initial_solution_ids = main_args.initial_solution_ids;
    if (main_args.has_goal) {
        parameters.has_goal = true;
        parameters.goal = global_cost_goal(local_scheme, main_args.goal);
    }
    return restarting_local_search(local_scheme, parameters).solution_pool;
}

template <typename LocalScheme>
SolutionPool<LocalScheme> run_iterated_local_search(
        const MainArgs& main_args,
        LocalScheme& local_scheme,
        const optimizationtools::Info& info)
{
    auto parameters = read_iterated_local_search_args<LocalScheme>(main_args.algorithm_argv);
    parameters.info = info;
    parameters.initial_solution_ids = main_args.initial_solution_ids;
    if (main_args.has_goal) {
        parameters.has_goal = true;
        parameters.goal = global_cost_goal(local_scheme, main_args.goal);
    }
    return iterated_local_search(local_scheme, parameters).solution_pool;
}

template <typename LocalScheme>
SolutionPool<LocalScheme> run_best_first_local_search(
        const MainArgs& main_args,
        LocalScheme& local_scheme,
        const optimizationtools::Info& info)
{
    auto parameters = read_best_first_local_search_args<LocalScheme>(main_args.algorithm_argv);
    parameters.info = info;
    parameters.initial_solution_ids = main_args.initial_solution_ids;
    if (main_args.has_goal) {
        parameters.has_goal = true;
        parameters.goal = global_cost_goal(local_scheme, main_args.goal);
    }
    return best_first_local_search(local_scheme, parameters).solution_pool;
}

template <typename LocalScheme>
SolutionPool<LocalScheme> run_genetic_local_search(
        const MainArgs& main_args,
        LocalScheme& local_scheme,
        const optimizationtools::Info& info)
{
    auto parameters = read_genetic_local_search_args<LocalScheme>(main_args.algorithm_argv);
    parameters.info = info;
    parameters.initial_solution_ids = main_args.initial_solution_ids;
    if (main_args.has_goal) {
        parameters.has_goal = true;
        parameters.goal = global_cost_goal(local_scheme, main_args.goal);
    }
    return genetic_local_search(local_scheme, parameters).solution_pool;
}

}

