#pragma once

#include "localsearchsolver/restarting_local_search.hpp"
#include "localsearchsolver/iterated_local_search.hpp"
#include "localsearchsolver/best_first_local_search.hpp"
#include "localsearchsolver/genetic_local_search.hpp"

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
        (",i", boost::program_options::value<std::vector<Counter>>(&parameters.initial_solution_ids)->multitoken(), "")
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
        (",i", boost::program_options::value<std::vector<Counter>>(&parameters.initial_solution_ids)->multitoken(), "")
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
        (",i", boost::program_options::value<std::vector<Counter>>(&parameters.initial_solution_ids)->multitoken(), "")
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
        (",i", boost::program_options::value<std::vector<Counter>>(&parameters.initial_solution_ids)->multitoken(), "")
        (",s", boost::program_options::value<Counter>(&parameters.maximum_size_of_the_population), "")
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
        const std::vector<std::string>& algorithm_args,
        LocalScheme& local_scheme,
        const optimizationtools::Info& info)
{
    std::vector<char*> algorithm_argv;
    for (Counter i = 0; i < (Counter)algorithm_args.size(); ++i)
        algorithm_argv.push_back(const_cast<char*>(algorithm_args[i].c_str()));

    auto parameters = read_restarting_local_search_args<LocalScheme>(algorithm_argv);
    parameters.info = info;
    if (parameters.initial_solution_ids.empty())
        parameters.initial_solution_ids = {0};
    return restarting_local_search(local_scheme, parameters).solution_pool;
}

template <typename LocalScheme>
SolutionPool<LocalScheme> run_iterated_local_search(
        const std::vector<std::string>& algorithm_args,
        LocalScheme& local_scheme,
        const optimizationtools::Info& info)
{
    std::vector<char*> algorithm_argv;
    for (Counter i = 0; i < (Counter)algorithm_args.size(); ++i)
        algorithm_argv.push_back(const_cast<char*>(algorithm_args[i].c_str()));

    auto parameters = read_iterated_local_search_args<LocalScheme>(algorithm_argv);
    parameters.info = info;
    if (parameters.initial_solution_ids.empty())
        parameters.initial_solution_ids = {0};
    return iterated_local_search(local_scheme, parameters).solution_pool;
}

template <typename LocalScheme>
SolutionPool<LocalScheme> run_best_first_local_search(
        const std::vector<std::string>& algorithm_args,
        LocalScheme& local_scheme,
        const optimizationtools::Info& info)
{
    std::vector<char*> algorithm_argv;
    for (Counter i = 0; i < (Counter)algorithm_args.size(); ++i)
        algorithm_argv.push_back(const_cast<char*>(algorithm_args[i].c_str()));

    auto parameters = read_best_first_local_search_args<LocalScheme>(algorithm_argv);
    parameters.info = info;
    if (parameters.initial_solution_ids.empty())
        parameters.initial_solution_ids = std::vector<Counter>(parameters.number_of_threads_2, 0);
    return best_first_local_search(local_scheme, parameters).solution_pool;
}

template <typename LocalScheme>
SolutionPool<LocalScheme> run_genetic_local_search(
        const std::vector<std::string>& algorithm_args,
        LocalScheme& local_scheme,
        const optimizationtools::Info& info)
{
    std::vector<char*> algorithm_argv;
    for (Counter i = 0; i < (Counter)algorithm_args.size(); ++i)
        algorithm_argv.push_back(const_cast<char*>(algorithm_args[i].c_str()));

    auto parameters = read_genetic_local_search_args<LocalScheme>(algorithm_argv);
    parameters.info = info;
    return genetic_local_search(local_scheme, parameters).solution_pool;
}

struct MainArgs
{
    std::string instance_path = "";
    std::string format = "";
    std::vector<std::string> algorithm_args;
    std::vector<char*> algorithm_argv;
    std::vector<std::string> local_scheme_args;
    std::vector<char*> local_scheme_argv;
    optimizationtools::Info info = optimizationtools::Info();
    bool print_instance = false;
    bool print_solution = false;
};

MainArgs read_args(int argc, char *argv[])
{
    MainArgs main_args;
    std::string output_path = "";
    std::string certificate_path = "";
    std::string algorithm = "iterative_beam_search";
    std::string local_scheme_parameters = "";
    double time_limit = std::numeric_limits<double>::infinity();

    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i", boost::program_options::value<std::string>(&main_args.instance_path)->required(), "set input path (required)")
        ("output,o", boost::program_options::value<std::string>(&output_path), "set JSON output path")
        ("certificate,c", boost::program_options::value<std::string>(&certificate_path), "set certificate path")
        ("format,f", boost::program_options::value<std::string>(&main_args.format), "set input file format (default: orlibrary)")
        ("algorithm,a", boost::program_options::value<std::string>(&algorithm), "set algorithm")
        ("local-scheme,b", boost::program_options::value<std::string>(&local_scheme_parameters), "set localscheme parameters")
        ("time-limit,t", boost::program_options::value<double>(&time_limit), "Time limit in seconds\n  ex: 3600")
        ("only-write-at-the-end,e", "Only write output and certificate files at the end")
        ("verbose,v", "")
        ("print-instance", "")
        ("print-solution", "")
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

    main_args.print_instance = (vm.count("print-instance"));
    main_args.print_solution = (vm.count("print-solution"));

    main_args.algorithm_args = boost::program_options::split_unix(algorithm);
    for (std::string& s: main_args.algorithm_args)
        main_args.algorithm_argv.push_back(const_cast<char*>(s.c_str()));

    main_args.local_scheme_args = boost::program_options::split_unix(local_scheme_parameters);
    std::string dummy = "dummy";
    main_args.local_scheme_argv.push_back(const_cast<char*>(dummy.c_str()));
    for (std::string& s: main_args.local_scheme_args)
        main_args.local_scheme_argv.push_back(const_cast<char*>(s.c_str()));

    main_args.info = optimizationtools::Info()
        .set_verbose(vm.count("verbose"))
        .set_time_limit(time_limit)
        .set_certificate_path(certificate_path)
        .set_json_output_path(output_path)
        .set_only_write_at_the_end(vm.count("only-write-at-the-end"))
        .set_only_write_at_the_end(true)
        .set_sigint_handler()
        ;

    return main_args;
}

}

