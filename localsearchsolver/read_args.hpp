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

}

