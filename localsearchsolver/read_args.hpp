#pragma once

#include "localsearchsolver/a_star_local_search.hpp"
#include "localsearchsolver/genetic_local_search.hpp"

#include <boost/program_options.hpp>

namespace localsearchsolver
{

template <typename LocalScheme>
inline AStarLocalSearchOptionalParameters<LocalScheme> read_astar_args(
        const std::vector<char*> argv)
{
    AStarLocalSearchOptionalParameters<LocalScheme> parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("seed,s", boost::program_options::value<Seed>(&parameters.seed), "")
        (",x", boost::program_options::value<Seed>(&parameters.thread_number_1), "")
        (",y", boost::program_options::value<Seed>(&parameters.thread_number_2), "")
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
inline GeneticLocalSearchOptionalParameters<LocalScheme> read_genetic_args(
        const std::vector<char*> argv)
{
    GeneticLocalSearchOptionalParameters<LocalScheme> parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("seed,s", boost::program_options::value<Seed>(&parameters.seed), "")
        (",t", boost::program_options::value<Seed>(&parameters.thread_number), "")
        (",i", boost::program_options::value<std::vector<Counter>>(&parameters.initial_solution_ids)->multitoken(), "")
        (",s", boost::program_options::value<Counter>(&parameters.population_size_max), "")
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

