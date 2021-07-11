#pragma once

#include "examples/knapsackwithconflicts.hpp"
#include "examples/multidimensionalmultiplechoiceknapsack.hpp"
#include "examples/quadraticassignment.hpp"
#include "examples/travellingsalesman.hpp"
#include "examples/schedulingwithsdsttwt.hpp"
#include "examples/permutationflowshopschedulingmakespan.hpp"
#include "examples/permutationflowshopschedulingtt.hpp"

#include <boost/program_options.hpp>

namespace localsearchsolver
{

inline knapsackwithconflicts::LocalScheme::Parameters read_knapsackwithconflicts_args(
        const std::vector<char*> argv)
{
    knapsackwithconflicts::LocalScheme::Parameters parameters;
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

inline multidimensionalmultiplechoiceknapsack::LocalScheme::Parameters read_multidimensionalmultiplechoiceknapsack_args(
        const std::vector<char*> argv)
{
    multidimensionalmultiplechoiceknapsack::LocalScheme::Parameters parameters;
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

inline quadraticassignment::LocalScheme::Parameters read_quadraticassignment_args(
        const std::vector<char*> argv)
{
    quadraticassignment::LocalScheme::Parameters parameters;
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

inline travellingsalesman::LocalScheme::Parameters read_travellingsalesman_args(
        const std::vector<char*> argv)
{
    travellingsalesman::LocalScheme::Parameters parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("oropt-size-max,o", boost::program_options::value<travellingsalesman::VertexPos>(&parameters.oropt_size_max), "")
        ("twoopt,t", boost::program_options::value<bool>(&parameters.twoopt), "")
        ("swap,s", boost::program_options::value<bool>(&parameters.swap), "")
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

inline schedulingwithsdsttwt::LocalScheme::Parameters read_schedulingwithsdsttwt_args(
        const std::vector<char*> argv)
{
    schedulingwithsdsttwt::LocalScheme::Parameters parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("bloc-size-max,b", boost::program_options::value<travellingsalesman::VertexPos>(&parameters.bloc_size_max), "")
        ("swap,s", boost::program_options::value<bool>(&parameters.swap), "")
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

inline permutationflowshopschedulingmakespan::LocalScheme::Parameters read_permutationflowshopschedulingmakespan_args(
        const std::vector<char*> argv)
{
    permutationflowshopschedulingmakespan::LocalScheme::Parameters parameters;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("bloc-size-max,b", boost::program_options::value<travellingsalesman::VertexPos>(&parameters.bloc_size_max), "")
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

inline permutationflowshopschedulingtt::LocalScheme::Parameters read_permutationflowshopschedulingtt_args(
        const std::vector<char*> argv)
{
    permutationflowshopschedulingtt::LocalScheme::Parameters parameters;
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

}
