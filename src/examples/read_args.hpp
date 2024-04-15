#pragma once

#include "localsearchsolver/multi_start_local_search.hpp"
#include "localsearchsolver/iterated_local_search.hpp"
#include "localsearchsolver/best_first_local_search.hpp"
#include "localsearchsolver/genetic_local_search.hpp"
#include "localsearchsolver/sequencing.hpp"

#include <boost/program_options.hpp>

namespace localsearchsolver
{

boost::program_options::options_description setup_args()
{
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i", boost::program_options::value<std::string>()->required(), "set input path (required)")
        ("format,f", boost::program_options::value<std::string>()->default_value(""), "set input file format")
        ("output,o", boost::program_options::value<std::string>()->default_value(""), "set JSON output path")
        ("certificate,c", boost::program_options::value<std::string>()->default_value(""), "set certificate path")
        ("algorithm,a", boost::program_options::value<std::string>(), "set algorithm")
        ("initial-solutions,", boost::program_options::value<std::vector<Counter>>()->multitoken(), "")
        ("time-limit,t", boost::program_options::value<double>(), "Time limit in seconds\n  ex: 3600")
        ("goal,g", boost::program_options::value<double>(), "set goal")
        ("only-write-at-the-end,e", "Only write output and certificate files at the end")
        ("verbosity-level,v", boost::program_options::value<int>(), "set verbosity level")
        ("print-checker", boost::program_options::value<int>()->default_value(1), "print checker")
        ("seed,s", boost::program_options::value<Seed>(), "")

        ("maximum-number-of-nodes", boost::program_options::value<Counter>(), "")
        ("maximum-number-of-restarts", boost::program_options::value<Counter>(), "")
        ("maximum-number-of-iterations", boost::program_options::value<Counter>(), "")
        ("minimum-number-of-perturbations", boost::program_options::value<Counter>(), "")
        ("number-of-threads", boost::program_options::value<Seed>(), "")
        ("number-of-threads-1", boost::program_options::value<Seed>(), "")
        ("number-of-threads-2", boost::program_options::value<Seed>(), "")
        ("maximum-size-of-the-population", boost::program_options::value<Counter>(), "")

        ("reverse", boost::program_options::value<bool>(), "")
        ;
    return desc;
}

template <typename LocalScheme>
void read_args(
        const LocalScheme& local_scheme,
        Parameters<LocalScheme>& parameters,
        const boost::program_options::variables_map& vm)
{
    parameters.timer.set_sigint_handler();
    parameters.messages_to_stdout = true;
    if (vm.count("time-limit"))
        parameters.timer.set_time_limit(vm["time-limit"].as<double>());
    if (vm.count("verbosity-level"))
        parameters.verbosity_level = vm["verbosity-level"].as<int>();
    if (vm.count("log"))
        parameters.log_path = vm["log"].as<std::string>();
    parameters.log_to_stderr = vm.count("log-to-stderr");

    parameters.has_goal = (vm.count("goal"));
    if (vm.count("goal"))
        parameters.goal = global_cost_goal(local_scheme, vm["goal"].as<double>());

    bool only_write_at_the_end = vm.count("only-write-at-the-end");
    if (!only_write_at_the_end) {
        std::string certificate_path = vm["certificate"].as<std::string>();
        std::string json_output_path = vm["output"].as<std::string>();
        parameters.new_solution_callback = [
            local_scheme,
            json_output_path,
            certificate_path](
                    const Output<LocalScheme>& output)
        {
            output.write_json_output(json_output_path);
            solution_write(
                    local_scheme,
                    output.solution_pool.best(),
                    certificate_path);
        };
    }
}

void setup_sequencing_args(
        boost::program_options::options_description& desc)
{
    desc.add_options()
        ("shift-block-maximum-length,", boost::program_options::value<sequencing::ElementPos>(), "")
        ("swap-block-maximum-length,", boost::program_options::value<sequencing::ElementPos>(), "")
        ("reverse,", boost::program_options::value<bool>(), "")
        ("shift-reverse-block-maximum-length,", boost::program_options::value<sequencing::ElementPos>(), "")
        ("double-bridge-number-of-perturbations,", boost::program_options::value<sequencing::ElementPos>(), "")
        ("ruin-and-recreate-number-of-perturbations,", boost::program_options::value<sequencing::ElementPos>(), "")
        ("ruin-and-recreate-number-of-elements-removed,", boost::program_options::value<sequencing::ElementPos>(), "")
        ("ox-weight,", boost::program_options::value<double>(), "")
        ("sjox-weight,", boost::program_options::value<double>(), "")
        ("sbox-weight,", boost::program_options::value<double>(), "")
        ;
}

template <typename SequencingScheme>
inline sequencing::Parameters read_sequencing_args(
        const boost::program_options::variables_map& vm)
{
    sequencing::Parameters parameters = SequencingScheme::sequencing_parameters();

    if (vm.count("shift-block-maximum-length"))
        parameters.shift_block_maximum_length = vm["shift-block-maximum-length"].as<sequencing::ElementPos>();
    if (vm.count("swap-block-maximum-length"))
        parameters.swap_block_maximum_length = vm["swap-block-maximum-length"].as<sequencing::ElementPos>();
    if (vm.count("reverse"))
        parameters.reverse = vm["reverse"].as<bool>();
    if (vm.count("shift-reverse-block-maximum-length"))
        parameters.shift_reverse_block_maximum_length = vm["shift-reverse-block-maximum-length"].as<sequencing::ElementPos>();

    if (vm.count("double-bridge-number-of-perturbations"))
        parameters.double_bridge_number_of_perturbations = vm["double-bridge-number-of-perturbations"].as<sequencing::ElementPos>();
    if (vm.count("ruin-and-recreate-number-of-perturbations"))
        parameters.ruin_and_recreate_number_of_perturbations = vm["ruin-and-recreate-number-of-perturbations"].as<sequencing::ElementPos>();
    if (vm.count("ruin-number-of-elements-removed"))
        parameters.ruin_number_of_elements_removed = vm["ruin-number-of-elements-removed"].as<sequencing::ElementPos>();

    if (vm.count("order-crossover-weight"))
        parameters.order_crossover_weight = vm["order-crossover-weight"].as<double>();

    return parameters;
}

template <typename LocalScheme>
void write_outputs(
        const Output<LocalScheme>& output,
        const boost::program_options::variables_map& vm)
{
    std::string json_output_path = vm["output"].as<std::string>();
    std::string certificate_path = vm["certificate"].as<std::string>();
    output.write_json_output(json_output_path);
    solution_write(
            output.solution_pool.local_scheme(),
            output.solution_pool.best(),
            certificate_path);
}

template <typename LocalScheme>
const Output<LocalScheme> run_multi_start_local_search(
        LocalScheme& local_scheme,
        const boost::program_options::variables_map& vm)
{
    MultiStartLocalSearchParameters<LocalScheme> parameters;
    read_args(local_scheme, parameters, vm);
    if (vm.count("maximum-number-of-restarts"))
        parameters.maximum_number_of_restarts = vm["maximum-number-of-restarts"].as<Counter>();
    const Output<LocalScheme> output = multi_start_local_search(local_scheme, parameters);
    write_outputs(output, vm);
    return output;
}

template <typename LocalScheme>
const Output<LocalScheme> run_iterated_local_search(
        LocalScheme& local_scheme,
        const boost::program_options::variables_map& vm)
{
    IteratedLocalSearchParameters<LocalScheme> parameters;
    read_args(local_scheme, parameters, vm);
    if (vm.count("maximum-number-of-restarts"))
        parameters.maximum_number_of_restarts = vm["maximum-number-of-restarts"].as<Counter>();
    if (vm.count("maximum-number-of-iterations"))
        parameters.maximum_number_of_iterations = vm["maximum-number-of-iterations"].as<Counter>();
    if (vm.count("minimum-number-of-perturbations"))
        parameters.minimum_number_of_perturbations = vm["minimum-number-of-perturbations"].as<Counter>();
    const Output<LocalScheme> output = iterated_local_search(local_scheme, parameters);
    write_outputs(output, vm);
    return output;
}

template <typename LocalScheme>
const Output<LocalScheme> run_best_first_local_search(
        LocalScheme& local_scheme,
        const boost::program_options::variables_map& vm)
{
    BestFirstLocalSearchParameters<LocalScheme> parameters;
    read_args(local_scheme, parameters, vm);
    if (vm.count("number-of-threads-1"))
        parameters.number_of_threads_1 = vm["number-of-threads-1"].as<Counter>();
    if (vm.count("number-of-threads-2"))
        parameters.number_of_threads_2 = vm["number-of-threads-2"].as<Counter>();
    if (vm.count("maximum-number-of-nodes"))
        parameters.maximum_number_of_nodes = vm["maximum-number-of-nodes"].as<Counter>();
    const Output<LocalScheme> output = best_first_local_search(local_scheme, parameters);
    write_outputs(output, vm);
    return output;
}

template <typename LocalScheme>
const Output<LocalScheme> run_genetic_local_search(
        LocalScheme& local_scheme,
        const boost::program_options::variables_map& vm)
{
    GeneticLocalSearchParameters<LocalScheme> parameters;
    read_args(local_scheme, parameters, vm);
    if (vm.count("number-of-threads"))
        parameters.number_of_threads = vm["number-of-threads"].as<Counter>();
    if (vm.count("maximum-number-of-iterations"))
        parameters.maximum_number_of_iterations = vm["maximum-number-of-iterations"].as<Counter>();
    if (vm.count("maximum-size-of-the-population"))
        parameters.maximum_size_of_the_population = vm["maximum-size-of-the-population"].as<Counter>();
    const Output<LocalScheme> output = genetic_local_search(local_scheme, parameters);
    write_outputs(output, vm);
    return output;
}

}
