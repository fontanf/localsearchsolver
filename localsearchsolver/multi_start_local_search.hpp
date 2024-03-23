#pragma once

#include "localsearchsolver/algorithm_formatter.hpp"

#include <random>

namespace localsearchsolver
{

template <typename LocalScheme>
struct MultiStartLocalSearchParameters: Parameters<LocalScheme>
{
    using Solution = typename LocalScheme::Solution;
    using GlobalCost = typename LocalScheme::GlobalCost;


    /** Maximum number of restarts. */
    Counter maximum_number_of_restarts = -1;


    virtual nlohmann::json to_json(
            const LocalScheme& local_scheme) const override
    {
        nlohmann::json json = Parameters<LocalScheme>::to_json(local_scheme);
        json.merge_patch({
            {"MaximumNumberOfRestarts", maximum_number_of_restarts}});
        return json;
    }

    virtual int format_width() const override { return 30; }

    virtual void format(
            std::ostream& os,
            const LocalScheme& local_scheme) const override
    {
        Parameters<LocalScheme>::format(os, local_scheme);
        int width = format_width();
        os
            << std::setw(width) << std::left << "Maximum number of restarts: " << maximum_number_of_restarts << std::endl
            ;
    }
};

template <typename LocalScheme>
struct MultiStartLocalSearchOutput: Output<LocalScheme>
{
    /** Constructor. */
    MultiStartLocalSearchOutput(
            const LocalScheme& local_scheme,
            Counter maximum_size_of_the_solution_pool):
        Output<LocalScheme>(local_scheme, maximum_size_of_the_solution_pool) { }


    /** Number of restarts. */
    Counter number_of_restarts = 0;


    virtual nlohmann::json to_json() const override
    {
        nlohmann::json json = Output<LocalScheme>::to_json();
        json.merge_patch({
            {"NumberOfRestarts", number_of_restarts}});
        return json;
    }

    virtual int format_width() const override { return 30; }

    virtual void format(std::ostream& os) const override
    {
        Output<LocalScheme>::format(os);
        int width = format_width();
        os
            << std::setw(width) << std::left << "Number of restarts: " << number_of_restarts << std::endl
            ;
    }

};

template <typename LocalScheme>
inline const MultiStartLocalSearchOutput<LocalScheme> multi_start_local_search(
        LocalScheme& local_scheme,
        const MultiStartLocalSearchParameters<LocalScheme>& parameters = {});

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// Template implementations //////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename LocalScheme>
inline const MultiStartLocalSearchOutput<LocalScheme> multi_start_local_search(
        LocalScheme& local_scheme,
        const MultiStartLocalSearchParameters<LocalScheme>& parameters)
{
    MultiStartLocalSearchOutput<LocalScheme> output(
            local_scheme,
            parameters.maximum_size_of_the_solution_pool);
    AlgorithmFormatter<LocalScheme> algorithm_formatter(local_scheme, parameters, output);
    algorithm_formatter.start("Multi-start local search");
    algorithm_formatter.print_header();

    std::mt19937_64 generator(parameters.seed);

    Counter number_of_initial_solutions
        = (Counter)parameters.initial_solution_ids.size()
        + (Counter)parameters.initial_solutions.size();
    for (output.number_of_restarts = 1;; output.number_of_restarts++) {
        // Check maximum number of restarts.
        if (parameters.maximum_number_of_restarts >= 0
                && output.number_of_restarts > parameters.maximum_number_of_restarts)
            break;

        // Check time.
        if (parameters.timer.needs_to_end())
            break;

        // Check goal.
        if (parameters.has_goal
                && output.solution_pool.size() > 0
                && !strictly_better(
                    local_scheme,
                    parameters.goal,
                    local_scheme.global_cost(output.solution_pool.best())))
            break;

        // Generate initial solution.
        Counter initial_solution_pos = (output.number_of_restarts - 1) % number_of_initial_solutions;
        auto solution = (initial_solution_pos < (Counter)parameters.initial_solution_ids.size())?
            local_scheme.initial_solution(parameters.initial_solution_ids[initial_solution_pos], generator):
            parameters.initial_solutions[initial_solution_pos - (Counter)parameters.initial_solution_ids.size()];
        // Local Search.
        local_scheme.local_search(solution, generator);

        // Check for a new best solution.
        if (output.solution_pool.size() == 0
                || strictly_better(
                    local_scheme,
                    local_scheme.global_cost(solution),
                    local_scheme.global_cost(output.solution_pool.worst()))) {
            std::stringstream ss;
            ss << "iteration " << output.number_of_restarts;
            algorithm_formatter.update_solution(solution, ss);
        }
    }

    algorithm_formatter.end();
    return output;
}

}
