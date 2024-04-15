#pragma once

#include "localsearchsolver/algorithm_formatter.hpp"

namespace localsearchsolver
{

template <typename LocalScheme>
struct IteratedLocalSearchParameters: Parameters<LocalScheme>
{
    using Solution = typename LocalScheme::Solution;
    using GlobalCost = typename LocalScheme::GlobalCost;

    /** Maximum number of iterations. */
    Counter maximum_number_of_iterations = -1;

    /** Maximum number of restarts. */
    Counter maximum_number_of_restarts = -1;

    /**
     * Minimum number of perturbation.
     *
     * In the literature, such a strategy is known as "Evolutionary Local
     * Search"
     *
     * See "Evolutionary Local Search for the Super-Peer Selection Problem and
     * the p-Hub Median Problem" (WolfPeter et Merz, 2007)
     * https://doi.org/10.1007/978-3-540-75514-2_1
     */
    Counter minimum_number_of_perturbations = 1;


    virtual nlohmann::json to_json(
            const LocalScheme& local_scheme) const override
    {
        nlohmann::json json = Parameters<LocalScheme>::to_json(local_scheme);
        json.merge_patch({
            {"MaximumNumberOfIterations", maximum_number_of_iterations},
            {"MaximumNumberOfRestarts", maximum_number_of_restarts},
            {"MinimumNumberOfPerturbations", minimum_number_of_perturbations}});
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
            << std::setw(width) << std::left << "Maximum number of iterations: " << maximum_number_of_iterations << std::endl
            << std::setw(width) << std::left << "Maximum number of restarts: " << maximum_number_of_restarts << std::endl
            << std::setw(width) << std::left << "Minimum number of perturbations: " << minimum_number_of_perturbations << std::endl
            ;
    }
};

template <typename LocalScheme>
struct IteratedLocalSearchOutput: Output<LocalScheme>
{
    /** Constructor. */
    IteratedLocalSearchOutput(
            const LocalScheme& local_scheme,
            Counter maximum_size_of_the_solution_pool):
        Output<LocalScheme>(local_scheme, maximum_size_of_the_solution_pool) { }


    /** Number of iterations. */
    Counter number_of_iterations = 0;

    /** Number of restarts. */
    Counter number_of_restarts = 0;


    virtual nlohmann::json to_json() const override
    {
        nlohmann::json json = Output<LocalScheme>::to_json();
        json.merge_patch({
            {"NumberOfIterations", number_of_iterations},
            {"NumberOfRestarts", number_of_restarts}});
        return json;
    }

    virtual int format_width() const override { return 30; }

    virtual void format(std::ostream& os) const override
    {
        Output<LocalScheme>::format(os);
        int width = format_width();
        os
            << std::setw(width) << std::left << "Number of iterations: " << number_of_iterations << std::endl
            << std::setw(width) << std::left << "Number of restarts: " << number_of_restarts << std::endl
            ;
    }
};

template <typename LocalScheme>
inline const IteratedLocalSearchOutput<LocalScheme> iterated_local_search(
        LocalScheme& local_scheme,
        const IteratedLocalSearchParameters<LocalScheme>& parameters = {});

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// Template implementations //////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename LocalScheme>
inline const IteratedLocalSearchOutput<LocalScheme> iterated_local_search(
        LocalScheme& local_scheme,
        const IteratedLocalSearchParameters<LocalScheme>& parameters)
{
    using Solution = typename LocalScheme::Solution;
    using Perturbation = typename LocalScheme::Perturbation;

    IteratedLocalSearchOutput<LocalScheme> output(
            local_scheme,
            parameters.maximum_size_of_the_solution_pool);
    AlgorithmFormatter<LocalScheme> algorithm_formatter(local_scheme, parameters, output);
    algorithm_formatter.start("Iterated local search");
    algorithm_formatter.print_header();

    auto move_compare = [](const Perturbation& move_1, const Perturbation& move_2) -> bool
    {
        return move_1.global_cost < move_2.global_cost;
    };

    std::mt19937_64 generator(parameters.seed);

    std::vector<Solution> initial_solutions;
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

        // Generate initial solutions.
        if (initial_solutions.empty()) {
            for (Counter initial_solution_pos = 0;
                    initial_solution_pos < number_of_initial_solutions;
                    initial_solution_pos++) {

                auto solution_tmp = (initial_solution_pos < (Counter)parameters.initial_solution_ids.size())?
                    local_scheme.initial_solution(parameters.initial_solution_ids[initial_solution_pos], generator):
                    parameters.initial_solutions[initial_solution_pos - (Counter)parameters.initial_solution_ids.size()];
                local_scheme.local_search(solution_tmp, generator);

                // Check for a new best solution.
                if (output.solution_pool.size() == 0
                        || strictly_better(
                            local_scheme,
                            local_scheme.global_cost(solution_tmp),
                            local_scheme.global_cost(output.solution_pool.worst()))) {
                    std::stringstream ss;
                    ss << "s" << output.number_of_restarts;
                    algorithm_formatter.update_solution(solution_tmp, ss);
                }

                initial_solutions.push_back(solution_tmp);
            }

            std::sort(
                    initial_solutions.begin(),
                    initial_solutions.end(),
                    [&local_scheme](const Solution& solution_1, const Solution& solution_2) -> bool
                    {
                        return local_scheme.global_cost(solution_1)
                            > local_scheme.global_cost(solution_2);
                    });
        }

        Solution solution_cur = initial_solutions.back();
        initial_solutions.pop_back();
        Counter perturbation_id = 0;
        std::vector<Perturbation> perturbations = local_scheme.perturbations(solution_cur, generator);
        // Sort moves.
        std::sort(perturbations.begin(), perturbations.end(), move_compare);
        Counter depth = 1;
        Solution solution_next = solution_cur;
        bool better_found = false;
        for (;; ++output.number_of_iterations) {

            // Check end.
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

            if (perturbation_id >= parameters.minimum_number_of_perturbations
                    && better_found) {
                solution_cur = solution_next;
                better_found = false;
                perturbation_id = 0;
                depth++;
                perturbations = local_scheme.perturbations(solution_cur, generator);
                // Sort moves.
                std::sort(perturbations.begin(), perturbations.end(), move_compare);
            }

            if (perturbation_id >= (Counter)perturbations.size())
                break;

            // Apply perturbation and local search.
            Solution solution_tmp = solution_cur;
            auto perturbation = perturbations[perturbation_id];
            local_scheme.apply_perturbation(solution_tmp, perturbation, generator);
            local_scheme.local_search(solution_tmp, generator, perturbation);

            // Check for a new best solution.
            if (output.solution_pool.size() == 0
                || strictly_better(
                    local_scheme,
                    local_scheme.global_cost(solution_tmp),
                    local_scheme.global_cost(output.solution_pool.worst()))) {
                std::stringstream ss;
                ss << "s" << output.number_of_restarts
                    << " d" << depth
                    << " c" << perturbation_id << "/" << perturbations.size()
                    << " i" << output.number_of_iterations;
                algorithm_formatter.update_solution(solution_tmp, ss);
            }

            if (local_scheme.global_cost(solution_next)
                    > local_scheme.global_cost(solution_tmp)) {
                solution_next = solution_tmp;
                better_found = true;
            }

            perturbation_id++;
        }
    }

    algorithm_formatter.end();
    return output;
}

}
