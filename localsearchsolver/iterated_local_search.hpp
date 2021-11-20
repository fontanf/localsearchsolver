#pragma once

#include "localsearchsolver/common.hpp"

#include <unordered_set>
#include <thread>

namespace localsearchsolver
{

template <typename LocalScheme>
using IteratedLocalSearchCallback = std::function<void(const typename LocalScheme::Solution&)>;

template <typename LocalScheme>
struct IteratedLocalSearchOptionalParameters
{
    typedef typename LocalScheme::Solution Solution;

    /** Maximum number of iterations. */
    Counter maximum_number_of_iterations = -1;
    /** Ids of generated initial solutions. */
    std::vector<Counter> initial_solution_ids = {0};
    /** User-provided initial solutions. */
    std::vector<Solution> initial_solutions;
    /** Maximum size of the solution pool. */
    Counter maximum_size_of_the_solution_pool = 1;
    /** Seed. */
    Seed seed = 0;
    /** Callback function called when a new best solution is found. */
    IteratedLocalSearchCallback<LocalScheme> new_solution_callback
        = [](const Solution& solution) { (void)solution; };
    /** Info structure. */
    optimizationtools::Info info;
};

template <typename LocalScheme>
struct IteratedLocalSearchOutput
{
    /** Constructor. */
    IteratedLocalSearchOutput(
            const LocalScheme& local_scheme,
            Counter maximum_size_of_the_solution_pool):
        solution_pool(local_scheme, maximum_size_of_the_solution_pool) { }

    /** Solution pool. */
    SolutionPool<LocalScheme> solution_pool;
    /** Number of iterations. */
    Counter number_of_iterations = 0;
};

template <typename LocalScheme>
inline IteratedLocalSearchOutput<LocalScheme> iterated_local_search(
        LocalScheme& local_scheme,
        IteratedLocalSearchOptionalParameters<LocalScheme> parameters = {});

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// Template implementations //////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename LocalScheme>
inline IteratedLocalSearchOutput<LocalScheme> iterated_local_search(
        LocalScheme& local_scheme,
        IteratedLocalSearchOptionalParameters<LocalScheme> parameters)
{
    typedef typename LocalScheme::Solution Solution;
    typedef typename LocalScheme::Move Move;

    // Initial display.
    VER(parameters.info,
               "=======================================" << std::endl
            << "          Local Search Solver          " << std::endl
            << "=======================================" << std::endl
            << std::endl
            << "Algorithm" << std::endl
            << "---------" << std::endl
            << "Iterated Local Search" << std::endl
            << std::endl
            << "Parameters" << std::endl
            << "----------" << std::endl
            << "Maximum number of iterations:    " << parameters.maximum_number_of_iterations << std::endl
            << "Seed:                            " << parameters.seed << std::endl
            << "Maximum size of the pool:        " << parameters.maximum_size_of_the_solution_pool << std::endl
            << "Time limit:                      " << parameters.info.time_limit << std::endl
            << std::endl
       );

    auto move_compare = [](const Move& move_1, const Move& move_2) -> bool
    {
        return move_1.global_cost < move_2.global_cost;
    };

    //std::cout << "iterated_local_search start" << std::endl;
    IteratedLocalSearchOutput<LocalScheme> output(
            local_scheme,
            parameters.maximum_size_of_the_solution_pool);
    output.solution_pool.display_init(parameters.info);

    std::mt19937_64 generator(parameters.seed);

    // Generate initial solutions.
    Counter number_of_initial_solutions
        = (Counter)parameters.initial_solution_ids.size()
        + (Counter)parameters.initial_solutions.size();
    Solution solution_cur;
    for (Counter initial_solution_pos = 0;
            initial_solution_pos < number_of_initial_solutions;
            initial_solution_pos++) {

        auto solution_tmp = (initial_solution_pos < (Counter)parameters.initial_solution_ids.size())?
            local_scheme.initial_solution(parameters.initial_solution_ids[initial_solution_pos], generator):
            parameters.initial_solutions[initial_solution_pos - (Counter)parameters.initial_solution_ids.size()];
        local_scheme.local_search(solution_tmp, generator);

        // Check for a new best solution.
        if (output.solution_pool.size() == 0
                || local_scheme.global_cost(output.solution_pool.worst())
                > local_scheme.global_cost(solution_tmp)) {
            std::stringstream ss;
            ss << "initial solution " << initial_solution_pos;
            output.solution_pool.add(solution_tmp, ss, parameters.info);
            output.solution_pool.display(ss, parameters.info);
            parameters.new_solution_callback(solution_tmp);
        }

        if (initial_solution_pos == 0
                || local_scheme.global_cost(solution_cur)
                > local_scheme.global_cost(solution_tmp)) {
            solution_cur = solution_tmp;
        }
    }

    Counter perturbation_id = 0;
    std::vector<Move> perturbations = local_scheme.perturbations(solution_cur, generator);
    auto global_cost_cur = local_scheme.global_cost(solution_cur);
    for (Move& move: perturbations)
        move.global_cost = update_move_cost(move.global_cost, global_cost_cur);
    // Sort moves.
    std::sort(perturbations.begin(), perturbations.end(), move_compare);
    for (output.number_of_iterations = 0;; ++output.number_of_iterations) {

        // Check end.
        if (parameters.info.needs_to_end())
            break;

        if (perturbation_id >= (Counter)perturbations.size())
            break;

        // Apply perturbation and local search.
        Solution solution_tmp = solution_cur;
        auto move = perturbations[perturbation_id];
        local_scheme.apply_move(solution_tmp, move);
        local_scheme.local_search(solution_tmp, generator, move);

        // Check for a new best solution.
        if (output.solution_pool.size() == 0
                || local_scheme.global_cost(output.solution_pool.worst())
                > local_scheme.global_cost(solution_tmp)) {
            std::stringstream ss;
            ss << "iteration " << output.number_of_iterations
                << " child " << perturbation_id;
            output.solution_pool.add(solution_tmp, ss, parameters.info);
            output.solution_pool.display(ss, parameters.info);
            parameters.new_solution_callback(solution_tmp);
        }

        if (local_scheme.global_cost(solution_cur)
                > local_scheme.global_cost(solution_tmp)) {
            solution_cur = solution_tmp;
            perturbations = local_scheme.perturbations(solution_cur, generator);
            perturbation_id = 0;
            auto global_cost_cur = local_scheme.global_cost(solution_cur);
            for (Move& move: perturbations)
                move.global_cost = update_move_cost(move.global_cost, global_cost_cur);
            // Sort moves.
            std::sort(perturbations.begin(), perturbations.end(), move_compare);
        } else {
            perturbation_id++;
        }
    }

    output.solution_pool.display_end(parameters.info);
    return output;
}

}

