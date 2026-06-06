#pragma once

#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>

#include "optimizationtools//utils//utils.hpp"

namespace localsearchsolver
{

using Counter = int64_t;
using Distance = int64_t;

template <typename Solution, typename Cost>
using PenalizedCostCallback = std::function<Cost(const Solution&)>;

template <typename Solution>
using DistanceCallback = std::function<Distance(const Solution&, const Solution&)>;

/*
 * Class to manage a population of solutions.
 *
 * This population is designed to handle varying penalties for constraints.
 *
 * It is based on:
 *
 * - "A Hybrid Genetic Algorithm for Multidepot and Periodic Vehicle Routing
 *   Problems" (Vidal et al., 2012)
 *   https://doi.org/10.1287/opre.1120.1048
 */
template <typename Solution, typename Cost>
class Population
{

public:

    struct Parameters
    {
        /** Minimum size of the population. */
        Counter minimum_size = 25;

        /**
         * Maximum size of the population.
         *
         * When adding a new solution to the population makes its size to go
         * above this value, solutions are removed from the population until it
         * reaches the size `minimum_size`.
         */
        Counter maximum_size = 25 + 40;

        /**
         * Parameter used when computing the diversity contribution of a
         * solution.
         *
         * The diversity contribution of a solution is its average distance to
         * its `number_of_closest_neighbors` closest neighbors.
         */
        Counter number_of_closest_neighbors = 3;

        /**
         * Number of elite solutions.
         */
        Counter number_of_elite_solutions = 8;
    };

    struct PopulationSolution
    {
        PopulationSolution(const Solution& solution):
            solution(solution) { }

        /** Solution. */
        Solution solution;

        Cost penalized_cost = Cost{};

        /** Penalized cost rank. */
        Counter penalized_cost_rank = -1;

        /** Diversity contribution. */
        double diversity = std::numeric_limits<double>::infinity();

        /** Diversity rank. */
        Counter diversity_rank = -1;

        /** Biased fitness. */
        double biased_fitness = 0;

        bool to_remove = false;
    };

    /** Constructor. */
    Population(
            PenalizedCostCallback<Solution, Cost> penalized_cost_callback,
            DistanceCallback<Solution> distance_callback,
            const Parameters& parameters):
        penalized_cost_callback_(penalized_cost_callback),
        distance_callback_(distance_callback),
        parameters_(parameters) {}

    /** Get population parameters. */
    const Parameters& parameters() const { return parameters_; }

    /** Get the size of the population. */
    Counter size() const { return solutions_.size(); }

    /** Get a solution of the population. */
    const PopulationSolution& solution(Counter solution_id) const { return solutions_[solution_id]; }

    /** Add a solution to the population. */
    void add(
            const Solution& solution,
            std::mt19937_64& generator);

    /** Get two parent solutions from a binary tournament. */
    std::pair<Solution, Solution> binary_tournament(
        std::mt19937_64& generator);

    /** Get one parent solution from a binary tournament. */
    Solution binary_tournament_single(
        std::mt19937_64& generator);

    /** Get the best solution of the population (lowest penalized cost). */
    const Solution& best_solution() const;

private:

    /*
     * Private methods
     */

    /** Compute penalized costs, ranks, distances, diversity, and biased fitness. */
    void compute_fitness(std::mt19937_64& generator);

    void survivor_selection(
        std::mt19937_64& generator);

    /*
     * Private attributes
     */

    /** Penalized cost callback. */
    PenalizedCostCallback<Solution, Cost> penalized_cost_callback_;

    /** Distance callback. */
    DistanceCallback<Solution> distance_callback_;

    /** Parameters. */
    Parameters parameters_;

    /** Solutions. */
    std::vector<PopulationSolution> solutions_;

    /** Pairwise distances between solutions (valid when fitness_valid_ is true). */
    std::vector<std::vector<Distance>> distances_;

    /** Whether biased fitness values are up to date. */
    bool fitness_valid_ = false;

};

}

template <typename Solution, typename Cost>
void localsearchsolver::Population<Solution, Cost>::add(
        const Solution& solution,
        std::mt19937_64& generator)
{
    PopulationSolution population_solution(solution);
    this->solutions_.push_back(population_solution);
    fitness_valid_ = false;

    if (this->size() > parameters_.maximum_size)
        this->survivor_selection(generator);
}

template <typename Solution, typename Cost>
void localsearchsolver::Population<Solution, Cost>::compute_fitness(
        std::mt19937_64& generator)
{
    // Compute the penalized cost of each solution.
    for (Counter solution_id = 0; solution_id < this->size(); ++solution_id) {
        PopulationSolution& solution = this->solutions_[solution_id];
        solution.penalized_cost = penalized_cost_callback_(solution.solution);
    }

    // Compute the penalized cost rank of each solution.
    std::vector<Counter> penalized_cost_order(this->size());
    std::iota(penalized_cost_order.begin(), penalized_cost_order.end(), 0);
    std::shuffle(penalized_cost_order.begin(), penalized_cost_order.end(), generator);
    std::sort(
            penalized_cost_order.begin(), penalized_cost_order.end(),
            [this](Counter solution_1_id, Counter solution_2_id) -> bool
            {
                return this->solutions_[solution_1_id].penalized_cost
                    < this->solutions_[solution_2_id].penalized_cost;
            });
    for (Counter pos = 0; pos < this->size(); ++pos)
        this->solutions_[penalized_cost_order[pos]].penalized_cost_rank = pos;

    // Compute the distances between each pair of solutions.
    distances_.assign(this->size(), std::vector<Distance>(this->size(), 0));
    for (Counter solution_1_id = 0; solution_1_id < this->size(); ++solution_1_id) {
        for (Counter solution_2_id = 0; solution_2_id < solution_1_id; ++solution_2_id) {
            Distance distance = distance_callback_(
                    this->solutions_[solution_1_id].solution,
                    this->solutions_[solution_2_id].solution);
            distances_[solution_1_id][solution_2_id] = distance;
            distances_[solution_2_id][solution_1_id] = distance;
        }
    }

    // Compute diversity contributions.
    for (Counter solution_id = 0; solution_id < this->size(); ++solution_id) {
        std::vector<Distance> neighbor_distances;
        neighbor_distances.reserve(this->size() - 1);
        for (Counter solution_2_id = 0; solution_2_id < this->size(); ++solution_2_id) {
            if (solution_2_id == solution_id)
                continue;
            neighbor_distances.push_back(distances_[solution_id][solution_2_id]);
        }
        Counter k = std::min(
                parameters_.number_of_closest_neighbors,
                (Counter)neighbor_distances.size());
        this->solutions_[solution_id].diversity = 0;
        if (k > 0) {
            std::nth_element(
                    neighbor_distances.begin(),
                    neighbor_distances.begin() + k - 1,
                    neighbor_distances.end());
            for (Counter pos = 0; pos < k; ++pos)
                this->solutions_[solution_id].diversity += neighbor_distances[pos];
            this->solutions_[solution_id].diversity /= k;
        }
    }

    // Compute the diversity rank of each solution.
    std::vector<Counter> diversity_order(this->size());
    std::iota(diversity_order.begin(), diversity_order.end(), 0);
    std::shuffle(diversity_order.begin(), diversity_order.end(), generator);
    std::sort(
            diversity_order.begin(), diversity_order.end(),
            [this](Counter solution_id_1, Counter solution_id_2) -> bool
            {
                return this->solutions_[solution_id_1].diversity
                    > this->solutions_[solution_id_2].diversity;
            });
    for (Counter pos = 0; pos < this->size(); ++pos)
        this->solutions_[diversity_order[pos]].diversity_rank = pos;

    // Compute the biased fitness of each solution.
    for (Counter solution_id = 0; solution_id < this->size(); ++solution_id) {
        PopulationSolution& solution = this->solutions_[solution_id];
        solution.biased_fitness = solution.penalized_cost_rank
            + (1.0 - (double)parameters_.number_of_elite_solutions / this->size())
            * solution.diversity_rank;
    }

    fitness_valid_ = true;
}

template <typename Solution, typename Cost>
void localsearchsolver::Population<Solution, Cost>::survivor_selection(
        std::mt19937_64& generator)
{
    compute_fitness(generator);

    Counter number_of_solutions_removed = 0;
    while ((Counter)this->solutions_.size() - number_of_solutions_removed > parameters_.minimum_size) {
        // Recompute diversity contributions for the remaining (non-removed) solutions.
        for (Counter solution_id = 0; solution_id < this->size(); ++solution_id) {
            if (this->solutions_[solution_id].to_remove)
                continue;
            std::vector<Distance> neighbor_distances;
            for (Counter solution_2_id = 0; solution_2_id < this->size(); ++solution_2_id) {
                if (solution_2_id == solution_id)
                    continue;
                if (this->solutions_[solution_2_id].to_remove)
                    continue;
                neighbor_distances.push_back(distances_[solution_id][solution_2_id]);
            }
            Counter k = std::min(
                    parameters_.number_of_closest_neighbors,
                    (Counter)neighbor_distances.size());
            this->solutions_[solution_id].diversity = 0;
            if (k > 0) {
                std::nth_element(
                        neighbor_distances.begin(),
                        neighbor_distances.begin() + k - 1,
                        neighbor_distances.end());
                for (Counter pos = 0; pos < k; ++pos)
                    this->solutions_[solution_id].diversity += neighbor_distances[pos];
                this->solutions_[solution_id].diversity /= k;
            }
        }

        // Recompute the diversity rank of each solution.
        std::vector<Counter> diversity_order(this->size());
        std::iota(diversity_order.begin(), diversity_order.end(), 0);
        std::shuffle(diversity_order.begin(), diversity_order.end(), generator);
        std::sort(
                diversity_order.begin(), diversity_order.end(),
                [this](Counter solution_id_1, Counter solution_id_2) -> bool
                {
                    return this->solutions_[solution_id_1].diversity
                        > this->solutions_[solution_id_2].diversity;
                });
        for (Counter pos = 0; pos < this->size(); ++pos)
            this->solutions_[diversity_order[pos]].diversity_rank = pos;

        // Recompute the biased fitness of each solution.
        for (Counter solution_id = 0; solution_id < this->size(); ++solution_id) {
            PopulationSolution& solution = this->solutions_[solution_id];
            solution.biased_fitness = solution.penalized_cost_rank
                + (1.0 - (double)parameters_.number_of_elite_solutions / this->size())
                * solution.diversity_rank;
        }

        // Remove the solution with the worst biased fitness.
        // Give priority to solutions which are at distance 0 from another solution (clone).
        Counter solution_worst_id = -1;
        bool is_worst_clone = false;
        double biased_fitness_worst = 0;
        for (Counter solution_id = 0; solution_id < this->size(); ++solution_id) {
            const PopulationSolution& solution = this->solutions_[solution_id];
            if (solution.to_remove)
                continue;
            bool is_clone = false;
            for (Counter solution_2_id = 0; solution_2_id < this->size(); ++solution_2_id) {
                if (solution_2_id == solution_id)
                    continue;
                const PopulationSolution& solution_2 = this->solutions_[solution_2_id];
                if (solution_2.to_remove)
                    continue;
                if (distances_[solution_id][solution_2_id] == 0) {
                    is_clone = true;
                    break;
                }
            }
            if (solution_worst_id == -1
                    || (!is_worst_clone && is_clone)
                    || (is_worst_clone == is_clone
                        && biased_fitness_worst < solution.biased_fitness)) {
                solution_worst_id = solution_id;
                is_worst_clone = is_clone;
                biased_fitness_worst = solution.biased_fitness;
            }
        }
        this->solutions_[solution_worst_id].to_remove = true;
        number_of_solutions_removed++;
    }

    // Remove solutions marked for removal.
    for (Counter solution_id = 0; solution_id < this->size();) {
        if (this->solutions_[solution_id].to_remove) {
            this->solutions_[solution_id] = this->solutions_.back();
            this->solutions_.pop_back();
        } else {
            solution_id++;
        }
    }

    fitness_valid_ = false;
}

template <typename Solution, typename Cost>
std::pair<Solution, Solution> localsearchsolver::Population<Solution, Cost>::binary_tournament(
        std::mt19937_64& generator)
{
    if (this->size() < 4)
        throw std::logic_error("binary_tournament requires at least 4 solutions");

    if (!fitness_valid_)
        compute_fitness(generator);

    // Draw 4 random solutions.
    auto solution_ids = optimizationtools::bob_floyd(
            (Counter)4, this->size(), generator);
    std::shuffle(solution_ids.begin(), solution_ids.end(), generator);

    Counter solution_id_1 = (this->solutions_[solution_ids[0]].biased_fitness
            < this->solutions_[solution_ids[1]].biased_fitness)?
        solution_ids[0]: solution_ids[1];

    Counter solution_id_2 = (this->solutions_[solution_ids[2]].biased_fitness
            < this->solutions_[solution_ids[3]].biased_fitness)?
        solution_ids[2]: solution_ids[3];

    return {
        this->solutions_[solution_id_1].solution,
        this->solutions_[solution_id_2].solution};
}

template <typename Solution, typename Cost>
Solution localsearchsolver::Population<Solution, Cost>::binary_tournament_single(
        std::mt19937_64& generator)
{
    if (solutions_.empty())
        throw std::logic_error("binary_tournament_single requires at least 1 solution");

    if (this->size() == 1)
        return this->solutions_[0].solution;

    if (!fitness_valid_)
        compute_fitness(generator);

    // Draw 2 random solutions.
    auto solution_ids = optimizationtools::bob_floyd(
            (Counter)2, this->size(), generator);
    std::shuffle(solution_ids.begin(), solution_ids.end(), generator);

    Counter solution_id_1 = (this->solutions_[solution_ids[0]].biased_fitness
            < this->solutions_[solution_ids[1]].biased_fitness)?
        solution_ids[0]: solution_ids[1];

    return this->solutions_[solution_id_1].solution;
}

template <typename Solution, typename Cost>
const Solution& localsearchsolver::Population<Solution, Cost>::best_solution() const
{
    if (solutions_.empty())
        throw std::logic_error("best_solution requires at least 1 solution");

    Counter best_id = 0;
    Cost best_cost = penalized_cost_callback_(this->solutions_[0].solution);
    for (Counter solution_id = 1; solution_id < this->size(); ++solution_id) {
        Cost cost = penalized_cost_callback_(this->solutions_[solution_id].solution);
        if (cost < best_cost) {
            best_cost = cost;
            best_id = solution_id;
        }
    }
    return this->solutions_[best_id].solution;
}
