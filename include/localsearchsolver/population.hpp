#pragma once

#include <cstdint>
#include <limits>
#include <vector>
#include <functional>
//#include <iostream>

#include "optimizationtools//utils//utils.hpp"

namespace localsearchsolver
{

using Counter = int64_t;
using Distance = int64_t;

template <typename Solution, typename Cost>
using PenalizedCostCallback = std::function<Cost(const Solution&)>;

template <typename Solution>
using DistanceCallback = std::function<int(const Solution&, const Solution& solution)>;

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
         * its `n_close` closest neighbors.
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

        Cost penalized_cost;

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

    /** Get of solution of the population. */
    const PopulationSolution& solution(Counter solution_id) { return solutions_[solution_id]; }

    /** Add a solution to the population. */
    void add(
            const Solution& solution,
            std::mt19937_64& generator);

    /** Get two parent solutions from a binary tournament. */
    std::pair<const Solution&, const Solution&> binary_tournament(
        std::mt19937_64& generator);

    /** Get one parent solutions from a binary tournament. */
    const Solution& binary_tournament_single(
        std::mt19937_64& generator);

private:

    /*
     * Private methods
     */

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

};

}

template <typename Solution, typename Cost>
void localsearchsolver::Population<Solution, Cost>::add(
        const Solution& solution,
        std::mt19937_64& generator)
{
    // Add the new solution to the population.
    PopulationSolution population_solution(solution);
    this->solutions_.push_back(population_solution);

    // If the size of the population goes above the maximum size allowed,
    // run the survivor selection.
    if (this->size() > parameters_.maximum_size)
        this->survivor_selection(generator);
}

template <typename Solution, typename Cost>
void localsearchsolver::Population<Solution, Cost>::survivor_selection(
        std::mt19937_64& generator)
{
    // Compute the penalized cost of each solution.
    //std::cout << "compute penalized costs..." << std::endl;
    for (Counter solution_id = 0; solution_id < this->size(); ++solution_id) {
        PopulationSolution& solution = this->solutions_[solution_id];
        solution.penalized_cost = penalized_cost_callback_(solution.solution);
    }

    // Compute the penalized cost rank of each solution.
    //std::cout << "compute penalized cost ranks..." << std::endl;
    std::vector<Counter> penalized_cost_ranks(this->size());
    std::iota(penalized_cost_ranks.begin(), penalized_cost_ranks.end(), 0);
    std::shuffle(penalized_cost_ranks.begin(), penalized_cost_ranks.end(), generator);
    sort(
            penalized_cost_ranks.begin(), penalized_cost_ranks.end(),
            [this](Counter solution_1_id, Counter solution_2_id) -> bool
            {
                const PopulationSolution& solution_1 = this->solutions_[solution_1_id];
                const PopulationSolution& solution_2 = this->solutions_[solution_2_id];
                return solution_1.penalized_cost < solution_2.penalized_cost;
            });
    for (Counter pos = 0; pos < this->size(); ++pos)
        this->solutions_[penalized_cost_ranks[pos]].penalized_cost_rank = pos;

    // Compute the distances between each pair of solutions.
    //std::cout << "compute distances..." << std::endl;
    std::vector<std::vector<Distance>> distances(
            this->size(), std::vector<Distance>(this->size(), 0));
    for (Counter solution_1_id = 0; solution_1_id < this->size(); ++solution_1_id) {
        const PopulationSolution& solution_1 = this->solutions_[solution_1_id];
        for (Counter solution_2_id = 0; solution_2_id < solution_1_id; ++solution_2_id) {
            const PopulationSolution& solution_2 = this->solutions_[solution_2_id];
            Distance distance = distance_callback_(solution_1.solution, solution_2.solution);
            distances[solution_1_id][solution_2_id] = distance;
            distances[solution_2_id][solution_1_id] = distance;
        }
    }

    Counter number_of_solutions_removed = 0;
    while (this->solutions_.size() - number_of_solutions_removed > parameters_.minimum_size) {
        // For each solution compute its closest neighbors and its diversity
        // contribution.
        //std::cout << "compute diversity contributions..." << std::endl;
        for (Counter solution_id = 0; solution_id < this->size(); ++solution_id) {
            std::vector<Distance> solution_distances;
            PopulationSolution& solution = this->solutions_[solution_id];
            for (Counter solution_2_id = 0; solution_2_id < this->size(); ++solution_2_id) {
                if (solution_2_id == solution_id)
                    continue;
                const PopulationSolution& solution_2 = this->solutions_[solution_2_id];
                if (solution_2.to_remove)
                    continue;
                Distance distance = distances[solution_id][solution_2_id];
                solution_distances.push_back(distance);
            }
            std::nth_element(
                    solution_distances.begin(),
                    solution_distances.begin() + parameters_.number_of_closest_neighbors - 1,
                    solution_distances.end());
            solution.diversity = 0;
            for (Counter pos = 0; pos < parameters_.number_of_closest_neighbors; ++pos)
                solution.diversity += solution_distances[pos];
            solution.diversity /= parameters_.number_of_closest_neighbors;

        }

        // Compute the diversity rank of each solution.
        //std::cout << "compute diversity contribution ranks..." << std::endl;
        std::vector<Counter> diversity_ranks(this->size());
        std::iota(diversity_ranks.begin(), diversity_ranks.end(), 0);
        std::shuffle(diversity_ranks.begin(), diversity_ranks.end(), generator);
        sort(
                diversity_ranks.begin(), diversity_ranks.end(),
                [this](Counter solution_id_1, Counter solution_id_2) -> bool
                {
                    return this->solutions_[solution_id_1].diversity
                        > this->solutions_[solution_id_2].diversity;
                });
        for (Counter pos = 0; pos < this->size(); ++pos)
            this->solutions_[diversity_ranks[pos]].diversity_rank = pos;

        // Compute the biased fitness of each solution.
        //std::cout << "compute biased fitnesses..." << std::endl;
        for (Counter solution_id = 0; solution_id < this->size(); ++solution_id) {
            PopulationSolution& solution = this->solutions_[solution_id];
            solution.biased_fitness = solution.penalized_cost_rank
                + (1 - (double)parameters_.number_of_elite_solutions / this->size())
                * solution.diversity_rank;
        }

        // Remove the solution with the worst biased fitness.
        // Give priority to solutions which are at distance 0 from another
        // solution of the population (clone).
        //std::cout << "find solution to remove..." << std::endl;
        Counter solution_worst_id = -1;
        bool is_worst_clone = false;
        double biased_fitness_worst = 0;
        for (Counter solution_id = 0; solution_id < this->size(); ++ solution_id) {
            const PopulationSolution& solution = this->solutions_[solution_id];
            if (solution.to_remove)
                continue;
            bool is_clone = false;
            for (Counter solution_2_id = 0; solution_2_id < this->size(); ++solution_2_id) {
                const PopulationSolution& solution_2 = this->solutions_[solution_2_id];
                if (solution_2.to_remove)
                    continue;
                if (distances[solution_id][solution_2_id] == 0) {
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
        //std::cout << "solution_worst_id " << solution_worst_id << std::endl;
        //std::cout << "penalized_cost " << solutions_[solution_worst_id].penalized_cost
        //    << " rank " << solutions_[solution_worst_id].penalized_cost_rank
        //    << " diversity " << solutions_[solution_worst_id].diversity
        //    << " rank " << solutions_[solution_worst_id].diversity_rank
        //    << " biased_fitness " << solutions_[solution_worst_id].biased_fitness
        //    << std::endl;
        this->solutions_[solution_worst_id].to_remove = true;
        number_of_solutions_removed++;
    }

    // Remove solution marked to be removed.
    //std::cout << "remove solutions..." << std::endl;
    for (Counter solution_id = 0; solution_id < this->size();) {
        PopulationSolution& solution = this->solutions_[solution_id];
        if (solution.to_remove) {
            this->solutions_[solution_id] = this->solutions_.back();
            this->solutions_.pop_back();
        } else {
            solution_id++;
        }
    }
}

template <typename Solution, typename Cost>
std::pair<const Solution&, const Solution&> localsearchsolver::Population<Solution, Cost>::binary_tournament(
        std::mt19937_64& generator)
{
    // Draw 4 random solutions.
    auto solution_ids = optimizationtools::bob_floyd(
            (Counter)4, this->size(), generator);
    std::shuffle(solution_ids.begin(), solution_ids.end(), generator);

    // Compute solution_id_1.
    Counter solution_id_1 = -1;
    if (this->solutions_[solution_ids[0]].biased_fitness
            < this->solutions_[solution_ids[1]].biased_fitness) {
        solution_id_1 = solution_ids[0];
    } else {
        solution_id_1 = solution_ids[1];
    }

    // Compute solution_id_2.
    Counter solution_id_2 = -1;
    if (this->solutions_[solution_ids[2]].biased_fitness
            < this->solutions_[solution_ids[3]].biased_fitness) {
        solution_id_2 = solution_ids[2];
    } else {
        solution_id_2 = solution_ids[3];
    }

    return {
        this->solutions_[solution_id_1].solution,
        this->solutions_[solution_id_2].solution};
}

template <typename Solution, typename Cost>
const Solution& localsearchsolver::Population<Solution, Cost>::binary_tournament_single(
        std::mt19937_64& generator)
{
    if (this->size() == 1)
        return this->solutions_[0].solution;

    // Draw 4 random solutions.
    auto solution_ids = optimizationtools::bob_floyd(
            (Counter)2, this->size(), generator);
    std::shuffle(solution_ids.begin(), solution_ids.end(), generator);

    // Compute solution_id_1.
    Counter solution_id_1 = -1;
    if (this->solutions_[solution_ids[0]].biased_fitness
            < this->solutions_[solution_ids[1]].biased_fitness) {
        solution_id_1 = solution_ids[0];
    } else {
        solution_id_1 = solution_ids[1];
    }

    return this->solutions_[solution_id_1].solution;
}
