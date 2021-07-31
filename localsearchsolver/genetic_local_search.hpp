#pragma once

#include "localsearchsolver/a_star_local_search.hpp"

namespace localsearchsolver
{

template <typename LocalScheme>
using GeneticLocalSearchCallback = std::function<void(const typename LocalScheme::Solution&)>;

template <typename LocalScheme>
struct GeneticLocalSearchOptionalParameters
{
    typedef typename LocalScheme::Solution Solution;

    /** Number of threads. */
    Counter thread_number = 1;
    /** Maximum number of genetic iterations. */
    Counter genetic_iteration_number_max = -1;
    /** Number of nodes of internal A* Local Search calls. */
    Counter local_search_node_number_max = 1000;
    /** Maximum size of the population. */
    Counter population_size_max = 10;
    /** Ids of generated initial solutions. */
    std::vector<Counter> initial_solution_ids = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    /** User-provided initial solutions. */
    std::vector<Solution> initial_solutions;
    /** Seed. */
    Seed seed = 0;
    /** Callback function called when a new best solution is found. */
    GeneticLocalSearchCallback<LocalScheme> new_solution_callback
        = [](const Solution& solution) { (void)solution; };

    optimizationtools::Info info;
};

template <typename LocalScheme>
struct GeneticLocalSearchOutput
{
    /** Constructor. */
    GeneticLocalSearchOutput(const LocalScheme& local_scheme):
        solution_pool(local_scheme, 1) { }

    /** Solution pool. */
    SolutionPool<LocalScheme> solution_pool;
    /** Number of genetic iterations. */
    Counter genetic_iteration_number = 0;
};

template <typename LocalScheme>
inline GeneticLocalSearchOutput<LocalScheme> genetic_local_search(
        LocalScheme& local_scheme,
        GeneticLocalSearchOptionalParameters<LocalScheme> parameters = {});

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// Template implementations //////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename LocalScheme>
class Population
{
    typedef typename LocalScheme::Solution Solution;
    typedef typename LocalScheme::GlobalCost GlobalCost;

    struct PopulationSolution
    {
        Solution solution;
        double diversity = std::numeric_limits<double>::infinity();
        Counter diversity_rank = -1;
        Counter global_cost_rank = -1;
        Counter rank = -1;
    };

public:

    /** Constructor. */
    Population(const LocalScheme& local_scheme, Counter population_size_max):
        local_scheme_(local_scheme),
        population_size_max_(population_size_max)
    { }

    /** Destructor. */
    virtual ~Population() { }

    /** Return the current number of solutions in the solution pool. */
    Counter size() const { return solutions_.size(); }

    /**
     * Return a pair of distinct solutions to use as parents for a crossover
     * operator.
     */
    std::pair<Solution, Solution> get_parents(
            std::mt19937_64& generator)
    {
        // Draw 4 random solutions.
        auto solution_ids = optimizationtools::bob_floyd(
                (Counter)4, size() - 1, generator);
        std::shuffle(solution_ids.begin(), solution_ids.end(), generator);

        // Compute solution_id_1.
        Counter solution_id_1 = -1;
        if (solutions_[solution_ids[0]].rank
                < solutions_[solution_ids[1]].rank) {
            solution_id_1 = solution_ids[0];
        } else {
            solution_id_1 = solution_ids[1];
        }

        // Compute solution_id_2.
        Counter solution_id_2 = -1;
        if (solutions_[solution_ids[2]].rank
                < solutions_[solution_ids[3]].rank) {
            solution_id_2 = solution_ids[2];
        } else {
            solution_id_2 = solution_ids[3];
        }

        return {
            solutions_[solution_id_1].solution,
            solutions_[solution_id_2].solution};
    }

    /**
     * Add a solution to the population.
     *
     * If the solution is already in the population, it is not added again.
     *
     * If the size of the population exceeds the maximum size of the
     * population, then the worst solution is removed from the population.
     */
    void add(const Solution& solution, std::mt19937_64& generator)
    {
        // Check if the solution is already in the population.
        for (Counter solution_id = 0; solution_id < size(); ++solution_id) {
            const Solution& solution_2 = solutions_[solution_id].solution;
            if (local_scheme_.distance(solution, solution_2) == 0) {
                std::cout << "New solution is identical to solution " << solution_id << "." << std::endl;
                return;
            }
        }

        // Add the new solution to the population.
        PopulationSolution population_solution;
        population_solution.solution = solution;
        std::cout << "Add new solution." << std::endl;
        solutions_.push_back(population_solution);
        update_scores(generator);

        // If the maximum population size is exceeded, remove the worst
        // solution.
        if (size() > population_size_max_) {
            Counter solution_id_worse = -1;
            for (Counter solution_id = 0; solution_id < size(); ++solution_id)
                if (solutions_[solution_id].rank == size() - 1)
                    solution_id_worse = solution_id;
            std::cout << "Remove solution " << solution_id_worse << "." << std::endl;
            solutions_[solution_id_worse] = solutions_.back();
            solutions_.pop_back();
            update_scores(generator);
        }
    }

    void print()
    {
        for (Counter solution_id = 0; solution_id < size(); ++solution_id) {
            const PopulationSolution& solution = solutions_[solution_id];
            std::cout << "Solution " << solution_id << ":"
                << " cost " << to_string(local_scheme_.global_cost(solution.solution))
                << "; cost rank " << solution.global_cost_rank
                << "; diversity: " << solution.diversity
                << "; diversity rank: " << solution.diversity_rank
                << "; rank: " << solution.rank
                << std::endl;
        }
    }

private:

    /*
     * Private methods.
     */

    /**
     * Update the scores of the solutions of the population.
     *
     * This method should be called each time the population is modified (when
     * a solution is added or removed) since it affects these values.
     */
    void update_scores(std::mt19937_64& generator)
    {
        // Compute diversities.
        for (Counter solution_id = 0; solution_id < size(); ++solution_id) {
            const Solution& solution = solutions_[solution_id].solution;
            solutions_[solution_id].diversity = std::numeric_limits<double>::infinity();
            for (Counter solution_id_2 = 0; solution_id_2 < size(); ++solution_id_2) {
                if (solution_id_2 == solution_id)
                    continue;
                const Solution& solution_2 = solutions_[solution_id_2].solution;
                double d = local_scheme_.distance(solution, solution_2);
                if (solutions_[solution_id].diversity > d)
                    solutions_[solution_id].diversity = d;
            }
        }

        // Compute diversity_ranks.
        std::vector<Counter> diversity_ranks(size());
        std::iota(diversity_ranks.begin(), diversity_ranks.end(), 0);
        std::shuffle(diversity_ranks.begin(), diversity_ranks.end(), generator);
        sort(
                diversity_ranks.begin(), diversity_ranks.end(),
                [this](Counter solution_id_1, Counter solution_id_2) -> bool
                {
                    return solutions_[solution_id_1].diversity
                            > solutions_[solution_id_2].diversity;
                });
        for (Counter pos = 0; pos < size(); ++pos)
            solutions_[diversity_ranks[pos]].diversity_rank = pos;

        // Compute global_cost_ranks.
        std::vector<Counter> global_cost_ranks(size());
        std::iota(global_cost_ranks.begin(), global_cost_ranks.end(), 0);
        std::shuffle(global_cost_ranks.begin(), global_cost_ranks.end(), generator);
        sort(
                global_cost_ranks.begin(), global_cost_ranks.end(),
                [this](Counter solution_id_1, Counter solution_id_2) -> bool
                {
                    return local_scheme_.global_cost(solutions_[solution_id_1].solution)
                            < local_scheme_.global_cost(solutions_[solution_id_2].solution);
                });
        for (Counter pos = 0; pos < size(); ++pos)
            solutions_[global_cost_ranks[pos]].global_cost_rank = pos;

        // Compute ranks.
        std::vector<Counter> ranks(size());
        std::iota(ranks.begin(), ranks.end(), 0);
        std::shuffle(ranks.begin(), ranks.end(), generator);
        sort(
                ranks.begin(), ranks.end(),
                [this](Counter solution_id_1, Counter solution_id_2) -> bool
                {
                    return solutions_[solution_id_1].diversity_rank
                            + solutions_[solution_id_1].global_cost_rank
                            < solutions_[solution_id_2].diversity_rank
                            + solutions_[solution_id_2].global_cost_rank;
                });
        for (Counter pos = 0; pos < size(); ++pos)
            solutions_[ranks[pos]].rank = pos;
        print();
    }

    /*
     * Private attributes.
     */

    /** LocalScheme. */
    const LocalScheme& local_scheme_;

    /** Maximum size of the population. */
    Counter population_size_max_ = -1;

    /** Set of solutions. */
    std::vector<PopulationSolution> solutions_;

};

template <typename LocalScheme>
struct GeneticLocalSearchData
{
    GeneticLocalSearchData(
            LocalScheme& local_scheme,
            GeneticLocalSearchOptionalParameters<LocalScheme>& parameters,
            GeneticLocalSearchOutput<LocalScheme>& output):
        local_scheme(local_scheme),
        parameters(parameters),
        population(local_scheme, parameters.population_size_max),
        output(output)
        {  }

    /** Local scheme. */
    LocalScheme& local_scheme;
    /** Genetic Local Search optional parameters. */
    GeneticLocalSearchOptionalParameters<LocalScheme>& parameters;
    /** Population. */
    Population<LocalScheme> population;
    /** Genetic Local Search output structure. */
    GeneticLocalSearchOutput<LocalScheme>& output;

    /** Position storing which is the next initial solution to generate/add. */
    Counter initial_solution_pos = 0;
    /** Mutex to use when manipulating the structure. */
    std::mutex mutex;
};

template <typename LocalScheme>
inline void genetic_local_search_worker(
        GeneticLocalSearchData<LocalScheme>& data,
        Counter thread_id)
{
    LocalScheme local_scheme(data.local_scheme);
    std::mt19937_64 generator(data.parameters.seed + thread_id);

    // Generate initial solutions.
    for (;;) {
        data.mutex.lock();
        // No more initial solutions to generate.
        if (data.initial_solution_pos >= (Counter)data.parameters.initial_solution_ids.size()) {
            data.mutex.unlock();
            break;
        }
        Counter initial_solution_pos = data.initial_solution_pos;
        data.initial_solution_pos++;
        data.mutex.unlock();

        // Generate initial solution.
        auto solution_initial = local_scheme.initial_solution(
                data.parameters.initial_solution_ids[initial_solution_pos], generator);
        // Run A* Local Search.
        AStarLocalSearchOptionalParameters<LocalScheme> parameters_astar;
        parameters_astar.initial_solution_ids = {};
        parameters_astar.initial_solutions.push_back(solution_initial);
        parameters_astar.node_number_max = data.parameters.local_search_node_number_max;
        auto output_astar = a_star_local_search(local_scheme, parameters_astar);
        auto solution = output_astar.solution_pool.best();
        //std::cout << to_string(local_scheme.global_cost(solution)) << std::endl;

        // Lock mutex since we will modify the shared structure.
        data.mutex.lock();
        // Add the new solution to the population.
        data.population.add(solution, generator);
        // Check for a new best solution.
        if (local_scheme.global_cost(data.output.solution_pool.worst())
                > local_scheme.global_cost(solution)) {
            std::stringstream ss;
            ss << "initial solution " << initial_solution_pos
                << " (thread " << thread_id << ")";
            auto res = data.output.solution_pool.add(solution, ss, data.parameters.info);
            if (res == 2) {
                data.output.solution_pool.display(ss, data.parameters.info);
                data.parameters.new_solution_callback(solution);
            }
        }
        // Unlock mutex.
        data.mutex.unlock();
    }

    for (;;) {

        // Check time.
        if (!data.parameters.info.check_time())
            break;

        data.mutex.lock();

        // Check genetic iteration limit.
        if (data.parameters.genetic_iteration_number_max != -1
                && data.output.genetic_iteration_number
                >= data.parameters.genetic_iteration_number_max) {
            data.mutex.unlock();
            return;
        }

        // We need at least two solutions.
        if (data.population.size() < 4) {
            data.mutex.unlock();
            continue;
        }

        Counter genetic_iteration_number = data.output.genetic_iteration_number;
        data.output.genetic_iteration_number++;

        // Draw parent solutions.
        auto p = data.population.get_parents(generator);
        const auto& solution_parent_1 = p.first;
        const auto& solution_parent_2 = p.second;

        data.mutex.unlock();

        // Crossover.
        auto solution_child = local_scheme.crossover(
                solution_parent_1, solution_parent_2, generator);

        // Run A* Local Search.
        AStarLocalSearchOptionalParameters<LocalScheme> parameters_astar;
        parameters_astar.initial_solution_ids = {};
        parameters_astar.initial_solutions.push_back(solution_child);
        parameters_astar.node_number_max = data.parameters.local_search_node_number_max;
        auto output_astar = a_star_local_search(local_scheme, parameters_astar);
        auto solution = output_astar.solution_pool.best();
        //std::cout << to_string(local_scheme.global_cost(solution)) << std::endl;

        // Lock mutex since we will modify the shared structure.
        data.mutex.lock();
        // Add the new solution to the population.
        data.population.add(solution, generator);
        // Check for a new best solution.
        if (local_scheme.global_cost(data.output.solution_pool.worst())
                > local_scheme.global_cost(solution)) {
            std::stringstream ss;
            ss << "iteration " << genetic_iteration_number
                << " (thread " << thread_id << ")";
            int res = data.output.solution_pool.add(solution, ss, data.parameters.info);
            if (res == 2) {
                data.output.solution_pool.display(ss, data.parameters.info);
                data.parameters.new_solution_callback(solution);
            }
        }
        // Unlock mutex.
        data.mutex.unlock();
    }

}

template <typename LocalScheme>
inline GeneticLocalSearchOutput<LocalScheme> genetic_local_search(
        LocalScheme& local_scheme,
        GeneticLocalSearchOptionalParameters<LocalScheme> parameters)
{
    GeneticLocalSearchOutput<LocalScheme> output(local_scheme);
    output.solution_pool.display_init(parameters.info);
    std::vector<std::thread> threads;
    GeneticLocalSearchData<LocalScheme> data(local_scheme, parameters, output);
    for (Counter thread_id = 0; thread_id < parameters.thread_number; ++thread_id) {
        threads.push_back(std::thread(
                    genetic_local_search_worker<LocalScheme>,
                    std::ref(data),
                    thread_id));
    }
    for (Counter thread_id = 0; thread_id < (Counter)threads.size(); ++thread_id)
        threads[thread_id].join();

    output.solution_pool.display_end(parameters.info);
    return output;
}

}

