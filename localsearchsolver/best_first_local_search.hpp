#pragma once

#include "localsearchsolver/algorithm_formatter.hpp"

#include <unordered_set>
#include <thread>
#include <random>

namespace localsearchsolver
{

template <typename LocalScheme>
struct BestFirstLocalSearchParameters: Parameters<LocalScheme>
{
    using Solution = typename LocalScheme::Solution;
    using GlobalCost = typename LocalScheme::GlobalCost;

    /** Maximum number of nodes. */
    Counter maximum_number_of_nodes = -1;

    /** Number of threads running the algorithm independently in parallel. */
    Counter number_of_threads_1 = 1;

    /** Number of threads running on the same node pool in parallel. */
    Counter number_of_threads_2 = 1;


    virtual nlohmann::json to_json(
            const LocalScheme& local_scheme) const override
    {
        nlohmann::json json = Parameters<LocalScheme>::to_json(local_scheme);
        json.merge_patch({
            {"MaximumNumberOfNodes", maximum_number_of_nodes}});
        return json;
    }

    virtual int format_width() const override { return 26; }

    virtual void format(
            std::ostream& os,
            const LocalScheme& local_scheme) const override
    {
        Parameters<LocalScheme>::format(os, local_scheme);
        int width = format_width();
        os
            << std::setw(width) << std::left << "Maximum number of nodes: " << maximum_number_of_nodes << std::endl
            ;
    }
};

template <typename LocalScheme>
struct BestFirstLocalSearchOutput: Output<LocalScheme>
{
    /** Constructor. */
    BestFirstLocalSearchOutput(
            const LocalScheme& local_scheme,
            Counter maximum_size_of_the_solution_pool):
        Output<LocalScheme>(local_scheme, maximum_size_of_the_solution_pool) { }


    /** Number of nodes. */
    Counter number_of_nodes = 0;


    virtual nlohmann::json to_json() const override
    {
        nlohmann::json json = Output<LocalScheme>::to_json();
        json.merge_patch({
            {"NumberOfNodes", number_of_nodes}});
        return json;
    }

    virtual int format_width() const override { return 18; }

    virtual void format(std::ostream& os) const override
    {
        Output<LocalScheme>::format(os);
        int width = format_width();
        os
            << std::setw(width) << std::left << "Number of nodes: " << number_of_nodes << std::endl
            ;
    }
};

template <typename LocalScheme>
inline const BestFirstLocalSearchOutput<LocalScheme> best_first_local_search(
        LocalScheme& local_scheme,
        const BestFirstLocalSearchParameters<LocalScheme>& parameters = {});

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// Template implementations //////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename LocalScheme>
struct BestFirstLocalSearchNode
{
    using CompactSolution = typename LocalScheme::CompactSolution;
    using Perturbation = typename LocalScheme::Perturbation;

    std::shared_ptr<CompactSolution> compact_solution;
    std::vector<Perturbation> perturbations;
    Counter id = -1;
    Counter depth;
    Counter next_child_pos = 0;
    Counter child_id;
};

template <typename LocalScheme>
struct BestFirstLocalSearchNodeComparator
{
    BestFirstLocalSearchNodeComparator(const LocalScheme& branching_scheme):
        branching_scheme(branching_scheme) {  }

    const LocalScheme& branching_scheme;

    bool operator()(
            const std::shared_ptr<BestFirstLocalSearchNode<LocalScheme>>& node_1,
            const std::shared_ptr<BestFirstLocalSearchNode<LocalScheme>>& node_2) const {
        if (node_1->perturbations.back().global_cost
            != node_2->perturbations.back().global_cost)
            return node_1->perturbations.back().global_cost
                < node_2->perturbations.back().global_cost;
        return node_1->id < node_2->id;
    }
};

template <typename LocalScheme>
using CompactSolutionSet = std::unordered_set<
        std::shared_ptr<typename LocalScheme::CompactSolution>,
        typename LocalScheme::CompactSolutionHasher&,
        typename LocalScheme::CompactSolutionHasher&>;

template <typename LocalScheme>
using NodeSet = std::multiset<
        std::shared_ptr<BestFirstLocalSearchNode<LocalScheme>>,
        const BestFirstLocalSearchNodeComparator<LocalScheme>&>;

template <typename LocalScheme>
using PerturbationMap = std::unordered_map<
        typename LocalScheme::Perturbation,
        Counter,
        typename LocalScheme::PerturbationHasher&,
        typename LocalScheme::PerturbationHasher&>;

template <typename LocalScheme>
struct BestFirstLocalSearchData
{
    using CompactSolutionHasher = typename LocalScheme::CompactSolutionHasher;
    using PerturbationHasher = typename LocalScheme::PerturbationHasher;
    using Perturbation = typename LocalScheme::Perturbation;

    BestFirstLocalSearchData(
            LocalScheme& local_scheme,
            const BestFirstLocalSearchParameters<LocalScheme>& parameters,
            AlgorithmFormatter<LocalScheme>& algorithm_formatter,
            BestFirstLocalSearchOutput<LocalScheme>& output):
        local_scheme(local_scheme),
        parameters(parameters),
        algorithm_formatter(algorithm_formatter),
        output(output),
        compact_solution_hasher(local_scheme.compact_solution_hasher()),
        history{0, compact_solution_hasher, compact_solution_hasher},
        q(local_scheme),
        perturbation_hasher(local_scheme.perturbation_hasher()),
        perturbation_nodes{0, perturbation_hasher, perturbation_hasher}
        {  }

    LocalScheme& local_scheme;
    const BestFirstLocalSearchParameters<LocalScheme>& parameters;
    AlgorithmFormatter<LocalScheme>& algorithm_formatter;
    BestFirstLocalSearchOutput<LocalScheme>& output;
    CompactSolutionHasher compact_solution_hasher;
    CompactSolutionSet<LocalScheme> history;
    NodeSet<LocalScheme> q;
    Counter number_of_working_threads = 0;
    Counter initial_solution_pos = 0;
    Counter number_of_perturbations = -1;

    PerturbationHasher perturbation_hasher;

    /**
     * perturbation_nodes[perturbation] is last node at which perturbation "perturbation" has been used as a
     * perturbation.
     */
    PerturbationMap<LocalScheme> perturbation_nodes;

    /**
     * tabu_nodes[number_of_nodes] is a vector a nodes which are currently tabu
     * and should be released at node "number_of_nodes".
     */
    std::unordered_map<Counter, std::vector<std::shared_ptr<BestFirstLocalSearchNode<LocalScheme>>>> tabu_nodes;

    std::mutex mutex;
};

template <typename LocalScheme>
std::vector<typename LocalScheme::Perturbation> update_children(
        LocalScheme& local_scheme,
        typename LocalScheme::Solution& solution,
        std::mt19937_64& generator,
        Counter next_child_pos = 0)
{
    using Perturbation = typename LocalScheme::Perturbation;

    auto perturbation_compare = [](const Perturbation& perturbation_1, const Perturbation& perturbation_2) -> bool
    {
        return perturbation_1.global_cost < perturbation_2.global_cost;
    };

    // Get perturbation perturbations.
    auto perturbations_0 = local_scheme.perturbations(solution, generator);
    if (next_child_pos >= (Counter)perturbations_0.size())
        return {};
    // Sort perturbations.
    std::sort(perturbations_0.begin(), perturbations_0.end(), perturbation_compare);
    std::vector<Perturbation> perturbations;
    for (auto it = perturbations_0.begin() + next_child_pos;
            it != perturbations_0.end() && it != perturbations_0.begin() + next_child_pos + 1024;
            ++it)
        perturbations.push_back(*it);
    std::reverse(perturbations.begin(), perturbations.end());
    return perturbations;
}

template <typename LocalScheme>
inline void best_first_local_search_worker(
        BestFirstLocalSearchData<LocalScheme>& data,
        Counter thread_id)
{
    using CompactSolution = typename LocalScheme::CompactSolution;
    //std::cout << "best_first_local_search_worker start" << std::endl;

    LocalScheme local_scheme_tmp(data.local_scheme);
    LocalScheme& local_scheme = (data.parameters.number_of_threads_1 == 1 && data.parameters.number_of_threads_2 == 1)?
        data.local_scheme: local_scheme_tmp;
    Seed seed = data.parameters.seed + thread_id;
    std::mt19937_64 generator(seed);

    // Generate initial solutions.
    for (;;) {

        // Check end.
        if (data.parameters.timer.needs_to_end())
            break;

        // Check goal.
        if (data.parameters.has_goal
                && data.output.solution_pool.size() > 0
                && !strictly_better(
                    local_scheme,
                    data.parameters.goal,
                    local_scheme.global_cost(data.output.solution_pool.best())))
            break;

        data.mutex.lock();
        // Store initial_solution_pos to use it after releasing the mutex.
        Counter initial_solution_pos = data.initial_solution_pos;
        // No more initial solutions to generate.
        if (initial_solution_pos
                >= (Counter)data.parameters.initial_solution_ids.size()
                + (Counter)data.parameters.initial_solutions.size()) {
            data.mutex.unlock();
            break;
        }
        data.initial_solution_pos++;
        data.number_of_working_threads++;
        data.mutex.unlock();

        auto solution = (initial_solution_pos < (Counter)data.parameters.initial_solution_ids.size())?
            local_scheme.initial_solution(data.parameters.initial_solution_ids[initial_solution_pos], generator):
            data.parameters.initial_solutions[initial_solution_pos - (Counter)data.parameters.initial_solution_ids.size()];
        local_scheme.local_search(solution, generator);

        // Check for a new best solution.
        if (data.output.solution_pool.size() == 0
                || strictly_better(
                    local_scheme,
                    local_scheme.global_cost(solution),
                    local_scheme.global_cost(data.output.solution_pool.worst()))) {
            std::stringstream ss;
            ss << "s" << initial_solution_pos
                << " (t" << thread_id << ")";
            data.algorithm_formatter.update_solution(solution, ss);
        }

        // Get perturbation perturbations.
        auto perturbations = update_children(local_scheme, solution, generator);
        auto compact_solution = std::shared_ptr<CompactSolution>(
                new CompactSolution(local_scheme.solution2compact(solution)));

        data.mutex.lock();
        if (data.number_of_perturbations == -1)
            data.number_of_perturbations = local_scheme.perturbations(solution, generator).size();
        if (data.history.find(compact_solution) == data.history.end()
                && perturbations.size() > 0) {
            BestFirstLocalSearchNode<LocalScheme> root;
            root.id = -1;
            root.compact_solution = compact_solution;
            root.perturbations = perturbations;
            root.child_id = 0;
            root.depth = 1;
            data.history.insert(compact_solution);
            data.q.insert(std::shared_ptr<BestFirstLocalSearchNode<LocalScheme>>(
                        new BestFirstLocalSearchNode<LocalScheme>(root)));
        }
        data.number_of_working_threads--;
        data.mutex.unlock();
    }

    for (;;) {

        // Check end.
        if (data.parameters.timer.needs_to_end())
            break;

        // Check goal.
        if (data.parameters.has_goal
                && data.output.solution_pool.size() > 0
                && !strictly_better(
                    local_scheme,
                    data.parameters.goal,
                    local_scheme.global_cost(data.output.solution_pool.best()))) {
            break;
        }

        data.mutex.lock();

        // Check for algorithm end.
        if (data.q.empty() && data.number_of_working_threads == 0) {
            data.mutex.unlock();
            break;
        }

        if (data.q.empty()) {
            data.mutex.unlock();
            continue;
        }

        // Check node limit.
        if (data.parameters.maximum_number_of_nodes != -1
                && data.output.number_of_nodes >= data.parameters.maximum_number_of_nodes) {
            data.mutex.unlock();
            return;
        }

        Counter node_id = data.output.number_of_nodes;
        data.output.number_of_nodes++;

        // Compute tabu tenure.
        Counter tabu_tenure = -1;
        if (node_id < data.number_of_perturbations) {
            tabu_tenure = data.number_of_perturbations / 10;
        } else if (node_id < 10 * data.number_of_perturbations) {
            tabu_tenure = data.number_of_perturbations / 2;
        } else {
            tabu_tenure = data.number_of_perturbations;
        }

        // Add back nodes which are not tabu anymore.
        if (data.tabu_nodes.find(node_id) != data.tabu_nodes.end()) {
            for (auto node: data.tabu_nodes[node_id])
                data.q.insert(node);
            data.tabu_nodes.erase(node_id);
        }

        // Check goal.
        if (data.parameters.has_goal
                && data.output.solution_pool.size() > 0
                && !strictly_better(
                    local_scheme,
                    data.parameters.goal,
                    local_scheme.global_cost(data.output.solution_pool.best()))) {
            break;
        }

        // Draw next node from the queue.
        auto node_cur = *data.q.begin();
        data.q.erase(data.q.begin());
        auto perturbation = node_cur->perturbations.back();
        // Check if the perturbation is tabu.
        while (data.perturbation_hasher.hashable(perturbation)
                && data.perturbation_nodes.find(perturbation) != data.perturbation_nodes.end()
                && data.perturbation_nodes[perturbation] > node_id - tabu_tenure) {
            // Add a new node for this perturbation to the tabu node set.
            BestFirstLocalSearchNode<LocalScheme> node_tmp;
            node_tmp.compact_solution = node_cur->compact_solution;
            node_tmp.perturbations = {perturbation};
            node_tmp.id = node_cur->id;
            node_tmp.depth = node_cur->depth;
            node_tmp.child_id = node_cur->child_id;
            node_tmp.next_child_pos = -1;
            if (data.tabu_nodes.find(data.perturbation_nodes[perturbation] + tabu_tenure) == data.tabu_nodes.end())
                data.tabu_nodes[data.perturbation_nodes[perturbation] + tabu_tenure] = {};
            data.tabu_nodes[data.perturbation_nodes[perturbation] + tabu_tenure].push_back(
                    std::shared_ptr<BestFirstLocalSearchNode<LocalScheme>>(
                        new BestFirstLocalSearchNode<LocalScheme>(node_tmp)));

            // Remove the perturbation from the father's children.
            node_cur->perturbations.pop_back();
            if (node_cur->next_child_pos != -1) {
                node_cur->next_child_pos++;
                if (node_cur->perturbations.empty()) {
                    auto solution_tmp = local_scheme.compact2solution(*node_cur->compact_solution);
                    node_cur->perturbations = update_children(
                            local_scheme,
                            solution_tmp,
                            generator,
                            node_cur->next_child_pos);
                }
            }
            if (!node_cur->perturbations.empty())
                data.q.insert(node_cur);

            node_cur = nullptr;
            if (data.q.empty())
                break;

            node_cur = *data.q.begin();
            data.q.erase(data.q.begin());
            perturbation = node_cur->perturbations.back();
        }
        if (node_cur == nullptr) {
            data.mutex.unlock();
            continue;
        }
        if (data.perturbation_hasher.hashable(perturbation))
            data.perturbation_nodes[perturbation] = node_id;
        data.number_of_working_threads++;
        //std::cout << "node " << node_id
        //    << " depth " << node_cur->depth
        //    << " cost " << to_string(local_scheme, perturbation.global_cost)
        //    << " thread " << thread_id
        //    << std::endl;
        data.mutex.unlock();

        auto solution = local_scheme.compact2solution(*node_cur->compact_solution);
        // Remove the perturbation from the father's children.
        node_cur->perturbations.pop_back();
        if (node_cur->next_child_pos != -1) {
            node_cur->next_child_pos++;
            if (node_cur->perturbations.empty())
                node_cur->perturbations = update_children(
                        local_scheme,
                        solution,
                        generator,
                        node_cur->next_child_pos);
        }
        // Apply perturbation and local search.
        local_scheme.apply_perturbation(solution, perturbation, generator);
        local_scheme.local_search(solution, generator, perturbation);

        // Check for a new best solution.
        if (data.output.solution_pool.size() == 0
                || strictly_better(
                    local_scheme,
                    local_scheme.global_cost(solution),
                    local_scheme.global_cost(data.output.solution_pool.worst()))) {
            std::stringstream ss;
            ss << "n" << node_id
                << " d" << node_cur->depth
                << " c" << node_cur->next_child_pos - 1
                << " (t" << thread_id << ")";
            data.algorithm_formatter.update_solution(solution, ss);
        }

        // Get perturbation perturbations.
        auto perturbations = update_children(local_scheme, solution, generator);
        auto compact_solution = std::shared_ptr<CompactSolution>(
                    new CompactSolution(local_scheme.solution2compact(solution)));

        data.mutex.lock();

        if (data.history.find(compact_solution) == data.history.end()
                && !perturbations.empty()) {
            BestFirstLocalSearchNode<LocalScheme> child;
            child.compact_solution = compact_solution;
            child.perturbations = perturbations;
            child.id = node_id;
            child.depth = node_cur->depth + 1;
            child.child_id = node_cur->next_child_pos - 1;
            data.history.insert(compact_solution);
            data.q.insert(std::shared_ptr<BestFirstLocalSearchNode<LocalScheme>>(
                        new BestFirstLocalSearchNode<LocalScheme>(child)));
        }

        if (!node_cur->perturbations.empty())
            data.q.insert(node_cur);

        data.number_of_working_threads--;
        data.mutex.unlock();
    }
}

template <typename LocalScheme>
inline const BestFirstLocalSearchOutput<LocalScheme> best_first_local_search(
        LocalScheme& local_scheme,
        const BestFirstLocalSearchParameters<LocalScheme>& parameters)
{
    BestFirstLocalSearchOutput<LocalScheme> output(
            local_scheme,
            parameters.maximum_size_of_the_solution_pool);
    AlgorithmFormatter<LocalScheme> algorithm_formatter(local_scheme, parameters, output);
    algorithm_formatter.start("Best first local search");
    algorithm_formatter.print_header();

    std::vector<std::thread> threads;
    std::vector<std::shared_ptr<BestFirstLocalSearchData<LocalScheme>>> datas;
    Counter thread_id = 0;
    for (Counter thread_id_1 = 0; thread_id_1 < parameters.number_of_threads_1; ++thread_id_1) {
        datas.push_back(std::shared_ptr<BestFirstLocalSearchData<LocalScheme>>(
                    new BestFirstLocalSearchData<LocalScheme>(
                        local_scheme,
                        parameters,
                        algorithm_formatter,
                        output)));
        for (Counter thread_id_2 = 0; thread_id_2 < parameters.number_of_threads_2; ++thread_id_2) {
            threads.push_back(std::thread(
                        best_first_local_search_worker<LocalScheme>,
                        std::ref(*datas[thread_id_1]),
                        thread_id));
            thread_id++;
        }
    }
    for (Counter thread_id = 0; thread_id < (Counter)threads.size(); ++thread_id)
        threads[thread_id].join();

    algorithm_formatter.end();
    return output;
}

}
