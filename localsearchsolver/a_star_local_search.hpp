#pragma once

#include "localsearchsolver/common.hpp"

#include <unordered_set>
#include <thread>

namespace localsearchsolver
{

template <typename LocalScheme>
using AStarLocalSearchCallback = std::function<void(const typename LocalScheme::Solution&)>;

template <typename LocalScheme>
struct AStarLocalSearchOptionalParameters
{
    typedef typename LocalScheme::Solution Solution;

    /** Maximum number of nodes. */
    Counter maximum_number_of_nodes = -1;
    /** Number of threads running the algorithm independently in parallel. */
    Counter number_of_threads_1 = 1;
    /** Number of threads running on the same node pool in parallel. */
    Counter number_of_threads_2 = 1;
    /** Ids of generated initial solutions. */
    std::vector<Counter> initial_solution_ids = {0};
    /** User-provided initial solutions. */
    std::vector<Solution> initial_solutions;
    /** Seed. */
    Seed seed = 0;
    /** Callback function called when a new best solution is found. */
    AStarLocalSearchCallback<LocalScheme> new_solution_callback
        = [](const Solution& solution) { (void)solution; };

    optimizationtools::Info info;
};

template <typename LocalScheme>
struct AStarLocalSearchOutput
{
    /** Constructor. */
    AStarLocalSearchOutput(
            const LocalScheme& local_scheme):
        solution_pool(local_scheme, 1) { }

    /** Solution pool. */
    SolutionPool<LocalScheme> solution_pool;
};

template <typename LocalScheme>
inline AStarLocalSearchOutput<LocalScheme> a_star_local_search(
        LocalScheme& local_scheme,
        AStarLocalSearchOptionalParameters<LocalScheme> parameters = {});

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// Template implementations //////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename LocalScheme>
struct AStarLocalSearchNode
{
    typedef typename LocalScheme::CompactSolution CompactSolution;
    typedef typename LocalScheme::Move Move;

    std::shared_ptr<CompactSolution> compact_solution;
    std::vector<Move> perturbations;
    Counter depth;
    Counter next_child_pos = 0;
    Counter child_id;
};

template <typename LocalScheme>
struct AStarLocalSearchNodeComparator
{
    AStarLocalSearchNodeComparator(const LocalScheme& branching_scheme):
        branching_scheme(branching_scheme) {  }

    const LocalScheme& branching_scheme;

    bool operator()(
            const std::shared_ptr<AStarLocalSearchNode<LocalScheme>>& node_1,
            const std::shared_ptr<AStarLocalSearchNode<LocalScheme>>& node_2) const {
        return node_1->perturbations.back().global_cost
            < node_2->perturbations.back().global_cost;
    }
};

template <typename LocalScheme>
using CompactSolutionSet = std::unordered_set<
        std::shared_ptr<typename LocalScheme::CompactSolution>,
        typename LocalScheme::CompactSolutionHasher&,
        typename LocalScheme::CompactSolutionHasher&>;

template <typename LocalScheme>
using NodeSet = std::multiset<
        std::shared_ptr<AStarLocalSearchNode<LocalScheme>>,
        const AStarLocalSearchNodeComparator<LocalScheme>&>;

template <typename LocalScheme>
using MoveMap = std::unordered_map<
        typename LocalScheme::Move,
        Counter,
        typename LocalScheme::MoveHasher&,
        typename LocalScheme::MoveHasher&>;

template <typename LocalScheme>
struct AStarLocalSearchData
{
    typedef typename LocalScheme::CompactSolutionHasher CompactSolutionHasher;
    typedef typename LocalScheme::MoveHasher MoveHasher;
    typedef typename LocalScheme::Move Move;

    AStarLocalSearchData(
            LocalScheme& local_scheme,
            AStarLocalSearchOptionalParameters<LocalScheme>& parameters,
            AStarLocalSearchOutput<LocalScheme>& output):
        local_scheme(local_scheme),
        parameters(parameters),
        output(output),
        compact_solution_hasher(local_scheme.compact_solution_hasher()),
        history{0, compact_solution_hasher, compact_solution_hasher},
        q(local_scheme),
        move_hasher(local_scheme.move_hasher()),
        move_nodes{0, move_hasher, move_hasher}
        {  }

    LocalScheme& local_scheme;
    AStarLocalSearchOptionalParameters<LocalScheme>& parameters;
    AStarLocalSearchOutput<LocalScheme>& output;
    CompactSolutionHasher compact_solution_hasher;
    CompactSolutionSet<LocalScheme> history;
    NodeSet<LocalScheme> q;
    Counter number_of_working_threads = 0;
    Counter initial_solution_pos = 0;
    Counter number_of_perturbations = -1;

    MoveHasher move_hasher;

    /**
     * move_nodes[move] is last node at which move "move" has been used as a
     * perturbation.
     */
    MoveMap<LocalScheme> move_nodes;

    /**
     * tabu_nodes[number_of_nodes] is a vector a nodes which are currently tabu
     * and should be released at node "number_of_nodes".
     */
    std::unordered_map<Counter, std::vector<std::shared_ptr<AStarLocalSearchNode<LocalScheme>>>> tabu_nodes;

    /** Number of nodes explored. */
    Counter number_of_nodes = 0;

    std::mutex mutex;
};

template <typename LocalScheme>
std::vector<typename LocalScheme::Move> update_children(
        LocalScheme& local_scheme,
        typename LocalScheme::Solution& solution,
        std::mt19937_64& generator,
        Counter next_child_pos = 0)
{
    typedef typename LocalScheme::Move Move;

    auto move_compare = [](const Move& move_1, const Move& move_2) -> bool
    {
        return move_1.global_cost < move_2.global_cost;
    };

    // Get perturbation moves.
    auto moves_0 = local_scheme.perturbations(solution, generator);
    if (next_child_pos >= (Counter)moves_0.size())
        return {};
    // Update move costs.
    auto global_cost_cur = local_scheme.global_cost(solution);
    for (Move& move: moves_0)
        move.global_cost = update_move_cost(move.global_cost, global_cost_cur);
    // Sort moves.
    std::sort(moves_0.begin(), moves_0.end(), move_compare);
    std::vector<Move> moves;
    for (auto it = moves_0.begin() + next_child_pos;
            it != moves_0.end() && it != moves_0.begin() + next_child_pos + 1024;
            ++it)
        moves.push_back(*it);
    std::reverse(moves.begin(), moves.end());
    return moves;
}

template <typename LocalScheme>
inline void a_star_local_search_worker(
        AStarLocalSearchData<LocalScheme>& data,
        Counter thread_id)
{
    typedef typename LocalScheme::CompactSolution CompactSolution;
    //std::cout << "a_star_local_search_worker start" << std::endl;

    LocalScheme local_scheme(data.local_scheme);
    Seed seed = data.parameters.seed + thread_id;
    std::mt19937_64 generator(seed);

    // Generate initial solutions.
    for (;;) {
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
                || local_scheme.global_cost(data.output.solution_pool.worst())
                > local_scheme.global_cost(solution)) {
            std::stringstream ss;
            ss << "initial solution " << initial_solution_pos
                << " (thread " << thread_id << ")";
            int res = data.output.solution_pool.add(solution, ss, data.parameters.info);
            if (res == 2) {
                data.output.solution_pool.display(ss, data.parameters.info);
                data.parameters.new_solution_callback(solution);
            }
        }

        // Get perturbation moves.
        auto moves = update_children(local_scheme, solution, generator);
        auto compact_solution = std::shared_ptr<CompactSolution>(
                new CompactSolution(local_scheme.solution2compact(solution)));

        data.mutex.lock();
        if (data.number_of_perturbations == -1)
            data.number_of_perturbations = local_scheme.perturbations(solution, generator).size();
        if (data.history.find(compact_solution) == data.history.end()
                && moves.size() > 0) {
            AStarLocalSearchNode<LocalScheme> root;
            root.compact_solution = compact_solution;
            root.perturbations = moves;
            root.child_id = 0;
            root.depth = 1;
            data.history.insert(compact_solution);
            data.q.insert(std::shared_ptr<AStarLocalSearchNode<LocalScheme>>(
                        new AStarLocalSearchNode<LocalScheme>(root)));
        }
        data.number_of_working_threads--;
        data.mutex.unlock();
    }

    for (;;) {

        // Check end.
        if (data.parameters.info.needs_to_end())
            break;

        data.mutex.lock();

        // Check node limit.
        if (data.parameters.maximum_number_of_nodes != -1
                && data.number_of_nodes >= data.parameters.maximum_number_of_nodes) {
            data.mutex.unlock();
            return;
        }

        Counter node_id = data.number_of_nodes;
        data.number_of_nodes++;

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

        // Check for algorithm end.
        if (data.q.empty() && data.number_of_working_threads == 0) {
            data.mutex.unlock();
            break;
        }

        if (data.q.empty()) {
            data.mutex.unlock();
            continue;
        }

        // Draw next node from the queue.
        auto node_cur = *data.q.begin();
        data.q.erase(data.q.begin());
        auto move = node_cur->perturbations.back();
        // Check if the perturbation is tabu.
        while (data.move_hasher.hashable(move)
                && data.move_nodes.find(move) != data.move_nodes.end()
                && data.move_nodes[move] > node_id - tabu_tenure) {
            // Add a new node for this perturbation to the tabu node set.
            AStarLocalSearchNode<LocalScheme> node_tmp;
            node_tmp.compact_solution = node_cur->compact_solution;
            node_tmp.perturbations = {move};
            node_tmp.depth = node_cur->depth;
            node_tmp.child_id = node_cur->child_id;
            node_tmp.next_child_pos = -1;
            if (data.tabu_nodes.find(data.move_nodes[move] + tabu_tenure) == data.tabu_nodes.end())
                data.tabu_nodes[data.move_nodes[move] + tabu_tenure] = {};
            data.tabu_nodes[data.move_nodes[move] + tabu_tenure].push_back(
                    std::shared_ptr<AStarLocalSearchNode<LocalScheme>>(
                        new AStarLocalSearchNode<LocalScheme>(node_tmp)));

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
            move = node_cur->perturbations.back();
        }
        if (node_cur == nullptr) {
            data.mutex.unlock();
            continue;
        }
        if (data.move_hasher.hashable(move))
            data.move_nodes[move] = node_id;
        data.number_of_working_threads++;
        //std::cout << "node " << node_id
        //    << " depth " << node_cur->depth
        //    << " cost " << to_string(move.global_cost)
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
        local_scheme.apply_move(solution, move);
        local_scheme.local_search(solution, generator, move);

        // Check for a new best solution.
        if (data.output.solution_pool.size() == 0
                || local_scheme.global_cost(data.output.solution_pool.worst())
                > local_scheme.global_cost(solution)) {
            std::stringstream ss;
            ss << "node " << node_id
                << " depth " << node_cur->depth
                << " child " << node_cur->next_child_pos - 1
                << " (thread " << thread_id << ")";
            int res = data.output.solution_pool.add(solution, ss, data.parameters.info);
            if (res == 2) {
                data.output.solution_pool.display(ss, data.parameters.info);
                data.parameters.new_solution_callback(solution);
            }
        }

        // Get perturbation moves.
        auto moves = update_children(local_scheme, solution, generator);
        auto compact_solution = std::shared_ptr<CompactSolution>(
                    new CompactSolution(local_scheme.solution2compact(solution)));

        data.mutex.lock();

        if (data.history.find(compact_solution) == data.history.end()
                && !moves.empty()) {
            AStarLocalSearchNode<LocalScheme> child;
            child.compact_solution = compact_solution;
            child.perturbations = moves;
            child.depth = node_cur->depth + 1;
            child.child_id = node_cur->next_child_pos - 1;
            data.history.insert(compact_solution);
            data.q.insert(std::shared_ptr<AStarLocalSearchNode<LocalScheme>>(
                        new AStarLocalSearchNode<LocalScheme>(child)));
        }

        if (!node_cur->perturbations.empty())
            data.q.insert(node_cur);

        data.number_of_working_threads--;
        data.mutex.unlock();
    }
}

template <typename LocalScheme>
inline AStarLocalSearchOutput<LocalScheme> a_star_local_search(
        LocalScheme& local_scheme,
        AStarLocalSearchOptionalParameters<LocalScheme> parameters)
{
    //std::cout << "a_star_local_search start" << std::endl;
    AStarLocalSearchOutput<LocalScheme> output(local_scheme);
    output.solution_pool.display_init(parameters.info);
    std::vector<std::thread> threads;
    std::vector<std::shared_ptr<AStarLocalSearchData<LocalScheme>>> datas;
    Counter thread_id = 0;
    for (Counter thread_id_1 = 0; thread_id_1 < parameters.number_of_threads_1; ++thread_id_1) {
        datas.push_back(std::shared_ptr<AStarLocalSearchData<LocalScheme>>(
                    new AStarLocalSearchData<LocalScheme>(local_scheme, parameters, output)));
        for (Counter thread_id_2 = 0; thread_id_2 < parameters.number_of_threads_2; ++thread_id_2) {
            threads.push_back(std::thread(
                        a_star_local_search_worker<LocalScheme>,
                        std::ref(*datas[thread_id_1]),
                        thread_id));
            thread_id++;
        }
    }
    for (Counter thread_id = 0; thread_id < (Counter)threads.size(); ++thread_id)
        threads[thread_id].join();

    output.solution_pool.display_end(parameters.info);
    return output;
}

}

