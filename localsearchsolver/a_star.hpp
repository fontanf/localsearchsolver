#pragma once

#include "localsearchsolver/common.hpp"

#include <unordered_set>
#include <thread>

namespace localsearchsolver
{

template <typename LocalScheme>
using AStarCallback = std::function<void(const typename LocalScheme::Solution&)>;

template <typename LocalScheme>
struct AStarOptionalParameters
{
    typedef typename LocalScheme::Solution Solution;

    Counter node_number_max = -1;
    Counter thread_number_1 = 1;
    Counter thread_number_2 = 1;
    std::vector<Counter> initial_solution_ids = {0};
    Seed seed = 0;

    AStarCallback<LocalScheme> new_solution_callback
        = [](const Solution& solution) { (void)solution; };

    optimizationtools::Info info;
};

template <typename LocalScheme>
struct AStarOutput
{
    AStarOutput(
            const LocalScheme& local_scheme):
        solution_pool(local_scheme, 1) { }

    SolutionPool<LocalScheme> solution_pool;
};

template <typename LocalScheme>
struct AStarNode
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
struct AStarNodeComparator
{
    AStarNodeComparator(const LocalScheme& branching_scheme):
        branching_scheme(branching_scheme) {  }

    const LocalScheme& branching_scheme;

    bool operator()(
            const std::shared_ptr<AStarNode<LocalScheme>>& node_1,
            const std::shared_ptr<AStarNode<LocalScheme>>& node_2) const {
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
        std::shared_ptr<AStarNode<LocalScheme>>,
        const AStarNodeComparator<LocalScheme>&>;

template <typename LocalScheme>
using MoveMap = std::unordered_map<
        typename LocalScheme::Move,
        Counter,
        typename LocalScheme::MoveHasher&,
        typename LocalScheme::MoveHasher&>;

template <typename LocalScheme>
struct AStarData
{
    typedef typename LocalScheme::CompactSolutionHasher CompactSolutionHasher;
    typedef typename LocalScheme::MoveHasher MoveHasher;
    typedef typename LocalScheme::Move Move;

    AStarData(
            LocalScheme& local_scheme,
            AStarOptionalParameters<LocalScheme>& parameters,
            AStarOutput<LocalScheme>& output):
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
    AStarOptionalParameters<LocalScheme>& parameters;
    AStarOutput<LocalScheme>& output;
    CompactSolutionHasher compact_solution_hasher;
    CompactSolutionSet<LocalScheme> history;
    NodeSet<LocalScheme> q;
    Counter working_thread_number = 0;
    Counter initial_solution_pos = 0;
    Counter node_number = 0;
    Counter perturbation_number = -1;

    MoveHasher move_hasher;

    /**
     * move_nodes[move] is last node at which move "move" has been used as a
     * perturbation.
     */
    MoveMap<LocalScheme> move_nodes;

    /**
     * tabu_nodes[node_number] is a vector a nodes which are currently tabu
     * and should be released at node "node_number".
     */
    std::unordered_map<Counter, std::vector<std::shared_ptr<AStarNode<LocalScheme>>>> tabu_nodes;

    std::mutex mutex;
};

template <typename LocalScheme>
std::vector<typename LocalScheme::Move> update_children(
        LocalScheme& local_scheme,
        typename LocalScheme::Solution& solution,
        Counter next_child_pos = 0)
{
    typedef typename LocalScheme::Move Move;

    auto move_compare = [](const Move& move_1, const Move& move_2) -> bool
    {
        return move_1.global_cost < move_2.global_cost;
    };

    // Get perturbation moves.
    auto moves_0 = local_scheme.perturbations(solution);
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
inline void a_star_worker(
        AStarData<LocalScheme>& data,
        Counter thread_id)
{
    typedef typename LocalScheme::CompactSolution CompactSolution;

    LocalScheme local_scheme(data.local_scheme);
    std::mt19937_64 generator(data.parameters.seed + thread_id);

    for (;;) {
        data.mutex.lock();
        // No more initial solutions to generate.
        if (data.initial_solution_pos >= (Counter)data.parameters.initial_solution_ids.size()) {
            data.mutex.unlock();
            break;
        }
        Counter initial_solution_id = data.parameters.initial_solution_ids[data.initial_solution_pos];
        data.initial_solution_pos++;
        data.working_thread_number++;
        data.mutex.unlock();

        auto solution = local_scheme.initial_solution(initial_solution_id, generator);
        local_scheme.local_search(solution, generator);

        // Check for a new best solution.
        if (local_scheme.global_cost(data.output.solution_pool.best())
                > local_scheme.global_cost(solution)) {
            std::stringstream ss;
            ss << "initial solution " << initial_solution_id
                << " (thread " << thread_id << ")";
            data.output.solution_pool.add(solution, ss, data.parameters.info);
            data.output.solution_pool.display(ss, data.parameters.info);
            data.parameters.new_solution_callback(solution);
        }

        // Get perturbation moves.
        auto moves = update_children(local_scheme, solution);
        auto compact_solution = std::shared_ptr<CompactSolution>(
                new CompactSolution(local_scheme.solution2compact(solution)));

        data.mutex.lock();
        if (data.perturbation_number == -1)
            data.perturbation_number = local_scheme.perturbations(solution).size();
        if (data.history.find(compact_solution) == data.history.end()
                && moves.size() > 0) {
            AStarNode<LocalScheme> root;
            root.compact_solution = compact_solution;
            root.perturbations = moves;
            root.child_id = 0;
            root.depth = 1;
            data.history.insert(compact_solution);
            data.q.insert(std::shared_ptr<AStarNode<LocalScheme>>(
                        new AStarNode<LocalScheme>(root)));
        }
        data.working_thread_number--;
        data.mutex.unlock();
    }

    Counter tabu_tenure = data.perturbation_number / 10;

    for (;;) {

        // Check time.
        if (!data.parameters.info.check_time())
            break;

        data.mutex.lock();

        // Check node limit.
        if (data.parameters.node_number_max != -1
                && data.node_number >= data.parameters.node_number_max) {
            data.mutex.unlock();
            return;
        }

        Counter node_id = data.node_number;
        data.node_number++;

        // Add back nodes which are not tabu anymore.
        if (data.tabu_nodes.find(data.node_number) != data.tabu_nodes.end()) {
            for (auto node: data.tabu_nodes[data.node_number])
                data.q.insert(node);
            data.tabu_nodes.erase(data.node_number);
        }

        if (data.q.empty() && data.working_thread_number == 0) {
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
        while (data.move_nodes.find(move) != data.move_nodes.end()
                && data.move_nodes[move] > data.node_number - tabu_tenure) {
            // Add a new node for this perturbation to the tabu node set.
            AStarNode<LocalScheme> node_tmp;
            node_tmp.compact_solution = node_cur->compact_solution;
            node_tmp.perturbations = {move};
            node_tmp.depth = node_cur->depth;
            node_tmp.child_id = node_cur->child_id;
            node_tmp.next_child_pos = -1;
            if (data.tabu_nodes.find(data.move_nodes[move] + tabu_tenure) == data.tabu_nodes.end())
                data.tabu_nodes[data.move_nodes[move] + tabu_tenure] = {};
            data.tabu_nodes[data.move_nodes[move] + tabu_tenure].push_back(
                    std::shared_ptr<AStarNode<LocalScheme>>(
                        new AStarNode<LocalScheme>(node_tmp)));

            // Remove the perturbation from the father's children.
            node_cur->perturbations.pop_back();
            if (node_cur->next_child_pos != -1) {
                node_cur->next_child_pos++;
                if (node_cur->perturbations.empty()) {
                    auto solution_tmp = local_scheme.compact2solution(*node_cur->compact_solution);
                    node_cur->perturbations = update_children(
                            local_scheme,
                            solution_tmp,
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
        data.move_nodes[move] = data.node_number;
        data.working_thread_number++;
        //std::cout << "node " << data.node_number
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
                        node_cur->next_child_pos);
        }
        // Apply perturbation and local search.
        local_scheme.apply_move(solution, move);
        local_scheme.local_search(solution, generator, move);

        // Check for a new best solution.
        if (local_scheme.global_cost(data.output.solution_pool.best())
                > local_scheme.global_cost(solution)) {
            std::stringstream ss;
            ss << "node " << node_id
                << " depth " << node_cur->depth
                << " child " << node_cur->next_child_pos - 1
                << " (thread " << thread_id << ")";
            data.output.solution_pool.add(solution, ss, data.parameters.info);
            data.output.solution_pool.display(ss, data.parameters.info);
            data.parameters.new_solution_callback(solution);
        }

        // Get perturbation moves.
        auto moves = update_children(local_scheme, solution);
        auto compact_solution = std::shared_ptr<CompactSolution>(
                    new CompactSolution(local_scheme.solution2compact(solution)));

        data.mutex.lock();

        if (data.history.find(compact_solution) == data.history.end()
                && !moves.empty()) {
            AStarNode<LocalScheme> child;
            child.compact_solution = compact_solution;
            child.perturbations = moves;
            child.depth = node_cur->depth + 1;
            child.child_id = node_cur->next_child_pos - 1;
            data.history.insert(compact_solution);
            data.q.insert(std::shared_ptr<AStarNode<LocalScheme>>(
                        new AStarNode<LocalScheme>(child)));
        }

        if (!node_cur->perturbations.empty())
            data.q.insert(node_cur);

        data.working_thread_number--;
        data.mutex.unlock();
    }
}

template <typename LocalScheme>
inline AStarOutput<LocalScheme> a_star(
        LocalScheme& local_scheme,
        AStarOptionalParameters<LocalScheme> parameters = {})
{
    AStarOutput<LocalScheme> output(local_scheme);
    output.solution_pool.display_init(parameters.info);
    std::vector<std::thread> threads;
    std::vector<std::shared_ptr<AStarData<LocalScheme>>> datas;
    Counter thread_id = 0;
    for (Counter thread_id_1 = 0; thread_id_1 < parameters.thread_number_1; ++thread_id_1) {
        datas.push_back(std::shared_ptr<AStarData<LocalScheme>>(
                    new AStarData<LocalScheme>(local_scheme, parameters, output)));
        for (Counter thread_id_2 = 0; thread_id_2 < parameters.thread_number_2; ++thread_id_2) {
            threads.push_back(std::thread(
                        a_star_worker<LocalScheme>,
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

