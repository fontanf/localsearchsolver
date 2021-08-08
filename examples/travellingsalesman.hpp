#pragma once

/**
 * Travelling Salesman Problem.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/travellingsalesman.hpp
 *
 * TODO
 *
 */

#include "localsearchsolver/common.hpp"

#include "orproblems/travellingsalesman.hpp"

namespace localsearchsolver
{

namespace travellingsalesman
{

using namespace orproblems::travellingsalesman;

class LocalScheme
{

public:

    /** Global cost: <Vertex number, Length>; */
    using GlobalCost = std::tuple<VertexId, Distance>;

    inline VertexId&       number_of_vertices(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Distance&                   length(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline VertexId  number_of_vertices(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Distance              length(const GlobalCost& global_cost) { return std::get<1>(global_cost); }

    static GlobalCost global_cost_worst()
    {
        return {
            std::numeric_limits<VertexId>::max(),
            std::numeric_limits<Distance>::max(),
        };
    }

    /*
     * Solutions.
     */

    /** Vector of size n + 1, starting and ending with the depot. */
    using CompactSolution = std::vector<VertexId>;

    struct CompactSolutionHasher
    {
        std::hash<VertexId> hasher;

        inline bool operator()(
                const std::shared_ptr<CompactSolution>& compact_solution_1,
                const std::shared_ptr<CompactSolution>& compact_solution_2) const
        {
            return *compact_solution_1 == *compact_solution_2;
        }

        inline std::size_t operator()(
                const std::shared_ptr<CompactSolution>& compact_solution) const
        {
            size_t hash = 0;
            for (VertexId j: *compact_solution)
                optimizationtools::hash_combine(hash, hasher(j));
            return hash;
        }
    };

    inline CompactSolutionHasher compact_solution_hasher() const { return CompactSolutionHasher(); }

    struct Solution
    {
        std::vector<VertexId> vertices;
        std::vector<VertexPos> positions;
        Distance length = 0;
    };

    CompactSolution solution2compact(const Solution& solution)
    {
        return solution.vertices;
    }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        auto solution = empty_solution();
        compute(solution, compact_solution);
        return solution;
    }

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        VertexPos oropt_size_max = 8;
        bool swap = true;
        bool twoopt = true;
        bool shuffle_neighborhood_order = true;
        Counter number_of_perturbations = 10;
    };

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters),
        positions1_(instance.number_of_vertices() - 1),
        positions2_(instance.number_of_vertices() - 1)
    {
        std::iota(positions1_.begin(), positions1_.end(), 0);
        std::iota(positions2_.begin(), positions2_.end(), 0);
    }

    LocalScheme(const LocalScheme& local_scheme):
        LocalScheme(local_scheme.instance_, local_scheme.parameters_) { }

    virtual ~LocalScheme() { }

    /*
     * Initial solutions.
     */

    inline Solution empty_solution() const
    {
        Solution solution;
        solution.positions.resize(instance_.number_of_vertices(), -1);
        return solution;
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator)
    {
        std::vector<VertexId> vertices(instance_.number_of_vertices() + 1);
        vertices[0] = 0;
        vertices[instance_.number_of_vertices()] = 0;
        std::iota(vertices.begin() + 1, vertices.end() - 1, 1);
        std::shuffle(vertices.begin() + 1, vertices.end() - 1, generator);
        return compact2solution(vertices);
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            -solution.vertices.size(),
            solution.length,
        };
    }

    /*
     * Local search.
     */

    struct Move
    {
        VertexPos pos_1;
        VertexPos pos_2;
        VertexPos pos_3;
        VertexPos pos_4;
        // New edges.
        VertexId j11;
        VertexId j12;
        VertexId j21;
        VertexId j22;
        VertexId j31;
        VertexId j32;
        VertexId j41;
        VertexId j42;
        GlobalCost global_cost;
    };

    static Move move_null()
    {
        return {
            -1, -1, -1, -1,
            -1, -1, -1, -1,
            -1, -1, -1, -1,
            global_cost_worst()};
    }

    struct MoveHasher
    {
        inline bool hashable(const Move&) const { return false; }
        inline bool operator()(const Move&, const Move&) const { return false; }
        inline std::size_t operator()(const Move&) const { return 0; }
    };

    inline MoveHasher move_hasher() const { return MoveHasher(); }

    inline std::vector<Move> perturbations(
            const Solution& solution,
            std::mt19937_64& generator)
    {
        std::vector<Move> moves;
        for (Counter perturbation = 0; perturbation < parameters_.number_of_perturbations; ++perturbation) {
            std::vector<VertexPos> edges = optimizationtools::bob_floyd<VertexPos>(
                    4, solution.vertices.size() - 1, generator);
            std::sort(edges.begin(), edges.end());
            Move move;
            move.pos_1 = edges[0];
            move.pos_2 = edges[1];
            move.pos_3 = edges[2];
            move.pos_4 = edges[3];
            move.j11 = solution.vertices[move.pos_1];
            move.j12 = solution.vertices[move.pos_3 + 1];
            move.j21 = solution.vertices[move.pos_4];
            move.j22 = solution.vertices[move.pos_2 + 1];
            move.j31 = solution.vertices[move.pos_3];
            move.j32 = solution.vertices[move.pos_1 + 1];
            move.j41 = solution.vertices[move.pos_2];
            move.j42 = solution.vertices[move.pos_4 + 1];
            assert(move.pos_1 >= 0);
            assert(move.pos_4 + 1 < (VertexPos)solution.vertices.size());
            move.global_cost = global_cost(solution);
            moves.push_back(move);
        }
        return moves;
    }

    inline void apply_move(Solution& solution, const Move& move)
    {
        //print(std::cout, solution);
        //std::cout << "pos_1 " << move.pos_1
        //    << " pos_2 " << move.pos_2
        //    << " pos_3 " << move.pos_3
        //    << " pos_4 " << move.pos_4
        //    << std::endl;
        //std::cout << "j11 " << move.j11 << " j12 " << move.j12 << std::endl;
        //std::cout << "j21 " << move.j21 << " j22 " << move.j22 << std::endl;
        //std::cout << "j31 " << move.j31 << " j32 " << move.j32 << std::endl;
        //std::cout << "j41 " << move.j41 << " j42 " << move.j42 << std::endl;
        std::vector<VertexId> vertices;
        for (VertexPos pos = 0; pos <= move.pos_1; ++pos)
            vertices.push_back(solution.vertices[pos]);
        for (VertexPos pos = move.pos_3 + 1; pos <= move.pos_4; ++pos)
            vertices.push_back(solution.vertices[pos]);
        for (VertexPos pos = move.pos_2 + 1; pos <= move.pos_3; ++pos)
            vertices.push_back(solution.vertices[pos]);
        for (VertexPos pos = move.pos_1 + 1; pos <= move.pos_2; ++pos)
            vertices.push_back(solution.vertices[pos]);
        for (VertexPos pos = move.pos_4 + 1; pos < (VertexPos)solution.vertices.size(); ++pos)
            vertices.push_back(solution.vertices[pos]);
        assert((VertexPos)vertices.size() <= instance_.number_of_vertices() + 1);
        compute(solution, vertices);
    }

    /*
     *  j11   j21              j11 - j21
     *      x          =>
     *  j22   j12              j22 - j12
     */
    struct MoveTwoOpt
    {
        VertexId j11 = -1;
        VertexId j12 = -1;
        VertexId j21 = -1;
        VertexId j22 = -1;
        GlobalCost cost_difference = global_cost_worst();
    };

    /*
     *  j11   j22   j13             j11 - j22 - j13
     *      x     x         =>
     *  j21   j12   j23             j21 - j12 - j23
     */
    struct MoveSwap
    {
        VertexId j11 = -1;
        VertexId j12 = -1;
        VertexId j13 = -1;
        VertexId j21 = -1;
        VertexId j22 = -1;
        VertexId j23 = -1;
        GlobalCost cost_difference = global_cost_worst();
    };

    /*
     *  j111                   j122    j111 ----------------- j122
     *      \                 /
     *       j112 - ... - j121      =>      j112 - ... - j121
     *                                    /                   \
     *  j21 ------------------ j22    j21                       j22
     */
    struct MoveOrOpt
    {
        VertexId j111 = -1;
        VertexId j112 = -1;
        VertexId j121 = -1;
        VertexId j122 = -1;
        VertexId j21 = -1;
        VertexId j22 = -1;
        bool reverse = false;
        GlobalCost cost_difference = global_cost_worst();
    };

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Move& perturbation = move_null())
    {
        //print(std::cout, solution);
        //std::cout << to_string(global_cost(solution)) << std::endl;

        // Compute neighborhoods.
        std::vector<Counter> neighborhoods;
        if (parameters_.twoopt)
            neighborhoods.push_back(-1);
        if (parameters_.swap)
            neighborhoods.push_back(0);
        for (VertexPos bloc_size = 1; bloc_size <= parameters_.oropt_size_max; ++bloc_size)
            neighborhoods.push_back(bloc_size);

        // Strucutres for the 2-opt neighborhood.
        // Edges added to the solution since the last neighborhood exploration.
        // This structure is updated after moves from other neighborhoods are
        // applied.
        std::vector<std::pair<VertexId, VertexId>> twoopt_new_edges;
        // Structure used to store whether or not edges of the solution are
        // new.
        std::vector<int8_t> twoopt_is_edge_new(instance_.number_of_vertices() + 1);
        // Vector containing the improving moves of the neighborhood.
        std::vector<MoveTwoOpt> twoopt_improving_moves;

        // Strucutres for the swap neighborhood.
        // Edges added to the solution since the last neighborhood exploration.
        // This structure is updated after moves from other neighborhoods are
        // applied.
        std::vector<std::pair<VertexId, VertexId>> swap_new_edges;
        // Structure to get whether a vertex has been modified or not.
        std::vector<uint8_t> swap_is_vertex_modified(instance_.number_of_vertices());
        // Vector containing the improving moves of the neighborhood.
        std::vector<MoveSwap> swap_improving_moves;

        // Strucutres for the Or-opt neighborhood.
        // Edges added to the solution since the last neighborhood exploration.
        // This structure is updated after moves from other neighborhoods are
        // applied.
        std::vector<std::vector<std::pair<VertexId, VertexId>>> oropt_new_edges(parameters_.oropt_size_max + 1);
        // Vector containing the improving moves of the neighborhood.
        std::vector<std::vector<MoveOrOpt>> oropt_improving_moves(parameters_.oropt_size_max + 1);
        // Structure used to store whether or not edges of the solution are
        // new.
        std::vector<int8_t> oropt_is_edge_new(instance_.number_of_vertices() + 1);

        // Initialize move structures.
        // If we call the local_search on a solution which has not been
        // perturbed, it is because it is not locally optimal. Therefore, all
        // edges are considered new.
        // Otherwise, only the edges modified by the perturbation are
        // considered new.
        if (perturbation.pos_1 == -1) {
            for (VertexPos pos = 0; pos < (VertexPos)solution.vertices.size() - 1; ++pos) {
                VertexId j1 = solution.vertices[pos];
                VertexId j2 = solution.vertices[pos + 1];
                twoopt_new_edges.push_back({j1, j2});
                swap_new_edges.push_back({j1, j2});
                for (VertexPos bloc_size = 1; bloc_size <= parameters_.oropt_size_max; ++bloc_size)
                    oropt_new_edges[bloc_size].push_back({j1, j2});
            }
        } else {
            twoopt_new_edges.push_back({perturbation.j11, perturbation.j12});
            twoopt_new_edges.push_back({perturbation.j21, perturbation.j22});
            twoopt_new_edges.push_back({perturbation.j31, perturbation.j32});
            twoopt_new_edges.push_back({perturbation.j41, perturbation.j42});
            swap_new_edges.push_back({perturbation.j11, perturbation.j12});
            swap_new_edges.push_back({perturbation.j21, perturbation.j22});
            swap_new_edges.push_back({perturbation.j31, perturbation.j32});
            swap_new_edges.push_back({perturbation.j41, perturbation.j42});
            for (VertexPos bloc_size = 1; bloc_size <= parameters_.oropt_size_max; ++bloc_size) {
                oropt_new_edges[bloc_size].push_back({perturbation.j11, perturbation.j12});
                oropt_new_edges[bloc_size].push_back({perturbation.j21, perturbation.j22});
                oropt_new_edges[bloc_size].push_back({perturbation.j31, perturbation.j32});
                oropt_new_edges[bloc_size].push_back({perturbation.j41, perturbation.j42});
            }
        }

        Counter it = 0;
        for (;; ++it) {
            GlobalCost c_cur = global_cost(solution);
            //std::cout << "it " << it
            //    << " c " << to_string(c_cur)
            //    << std::endl;
            //print(std::cout, solution);

            if (parameters_.shuffle_neighborhood_order)
                std::shuffle(neighborhoods.begin(), neighborhoods.end(), generator);
            bool improved = false;
            // Loop through neighborhoods.
            for (Counter neighborhood: neighborhoods) {
                switch (neighborhood) {
                case -1: { // 2-opt neighborhood.
                    //std::cout << "Explore 2-opt neighborhood." << std::endl;

                    // Remove obsolete moves.
                    for (auto it = twoopt_improving_moves.begin();
                            it != twoopt_improving_moves.end();) {
                        VertexPos p11 = solution.positions[it->j11];
                        VertexPos p21 = solution.positions[it->j21];
                        if (solution.vertices[p11 + 1] != it->j12
                                || solution.vertices[p21 + 1] != it->j22) {
                            *it = *std::prev(twoopt_improving_moves.end());
                            twoopt_improving_moves.pop_back();
                        } else {
                            ++it;
                        }
                    }

                    // Evaluate new moves.
                    // For each new edge, we add a 2-opt exchange between it
                    // and all other edges.
                    // Update twoopt_is_edge_new.
                    std::fill(twoopt_is_edge_new.begin(), twoopt_is_edge_new.end(), 0);
                    for (auto edge: twoopt_new_edges) {
                        VertexId j1 = edge.first;
                        VertexId j2 = edge.second;
                        VertexPos p1 = solution.positions[j1];
                        // Even if it had been added since the previous
                        // 2-opt neighborhood exploration, an edge might have
                        // already been removed from the solution.
                        if (solution.vertices[p1 + 1] != j2)
                            continue;
                        twoopt_is_edge_new[p1] = 1;
                    }
                    for (VertexPos pos_1 = 0; pos_1 < instance_.number_of_vertices(); ++pos_1) {
                        // If the edge is not new, continue.
                        if (twoopt_is_edge_new[pos_1] == 0)
                            continue;
                        VertexId j11 = solution.vertices[pos_1];
                        VertexId j12 = solution.vertices[pos_1 + 1];
                        for (VertexPos pos_2 = 0; pos_2 < instance_.number_of_vertices(); ++pos_2) {
                            if (pos_1 == pos_2)
                                continue;
                            VertexId j21 = solution.vertices[pos_2];
                            VertexId j22 = solution.vertices[pos_2 + 1];
                            Distance length_difference =
                                - instance_.distance(j11, j12)
                                - instance_.distance(j21, j22)
                                + instance_.distance(j11, j21)
                                + instance_.distance(j22, j12);
                            if (length_difference >= 0)
                                continue;
                            MoveTwoOpt move;
                            move.j11 = j11;
                            move.j12 = j12;
                            move.j21 = j21;
                            move.j22 = j22;
                            move.cost_difference = {0, length_difference};
                            twoopt_improving_moves.push_back(move);
                        }
                    }
                    twoopt_new_edges.clear();

                    // If there is no improving move, then stop here.
                    if (twoopt_improving_moves.empty())
                        break;

                    // Otherwise, look for a pareto improving move.
                    //std::shuffle(twoopt_improving_moves.begin(), twoopt_improving_moves.end(), generator);
                    auto it_best = twoopt_improving_moves.begin();
                    for (auto it = twoopt_improving_moves.begin() + 1;
                            it != twoopt_improving_moves.end(); ++it) {
                        if (it->cost_difference >= it_best->cost_difference
                                || !dominates(it->cost_difference, it_best->cost_difference))
                            continue;
                        it_best = it;
                    }

                    // Apply best move.
                    improved = true;
                    VertexPos p11 = solution.positions[it_best->j11];
                    VertexPos p12 = solution.positions[it_best->j12];
                    VertexPos p21 = solution.positions[it_best->j21];
                    VertexPos p22 = solution.positions[it_best->j22];
                    if (p12 == 0)
                        p12 = solution.vertices.size() - 1;
                    if (p22 == 0)
                        p22 = solution.vertices.size() - 1;
                    //std::cout << "j11 " << it_best->j11
                    //    << " j12 " << it_best->j12
                    //    << " j21 " << it_best->j21
                    //    << " j22 " << it_best->j22
                    //    << std::endl;
                    //std::cout << "p11 " << p11
                    //    << " p12 " << p12
                    //    << " p21 " << p21
                    //    << " p22 " << p22
                    //    << std::endl;
                    //std::cout << "solution:";
                    //for (VertexId j: solution.vertices)
                    //    std::cout << " " << j;
                    //std::cout << std::endl;
                    std::vector<VertexId> vertices;
                    if (p11 < p21) {
                        for (VertexPos pos = 0; pos <= p11; ++pos)
                            vertices.push_back(solution.vertices[pos]);
                        for (VertexPos pos = p21; pos >= p12; --pos)
                            vertices.push_back(solution.vertices[pos]);
                        for (VertexPos pos = p22; pos < (VertexPos)solution.vertices.size(); ++pos)
                            vertices.push_back(solution.vertices[pos]);
                    } else {
                        for (VertexPos pos = 0; pos <= p21; ++pos)
                            vertices.push_back(solution.vertices[pos]);
                        for (VertexPos pos = p11; pos >= p22; --pos)
                            vertices.push_back(solution.vertices[pos]);
                        for (VertexPos pos = p12; pos < (VertexPos)solution.vertices.size(); ++pos)
                            vertices.push_back(solution.vertices[pos]);
                    }
                    //std::cout << "solution:";
                    //for (VertexId j: vertices)
                    //    std::cout << " " << j;
                    //std::cout << std::endl;
                    compute(solution, vertices);
                    if (global_cost(solution) != c_cur + it_best->cost_difference) {
                        throw std::logic_error("2-opt. Costs do not match:\n"
                                "* Current cost: " + to_string(c_cur) + "\n"
                                + "* Move cost difference: " + to_string(it_best->cost_difference) + "\n"
                                + "* Expected new cost: " + to_string(c_cur + it_best->cost_difference) + "\n"
                                + "* Actual new cost: " + to_string(global_cost(solution)) + "\n");
                    }

                    // Update move structures.
                    twoopt_new_edges.push_back({it_best->j11, it_best->j21});
                    twoopt_new_edges.push_back({it_best->j22, it_best->j12});
                    swap_new_edges.push_back({it_best->j11, it_best->j21});
                    swap_new_edges.push_back({it_best->j22, it_best->j12});
                    for (VertexPos bloc_size = 1; bloc_size <= parameters_.oropt_size_max; ++bloc_size) {
                        oropt_new_edges[bloc_size].push_back({it_best->j11, it_best->j21});
                        oropt_new_edges[bloc_size].push_back({it_best->j22, it_best->j12});
                    }

                    //std::cout << "Improve with 2-opt." << std::endl;
                    break;
                } case 0: { // Swap neighborhood.
                    //std::cout << "Explore swap neighborhood." << std::endl;

                    // Remove obsolete moves.
                    for (auto it = swap_improving_moves.begin();
                            it != swap_improving_moves.end();) {
                        VertexPos p11 = solution.positions[it->j11];
                        VertexPos p12 = solution.positions[it->j12];
                        VertexPos p21 = solution.positions[it->j21];
                        VertexPos p22 = solution.positions[it->j22];
                        if (p11 + 2 >= (VertexPos)solution.vertices.size()
                                || p21 + 2 >= (VertexPos)solution.vertices.size()
                                || solution.vertices[p11 + 1] != it->j12
                                || solution.vertices[p12 + 1] != it->j13
                                || solution.vertices[p21 + 1] != it->j22
                                || solution.vertices[p22 + 1] != it->j23) {
                            *it = *std::prev(swap_improving_moves.end());
                            swap_improving_moves.pop_back();
                        } else {
                            ++it;
                        }
                    }

                    // Evaluate new moves.
                    // For each new edge, we add a swap move between each
                    // vertex of the edge and all other vertices.
                    std::fill(swap_is_vertex_modified.begin(), swap_is_vertex_modified.end(), 0);
                    for (auto edge: twoopt_new_edges) {
                        VertexId j1 = edge.first;
                        VertexId j2 = edge.second;
                        VertexPos p1 = solution.positions[j1];
                        assert(p1 + 1 < (VertexPos)solution.vertices.size());
                        // Even if it had been added since the previous
                        // swap neighborhood exploration, an edge might have
                        // already been removed from the solution.
                        if (solution.vertices[p1 + 1] != j2)
                            continue;
                        swap_is_vertex_modified[j1] = 1;
                        swap_is_vertex_modified[j2] = 1;
                    }
                    for (VertexId j12 = 1; j12 < (VertexPos)solution.vertices.size() - 1; ++j12) {
                        // If the vertex has not been modified, continue.
                        if (swap_is_vertex_modified[j12] == 0)
                            continue;
                        if (j12 == 0)
                            continue;
                        VertexPos pos_1 = solution.positions[j12];
                        VertexId j11 = solution.vertices[pos_1 - 1];
                        VertexId j13 = solution.vertices[pos_1 + 1];
                        for (VertexId j22 = 1; j22 < (VertexPos)solution.vertices.size() - 1; ++j22) {
                            VertexPos pos_2 = solution.positions[j22];
                            // We do not consider swaps of consecutive vertices
                            // since they are included in shifts.
                            if (pos_2 == pos_1
                                    || pos_2 == pos_1 - 1
                                    || pos_2 == pos_1 + 1)
                                continue;
                            VertexId j21 = solution.vertices[pos_2 - 1];
                            VertexId j23 = solution.vertices[pos_2 + 1];
                            Distance length_difference =
                                - instance_.distance(j11, j12)
                                - instance_.distance(j12, j13)
                                - instance_.distance(j21, j22)
                                - instance_.distance(j22, j23)
                                + instance_.distance(j11, j22)
                                + instance_.distance(j22, j13)
                                + instance_.distance(j21, j12)
                                + instance_.distance(j12, j23);
                            if (length_difference >= 0)
                                continue;
                            MoveSwap move;
                            move.j11 = j11;
                            move.j12 = j12;
                            move.j13 = j13;
                            move.j21 = j21;
                            move.j22 = j22;
                            move.j23 = j23;
                            move.cost_difference = {0, length_difference};
                            swap_improving_moves.push_back(move);
                        }
                    }
                    swap_new_edges.clear();

                    // If there is no improving move, then stop here.
                    if (swap_improving_moves.empty())
                        break;

                    // Otherwise, look for a pareto improving move.
                    //std::shuffle(swap_improving_moves.begin(), swap_improving_moves.end(), generator);
                    auto it_best = swap_improving_moves.begin();
                    for (auto it = swap_improving_moves.begin() + 1;
                            it != swap_improving_moves.end(); ++it) {
                        if (it->cost_difference >= it_best->cost_difference
                                || !dominates(it->cost_difference, it_best->cost_difference))
                            continue;
                        it_best = it;
                    }

                    // Apply best move.
                    improved = true;
                    std::vector<VertexId> vertices;
                    //VertexPos p11 = solution.positions[it_best->j11];
                    VertexPos p12 = solution.positions[it_best->j12];
                    //VertexPos p13 = solution.positions[it_best->j13];
                    //VertexPos p21 = solution.positions[it_best->j21];
                    VertexPos p22 = solution.positions[it_best->j22];
                    //VertexPos p23 = solution.positions[it_best->j23];
                    //std::cout << "j11 " << it_best->j11
                    //    << " j12 " << it_best->j12
                    //    << " j13 " << it_best->j13
                    //    << " j21 " << it_best->j21
                    //    << " j22 " << it_best->j22
                    //    << " j23 " << it_best->j23
                    //    << std::endl;
                    //std::cout << "p11 " << p11
                    //    << " p12 " << p12
                    //    << " p13 " << p13
                    //    << " p21 " << p21
                    //    << " p22 " << p22
                    //    << " p23 " << p23
                    //    << std::endl;
                    for (VertexPos pos = 0; pos < (VertexPos)solution.vertices.size(); ++pos) {
                        if (pos == p12) {
                            vertices.push_back(it_best->j22);
                        } else if (pos == p22) {
                            vertices.push_back(it_best->j12);
                        } else {
                            vertices.push_back(solution.vertices[pos]);
                        }
                    }
                    compute(solution, vertices);
                    if (global_cost(solution) != c_cur + it_best->cost_difference) {
                        throw std::logic_error("Swap. Costs do not match:\n"
                                "* Current cost: " + to_string(c_cur) + "\n"
                                + "* Move cost difference: " + to_string(it_best->cost_difference) + "\n"
                                + "* Expected new cost: " + to_string(c_cur + it_best->cost_difference) + "\n"
                                + "* Actual new cost: " + to_string(global_cost(solution)) + "\n");
                    }

                    // Update move structures.
                    twoopt_new_edges.push_back({it_best->j11, it_best->j22});
                    twoopt_new_edges.push_back({it_best->j22, it_best->j13});
                    twoopt_new_edges.push_back({it_best->j21, it_best->j12});
                    twoopt_new_edges.push_back({it_best->j12, it_best->j23});
                    swap_new_edges.push_back({it_best->j11, it_best->j22});
                    swap_new_edges.push_back({it_best->j22, it_best->j13});
                    swap_new_edges.push_back({it_best->j21, it_best->j12});
                    swap_new_edges.push_back({it_best->j12, it_best->j23});
                    for (VertexPos bloc_size = 1; bloc_size <= parameters_.oropt_size_max; ++bloc_size) {
                        oropt_new_edges[bloc_size].push_back({it_best->j11, it_best->j22});
                        oropt_new_edges[bloc_size].push_back({it_best->j22, it_best->j13});
                        oropt_new_edges[bloc_size].push_back({it_best->j21, it_best->j12});
                        oropt_new_edges[bloc_size].push_back({it_best->j12, it_best->j23});
                    }

                    //std::cout << "Improve with swap." << std::endl;
                    break;
                } default: { // Or-opt neighborhood.
                    VertexPos bloc_size = neighborhood;
                    //std::cout << "Explore or-opt neighborhood " << bloc_size << "." << std::endl;

                    // Remove obsolete moves.
                    for (auto it = oropt_improving_moves[bloc_size].begin();
                            it != oropt_improving_moves[bloc_size].end();) {
                        VertexPos p111 = solution.positions[it->j111];
                        VertexPos p121 = solution.positions[it->j121];
                        VertexPos p122 = solution.positions[it->j122];
                        if (p122 == 0)
                            p122 = solution.vertices.size() - 1;
                        VertexPos p21 = solution.positions[it->j21];
                        VertexPos p22 = solution.positions[it->j22];
                        if (p22 == 0)
                            p22 = solution.vertices.size() - 1;
                        if (solution.vertices[p111 + 1] != it->j112
                                || solution.vertices[p121 + 1] != it->j122
                                || solution.vertices[p21 + 1] != it->j22
                                || (p22 > p111 && p21 < p122)
                                || p111 >= p122) {
                            *it = *std::prev(oropt_improving_moves[bloc_size].end());
                            oropt_improving_moves[bloc_size].pop_back();
                        } else {
                            ++it;
                        }
                    }

                    // Evaluate new moves.
                    // For each new edge, we add or-opt moves with the edge +
                    // bloc_size / edge - bloc_size and all other edges, and
                    // or-opt moves with all other edge and edge + bloc_size.
                    std::fill(oropt_is_edge_new.begin(), oropt_is_edge_new.end(), 0);
                    for (auto edge: oropt_new_edges[bloc_size]) {
                        VertexId j1 = edge.first;
                        VertexId j2 = edge.second;
                        VertexPos p1 = solution.positions[j1];
                        assert(p1 + 1 < (VertexPos)solution.vertices.size());
                        // Even if it had been added since the previous
                        // 2-opt neighborhood exploration, an edge might have
                        // already been removed from the solution.
                        if (solution.vertices[p1 + 1] != j2)
                            continue;
                        oropt_is_edge_new[p1] = 1;
                    }
                    for (VertexPos pos_1 = 0; pos_1 < instance_.number_of_vertices(); ++pos_1) {
                        // If the edge is not new, continue.
                        if (oropt_is_edge_new[pos_1] == 0)
                            continue;
                        VertexId ja1 = solution.vertices[pos_1];
                        VertexId ja2 = solution.vertices[pos_1 + 1];
                        for (VertexPos pos_2 = 0; pos_2 < instance_.number_of_vertices(); ++pos_2) {
                            if (pos_1 == pos_2)
                                continue;
                            VertexId jb1 = solution.vertices[pos_2];
                            VertexId jb2 = solution.vertices[pos_2 + 1];
                            // ja1 = j111
                            // jb1 = j21
                            if (pos_1 + bloc_size + 1 < (VertexPos)solution.vertices.size()) {
                                VertexId j111 = ja1;
                                VertexId j112 = ja2;
                                VertexId j121 = solution.vertices[pos_1 + bloc_size];
                                VertexId j122 = solution.vertices[pos_1 + bloc_size + 1];
                                VertexId j21 = jb1;
                                VertexId j22 = jb2;

                                VertexPos p111 = solution.positions[j111];
                                VertexPos p122 = solution.positions[j122];
                                if (p122 == 0)
                                    p122 = solution.vertices.size() - 1;
                                VertexPos p21 = solution.positions[j21];
                                VertexPos p22 = solution.positions[j22];
                                if (p22 == 0)
                                    p22 = solution.vertices.size() - 1;
                                if (p22 <= p111 || p21 >= p122) {
                                    if (p111 >= p122)
                                        throw std::logic_error(
                                                "ja1 = j111, jb1 = j21."
                                                " In an or-opt move,"
                                                " j122 must be after j111.");
                                    Distance length_difference =
                                        - instance_.distance(j111, j112)
                                        - instance_.distance(j121, j122)
                                        - instance_.distance(j21, j22)
                                        + instance_.distance(j111, j122)
                                        + instance_.distance(j21, j112)
                                        + instance_.distance(j121, j22);
                                    if (length_difference < 0) {
                                        MoveOrOpt move;
                                        move.j111 = j111;
                                        move.j112 = j112;
                                        move.j121 = j121;
                                        move.j122 = j122;
                                        move.j21 = j21;
                                        move.j22 = j22;
                                        move.cost_difference = {0, length_difference};
                                        oropt_improving_moves[bloc_size].push_back(move);
                                    }
                                }
                            }
                            // ja1 = j111
                            // jb1 = j21
                            // reverse
                            if (pos_1 + bloc_size + 1 < (VertexPos)solution.vertices.size()) {
                                VertexId j111 = ja1;
                                VertexId j112 = ja2;
                                VertexId j121 = solution.vertices[pos_1 + bloc_size];
                                VertexId j122 = solution.vertices[pos_1 + bloc_size + 1];
                                VertexId j21 = jb1;
                                VertexId j22 = jb2;

                                VertexPos p111 = solution.positions[j111];
                                VertexPos p122 = solution.positions[j122];
                                if (p122 == 0)
                                    p122 = solution.vertices.size() - 1;
                                VertexPos p21 = solution.positions[j21];
                                VertexPos p22 = solution.positions[j22];
                                if (p22 == 0)
                                    p22 = solution.vertices.size() - 1;
                                if (p22 <= p111 || p21 >= p122) {
                                    if (p111 >= p122)
                                        throw std::logic_error(
                                                "ja1 = j111, jb1 = j21, reverse."
                                                " In an or-opt move,"
                                                " j122 must be after j111.");
                                    Distance length_difference =
                                        - instance_.distance(j111, j112)
                                        - instance_.distance(j121, j122)
                                        - instance_.distance(j21, j22)
                                        + instance_.distance(j111, j122)
                                        + instance_.distance(j21, j121)
                                        + instance_.distance(j112, j22);
                                    if (length_difference < 0) {
                                        MoveOrOpt move;
                                        move.j111 = j111;
                                        move.j112 = j112;
                                        move.j121 = j121;
                                        move.j122 = j122;
                                        move.j21 = j21;
                                        move.j22 = j22;
                                        move.reverse = true;
                                        move.cost_difference = {0, length_difference};
                                        oropt_improving_moves[bloc_size].push_back(move);
                                    }
                                }
                            }
                            // ja1 = j121
                            // jb1 = j21
                            if (pos_1 - bloc_size >= 0) {
                                VertexId j111 = solution.vertices[pos_1 - bloc_size];
                                VertexId j112 = solution.vertices[pos_1 - bloc_size + 1];
                                VertexId j121 = ja1;
                                VertexId j122 = ja2;
                                VertexId j21 = jb1;
                                VertexId j22 = jb2;

                                VertexPos p111 = solution.positions[j111];
                                VertexPos p122 = solution.positions[j122];
                                if (p122 == 0)
                                    p122 = solution.vertices.size() - 1;
                                VertexPos p21 = solution.positions[j21];
                                VertexPos p22 = solution.positions[j22];
                                if (p22 == 0)
                                    p22 = solution.vertices.size() - 1;
                                if (p22 <= p111 || p21 >= p122) {
                                    if (j122 != 0 && p111 >= p122)
                                        throw std::logic_error(
                                                "ja1 = j121, jb1 = j21."
                                                " In an or-opt move,"
                                                " j122 must be after j111.");
                                    if (p111 >= p122)
                                        throw std::logic_error(
                                                "ja1 = j121, jb1 = j21."
                                                " In an or-opt move with bloc size 1,"
                                                " j112 must be equal to j121.");
                                    Distance length_difference =
                                        - instance_.distance(j111, j112)
                                        - instance_.distance(j121, j122)
                                        - instance_.distance(j21, j22)
                                        + instance_.distance(j111, j122)
                                        + instance_.distance(j21, j112)
                                        + instance_.distance(j121, j22);
                                    if (length_difference < 0) {
                                        MoveOrOpt move;
                                        move.j111 = j111;
                                        move.j112 = j112;
                                        move.j121 = j121;
                                        move.j122 = j122;
                                        move.j21 = j21;
                                        move.j22 = j22;
                                        move.cost_difference = {0, length_difference};
                                        oropt_improving_moves[bloc_size].push_back(move);
                                    }
                                }
                            }
                            // ja1 = j121
                            // jb1 = j21
                            // reverse
                            if (pos_1 - bloc_size >= 0) {
                                VertexId j111 = solution.vertices[pos_1 - bloc_size];
                                VertexId j112 = solution.vertices[pos_1 - bloc_size + 1];
                                VertexId j121 = ja1;
                                VertexId j122 = ja2;
                                VertexId j21 = jb1;
                                VertexId j22 = jb2;

                                VertexPos p111 = solution.positions[j111];
                                VertexPos p122 = solution.positions[j122];
                                if (p122 == 0)
                                    p122 = solution.vertices.size() - 1;
                                VertexPos p21 = solution.positions[j21];
                                VertexPos p22 = solution.positions[j22];
                                if (p22 == 0)
                                    p22 = solution.vertices.size() - 1;
                                if (p22 <= p111 || p21 >= p122) {
                                    if (j122 != 0 && p111 >= p122)
                                        throw std::logic_error(
                                                "ja1 = j121, jb1 = j21, reverse."
                                                " In an or-opt move,"
                                                " j122 must be after j111.");
                                    if (p111 >= p122)
                                        throw std::logic_error(
                                                "ja1 = j121, jb1 = j21, reverse."
                                                " In an or-opt move with bloc size 1,"
                                                " j112 must be equal to j121.");
                                    Distance length_difference =
                                        - instance_.distance(j111, j112)
                                        - instance_.distance(j121, j122)
                                        - instance_.distance(j21, j22)
                                        + instance_.distance(j111, j122)
                                        + instance_.distance(j21, j121)
                                        + instance_.distance(j112, j22);
                                    if (length_difference < 0) {
                                        MoveOrOpt move;
                                        move.j111 = j111;
                                        move.j112 = j112;
                                        move.j121 = j121;
                                        move.j122 = j122;
                                        move.j21 = j21;
                                        move.j22 = j22;
                                        move.reverse = true;
                                        move.cost_difference = {0, length_difference};
                                        oropt_improving_moves[bloc_size].push_back(move);
                                    }
                                }
                            }
                            // ja1 = j21
                            // jb1 = j111
                            if (pos_2 + bloc_size + 1 < (VertexPos)solution.vertices.size()) {
                                VertexId j111 = jb1;
                                VertexId j112 = jb2;
                                VertexId j121 = solution.vertices[pos_2 + bloc_size];
                                VertexId j122 = solution.vertices[pos_2 + bloc_size + 1];
                                VertexId j21 = ja1;
                                VertexId j22 = ja2;

                                VertexPos p111 = solution.positions[j111];
                                VertexPos p122 = solution.positions[j122];
                                if (p122 == 0)
                                    p122 = solution.vertices.size() - 1;
                                VertexPos p21 = solution.positions[j21];
                                VertexPos p22 = solution.positions[j22];
                                if (p22 == 0)
                                    p22 = solution.vertices.size() - 1;
                                if (p22 <= p111 || p21 >= p122) {
                                    if (p111 >= p122)
                                        throw std::logic_error(
                                                "ja1 = j21, jb1 = j111."
                                                " In an or-opt move,"
                                                " j122 must be after j111.");
                                    Distance length_difference =
                                        - instance_.distance(j111, j112)
                                        - instance_.distance(j121, j122)
                                        - instance_.distance(j21, j22)
                                        + instance_.distance(j111, j122)
                                        + instance_.distance(j21, j112)
                                        + instance_.distance(j121, j22);
                                    if (length_difference < 0) {
                                        MoveOrOpt move;
                                        move.j111 = j111;
                                        move.j112 = j112;
                                        move.j121 = j121;
                                        move.j122 = j122;
                                        move.j21 = j21;
                                        move.j22 = j22;
                                        move.cost_difference = {0, length_difference};
                                        oropt_improving_moves[bloc_size].push_back(move);
                                    }
                                }
                            }
                            // ja1 = j21
                            // jb1 = j111
                            // reverse
                            if (pos_2 + bloc_size + 1 < (VertexPos)solution.vertices.size()) {
                                VertexId j111 = jb1;
                                VertexId j112 = jb2;
                                VertexId j121 = solution.vertices[pos_2 + bloc_size];
                                VertexId j122 = solution.vertices[pos_2 + bloc_size + 1];
                                VertexId j21 = ja1;
                                VertexId j22 = ja2;

                                VertexPos p111 = solution.positions[j111];
                                VertexPos p122 = solution.positions[j122];
                                if (p122 == 0)
                                    p122 = solution.vertices.size() - 1;
                                VertexPos p21 = solution.positions[j21];
                                VertexPos p22 = solution.positions[j22];
                                if (p22 == 0)
                                    p22 = solution.vertices.size() - 1;
                                if (p22 <= p111 || p21 >= p122) {
                                    if (p111 >= p122)
                                        throw std::logic_error(
                                                "ja1 = j21, jb1 = j111, reverse."
                                                " In an or-opt move,"
                                                " j122 must be after j111.");
                                    Distance length_difference =
                                        - instance_.distance(j111, j112)
                                        - instance_.distance(j121, j122)
                                        - instance_.distance(j21, j22)
                                        + instance_.distance(j111, j122)
                                        + instance_.distance(j21, j121)
                                        + instance_.distance(j112, j22);
                                    if (length_difference < 0) {
                                        MoveOrOpt move;
                                        move.j111 = j111;
                                        move.j112 = j112;
                                        move.j121 = j121;
                                        move.j122 = j122;
                                        move.j21 = j21;
                                        move.j22 = j22;
                                        move.reverse = true;
                                        move.cost_difference = {0, length_difference};
                                        oropt_improving_moves[bloc_size].push_back(move);
                                    }
                                }
                            }
                        }
                    }
                    oropt_new_edges[bloc_size].clear();

                    // If there is no improving move, then stop here.
                    if (oropt_improving_moves[bloc_size].empty())
                        break;

                    // Otherwise, look for a pareto improving move.
                    //std::shuffle(
                    //        oropt_improving_moves[bloc_size].begin(),
                    //        oropt_improving_moves[bloc_size].end(),
                    //        generator);
                    auto it_best = oropt_improving_moves[bloc_size].begin();
                    for (auto it = oropt_improving_moves[bloc_size].begin() + 1;
                            it != oropt_improving_moves[bloc_size].end(); ++it) {
                        if (it->cost_difference >= it_best->cost_difference
                                || !dominates(it->cost_difference, it_best->cost_difference))
                            continue;
                        it_best = it;
                    }

                    // Apply best move.
                    improved = true;
                    std::vector<VertexId> vertices;
                    VertexPos p111 = solution.positions[it_best->j111];
                    VertexPos p112 = solution.positions[it_best->j112];
                    VertexPos p121 = solution.positions[it_best->j121];
                    VertexPos p122 = solution.positions[it_best->j122];
                    if (p122 == 0)
                        p122 = solution.vertices.size() - 1;
                    VertexPos p21 = solution.positions[it_best->j21];
                    VertexPos p22 = solution.positions[it_best->j22];
                    if (p22 == 0)
                        p22 = solution.vertices.size() - 1;
                    //std::cout << "j111 " << it_best->j111
                    //    << " j112 " << it_best->j112
                    //    << " j121 " << it_best->j121
                    //    << " j122 " << it_best->j122
                    //    << " j21 " << it_best->j21
                    //    << " j22 " << it_best->j22
                    //    << std::endl;
                    //std::cout << "p111 " << p111
                    //    << " p112 " << p112
                    //    << " p121 " << p121
                    //    << " p122 " << p122
                    //    << " p21 " << p21
                    //    << " p22 " << p22
                    //    << std::endl;
                    //std::cout << "solution:";
                    //for (VertexId j: solution.vertices)
                    //    std::cout << " " << j;
                    //std::cout << std::endl;
                    if (p121 > p21) {
                        for (VertexPos p = 0; p <= p21; ++p)
                            vertices.push_back(solution.vertices[p]);
                        if (!it_best->reverse) {
                            for (VertexPos p = p112; p <= p121; ++p)
                                vertices.push_back(solution.vertices[p]);
                        } else {
                            for (VertexPos p = p121; p >= p112; --p)
                                vertices.push_back(solution.vertices[p]);
                        }
                        for (VertexPos p = p22; p <= p111; ++p)
                            vertices.push_back(solution.vertices[p]);
                        for (VertexPos p = p122; p < (VertexPos)solution.vertices.size(); ++p)
                            vertices.push_back(solution.vertices[p]);
                    } else {
                        for (VertexPos p = 0; p <= p111; ++p)
                            vertices.push_back(solution.vertices[p]);
                        for (VertexPos p = p122; p <= p21; ++p)
                            vertices.push_back(solution.vertices[p]);
                        if (!it_best->reverse) {
                            for (VertexPos p = p112; p <= p121; ++p)
                                vertices.push_back(solution.vertices[p]);
                        } else {
                            for (VertexPos p = p121; p >= p112; --p)
                                vertices.push_back(solution.vertices[p]);
                        }
                        for (VertexPos p = p22; p < (VertexPos)solution.vertices.size(); ++p)
                            vertices.push_back(solution.vertices[p]);
                    }
                    //std::cout << "solution:";
                    //for (VertexId j: vertices)
                    //    std::cout << " " << j;
                    //std::cout << std::endl;
                    compute(solution, vertices);
                    if (global_cost(solution) != c_cur + it_best->cost_difference) {
                        throw std::logic_error("Or-opt " + std::to_string(bloc_size) + ". Costs do not match:\n"
                                "* Current cost: " + to_string(c_cur) + "\n"
                                + "* Move cost difference: " + to_string(it_best->cost_difference) + "\n"
                                + "* Expected new cost: " + to_string(c_cur + it_best->cost_difference) + "\n"
                                + "* Actual new cost: " + to_string(global_cost(solution)) + "\n");
                    }

                    // Update move structures.
                    if (!it_best->reverse) {
                        twoopt_new_edges.push_back({it_best->j111, it_best->j122});
                        twoopt_new_edges.push_back({it_best->j21, it_best->j112});
                        twoopt_new_edges.push_back({it_best->j121, it_best->j22});
                        swap_new_edges.push_back({it_best->j111, it_best->j122});
                        swap_new_edges.push_back({it_best->j21, it_best->j112});
                        swap_new_edges.push_back({it_best->j121, it_best->j22});
                        for (VertexPos bloc_size = 1; bloc_size <= parameters_.oropt_size_max; ++bloc_size) {
                            oropt_new_edges[bloc_size].push_back({it_best->j111, it_best->j122});
                            oropt_new_edges[bloc_size].push_back({it_best->j21, it_best->j112});
                            oropt_new_edges[bloc_size].push_back({it_best->j121, it_best->j22});
                        }
                    } else {
                        twoopt_new_edges.push_back({it_best->j111, it_best->j122});
                        twoopt_new_edges.push_back({it_best->j21, it_best->j121});
                        twoopt_new_edges.push_back({it_best->j112, it_best->j22});
                        swap_new_edges.push_back({it_best->j111, it_best->j122});
                        swap_new_edges.push_back({it_best->j21, it_best->j121});
                        swap_new_edges.push_back({it_best->j112, it_best->j22});
                        for (VertexPos bloc_size = 1; bloc_size <= parameters_.oropt_size_max; ++bloc_size) {
                            oropt_new_edges[bloc_size].push_back({it_best->j111, it_best->j122});
                            oropt_new_edges[bloc_size].push_back({it_best->j21, it_best->j121});
                            oropt_new_edges[bloc_size].push_back({it_best->j112, it_best->j22});
                        }
                    }

                    //std::cout << "Improve with or-opt " << bloc_size << "." << std::endl;
                    break;
                }
                }
                if (improved)
                    break;
            }
            if (!improved)
                break;
        }
        //print(std::cout, solution);
    }

    /*
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution)
    {
        os << "vertices:";
        for (VertexId j: solution.vertices)
            os << " " << j;
        os << std::endl;
        os << "length: " << solution.length << std::endl;
        return os;
    }

    inline void write(
            const Solution& solution,
            std::string filepath) const
    {
        if (filepath.empty())
            return;
        std::ofstream cert(filepath);
        if (!cert.good()) {
            std::cerr << "\033[31m" << "ERROR, unable to open file \"" << filepath << "\"" << "\033[0m" << std::endl;
            return;
        }

        for (VertexId j: solution.vertices)
            cert << j << " ";
    }

private:

    /*
     * Manipulate solutions.
     */

    inline void compute(
            Solution& solution,
            const std::vector<VertexId>& vertices)
    {
        if (vertices.front() != 0)
            throw std::logic_error("A tour must start with vertex 0.");
        if (vertices.back() != 0)
            throw std::logic_error("A tour must end with vertex 0.");
        solution.vertices = vertices;
        solution.length = 0;
        std::fill(solution.positions.begin(), solution.positions.end(), -1);
        for (VertexPos pos = 0; pos < (VertexPos)vertices.size() - 1; ++pos) {
            solution.length += instance_.distance(vertices[pos], vertices[pos + 1]);
            if (solution.positions[vertices[pos]] != -1)
                throw std::logic_error(
                        "Solution contains vertex "
                        + std::to_string(vertices[pos])
                        + " more than once.");
            solution.positions[vertices[pos]] = pos;
        }
    }

    /*
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

    std::vector<VertexPos> positions1_;
    std::vector<VertexPos> positions2_;

};

}

}

