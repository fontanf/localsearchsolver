#pragma once

#include "localsearchsolver/common.hpp"

#include "external/treesearchsolver/examples/travellingsalesman.hpp"

namespace localsearchsolver
{

namespace travellingsalesman
{

using VertexId = treesearchsolver::travellingsalesman::VertexId;
using VertexPos = treesearchsolver::travellingsalesman::VertexPos;
using Distance = treesearchsolver::travellingsalesman::Distance;
using Instance = treesearchsolver::travellingsalesman::Instance;

class LocalScheme
{

public:

    /** Global cost: <Vertex number, Length>; */
    using GlobalCost = std::tuple<VertexId, Distance>;

    inline VertexId&       vertex_number(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Distance&              length(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline VertexId  vertex_number(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Distance         length(const GlobalCost& global_cost) { return std::get<1>(global_cost); }

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
        VertexPos oropt_size_max = 3;
        bool swap = true;
        bool twoopt = true;
        bool shuffle_neighborhood_order = true;
        Counter perturbation_number = 10;
    };

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters),
        positions1_(instance.vertex_number() - 1),
        positions2_(instance.vertex_number() - 1)
    {
        VertexId n = instance_.vertex_number();
        std::iota(positions1_.begin(), positions1_.end(), 0);
        std::iota(positions2_.begin(), positions2_.end(), 0);
        // We do not consider swaps of consecutive vertices since they are
        // included in shifts.
        for (VertexPos pos_1 = 0; pos_1 < n - 1; ++pos_1)
            for (VertexPos pos_2 = pos_1 + 2; pos_2 < n - 1; ++pos_2)
                pairs_swap_.push_back({pos_1, pos_2});
        for (VertexPos pos_1 = 0; pos_1 < n; ++pos_1)
            for (VertexPos pos_2 = pos_1 + 1; pos_2 < instance_.vertex_number(); ++pos_2)
                pairs_twoopt_.push_back({pos_1, pos_2});
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
        return solution;
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator)
    {
        std::vector<VertexId> vertices(instance_.vertex_number() - 1);
        std::iota(vertices.begin(), vertices.end(), 1);
        std::shuffle(vertices.begin(), vertices.end(), generator);
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
        GlobalCost global_cost;
    };

    static Move move_null() { return {-1, -1, -1, -1, global_cost_worst()}; }

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
        for (Counter perturbation = 0; perturbation < parameters_.perturbation_number; ++perturbation) {
            std::vector<VertexPos> edges = optimizationtools::bob_floyd<VertexPos>(
                    4, solution.vertices.size() + 1, generator);
            std::sort(edges.begin(), edges.end());
            Move move;
            move.pos_1 = edges[0];
            move.pos_2 = edges[1];
            move.pos_3 = edges[2];
            move.pos_4 = edges[3];
            assert(move.pos_1 >= 0);
            assert(move.pos_4 <= (VertexPos)solution.vertices.size());
            move.global_cost = global_cost(solution);
            moves.push_back(move);
        }
        return moves;
    }

    inline void apply_move(Solution& solution, const Move& move)
    {
        std::vector<VertexId> vertices;
        for (VertexPos pos = 0; pos < move.pos_1; ++pos)
            vertices.push_back(solution.vertices[pos]);
        for (VertexPos pos = move.pos_3; pos < move.pos_4; ++pos)
            vertices.push_back(solution.vertices[pos]);
        for (VertexPos pos = move.pos_2; pos < move.pos_3; ++pos)
            vertices.push_back(solution.vertices[pos]);
        for (VertexPos pos = move.pos_1; pos < move.pos_2; ++pos)
            vertices.push_back(solution.vertices[pos]);
        for (VertexPos pos = move.pos_4; pos < (VertexPos)solution.vertices.size(); ++pos)
            vertices.push_back(solution.vertices[pos]);
        assert((VertexPos)vertices.size() <= instance_.vertex_number());
        compute(solution, vertices);
    }

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Move& = move_null())
    {
        //print(std::cout, solution);
        //std::cout << to_string(global_cost(solution)) << std::endl;

        VertexId n = instance_.vertex_number();
        Counter it = 0;
        std::vector<Counter> neighborhoods;
        if (parameters_.twoopt)
            neighborhoods.push_back(-1);
        if (parameters_.swap)
            neighborhoods.push_back(0);
        for (VertexPos bloc_size = 1; bloc_size <= parameters_.oropt_size_max; ++bloc_size)
            neighborhoods.push_back(bloc_size);
        for (;; ++it) {
            //std::cout << "it " << it
            //    << " c " << to_string(global_cost(solution))
            //    << std::endl;
            //print(std::cout, solution);

            if (parameters_.shuffle_neighborhood_order)
                std::shuffle(neighborhoods.begin(), neighborhoods.end(), generator);
            bool improved = false;
            // Loop through neighborhoods.
            for (Counter neighborhood: neighborhoods) {
                switch (neighborhood) {
                case -1: { // 2-opt neighborhood.
                    std::shuffle(pairs_twoopt_.begin(), pairs_twoopt_.end(), generator);
                    VertexPos pos_1_best = -1;
                    VertexPos pos_2_best = -1;
                    Distance length_best = solution.length;
                    for (auto pair: pairs_twoopt_) {
                        VertexPos pos_1 = pair.first;
                        VertexPos pos_2 = pair.second;
                        VertexId j1_prev = (pos_1 == 0)? 0: solution.vertices[pos_1 - 1];
                        VertexId j1_next = (pos_1 == n - 1)? 0: solution.vertices[pos_1];
                        VertexId j2_prev = (pos_2 == 0)? 0: solution.vertices[pos_2 - 1];
                        VertexId j2_next = (pos_2 == n - 1)? 0: solution.vertices[pos_2];
                        Distance length_new = solution.length
                            - instance_.distance(j1_prev, j1_next)
                            - instance_.distance(j2_prev, j2_next)
                            + instance_.distance(j1_prev, j2_prev)
                            + instance_.distance(j1_next, j2_next);
                        if (length_new >= length_best)
                            continue;
                        pos_1_best = pos_1;
                        pos_2_best = pos_2;
                        length_best = length_new;
                    }
                    if (pos_1_best != -1) {
                        improved = true;
                        // Apply best move.
                        std::vector<VertexId> vertices;
                        for (VertexPos pos = 0; pos < pos_1_best; ++pos)
                            vertices.push_back(solution.vertices[pos]);
                        for (VertexPos pos = pos_2_best - 1; pos >= pos_1_best; --pos)
                            vertices.push_back(solution.vertices[pos]);
                        for (VertexPos pos = pos_2_best; pos < n - 1; ++pos)
                            vertices.push_back(solution.vertices[pos]);
                        compute(solution, vertices);
                        assert(solution.length == length_best);
                    }
                    break;
                } case 0: { // Swap neighborhood.
                    std::shuffle(pairs_swap_.begin(), pairs_swap_.end(), generator);
                    VertexPos pos_1_best = -1;
                    VertexPos pos_2_best = -1;
                    Distance length_best = solution.length;
                    for (auto pair: pairs_swap_) {
                        VertexPos pos_1 = pair.first;
                        VertexPos pos_2 = pair.second;
                        VertexId j1_prev = (pos_1 == 0)? 0: solution.vertices[pos_1 - 1];
                        VertexId j1 = solution.vertices[pos_1];
                        VertexId j1_next = (pos_1 + 1 == n - 1)? 0: solution.vertices[pos_1 + 1];
                        VertexId j2_prev = (pos_2 == 0)? 0: solution.vertices[pos_2 - 1];
                        VertexId j2 = solution.vertices[pos_2];
                        VertexId j2_next = (pos_2 + 1 == n - 1)? 0: solution.vertices[pos_2 + 1];
                        Distance length_new = solution.length
                            - instance_.distance(j1_prev, j1)
                            - instance_.distance(j1, j1_next)
                            - instance_.distance(j2_prev, j2)
                            - instance_.distance(j2, j2_next)
                            + instance_.distance(j1_prev, j2)
                            + instance_.distance(j2, j1_next)
                            + instance_.distance(j2_prev, j1)
                            + instance_.distance(j1, j2_next);
                        if (length_new >= length_best)
                            continue;
                        pos_1_best = pos_1;
                        pos_2_best = pos_2;
                        length_best = length_new;
                    }
                    if (pos_1_best != -1) {
                        improved = true;
                        // Apply best move.
                        std::vector<VertexId> vertices;
                        for (VertexPos pos = 0; pos < n - 1; ++pos) {
                            if (pos == pos_1_best) {
                                vertices.push_back(solution.vertices[pos_2_best]);
                            } else if (pos == pos_2_best) {
                                vertices.push_back(solution.vertices[pos_1_best]);
                            } else {
                                vertices.push_back(solution.vertices[pos]);
                            }
                        }
                        compute(solution, vertices);
                        if (solution.length != length_best) {
                            std::cout << "swap"
                                << " pos_1_best " << pos_1_best
                                << " pos_2_best " << pos_2_best
                                << " length_best " << length_best
                                << " solution.length " << solution.length
                                << std::endl;
                            print(std::cout, solution);
                        }
                        assert(solution.length == length_best);
                    }
                    break;
                } default: { // Or-opt neighborhood.
                    VertexPos bloc_size = neighborhood;
                    std::shuffle(positions1_.begin(), positions1_.end(), generator);
                    std::shuffle(positions2_.begin(), positions2_.end(), generator);
                    VertexPos pos_best = -1;
                    VertexPos pos_new_best = -1;
                    Distance length_best = solution.length;
                    for (VertexPos pos: positions1_) {
                        if (pos > (VertexPos)solution.vertices.size() - bloc_size)
                            continue;
                        for (VertexPos pos_new: positions2_) {
                            if (pos == pos_new || pos_new > (VertexPos)solution.vertices.size() - 1 - bloc_size)
                                continue;
                            VertexId j_start_prev = (pos == 0)? 0: solution.vertices[pos - 1];
                            VertexId j_start = solution.vertices[pos];
                            VertexId j_end = solution.vertices[pos + bloc_size - 1];
                            VertexId j_end_next = (pos + bloc_size == n - 1)? 0: solution.vertices[pos + bloc_size];
                            VertexId j_new_prev = -1;
                            if (pos_new == 0) {
                                j_new_prev = 0;
                            } else if (pos_new < pos) {
                                j_new_prev = solution.vertices[pos_new - 1];
                            } else {
                                j_new_prev = solution.vertices[pos_new - 1 + bloc_size];
                            }
                            VertexId j_new_next = -1;
                            if (pos_new == n - 1) {
                                j_new_next = 0;
                            } else if (pos_new < pos) {
                                j_new_next = solution.vertices[pos_new];
                            } else {
                                j_new_next = solution.vertices[pos_new + bloc_size];
                            }
                            Distance length_new = solution.length
                                - instance_.distance(j_start_prev, j_start)
                                - instance_.distance(j_end, j_end_next)
                                - instance_.distance(j_new_prev, j_new_next)
                                + instance_.distance(j_start_prev, j_end_next)
                                + instance_.distance(j_new_prev, j_start)
                                + instance_.distance(j_end, j_new_next);
                            if (length_new >= length_best)
                                continue;
                            pos_best = pos;
                            pos_new_best = pos_new;
                            length_best = length_new;
                        }
                    }
                    if (pos_best != -1) {
                        improved = true;
                        // Apply best move.
                        std::vector<VertexId> vertices;
                        if (pos_best > pos_new_best) {
                            for (VertexPos p = 0; p < pos_new_best; ++p)
                                vertices.push_back(solution.vertices[p]);
                            for (VertexPos p = pos_best; p < pos_best + bloc_size; ++p)
                                vertices.push_back(solution.vertices[p]);
                            for (VertexPos p = pos_new_best; p < pos_best; ++p)
                                vertices.push_back(solution.vertices[p]);
                            for (VertexPos p = pos_best + bloc_size; p < (VertexPos)solution.vertices.size(); ++p)
                                vertices.push_back(solution.vertices[p]);
                        } else {
                            for (VertexPos p = 0; p < pos_best; ++p)
                                vertices.push_back(solution.vertices[p]);
                            for (VertexPos p = pos_best + bloc_size; p < pos_new_best + bloc_size; ++p)
                                vertices.push_back(solution.vertices[p]);
                            for (VertexPos p = pos_best; p < pos_best + bloc_size; ++p)
                                vertices.push_back(solution.vertices[p]);
                            for (VertexPos p = pos_new_best + bloc_size; p < (VertexPos)solution.vertices.size(); ++p)
                                vertices.push_back(solution.vertices[p]);
                        }
                        //VertexId j_new_prev = -1;
                        //if (pos_new_best == 0) {
                        //    j_new_prev = 0;
                        //} else if (pos_new_best < pos_best) {
                        //    j_new_prev = solution.vertices[pos_new_best - 1];
                        //} else {
                        //    j_new_prev = solution.vertices[pos_new_best - 1 + bloc_size];
                        //}
                        //VertexId j_new_next = -1;
                        //if (pos_new_best == n - 1) {
                        //    j_new_next = 0;
                        //} else if (pos_new_best < pos_best) {
                        //    j_new_next = solution.vertices[pos_new_best];
                        //} else {
                        //    j_new_next = solution.vertices[pos_new_best + bloc_size];
                        //}
                        //std::cout
                        //        << " j_start_prev " << ((pos_best == 0)? 0: solution.vertices[pos_best - 1])
                        //        << " j_start " << solution.vertices[pos_best]
                        //        << " j_end " << solution.vertices[pos_best + bloc_size - 1]
                        //        << " j_end_next " << ((pos_best + bloc_size == n - 1)? 0: solution.vertices[pos_best + bloc_size])
                        //        << " j_new_prev " << j_new_prev
                        //        << " j_new_next " << j_new_next
                        //        << std::endl;
                        compute(solution, vertices);
                        if (solution.length != length_best) {
                            std::cout << "or-opt"
                                << " bloc_size " << bloc_size
                                << " pos_best " << pos_best
                                << " pos_new_best " << pos_new_best
                                << " length_best " << length_best
                                << " solution.length " << solution.length
                                << std::endl;
                            print(std::cout, solution);
                        }
                        assert(solution.length == length_best);
                    }
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
        solution.vertices = vertices;
        solution.length = 0;
        VertexId j_prev = 0;
        for (VertexId j: solution.vertices) {
            solution.length += instance_.distance(j_prev, j);
            j_prev = j;
        }
        solution.length += instance_.distance(j_prev, 0);
    }

    /*
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

    std::vector<VertexPos> positions1_;
    std::vector<VertexPos> positions2_;
    std::vector<std::pair<VertexPos, VertexPos>> pairs_swap_;
    std::vector<std::pair<VertexPos, VertexPos>> pairs_twoopt_;

};

}

}

