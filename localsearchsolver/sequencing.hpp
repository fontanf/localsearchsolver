/**
 * Sequencing problems.
 */

#pragma once

#include "localsearchsolver/common.hpp"

#include "optimizationtools/containers/indexed_set.hpp"

namespace localsearchsolver
{

namespace sequencing
{

using SequenceId = int64_t;
using SequencePos = int64_t;
using ElementId = int64_t;
using ElementPos = int64_t;
using Mode = int64_t;

enum class Perturbations
{
    None,
    DoubleBridge,
    AdjacentSwaps,
    RuinAndRecreate,
    ForceAdd,
};

enum class Neighborhoods
{
    Shift,
    Swap,
    Reverse,
    ShiftReverse,

    Add,
    Remove,
    Replace,

    SwapTails,
    Split,
    InterShift,
    InterSwap,
    InterShiftReverse,
    InterSwapStar,

    ShiftChangeMode,
    ModeSwap,
    SwapWithModes,
};

std::string neighborhood2string(
        Neighborhoods neighborhood,
        ElementPos k1,
        ElementPos k2)
{
    switch (neighborhood) {
    case Neighborhoods::Shift:
        return std::to_string(k1) + "-shift";
    case Neighborhoods::Swap:
        return "(" + std::to_string(k1) + "," + std::to_string(k2) + ")" + "-swap";
    case Neighborhoods::Reverse:
        return "Reverse";
    case Neighborhoods::ShiftReverse:
        return std::to_string(k1) + "-shift-reverse";
    case Neighborhoods::Add:
        return "Add";
    case Neighborhoods::Remove:
        return "Remove";
    case Neighborhoods::Replace:
        return "Replace";
    case Neighborhoods::SwapTails:
        return "Inter-two-opt";
    case Neighborhoods::Split:
        return "Inter-two-opt-reverse";
    case Neighborhoods::InterShift:
        return std::to_string(k1) + "-inter-shift";
    case Neighborhoods::InterSwap:
        return "(" + std::to_string(k1) + "," + std::to_string(k2) + ")" + "-inter-swap";
    case Neighborhoods::InterShiftReverse:
        return std::to_string(k1) + "-inter-shift-reverse";
    case Neighborhoods::InterSwapStar:
        return "Inter-swap-star";
    case Neighborhoods::ShiftChangeMode:
        return "Shift-change-mode";
    case Neighborhoods::ModeSwap:
        return "Mode-swap";
    case Neighborhoods::SwapWithModes:
        return "Swap-with-modes";
    }
    return "";
}

struct Parameters
{

    bool shuffle_neighborhood_order = true;

    bool linking_constraints = false;

    /*
     * Neighborhoods - Intra.
     */

    ElementPos shift_block_maximum_length = 3;
    ElementPos shift_maximum_distance = 1024;
    ElementPos swap_maximum_distance = 1024;
    ElementPos swap_block_maximum_length = 2;
    bool reverse = false;
    ElementPos reverse_maximum_length = 1024;
    ElementPos shift_reverse_block_maximum_length = 0;

    /*
     * Neighborhoods - Inter.
     */

    bool swap_tails = false;
    bool split = false;
    ElementPos inter_shift_block_maximum_length = 0;
    ElementPos inter_swap_block_maximum_length = 0;
    ElementPos inter_shift_reverse_block_maximum_length = 0;
    bool inter_swap_star = false;

    /*
     * Neighborhoods - Sub-sequence.
     */

    bool add_remove = false;
    bool replace = false;

    /*
     * Neighborhoods - Modes.
     */

    bool shift_change_mode = false;
    bool mode_swap = false;
    bool swap_with_modes = false;

    /*
     * Perturbations.
     */

    /**
     * Perturbation "double-bridge".
     *
     * Thsi corresponds to a swap of two blocs of elements.
     *
     * This is a classical perturbation introduced for the Traveling Salesman
     * Problem and used in several variants.
     *
     * See:
     * - "An Iterated Local Search heuristic for the single machine total
     *   weighted tardiness scheduling problem with sequence-dependent setup
     *   times" (Subramanian et al. 2014)
     *   https://doi.org/10.1080/00207543.2014.883472
     */
    Counter double_bridge_number_of_perturbations = 0;

    /**
     * Perturbation "random adjacent swaps".
     *
     * See:
     * - "Iterated-greedy-based algorithms with beam search initialization for
     *   the permutation flowshop to minimise total tardiness"
     *   (Fernandez-Viagas et al, 2018)
     *   https://doi.org/10.1016/j.eswa.2017.10.050
     */
    Counter adjacent_swaps_number_of_perturbations = 0;

    /**
     * Perturbation "random adjacent swaps", number of random adjacent swaps
     * performed.
     *
     * The default value is chosen according to the value that yielded the best
     * results among the tested values in the computational experiments of the
     * article cited above.
     */
    Counter adjacent_swaps_number_of_swaps = 4;

    /**
     * Perturbation "ruin-and-recreate".
     */
    Counter ruin_and_recreate_number_of_perturbations = 0;

    ElementPos ruin_number_of_elements_removed = 4;

    double ruin_random_weight = 1.0;

    double ruin_nearest_weight = 0.0;

    /**
     * Ruin method "adjacent string removal".
     *
     * See:
     * - "Slack Induction by String Removals for Vehicle Routing Problems"
     *   (Christiaens et Greet Vanden Berghe, 2020)
     *   https://doi.org/10.1287/trsc.2019.0914
     */
    double ruin_adjacent_string_removal_weight = 0.0;

    ElementPos ruin_adjacent_string_removal_maximum_string_cardinality = 10;

    double ruin_adjacent_string_removal_split_rate = 0.5;

    double ruin_adjacent_string_removal_beta = 0.01;

    double recreate_random_weight = 0.0;

    double recreate_best_weight = 0.0;

    /**
     * Perturbation "force-add".
     */
    bool force_add = false;

    /*
     * Crossovers.
     */

    double crossover_ox_weight = 0;

    double crossover_sjox_weight = 0;

    double crossover_sbox_weight = 0;

    double crossover_srex1_weight = 0;

    double crossover_srex2_weight = 0;

};

template <typename SequencingScheme>
class LocalScheme
{

public:

    using GlobalCost = typename SequencingScheme::GlobalCost;

    std::string to_string(const GlobalCost& gc) const
    {
        return localsearchsolver::to_string(sequencing_scheme_, gc);
    }

    using SequenceData = typename SequencingScheme::SequenceData;

    struct SolutionElement
    {
        SequenceId sequence_id = -1;
        ElementPos pos = -1;
        Mode mode = -1;
    };

    struct SequenceElement
    {
        ElementId element_id;
        Mode mode;

        bool operator==(const SequenceElement& rhs) const
        {
            return element_id == rhs.element_id && mode == rhs.mode;
        }
    };

    struct Sequence
    {
        SequenceId sequence_id = -1;
        std::vector<SequenceElement> elements;
        SequenceData data;
    };

    struct SubSequence
    {
        SubSequence(
                const Sequence& sequence,
                ElementPos pos_1,
                ElementPos pos_2,
                bool reverse = false):
            sequence(&sequence),
            pos_1(pos_1),
            pos_2(pos_2),
            reverse(reverse) { }

        SubSequence(
                ElementId j,
                Mode mode):
            j(j),
            mode(mode) { }

        const Sequence* sequence = nullptr;
        ElementPos pos_1 = -1;
        ElementPos pos_2 = -1;
        bool reverse = false;
        ElementId j = -1;
        Mode mode = -1;
    };

    struct Move
    {
        Neighborhoods type;
        ElementPos k1 = -1;
        ElementPos k2 = -1;

        SequenceId sequence_id_1 = -1;
        SequenceId sequence_id_2 = -1;
        ElementId element_id = -1;
        ElementPos pos_1 = -1;
        ElementPos pos_2 = -1;
        ElementPos pos_3 = -1;
        ElementPos pos_4 = -1;
        Mode mode = -1;

        GlobalCost global_cost;
    };

    struct Neighborhood
    {
        std::vector<Move> improving_moves = {};
        // modified_sequences[i] == true iff sequence i has changed since last
        // neighborhood exploration.
        std::vector<bool> modified_sequences;
        Counter number_of_explorations = 0;
        Counter number_of_successes = 0;
        double time = 0.0;
    };


    struct Solution
    {
        std::vector<Sequence> sequences;
        GlobalCost global_cost;

        std::vector<bool> modified_sequences;
    };

    using CompactSolution = std::vector<std::vector<SequenceElement>>;

    struct CompactSolutionHasher
    {
        std::hash<ElementId> hasher_j;
        std::hash<Mode> hasher_mode;

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
            for (const auto& sequence: *compact_solution) {
                size_t hash_tmp = 0;
                for (const SequenceElement& se: sequence) {
                    optimizationtools::hash_combine(hash_tmp, hasher_j(se.element_id));
                    optimizationtools::hash_combine(hash_tmp, hasher_mode(se.mode));
                }
                optimizationtools::hash_combine(hash, hash_tmp);
            }
            return hash;
        }
    };

    inline CompactSolutionHasher compact_solution_hasher() const { return CompactSolutionHasher(); }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        auto solution = empty_solution();
        for (SequenceId i = 0; i < number_of_sequences(); ++i)
            for (const SequenceElement& se: compact_solution[i])
                append(solution.sequences[i], se);
        compute_global_cost(solution);
        return solution;
    }

    CompactSolution solution2compact(const Solution& solution)
    {
        CompactSolution compact_solution(number_of_sequences());
        for (SequenceId i = 0; i < number_of_sequences(); ++i)
            compact_solution[i] = solution.sequences[i].elements;
        return compact_solution;
    }

    /*
     * Constructors and destructor.
     */

    LocalScheme(
            SequencingScheme& sequencing_scheme,
            Parameters parameters = Parameters()):
        sequencing_scheme_(sequencing_scheme),
        parameters_(parameters)
    {
        SequencePos m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();

        // Initialize temporary structures.
        sequence_datas_cur_1_ = std::vector<std::vector<SequenceData>>(m);
        sequence_datas_cur_2_ = std::vector<std::vector<std::vector<SequenceData>>>(m);
        global_costs_cur_ = std::vector<GlobalCost>(m);
        if (parameters_.linking_constraints) {
            if (m > 1)
                partial_global_costs_cur_1_ = std::vector<GlobalCost>(m);
            if (m > 2)
                partial_global_costs_cur_2_ = std::vector<std::vector<GlobalCost>>(
                        m, std::vector<GlobalCost>(m));
        }
        elements_cur_ = std::vector<SolutionElement>(n);
        solution_cur_ = empty_solution();
        solution_tmp_ = empty_solution();

        compute_sorted_neighbors();

        // Initialize neighborhoods_.

        neighborhoods_[int(Neighborhoods::Shift)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.shift_block_maximum_length + 1,
                    {{Neighborhood()}});
        neighborhoods_[int(Neighborhoods::Swap)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.swap_block_maximum_length + 1);
        for (ElementPos k2 = 1; k2 <= parameters_.swap_block_maximum_length; ++k2)
            neighborhoods_[int(Neighborhoods::Swap)][k2]
                = std::vector<Neighborhood>(
                        parameters_.swap_block_maximum_length + 1,
                        {Neighborhood()});
        neighborhoods_[int(Neighborhoods::Reverse)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::ShiftReverse)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.shift_reverse_block_maximum_length + 1,
                    {{Neighborhood()}});

        neighborhoods_[int(Neighborhoods::Add)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::Remove)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::Replace)] = {{{Neighborhood()}}};

        neighborhoods_[int(Neighborhoods::SwapTails)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::Split)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::InterShift)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.inter_shift_block_maximum_length + 1,
                    {{Neighborhood()}});
        neighborhoods_[int(Neighborhoods::InterSwap)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.swap_block_maximum_length + 1);
        for (ElementPos k2 = 1; k2 <= parameters_.swap_block_maximum_length; ++k2)
            neighborhoods_[int(Neighborhoods::InterSwap)][k2]
                = std::vector<Neighborhood>(
                        parameters_.swap_block_maximum_length + 1,
                        {Neighborhood()});
        neighborhoods_[int(Neighborhoods::InterShiftReverse)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.inter_shift_reverse_block_maximum_length + 1,
                    {{Neighborhood()}});
        neighborhoods_[int(Neighborhoods::InterSwapStar)] = {{{Neighborhood()}}};

        neighborhoods_[int(Neighborhoods::ShiftChangeMode)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::ModeSwap)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::SwapWithModes)] = {{{Neighborhood()}}};

        for (int a = 0; a < (int)neighborhoods_.size(); ++a) {
            for (int k1 = 0; k1 < (int)neighborhoods_[a].size(); ++k1) {
                for (int k2 = 0; k2 < (int)neighborhoods_[a][k1].size(); ++k2) {
                    Neighborhood& neighborhood = neighborhoods_[a][k1][k2];
                    neighborhood.modified_sequences = std::vector<bool>(m, true);
                }
            }
        }

        if (!parameters_.linking_constraints
                && parameters_.inter_shift_block_maximum_length >= 1) {
            inter_shift_1_best_positions_
                = std::vector<std::vector<std::pair<ElementPos, GlobalCost>>>(
                        m, std::vector<std::pair<ElementPos, GlobalCost>>(
                            n, {-2, GlobalCost()}));
        }

        if (parameters_.inter_swap_star) {
            inter_swap_star_best_positions_
                = std::vector<std::vector<std::vector<ElementPos>>>(
                        m, std::vector<std::vector<ElementPos>>(
                            n, std::vector<ElementPos>(3, -2)));
        }
    }

    LocalScheme(const LocalScheme& sequencing_scheme):
        LocalScheme(sequencing_scheme.sequencing_scheme_, sequencing_scheme.parameters_) { }

    virtual ~LocalScheme() { }

    inline Sequence empty_sequence(SequenceId sequence_id) const
    {
        Sequence sequence;
        sequence.sequence_id = sequence_id;
        sequence.data = empty_sequence_data(sequence_id);
        return sequence;
    }

    inline Solution empty_solution() const
    {
        Solution solution;
        for (SequenceId i = 0; i < number_of_sequences(); ++i)
            solution.sequences.push_back(empty_sequence(i));
        return solution;
    }

    /*
     * Initial solution.
     */

    Solution split(
            const std::vector<SequenceElement>& elements,
            std::mt19937_64& generator)
    {
        SequencePos m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        ElementPos seq_size = elements.size();

        if (m == 1) {
            Solution solution = empty_solution();
            for (ElementPos pos = 0; pos < seq_size; ++pos)
                append(solution.sequences[0], elements[pos]);
        }

        //for (ElementPos pos = 0; pos < seq_size; ++pos)
        //    std::cout << " " << elements[pos].element_id;
        //std::cout << std::endl;

        std::uniform_int_distribution<SequenceId> d_i(0, m - 1);
        SequenceId i = d_i(generator);

        // Compute edge costs.
        std::vector<std::vector<GlobalCost>> edges(n + 1);
        for (ElementPos pos = 0; pos <= n; ++pos)
            edges[pos] = std::vector<GlobalCost>(pos);
        for (ElementPos pos_1 = 0; pos_1 < seq_size; ++pos_1) {
            Sequence sequence = empty_sequence(i);
            for (ElementPos pos_2 = pos_1 + 1; pos_2 <= seq_size; ++pos_2) {
                append(sequence, elements[pos_2 - 1]);
                edges[pos_2][pos_1] = sequencing_scheme_.global_cost(sequence.data);
            }
        }

        // Run Bellman algorithm.
        std::vector<std::vector<ElementPos>> prev(
                m + 1, std::vector<ElementPos>(n + 1, -1));
        std::vector<std::vector<int8_t>> distance_init(
                m + 1, std::vector<int8_t>(n + 1, false));
        std::vector<std::vector<GlobalCost>> distance(
                m + 1, std::vector<GlobalCost>(n + 1, GlobalCost()));
        distance[0][0] = sequencing_scheme_.global_cost(empty_sequence_data(i));
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            for (ElementPos pos_1 = 0; pos_1 < seq_size; ++pos_1) {
                for (ElementPos pos_2 = pos_1 + 1; pos_2 <= seq_size; ++pos_2) {
                    GlobalCost d = merge(
                            distance[i][pos_1],
                            edges[pos_2][pos_1]);
                    if (!distance_init[i + 1][pos_2]
                            || strictly_better(d, distance[i + 1][pos_2])) {
                        //std::cout << "Update " << i
                        //    << " pos_2 " << pos_2
                        //    << " from pos_1 " << pos_1
                        //    << " dv " << to_string(distance[i + 1][pos_2])
                        //    << " du " << to_string(distance[i][pos_1])
                        //    << " uv " << to_string(edges[pos_2][pos_1])
                        //    << " du+uv " << to_string(distance[i][pos_1] + edges[pos_2][pos_1])
                        //    << std::endl;
                        distance_init[i + 1][pos_2] = true;
                        distance[i + 1][pos_2] = d;
                        prev[i + 1][pos_2] = pos_1;
                    }
                }
            }
            //ElementPos pos = seq_size;
            //SequenceId i_cur = i + 1;
            //while (pos > 0) {
            //    pos = prev[i_cur][pos];
            //    i_cur--;
            //    std::cout << " " << pos;
            //}
            //std::cout << std::endl;
        }

        // Retrieve solution.
        //std::cout << "Retrieve solution" << std::endl;
        ElementPos pos = seq_size;
        ElementPos pos_prev = seq_size;
        Solution solution = empty_solution();
        for (SequenceId i_cur = 0; pos != 0; ++i_cur) {
            //std::cout << "i_cur " << i_cur
            //    << " / " << m
            //    << " pos " << pos
            //    << " pos_prev " << pos_prev
            //    << " prev " << prev[m - i_cur][pos]
            //    << std::endl;
            pos_prev = pos;
            pos = prev[m - i_cur][pos];
            for (ElementPos p = pos; p < pos_prev; ++p) {
                //std::cout
                //    << "i_cur " << i_cur << " / " << m
                //    << " p " << p << " / " << elements.size()
                //    << std::endl;
                append(solution.sequences[i_cur], elements[p]);
                //std::cout << to_string(sequencing_scheme_, sequencing_scheme_.global_cost(solution.sequences[i_cur].data)) << std::endl;
            }
        }
        compute_global_cost(solution);

        //print(std::cout, solution);

        return solution;
    }

    template<typename, typename T>
    struct HasInitialSolutionMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasInitialSolutionMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().initial_solution(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    Solution initial_solution(
            Counter initial_solution_id,
            std::mt19937_64& generator,
            std::false_type)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();

        switch (initial_solution_id) {
        case 1: {
            // Fix a random order, find the best position.
            std::vector<ElementId> elements(n);
            std::iota(elements.begin(), elements.end(), 0);
            std::shuffle(elements.begin(), elements.end(), generator);
            Solution solution = empty_solution();
            compute_temporary_structures(solution);
            for (ElementId j: elements) {
                auto improving_moves = explore_add(solution, j);
                if (!improving_moves.empty()) {
                    std::shuffle(
                            improving_moves.begin(),
                            improving_moves.end(),
                            generator);
                    Move move_best;
                    for (const Move move: improving_moves)
                        if (move_best.sequence_id_1 == -1 || dominates(
                                    move.global_cost,
                                    move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    compute_temporary_structures(solution, move_best.sequence_id_1, move_best.sequence_id_2);
                }
            }
            return solution;
        } case 2: {
            // Insert at the end, find the best element.
            std::vector<uint8_t> contains(sequencing_scheme_.number_of_elements(), 0);
            Solution solution = empty_solution();
            compute_temporary_structures(solution);
            for (;;) {
                auto improving_moves = explore_add_end(solution);
                if (!improving_moves.empty()) {
                    std::shuffle(
                            improving_moves.begin(),
                            improving_moves.end(),
                            generator);
                    Move move_best;
                    for (const Move move: improving_moves)
                        if (move_best.sequence_id_1 == -1 || dominates(
                                    move.global_cost,
                                    move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    compute_temporary_structures(solution, move_best.sequence_id_1, move_best.sequence_id_2);
                } else {
                    break;
                }
            }
            return solution;
        } case 3: {
            // Find the best element x position.
            // Warning, this one can be expensive.
            std::vector<uint8_t> contains(sequencing_scheme_.number_of_elements(), 0);
            Solution solution = empty_solution();
            compute_temporary_structures(solution);
            for (;;) {
                auto improving_moves = explore_add_2(solution);
                if (!improving_moves.empty()) {
                    std::shuffle(
                            improving_moves.begin(),
                            improving_moves.end(),
                            generator);
                    Move move_best;
                    for (const Move move: improving_moves)
                        if (move_best.sequence_id_1 == -1 || dominates(
                                    move.global_cost,
                                    move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    compute_temporary_structures(solution, move_best.sequence_id_1, move_best.sequence_id_2);
                } else {
                    break;
                }
            }
            return solution;
        } case 4: {
            // Random sequence splitted.
            std::vector<SequenceElement> elements(n);
            for (ElementId element_id = 0; element_id < n; ++element_id) {
                std::uniform_int_distribution<SequenceId> d_mode(0, number_of_modes(element_id) - 1);
                Mode mode = d_mode(generator);
                elements[element_id] = {element_id, mode};
            }
            std::shuffle(elements.begin(), elements.end(), generator);
            Solution solution = split(elements, generator);
            return solution;
        } default: {
            // Random permutation.
            std::vector<ElementId> elements(n);
            std::iota(elements.begin(), elements.end(), 0);
            std::shuffle(elements.begin(), elements.end(), generator);
            std::uniform_int_distribution<SequenceId> di(0, m - 1);
            Solution solution = empty_solution();
            for (ElementId j: elements) {
                std::uniform_int_distribution<SequenceId> d_mode(0, number_of_modes(j) - 1);
                Mode mode = d_mode(generator);
                SequenceId i = di(generator);
                append(solution.sequences[i], {j, mode});
            }
            compute_global_cost(solution);
            return solution;
        }
        }
    }

    Solution initial_solution(
            Counter initial_solution_id,
            std::mt19937_64& generator,
            std::true_type)
    {
        if (initial_solution_id >= 4)
            return sequencing_scheme_.initial_solution(initial_solution_id, generator);
        return initial_solution(initial_solution_id, generator, false);
    }

    inline Solution initial_solution(
            Counter initial_solution_id,
            std::mt19937_64& generator)
    {
        //std::cout << "initial_solution " << initial_solution_id << std::endl;
        auto begin = std::chrono::steady_clock::now();
        Solution solution = initial_solution(
                initial_solution_id,
                generator,
                std::integral_constant<
                    bool,
                    HasInitialSolutionMethod<SequencingScheme,
                    Solution(
                        Counter,
                        std::mt19937_64&)>::value>());
        auto end = std::chrono::steady_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
        initial_solution_time_ += time_span.count();
        number_of_initial_solution_calls_++;
        //std::cout << to_string(solution.global_cost) << std::endl;
        return solution;
    }

    /*
     * Crossovers.
     */

    inline Solution crossover(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        auto begin = std::chrono::steady_clock::now();
        std::discrete_distribution<Counter> d_crossover({
                parameters_.crossover_ox_weight,
                parameters_.crossover_sjox_weight,
                parameters_.crossover_sbox_weight,
                parameters_.crossover_srex1_weight,
                parameters_.crossover_srex2_weight,
                });
        Counter x = d_crossover(generator);
        Solution solution_1;
        Solution solution_2;
        switch (x) {
        case 0: {
            solution_1 = crossover_ox(solution_parent_1, solution_parent_2, generator);
            solution_2 = crossover_ox(solution_parent_2, solution_parent_1, generator);
            break;
        } case 1: {
            solution_1 = crossover_sjox(solution_parent_1, solution_parent_2, generator);
            solution_2 = crossover_sjox(solution_parent_2, solution_parent_1, generator);
            break;
        } case 2: {
            solution_1 = crossover_sbox(solution_parent_1, solution_parent_2, generator);
            solution_2 = crossover_sbox(solution_parent_2, solution_parent_1, generator);
            break;
        } case 3: {
            solution_1 = crossover_srex1(solution_parent_1, solution_parent_2, generator);
            solution_2 = crossover_srex1(solution_parent_2, solution_parent_1, generator);
            break;
        } case 4: {
            solution_1 = crossover_srex2(solution_parent_1, solution_parent_2, generator);
            solution_2 = crossover_srex2(solution_parent_2, solution_parent_1, generator);
            break;
        } default: {
            solution_1 = crossover_ox(solution_parent_1, solution_parent_2, generator);
            solution_2 = crossover_ox(solution_parent_2, solution_parent_1, generator);
            break;
        }
        }
        Solution solution;
        if (dominates(global_cost(solution_1), global_cost(solution_2))) {
            solution = solution_1;
        } else if (dominates(global_cost(solution_2), global_cost(solution_1))) {
            solution = solution_2;
        } else {
            solution = solution_1;
        }

        SequenceId m = number_of_sequences();
        solution.modified_sequences = std::vector<bool>(m, true);
        SequencePos number_of_non_empty_sequences = 0;
        SequencePos number_of_modified_sequences = m;
        for (SequenceId i1 = 0; i1 < m; ++i1) {
            if (solution.sequences[i1].elements.size() > 0)
                number_of_non_empty_sequences++;
            for (SequenceId i2 = 0; i2 < m; ++i2) {
                if (solution.sequences[i1].elements
                        == solution_parent_1.sequences[i2].elements
                        && equals(
                            sequencing_scheme_.global_cost(solution.sequences[i1].data),
                            sequencing_scheme_.global_cost(solution_parent_1.sequences[i2].data))) {
                    solution.modified_sequences[i1] = false;
                    number_of_modified_sequences--;
                    break;
                }
                if (solution.sequences[i1].elements
                        == solution_parent_2.sequences[i2].elements
                        && equals(
                            sequencing_scheme_.global_cost(solution.sequences[i1].data),
                            sequencing_scheme_.global_cost(solution_parent_2.sequences[i2].data))) {
                    solution.modified_sequences[i1] = false;
                    number_of_modified_sequences--;
                    break;
                }
            }
        }
        (void)number_of_modified_sequences;
        (void)number_of_non_empty_sequences;
        //std::cout << "number_of_modified_sequences " << number_of_modified_sequences << std::endl;
        //std::cout << "number_of_non_empty_sequences " << number_of_non_empty_sequences << std::endl;

        auto end = std::chrono::steady_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
        crossover_time_ += time_span.count();
        number_of_crossover_calls_++;
        return solution;
    }

    /**
     * Generate a new solution from two parent solutions using the OX crossover
     * operator.
     *
     * The OX crossover operator consists in selecting a random substring from
     * the first parent, copying this substring into the child while leaving
     * the rest of the positions empty, and finally completing the child’s
     * missing positions circularly with the visits from the second parent,
     * starting from the end cutting point.
     *
     * References:
     * - "A simple and effective hybrid genetic search for the job sequencing
     *   and tool switching problem" (Mecler et al., 2021)
     *   https://doi.org/10.1016/j.cor.2020.105153
     */
    inline Solution crossover_ox(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        SequenceId m = number_of_sequences();

        std::vector<SequenceElement> elements_parent_1;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            elements_parent_1.insert(
                    elements_parent_1.begin(),
                    solution_parent_1.sequences[sequence_id].elements.begin(),
                    solution_parent_1.sequences[sequence_id].elements.end());
        }
        std::vector<SequenceElement> elements_parent_2;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            elements_parent_2.insert(
                    elements_parent_2.begin(),
                    solution_parent_2.sequences[sequence_id].elements.begin(),
                    solution_parent_2.sequences[sequence_id].elements.end());
        }

        ElementPos seq_1_size = elements_parent_1.size();

        std::vector<ElementPos> edges = optimizationtools::bob_floyd<ElementPos>(
                2, seq_1_size + 1, generator);
        std::sort(edges.begin(), edges.end());

        ElementPos pos_1 = edges[0];
        ElementPos pos_2 = edges[1];

        std::vector<uint8_t> in_substring(seq_1_size, false);
        for (ElementPos pos = pos_1; pos < pos_2; ++pos) {
            ElementId j = elements_parent_1[pos].element_id;
            in_substring[j] = true;
        }

        std::vector<SequenceElement> elements;

        for (ElementPos pos = 0; pos < seq_1_size; ++pos) {
            if ((ElementPos)elements.size() == pos_1)
                for (ElementPos p = pos_1; p < pos_2; ++p)
                    elements.push_back(elements_parent_1[p]);
            ElementId j = elements_parent_2[pos].element_id;
            if (in_substring[j])
                continue;
            elements.push_back(elements_parent_2[pos]);
        }
        if ((ElementPos)elements.size() == pos_1)
            for (ElementPos p = pos_1; p < pos_2; ++p)
                elements.push_back(elements_parent_1[p]);

        return split(elements, generator);
    }

    /**
     * Generate a new solution from two parent solutions using the SJOX
     * crossover operator.
     *
     * References:
     * - "Two new robust genetic algorithms for the flowshop scheduling
     *   problem" (Ruiz et al., 2006)
     *   https://doi.org/10.1016/j.omega.2004.12.006
     */
    inline Solution crossover_sjox(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        SequenceId m = number_of_sequences();

        std::vector<SequenceElement> elements_parent_1;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            elements_parent_1.insert(
                    elements_parent_1.begin(),
                    solution_parent_1.sequences[sequence_id].elements.begin(),
                    solution_parent_1.sequences[sequence_id].elements.end());
        }
        std::vector<SequenceElement> elements_parent_2;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            elements_parent_2.insert(
                    elements_parent_2.begin(),
                    solution_parent_2.sequences[sequence_id].elements.begin(),
                    solution_parent_2.sequences[sequence_id].elements.end());
        }

        ElementPos seq_1_size = elements_parent_1.size();
        std::vector<ElementPos> positions(seq_1_size, -1);

        std::vector<SequenceElement> elements;

        // Add elements from parent_1 up to a given cut point.
        std::uniform_int_distribution<ElementPos> d_point(1, seq_1_size);
        ElementPos pos_0 = d_point(generator);
        for (ElementPos pos = 0; pos < pos_0; ++pos) {
            ElementId j = elements_parent_1[pos].element_id;
            positions[j] = pos;
            elements.push_back(elements_parent_1[pos]);
        }

        // Add elements from parent_2 keeping the relative order.
        for (ElementPos pos = 0; pos < seq_1_size; ++pos) {
            // Add elements which have the same positions in both parents.
            for (;;) {
                ElementPos p = elements.size();
                if (p == seq_1_size)
                    break;
                ElementId j1 = elements_parent_1[p].element_id;
                ElementId j2 = elements_parent_2[p].element_id;
                if (j1 == j2) {
                    positions[j1] = p;
                    elements.push_back(elements_parent_1[p]);
                    continue;
                }
                break;
            }
            ElementId j = elements_parent_2[pos].element_id;
            if (positions[j] != -1)
                continue;
            positions[j] = elements.size();
            elements.push_back(elements_parent_2[pos]);
        }

        return split(elements, generator);
    }

    /**
     * Generate a new solution from two parent solutions using the SBOX
     * crossover operator.
     *
     * References:
     * - "Two new robust genetic algorithms for the flowshop scheduling
     *   problem" (Ruiz et al., 2006)
     *   https://doi.org/10.1016/j.omega.2004.12.006
     */
    inline Solution crossover_sbox(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        SequenceId m = number_of_sequences();

        std::vector<SequenceElement> elements_parent_1;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            elements_parent_1.insert(
                    elements_parent_1.begin(),
                    solution_parent_1.sequences[sequence_id].elements.begin(),
                    solution_parent_1.sequences[sequence_id].elements.end());
        }
        std::vector<SequenceElement> elements_parent_2;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            elements_parent_2.insert(
                    elements_parent_2.begin(),
                    solution_parent_2.sequences[sequence_id].elements.begin(),
                    solution_parent_2.sequences[sequence_id].elements.end());
        }

        ElementPos seq_1_size = elements_parent_1.size();
        std::vector<ElementPos> positions(seq_1_size, -1);

        std::vector<SequenceElement> elements;

        // Add elements from parent_1 up to a given cut point.
        std::uniform_int_distribution<ElementPos> d_point(1, seq_1_size);
        ElementPos pos_0 = d_point(generator);
        for (ElementPos pos = 0; pos < pos_0; ++pos) {
            ElementId j = elements_parent_1[pos].element_id;
            positions[j] = pos;
            elements.push_back(elements_parent_1[pos]);
        }

        // Add elements from parent_2 keeping the relative order.
        for (ElementPos pos = 0; pos < seq_1_size; ++pos) {
            // Add elements which have the same positions in both parents.
            for (;;) {
                ElementPos p = elements.size();
                if (p == seq_1_size)
                    break;
                ElementId j1 = elements_parent_1[p].element_id;
                ElementId j2 = elements_parent_2[p].element_id;
                if (p <= seq_1_size - 1) {
                    ElementId j1_next = elements_parent_1[p + 1].element_id;
                    ElementId j2_next = elements_parent_2[p + 1].element_id;
                    if (j1 == j2 && j1_next == j2_next) {
                        positions[j1] = p;
                        elements.push_back(elements_parent_1[p + 1]);
                        continue;
                    }
                }
                if (p >= 1) {
                    ElementId j1_prev = elements_parent_1[p - 1].element_id;
                    ElementId j2_prev = elements_parent_2[p - 1].element_id;
                    if (j1 == j2 && j1_prev == j2_prev) {
                        positions[j1] = p;
                        elements.push_back(elements_parent_1[p - 1]);
                        continue;
                    }
                }
                break;
            }
            ElementId j = elements_parent_2[pos].element_id;
            if (positions[j] != -1)
                continue;
            positions[j] = elements.size();
            elements.push_back(elements_parent_2[pos]);
        }

        return split(elements, generator);
    }

    /**
     * Crossover SREX 1
     *
     * The idea is to keep some routes from the first parent solution and to
     * fill the rest of the child solution with routes or partial routes from
     * the second parent solution.
     *
     * This variant kind of assumes that the sequences are homogeneous.
     */
    inline Solution crossover_srex1(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        //std::cout << "crossover_srex1 start" << std::endl;
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();

        // Compute the number of sequences to remove.
        SequenceId m1 = 0;
        SequenceId m2 = 0;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!solution_parent_1.sequences[sequence_id].elements.empty())
                m1++;
            if (!solution_parent_2.sequences[sequence_id].elements.empty())
                m2++;
        }
        std::uniform_int_distribution<SequencePos> d_m(1, std::min(m1, m2));
        SequencePos number_of_sequences_to_remove = (std::min(m1, m2) <= 1)? 1: d_m(generator);
        //std::cout << "m1 " << m1
        //    << " m2 " << m2
        //    << " number_of_sequences_to_remove " << number_of_sequences_to_remove
        //    << std::endl;

        // Child solution.
        Solution solution = solution_parent_1;

        // Get elements from the first parent solution and store the sequences
        // they belong to.
        std::vector<SequenceId> sequences(n, -1);
        optimizationtools::IndexedSet elts(n);
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            const auto& elements = solution_parent_1.sequences[sequence_id].elements;
            ElementPos seq_size = (ElementPos)elements.size();
            for (ElementPos pos = 0; pos < seq_size; ++pos) {
                ElementId j = elements[pos].element_id;
                elts.add(j);
                sequences[j] = sequence_id;
            }
        }

        if (sorted_successors_.empty()) {

            // If the 'distance' method has not been provided, remove random
            // sequences.
            std::vector<SequenceId> seqs_1;
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (!solution_parent_1.sequences[sequence_id].elements.empty())
                    seqs_1.push_back(sequence_id);
            std::shuffle(seqs_1.begin(), seqs_1.end(), generator);
            for (auto it = seqs_1.begin();
                    it != seqs_1.begin() + number_of_sequences_to_remove; ++it) {
                SequenceId i = *it;
                // Remove the sequence containing j2 from the child solution.
                for (const auto& se: solution.sequences[i].elements)
                    elts.remove(se.element_id);
                solution.sequences[i] = empty_sequence(i);
            }

        } else {

            // Draw seed element.
            std::uniform_int_distribution<ElementPos> d(0, elts.size() - 1);
            ElementPos pos = d(generator);
            ElementId j = *(elts.begin() + pos);

            // Remove sequence containing the seed element from the child solution..
            SequencePos number_of_removed_sequences = 1;
            SequenceId i_j = sequences[j];
            for (const auto& se: solution.sequences[i_j].elements)
                elts.remove(se.element_id);
            solution.sequences[i_j] = empty_sequence(i_j);

            // Loop through closest elements from the seed element.
            for (ElementPos p = 0;
                    number_of_removed_sequences < number_of_sequences_to_remove; ++p) {
                ElementId j2 = (p % 2 == 0)?
                    sorted_successors_[j][p / 2]:
                    sorted_predecessors_[j][p / 2];
                SequenceId i_j2 = sequences[j2];
                // If j2 still belongs to the child solution.
                if (elts.contains(j2)) {
                    // Remove the sequence containing j2 from the child solution.
                    for (const auto& se: solution.sequences[i_j2].elements)
                        elts.remove(se.element_id);
                    solution.sequences[i_j2] = empty_sequence(i_j2);
                    number_of_removed_sequences++;
                }
            }

        }

        // Add sequences from the second parent solution.
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!solution.sequences[sequence_id].elements.empty())
                continue;
            // Find the sequence from the second parent with the highest number
            // of unscheduled elements.
            SequenceId i2_best = -1;
            ElementPos seq_size_best = -1;
            for (SequenceId i2 = 0; i2 < m; ++i2) {
                ElementPos seq_size = 0;
                for (const auto& se: solution_parent_2.sequences[i2].elements)
                    if (!elts.contains(se.element_id))
                        seq_size++;
                if (seq_size <= (SequencePos)solution_parent_2.sequences[i2].elements.size() / 2)
                    continue;
                if (i2_best == -1 || seq_size_best < seq_size) {
                    i2_best = i2;
                    seq_size_best = seq_size;
                }
            }
            if (i2_best == -1)
                break;
            //std::cout << "i2_best " << i2_best << " seq_size_best " << seq_size_best << std::endl;

            solution.sequences[sequence_id] = empty_sequence(sequence_id);
            for (const auto& se: solution_parent_2.sequences[i2_best].elements) {
                if (!elts.contains(se.element_id)) {
                    append(solution.sequences[sequence_id], se);
                    elts.add(se.element_id);
                }
            }
            //std::cout << "i " << i
            //    << " gc " << to_string(global_cost(solution.sequences[i]))
            //    << std::endl;
        }
        //std::cout << "elts.size() " << elts.size() << " / " << n << std::endl;

        // Add remaining elements at random positions.
        compute_temporary_structures(solution);
        for (ElementId element_id = 0; element_id < n; ++element_id) {
            if (elts.contains(element_id))
                continue;
            auto improving_moves = explore_add(solution, element_id);
            if (!improving_moves.empty()) {
                std::shuffle(
                        improving_moves.begin(),
                        improving_moves.end(),
                        generator);
                Move move_best;
                for (const Move move: improving_moves)
                    if (move_best.sequence_id_1 == -1 || dominates(
                                move.global_cost,
                                move_best.global_cost))
                        move_best = move;
                apply_move(solution, move_best);
                compute_temporary_structures(solution, move_best.sequence_id_1, move_best.sequence_id_2);
            }
        }

        compute_global_cost(solution);
        //std::cout << "crossover_srex1 end" << std::endl;
        return solution;
    }

    /**
     * Crossover SREX 2.
     *
     * This variant kind of assumes that the sequences are heterogeneous.
     */
    inline Solution crossover_srex2(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();

        // Compute the number of sequences to remove.
        SequenceId m1 = 0;
        SequenceId m2 = 0;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!solution_parent_1.sequences[sequence_id].elements.empty())
                m1++;
            if (!solution_parent_2.sequences[sequence_id].elements.empty())
                m2++;
        }
        std::uniform_int_distribution<SequencePos> d_m(1, std::min(m1, m2));
        SequencePos number_of_sequences_to_remove = (std::min(m1, m2) <= 1)? 1: d_m(generator);
        //std::cout << "m1 " << m1
        //    << " m2 " << m2
        //    << " number_of_sequences_to_remove " << number_of_sequences_to_remove
        //    << std::endl;

        // Child solution.
        Solution solution = solution_parent_1;

        // Get elements from the first parent solution and store the sequences
        // they belong to.
        std::vector<SequenceId> sequences(n, -1);
        optimizationtools::IndexedSet elts(n);
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            const auto& elements = solution_parent_1.sequences[sequence_id].elements;
            ElementPos seq_size = (ElementPos)elements.size();
            for (ElementPos pos = 0; pos < seq_size; ++pos) {
                ElementId j = elements[pos].element_id;
                elts.add(j);
                sequences[j] = sequence_id;
            }
        }

        if (sorted_successors_.empty()) {

            // If the 'distance' method has not been provided, remove random
            // sequences.
            std::vector<SequenceId> seqs_1;
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (!solution_parent_1.sequences[sequence_id].elements.empty())
                    seqs_1.push_back(sequence_id);
            std::shuffle(seqs_1.begin(), seqs_1.end(), generator);
            for (auto it = seqs_1.begin();
                    it != seqs_1.begin() + number_of_sequences_to_remove; ++it) {
                SequenceId i = *it;
                // Remove the sequence containing j2 from the child solution.
                for (const auto& se: solution.sequences[i].elements)
                    elts.remove(se.element_id);
                solution.sequences[i] = empty_sequence(i);
            }

        } else {

            // Draw seed element.
            std::uniform_int_distribution<ElementPos> d(0, elts.size() - 1);
            ElementPos pos = d(generator);
            ElementId j = *(elts.begin() + pos);

            // Remove sequence containing the seed element from the child solution..
            SequencePos number_of_removed_sequences = 1;
            SequenceId i_j = sequences[j];
            for (const auto& se: solution.sequences[i_j].elements)
                elts.remove(se.element_id);
            solution.sequences[i_j] = empty_sequence(i_j);

            // Loop through closest elements from the seed element.
            for (ElementPos p = 0;
                    number_of_removed_sequences < number_of_sequences_to_remove; ++p) {
                ElementId j2 = (p % 2 == 0)?
                    sorted_successors_[j][p / 2]:
                    sorted_predecessors_[j][p / 2];
                SequenceId i_j2 = sequences[j2];
                // If j2 still belongs to the child solution.
                if (elts.contains(j2)) {
                    // Remove the sequence containing j2 from the child solution.
                    for (const auto& se: solution.sequences[i_j2].elements)
                        elts.remove(se.element_id);
                    solution.sequences[i_j2] = empty_sequence(i_j2);
                    number_of_removed_sequences++;
                }
            }

        }

        // Add sequences from the second parent.
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!solution.sequences[sequence_id].elements.empty())
                continue;
            solution.sequences[sequence_id] = empty_sequence(sequence_id);
            for (const auto& se: solution_parent_2.sequences[sequence_id].elements) {
                if (!elts.contains(se.element_id)) {
                    append(solution.sequences[sequence_id], se);
                    elts.add(se.element_id);
                }
            }
        }

        // Add remaining elements at random positions.
        for (ElementId element_id = 0; element_id < n; ++element_id) {
            if (elts.contains(element_id))
                continue;
            auto improving_moves = explore_add(solution, element_id);
            if (!improving_moves.empty()) {
                std::shuffle(
                        improving_moves.begin(),
                        improving_moves.end(),
                        generator);
                Move move_best;
                for (const Move move: improving_moves)
                    if (move_best.sequence_id_1 == -1 || dominates(
                                move.global_cost,
                                move_best.global_cost))
                        move_best = move;
                apply_move(solution, move_best);
                compute_temporary_structures(solution, move_best.sequence_id_1, move_best.sequence_id_2);
            }
        }

        compute_global_cost(solution);
        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return solution.global_cost;
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return sequencing_scheme_.bound(sequence_data);
    }

    inline GlobalCost bound(const Sequence& sequence) const
    {
        return sequencing_scheme_.bound(sequence.data);
    }

    GlobalCost global_cost_goal(double value) const
    {
        return localsearchsolver::global_cost_goal(sequencing_scheme_, value);
    }

    inline bool strictly_better(
            const GlobalCost& gc1,
            const GlobalCost& gc2) const
    {
        return localsearchsolver::strictly_better(sequencing_scheme_, gc1, gc2);
    }

    inline bool dominates(
            const GlobalCost& gc1,
            const GlobalCost& gc2) const
    {
        return localsearchsolver::dominates(sequencing_scheme_, gc1, gc2);
    }

    inline bool equals(
            const GlobalCost& gc1,
            const GlobalCost& gc2) const
    {
        return !strictly_better(gc1, gc2)
            && !strictly_better(gc2, gc1);
    }

    /**
     * Return the distance between two solutions.
     *
     * The distance between two solutions (represented as element permutations)
     * is defined as the number of different edges between them (broken-pairs
     * distance).
     *
     * References:
     * - "Two memetic algorithms for heterogeneous fleet vehicle routing
     *   problems" (Prins, 2009)
     *   https://doi.org/10.1016/j.engappai.2008.10.006
     * - "A simple and effective hybrid genetic search for the job sequencing
     *   and tool switching problem" (Mecler et al., 2021)
     *   https://doi.org/10.1016/j.cor.2020.105153
     */
    inline ElementPos distance(
            const Solution& solution_1,
            const Solution& solution_2) const
    {
        SequencePos m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();

        // For each element, get its next element in the first solution.
        std::vector<ElementId> next_1(n, -1);
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            const auto& elements = solution_1.sequences[sequence_id].elements;
            ElementPos seq_size = elements.size();
            for (ElementPos pos = 0; pos < seq_size - 1; ++pos) {
                ElementId j = elements[pos].element_id;
                ElementId j_next = elements[pos + 1].element_id;
                next_1[j] = j_next;
            }
        }

        // Compare with second solution.
        ElementPos d = 0;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            const auto& elements = solution_2.sequences[sequence_id].elements;
            ElementPos seq_size = elements.size();
            for (ElementPos pos = 0; pos < seq_size - 1; ++pos) {
                ElementId j = elements[pos].element_id;
                ElementId j_next = elements[pos + 1].element_id;
                if (j_next != next_1[j])
                    d++;
            }
        }
        return d;
    }

    /*
     * Local search.
     */

    struct Perturbation
    {
        Perturbation(): type(Perturbations::None) { }

        /** Type of perturbation. */
        Perturbations type;
        /** ForceAdd: element to add. */
        ElementId force_add_j = -1;
        /** Global cost of the perturbation. */
        GlobalCost global_cost;
    };

    struct PerturbationHasher
    {
        std::hash<ElementPos> hasher;

        inline bool hashable(const Perturbation& perturbation) const
        {
            if (perturbation.type == Perturbations::ForceAdd)
                return true;
            return false;
        }

        inline bool operator()(
                const Perturbation& perturbation_1,
                const Perturbation& perturbation_2) const
        {
            if (perturbation_1.type == Perturbations::ForceAdd
                    && perturbation_2.type == Perturbations::ForceAdd
                    && perturbation_1.force_add_j == perturbation_2.force_add_j)
                return true;
            return false;
        }

        inline std::size_t operator()(
                const Perturbation& perturbation) const
        {
            return hasher(perturbation.force_add_j);
        }
    };

    inline PerturbationHasher perturbation_hasher() const { return PerturbationHasher(); }

    /*
     * Perturbations.
     */

    inline std::vector<Perturbation> perturbations(
            const Solution& solution,
            std::mt19937_64& generator)
    {
        std::vector<Perturbation> perturbations;

        // Double-bridge.
        for (Counter perturbation_id = 0;
                perturbation_id < parameters_.double_bridge_number_of_perturbations;
                ++perturbation_id) {
            Perturbation perturbation;
            perturbation.type = Perturbations::DoubleBridge;
            perturbation.global_cost = global_cost(solution);
            perturbations.push_back(perturbation);
        }

        // Adjacent swaps.
        for (Counter perturbation_id = 0;
                perturbation_id < parameters_.adjacent_swaps_number_of_perturbations;
                ++perturbation_id) {
            Perturbation perturbation;
            perturbation.type = Perturbations::AdjacentSwaps;
            perturbation.global_cost = global_cost(solution);
            perturbations.push_back(perturbation);
        }

        // Ruin-and-recreate.
        for (Counter perturbation_id = 0;
                perturbation_id < parameters_.ruin_and_recreate_number_of_perturbations;
                ++perturbation_id) {
            Perturbation perturbation;
            perturbation.type = Perturbations::RuinAndRecreate;
            perturbation.global_cost = global_cost(solution);
            perturbations.push_back(perturbation);
        }

        // Force-add.
        if (parameters_.force_add) {
            std::vector<uint8_t> contains(sequencing_scheme_.number_of_elements(), 0);
            for (SequenceId i = 0; i < number_of_sequences(); ++i)
                for (SequenceElement se: solution.sequences[i].elements)
                    contains[se.element_id] = 1;
            for (ElementId j = 0; j < sequencing_scheme_.number_of_elements(); ++j) {
                if (contains[j])
                    continue;
                Perturbation perturbation;
                perturbation.type = Perturbations::ForceAdd;
                perturbation.force_add_j = j;
                perturbation.global_cost = global_cost(solution);
                perturbations.push_back(perturbation);
            }
        }

        std::shuffle(perturbations.begin(), perturbations.end(), generator);
        return perturbations;
    }

    inline void apply_perturbation(
            Solution& solution,
            Perturbation& perturbation,
            std::mt19937_64& generator)
    {
        SequenceId m = number_of_sequences();

        switch (perturbation.type) {
        case Perturbations::None: {
            break;

        } case Perturbations::DoubleBridge: {
            std::uniform_int_distribution<SequenceId> distribution(0, m - 1);
            SequenceId i = distribution(generator);
            auto positions = optimizationtools::bob_floyd<ElementPos>(
                    4,
                    solution.sequences[i].elements.size() + 1,
                    generator);
            std::sort(positions.begin(), positions.end());

            const auto& elements = solution.sequences[i].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[i];
            sequence_tmp = empty_sequence(i);
            for (ElementPos pos = 0; pos < positions[0]; ++pos)
                append(sequence_tmp, elements[pos]);
            for (ElementPos pos = positions[2]; pos < positions[3]; ++pos)
                append(sequence_tmp, elements[pos]);
            for (ElementPos pos = positions[1]; pos < positions[2]; ++pos)
                append(sequence_tmp, elements[pos]);
            for (ElementPos pos = positions[0]; pos < positions[1]; ++pos)
                append(sequence_tmp, elements[pos]);
            for (ElementPos pos = positions[3]; pos < seq_size; ++pos)
                append(sequence_tmp, elements[pos]);
            compute_global_cost(solution_tmp_);
            solution = solution_tmp_;
            break;

        } case Perturbations::AdjacentSwaps: {
            // TODO
            break;

        } case Perturbations::RuinAndRecreate: {

            std::discrete_distribution<Counter> d_ruin({
                    parameters_.ruin_random_weight,
                    parameters_.ruin_nearest_weight,
                    parameters_.ruin_adjacent_string_removal_weight,
                    });
            Counter x_ruin = d_ruin(generator);

            std::discrete_distribution<Counter> d_recreate({
                    parameters_.recreate_random_weight,
                    parameters_.recreate_best_weight,
                    });
            Counter x_recreate = d_recreate(generator);

            //std::cout << "RuinAndRecreate" << std::endl;
            //std::cout << to_string(solution.global_cost) << std::endl;
            ElementId n = sequencing_scheme_.number_of_elements();
            std::vector<SequenceId> sequences(n, -1);
            std::vector<Mode> modes(n, -1);
            std::vector<ElementPos> positions(n, -1);
            optimizationtools::IndexedSet elts(n);
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
                const auto& elements = solution.sequences[sequence_id].elements;
                ElementPos seq_size = (ElementPos)elements.size();
                for (ElementPos pos = 0; pos < seq_size; ++pos) {
                    ElementId j = elements[pos].element_id;
                    elts.add(j);
                    sequences[j] = sequence_id;
                    modes[j] = elements[pos].mode;
                    positions[j] = pos;
                }
            }

            Solution solution_cur_ = empty_solution();
            solution_cur_.modified_sequences = std::vector<bool>(m, false);

            ElementPos number_of_elements_to_remove = std::min(
                    parameters_.ruin_number_of_elements_removed,
                    elts.size());
            ElementPos number_of_elements_removed = 0;
            switch (x_ruin) {
            case 0: {  // Random.
                elts.shuffle_in(generator);
                while (number_of_elements_removed < number_of_elements_to_remove) {
                    ElementId j = *(elts.begin());
                    elts.remove(j);
                    solution_cur_.modified_sequences[sequences[j]] = true;
                    number_of_elements_removed++;
                }
                break;
            } case 1: {  // Nearest.
                std::uniform_int_distribution<ElementPos> d(0, elts.size() - 1);
                ElementPos pos = d(generator);
                ElementId j = *(elts.begin() + pos);
                elts.remove(j);
                solution_cur_.modified_sequences[sequences[j]] = true;
                number_of_elements_removed++;
                for (ElementPos p = 0;
                        number_of_elements_removed < number_of_elements_to_remove; ++p) {
                    ElementId j2 = (p % 2 == 0)?
                        sorted_successors_[j][p / 2]:
                        sorted_predecessors_[j][p / 2];
                    if (elts.contains(j2)) {
                        elts.remove(j2);
                        solution_cur_.modified_sequences[sequences[j2]] = true;
                        number_of_elements_removed++;
                    }
                }
                break;
            } case 2: {  // Adjacent string removal.
                // Compute average sequence size.
                double average_sequence_size = 0;
                SequenceId number_of_non_empty_sequences = 0;
                for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
                    ElementPos seq_size = solution.sequences[sequence_id].elements.size();
                    average_sequence_size += seq_size;
                    if (seq_size > 0)
                        number_of_non_empty_sequences++;
                }
                average_sequence_size /= number_of_non_empty_sequences;
                // Compute maximum string cardinality.
                double maximum_string_cardinality = std::min(
                        (double)parameters_.ruin_adjacent_string_removal_maximum_string_cardinality,
                        average_sequence_size);
                //std::cout << "maximum_string_cardinality "
                //    << maximum_string_cardinality
                //    << std::endl;
                // Compute maximum number of strings to remove.
                double maximum_number_of_strings_to_remove
                    = 4.0 * (double)parameters_.ruin_number_of_elements_removed
                    / (1 + maximum_string_cardinality) - 1;
                // Compute number of strings to remove.
                std::uniform_real_distribution<double> d_str(1, maximum_number_of_strings_to_remove + 1);
                SequenceId number_of_strings_to_remove = d_str(generator);
                //std::cout << "number_of_strings_to_remove "
                //    << number_of_strings_to_remove
                //    << " / " << maximum_number_of_strings_to_remove
                //    << std::endl;
                // Draw seed element.
                std::uniform_int_distribution<ElementPos> d(0, elts.size() - 1);
                ElementPos pos = d(generator);
                ElementId j = *(elts.begin() + pos);
                std::vector<int> ruined_sequences(m, 0);
                SequencePos number_of_sequences_ruined = 0;
                for (ElementPos p = 0;
                        number_of_sequences_ruined < number_of_strings_to_remove; ++p) {
                    ElementId j2 = (p % 2 == 0)?
                        sorted_successors_[j][p / 2]:
                        sorted_predecessors_[j][p / 2];
                    SequenceId i_j2 = sequences[j2];
                    if (elts.contains(j2) && ruined_sequences[i_j2] == 0) {
                        // Compute the cardinality of the string to remove.
                        ElementPos pos_j2 = positions[j2];
                        double seq_size = solution.sequences[i_j2].elements.size();
                        //std::cout << "sequence " << i_j2
                        //    << " size " << seq_size
                        //    << " j2 " << j2
                        //    << " pos_j2 " << pos_j2
                        //    << std::endl;
                        std::uniform_real_distribution<double> d_card(
                                1.0,
                                std::min(seq_size, maximum_string_cardinality));
                        ElementPos string_cardinality = d_card(generator);
                        //std::cout << "string_cardinality "
                        //    << string_cardinality
                        //    << " / " << maximum_string_cardinality
                        //    << std::endl;
                        // Choose between the 'string' procedure and the 'split
                        // string' procedure.
                        std::uniform_real_distribution<double> d01(0, 1);
                        if (string_cardinality == seq_size
                                || d01(generator) >= parameters_.ruin_adjacent_string_removal_split_rate) {
                            // String procedure.
                            // Draw start of string to remove.
                            std::uniform_int_distribution<ElementPos> d_start(
                                    std::max(
                                        (ElementPos)0,
                                        pos_j2 - string_cardinality + 1),
                                    std::min(
                                        (ElementPos)seq_size - string_cardinality,
                                        pos_j2));
                            ElementPos start = d_start(generator);
                            //std::cout << "start " << start << std::endl;
                            // Remove string.
                            for (ElementPos pos_j3 = start;
                                    pos_j3 < start + string_cardinality; ++pos_j3) {
                                ElementId j3 = solution.sequences[i_j2].elements[pos_j3].element_id;
                                elts.remove(j3);
                            }
                            if (elts.contains(j2)) {
                                throw std::runtime_error("");
                            }
                        } else {
                            // Split string procedure.
                            // Compute m.
                            ElementPos m = 1;
                            while (string_cardinality + m < seq_size
                                    && d01(generator) > parameters_.ruin_adjacent_string_removal_beta) {
                                m++;
                            }
                            string_cardinality += m;
                            //std::cout
                            //    << "m " << m
                            //    << " string_cardinality " << string_cardinality
                            //    << std::endl;
                            // Draw start of string to remove.
                            ElementPos d_start_min = std::max(
                                        (ElementPos)0,
                                        pos_j2 - string_cardinality + 1);
                            ElementPos d_start_max = std::min(
                                        (ElementPos)seq_size - string_cardinality,
                                        pos_j2);
                            //std::cout << "d_start_min " << d_start_min << " d_start_max " << d_start_max << std::endl;
                            std::uniform_int_distribution<ElementPos> d_start(
                                    d_start_min, d_start_max);
                            ElementPos start = d_start(generator);
                            //std::cout << "start " << start << std::endl;
                            // Draw start of the substring to remove.
                            std::uniform_int_distribution<ElementPos> d_sub(
                                    start, start + string_cardinality - m);
                            ElementPos start_sub = d_sub(generator);
                            //std::cout << "start_sub " << start << std::endl;
                            // Remove splitted string.
                            for (ElementPos pos_j3 = start;
                                    pos_j3 < start_sub; ++pos_j3) {
                                ElementId j3 = solution.sequences[i_j2].elements[pos_j3].element_id;
                                elts.remove(j3);
                            }
                            for (ElementPos pos_j3 = start_sub + m;
                                    pos_j3 < start + string_cardinality; ++pos_j3) {
                                ElementId j3 = solution.sequences[i_j2].elements[pos_j3].element_id;
                                elts.remove(j3);
                            }
                        }
                        // Update structures.
                        ruined_sequences[i_j2] = 1;
                        number_of_sequences_ruined++;
                        solution_cur_.modified_sequences[i_j2] = true;
                    }
                }
                break;
            } default: {
                break;
            }
            }

            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
                const auto& elements = solution.sequences[sequence_id].elements;
                SequencePos seq_size = elements.size();
                Sequence& sequence_cur = solution_cur_.sequences[sequence_id];
                for (ElementPos pos = 0; pos < seq_size; ++pos) {
                    ElementId j = elements[pos].element_id;
                    if (elts.contains(j))
                        append(sequence_cur, elements[pos]);
                }
            }

            // Add back removed elements.
            compute_global_cost(solution_cur_);
            compute_temporary_structures(solution_cur_);
            elts.shuffle_out(generator);
            switch (x_recreate) {
            case 0: {  // Random.
                for (auto it = elts.out_begin(); it != elts.out_end(); ++it) {
                    ElementId j = *it;
                    // Draw sequence.
                    std::uniform_int_distribution<SequenceId> d_i(0, m - 1);
                    SequenceId i = d_i(generator);
                    // Draw position in sequence.
                    const auto& elements = solution_cur_.sequences[i].elements;
                    SequencePos seq_size = elements.size();
                    std::uniform_int_distribution<SequenceId> d_pos(0, seq_size);
                    ElementPos pos = d_pos(generator);
                    // Draw mode.
                    std::uniform_int_distribution<SequenceId> d_mode(0, number_of_modes(j) - 1);
                    Mode mode = d_mode(generator);
                    // Apply move.
                    Move move;
                    move.type = Neighborhoods::Add;
                    move.sequence_id_1 = i;
                    move.element_id = j;
                    move.mode = mode;
                    move.pos_1 = pos;
                    apply_move(solution_cur_, move);
                    solution_cur_.modified_sequences[i] = true;
                }
                break;
            } case 1: {  // Best.
                for (auto it = elts.out_begin(); it != elts.out_end(); ++it) {
                    ElementId j = *it;
                    //std::cout << to_string(solution_cur_.global_cost) << std::endl;
                    ElementId j_prec_old = -1;
                    if (positions[j] > 0) {
                        const auto& elements = solution.sequences[sequences[j]].elements;
                        j_prec_old = elements[positions[j] - 1].element_id;
                    }
                    auto improving_moves = explore_add(
                            solution_cur_,
                            j,
                            sequences[j],
                            j_prec_old,
                            modes[j]);
                    //std::cout << "improving_moves.size() " << improving_moves.size() << std::endl;
                    if (!improving_moves.empty()) {
                        std::shuffle(
                                improving_moves.begin(),
                                improving_moves.end(),
                                generator);
                        Move move_best;
                        for (const Move move: improving_moves)
                            if (move_best.sequence_id_1 == -1 || dominates(
                                        move.global_cost,
                                        move_best.global_cost))
                                move_best = move;
                        apply_move(solution_cur_, move_best);
                        solution_cur_.modified_sequences[move_best.sequence_id_1] = true;
                        compute_temporary_structures(solution_cur_, move_best.sequence_id_1, move_best.sequence_id_2);
                    }
                }
                break;
            } default: {
                break;
            }
            }
            solution = solution_cur_;
            //std::cout << to_string(solution.global_cost) << std::endl;
            break;

        } case Perturbations::ForceAdd: {
            ElementId j = perturbation.force_add_j;
            // Draw sequence.
            std::uniform_int_distribution<SequenceId> d_i(0, m - 1);
            SequenceId i = d_i(generator);
            // Draw position in sequence.
            const auto& elements = solution.sequences[i].elements;
            SequencePos seq_size = elements.size();
            std::uniform_int_distribution<SequenceId> d_pos(0, seq_size);
            ElementPos pos = d_pos(generator);
            // Draw mode.
            std::uniform_int_distribution<SequenceId> d_mode(0, number_of_modes(j) - 1);
            Mode mode = d_mode(generator);
            // Apply move.
            Move move;
            move.type = Neighborhoods::Add;
            move.sequence_id_1 = i;
            move.element_id = j;
            move.mode = mode;
            move.pos_1 = pos;
            apply_move(solution, move);
            break;
        }
        }
    }

    /*
     * Local Search.
     */

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Perturbation& perturbation = Perturbation())
    {
        //std::cout << "local_search" << std::endl;
        //if (tabu.element_id != -1)
        //    std::cout << "j " << tabu.element_id << " j_prev " << tabu.j_prev
        //        << std::endl;
        //print(std::cout, solution);
        //std::cout << to_string(global_cost(solution)) << std::endl;

        auto local_search_begin = std::chrono::steady_clock::now();
        SequencePos m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();

        // Get neighborhoods.
        std::vector<std::tuple<Neighborhoods, ElementPos, ElementPos>> neighborhoods;
        for (ElementPos block_size = 1;
                block_size <= parameters_.shift_block_maximum_length;
                ++block_size) {
            neighborhoods.push_back({Neighborhoods::Shift, block_size, 0});
        }
        for (ElementPos block_size_1 = 1;
                block_size_1 <= parameters_.swap_block_maximum_length;
                ++block_size_1) {
            for (ElementPos block_size_2 = 1;
                    block_size_2 <= block_size_1;
                    ++block_size_2) {
                neighborhoods.push_back({
                        Neighborhoods::Swap, block_size_1, block_size_2});
            }
        }
        if (parameters_.reverse)
            neighborhoods.push_back({Neighborhoods::Reverse, 0, 0});
        for (ElementPos block_size = 2;
                block_size <= parameters_.shift_reverse_block_maximum_length;
                ++block_size) {
            neighborhoods.push_back({
                    Neighborhoods::ShiftReverse, block_size, 0});
        }

        if (parameters_.add_remove) {
            neighborhoods.push_back({Neighborhoods::Add, 0, 0});
            neighborhoods.push_back({Neighborhoods::Remove, 0, 0});
        }
        if (parameters_.replace)
            neighborhoods.push_back({Neighborhoods::Replace, 0, 0});

        if (m > 1) {
            if (parameters_.swap_tails)
                neighborhoods.push_back({Neighborhoods::SwapTails, 0, 0});
            if (parameters_.split)
                neighborhoods.push_back({Neighborhoods::Split, 0, 0});
            for (ElementPos block_size = 1;
                    block_size <= parameters_.inter_shift_block_maximum_length;
                    ++block_size) {
                neighborhoods.push_back({
                        Neighborhoods::InterShift, block_size, 0});
            }
            for (ElementPos block_size_1 = 1;
                    block_size_1 <= parameters_.inter_swap_block_maximum_length;
                    ++block_size_1) {
                for (ElementPos block_size_2 = 1;
                        block_size_2 <= block_size_1;
                        ++block_size_2) {
                    neighborhoods.push_back({
                            Neighborhoods::InterSwap,
                            block_size_1,
                            block_size_2});
                }
            }
            for (ElementPos block_size = 2;
                    block_size <= parameters_.inter_shift_reverse_block_maximum_length;
                    ++block_size) {
                neighborhoods.push_back({
                        Neighborhoods::InterShiftReverse, block_size, 0});
            }
            if (parameters_.inter_swap_star)
                neighborhoods.push_back({Neighborhoods::InterSwapStar, 0, 0});
        }

        if (maximum_number_of_modes_ > 1) {
            if (parameters_.shift_change_mode)
                neighborhoods.push_back({Neighborhoods::ShiftChangeMode, 0, 0});
            if (parameters_.mode_swap)
                neighborhoods.push_back({Neighborhoods::ModeSwap, 0, 0});
            if (parameters_.swap_with_modes)
                neighborhoods.push_back({Neighborhoods::SwapWithModes, 0, 0});
        }

        for (int a = 0; a < (int)neighborhoods_.size(); ++a) {
            for (int k1 = 0; k1 < (int)neighborhoods_[a].size(); ++k1) {
                for (int k2 = 0; k2 < (int)neighborhoods_[a][k1].size(); ++k2) {
                    std::string s = neighborhood2string((Neighborhoods)a, k1, k2);
                    Neighborhood& neighborhood = neighborhoods_[a][k1][k2];
                    if (solution.modified_sequences.empty()) {
                        std::fill(
                                neighborhood.modified_sequences.begin(),
                                neighborhood.modified_sequences.end(),
                                true);
                    } else {
                        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                            neighborhood.modified_sequences[sequence_id] = solution.modified_sequences[sequence_id];
                    }
                }
            }
        }

        compute_temporary_structures(solution);

        if (parameters_.inter_shift_block_maximum_length >= 1) {
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
                std::fill(
                        inter_shift_1_best_positions_[sequence_id].begin(),
                        inter_shift_1_best_positions_[sequence_id].end(),
                        std::make_pair(-2, GlobalCost()));
            }
        }

        if (parameters_.inter_swap_star) {
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
                for (ElementId element_id = 0; element_id < n; ++element_id) {
                    std::fill(
                            inter_swap_star_best_positions_[sequence_id][element_id].begin(),
                            inter_swap_star_best_positions_[sequence_id][element_id].end(),
                            -2);
                }
            }
        }

        Counter it = 0;
        for (;; ++it) {
            //std::cout << "it " << it
            //    << " c " << to_string(global_cost(solution))
            //    << std::endl;
            //print(std::cout, solution);

            GlobalCost gc_old = global_cost(solution);

            std::shuffle(neighborhoods.begin(), neighborhoods.end(), generator);

            bool improved = false;
            // Loop through neighborhoods.
            for (auto neighborhood_id: neighborhoods) {
                auto begin = std::chrono::steady_clock::now();
                Neighborhoods type = std::get<0>(neighborhood_id);
                ElementPos k1 = std::get<1>(neighborhood_id);
                ElementPos k2 = std::get<2>(neighborhood_id);
                Neighborhood& neighborhood = neighborhoods_[(int)type][k1][k2];

                //std::cout << neighborhood2string(type, k1, k2) << std::endl;
                //for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                //    std::cout << " " << neighborhood.modified_sequences[i];
                //std::cout << std::endl;

                // Remove moves which have changed from improving_moves.
                for (auto it = neighborhood.improving_moves.begin();
                        it != neighborhood.improving_moves.end();) {
                    if (neighborhood.modified_sequences[it->sequence_id_1]
                            || (it->sequence_id_2 != -1
                                && neighborhood.modified_sequences[it->sequence_id_2])) {
                        *it = neighborhood.improving_moves.back();
                        neighborhood.improving_moves.pop_back();
                    } else {
                        it++;
                    }
                }

                switch (std::get<0>(neighborhood_id)) {
                case Neighborhoods::Shift:
                    explore_shift(solution, k1);
                    break;
                case Neighborhoods::Swap:
                    explore_swap(solution, k1, k2);
                    break;
                case Neighborhoods::Reverse:
                    explore_reverse(solution);
                    break;
                case Neighborhoods::ShiftReverse:
                    explore_shift(solution, k1, true);
                    break;
                case Neighborhoods::Add:
                    explore_add(solution);
                    break;
                case Neighborhoods::Remove:
                    explore_remove(solution, perturbation);
                    break;
                case Neighborhoods::Replace:
                    explore_replace(solution, perturbation);
                    break;
                case Neighborhoods::ShiftChangeMode:
                    explore_shift_change_mode(solution);
                    break;
                case Neighborhoods::ModeSwap:
                    explore_mode_swap(solution);
                    break;
                case Neighborhoods::SwapWithModes:
                    explore_swap_with_modes(solution);
                    break;
                case Neighborhoods::SwapTails:
                    explore_swap_tails(solution);
                    break;
                case Neighborhoods::Split:
                    explore_split(solution);
                    break;
                case Neighborhoods::InterShift:
                    explore_inter_shift(solution, k1);
                    //if (!parameters_.linking_constraints && k1 == 1) {
                    //    explore_inter_shift_1(solution);
                    //} else {
                    //    explore_inter_shift(solution, k1);
                    //}
                    break;
                case Neighborhoods::InterSwap:
                    explore_inter_swap(solution, k1, k2);
                    break;
                case Neighborhoods::InterShiftReverse:
                    explore_inter_shift(solution, k1, true);
                    break;
                case Neighborhoods::InterSwapStar:
                    explore_inter_swap_star(solution);
                    break;
                }

                // Update neighborhood.
                std::fill(
                        neighborhood.modified_sequences.begin(),
                        neighborhood.modified_sequences.end(),
                        false);
                neighborhood.number_of_explorations++;

                if (!neighborhood.improving_moves.empty()) {
                    //std::cout << neighborhood2string(type, k1, k2) << std::endl;
                    improved = true;
                    neighborhood.number_of_successes++;
                    std::shuffle(
                            neighborhood.improving_moves.begin(),
                            neighborhood.improving_moves.end(),
                            generator);
                    Move move_best;
                    for (const Move move: neighborhood.improving_moves)
                        if (move_best.sequence_id_1 == -1 || dominates(
                                    move.global_cost,
                                    move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    // Check new current solution cost.
                    if (!equals(
                                diff(global_cost(solution), gc_old),
                                move_best.global_cost)) {
                        std::cout
                            << "k1 " << move_best.k1
                            << " k2 " << move_best.k2
                            << " i1 " << move_best.sequence_id_1
                            << " i2 " << move_best.sequence_id_2
                            << " pos_1 " << move_best.pos_1
                            << " pos_2 " << move_best.pos_2
                            << " pos_3 " << move_best.pos_3
                            << " pos_4 " << move_best.pos_4
                            << std::endl;
                        throw std::logic_error(
                                neighborhood2string(type, k1, k2)
                                + ", costs do not match:\n"
                                + "* Old cost: "
                                + to_string(gc_old) + "\n"
                                + "* Move cost: "
                                + to_string(move_best.global_cost) + "\n"
                                + "* New cost: "
                                + to_string(global_cost(solution)) + "\n");
                    }
                    compute_temporary_structures(solution, move_best.sequence_id_1, move_best.sequence_id_2);

                    // Update modified sequences.
                    for (int a = 0; a < (int)neighborhoods_.size(); ++a) {
                        for (int k1 = 0; k1 < (int)neighborhoods_[a].size(); ++k1) {
                            for (int k2 = 0; k2 < (int)neighborhoods_[a][k1].size(); ++k2) {
                                Neighborhood& neighborhood_2 = neighborhoods_[a][k1][k2];
                                neighborhood_2.modified_sequences[move_best.sequence_id_1] = true;
                                if (move_best.sequence_id_2 != -1)
                                    neighborhood_2.modified_sequences[move_best.sequence_id_2] = true;
                            }
                        }
                    }
                }
                auto end = std::chrono::steady_clock::now();
                auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
                neighborhood.time += time_span.count();

                if (improved)
                    break;
            }

            if (!improved)
                break;
        }

        //print(std::cout, solution);
        //std::cout << to_string(global_cost(solution)) << std::endl;
        auto local_search_end = std::chrono::steady_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(local_search_end - local_search_begin);
        local_search_time_ += time_span.count();
        number_of_local_search_calls_++;
        number_of_local_search_iterations_ += it;
    }

    /*
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution) const
    {
        for (SequenceId i = 0; i < number_of_sequences(); ++i) {
            os << "Sequence " << i << ":";
            for (SequenceElement se: solution.sequences[i].elements) {
                if (maximum_number_of_modes_ == 1) {
                    os << " " << se.element_id;
                } else {
                    os << " " << se.element_id << " " << se.mode;
                }
            }
            os << std::endl;
            os << "    Cost: " << to_string(sequencing_scheme_.global_cost(solution.sequences[i].data)) << std::endl;
        }
        os << "Total cost: " << to_string(solution.global_cost) << std::endl;
        return os;
    }

    inline void write(
            const Solution& solution,
            std::string certificate_path,
            ElementPos offset = 0) const
    {
        if (certificate_path.empty())
            return;
        std::ofstream file(certificate_path);
        if (!file.good()) {
            throw std::runtime_error(
                    "Unable to open file \"" + certificate_path + "\".");
        }

        SequenceId m = number_of_sequences();
        if (m == 1) {
            for (SequenceElement se: solution.sequences[0].elements) {
                if (maximum_number_of_modes_ == 1) {
                    file << se.element_id + offset << " ";
                } else {
                    file << se.element_id + offset << " " << se.mode << " ";
                }
            }
        } else {
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
                file << solution.sequences[sequence_id].elements.size() << std::endl;
                for (SequenceElement se: solution.sequences[sequence_id].elements) {
                    if (maximum_number_of_modes_ == 1) {
                        file << se.element_id + offset << " ";
                    } else {
                        file << se.element_id + offset << " " << se.mode << " ";
                    }
                }
                file << std::endl;
            }
        }
    }

    void print_parameters(
            optimizationtools::Info& info) const
    {
        SequencePos m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        info.os()
            << "Number of sequences:                           " << m << std::endl
            << "Number of elements:                            " << n << std::endl
            << "Maximum number of modes:                       " << maximum_number_of_modes_ << std::endl
            << "Neighborhoods" << std::endl
            << "    Shift" << std::endl
            << "        Block maximum length:                  " << parameters_.shift_block_maximum_length << std::endl
            << "    Swap" << std::endl
            << "        Block maximum length:                  " << parameters_.swap_block_maximum_length << std::endl
            << "    Reverse:                                   " << parameters_.reverse << std::endl
            << "    Shift-reverse" << std::endl
            << "        Block maximum length:                  " << parameters_.shift_reverse_block_maximum_length << std::endl
            << "    Add/Remove:                                " << parameters_.add_remove << std::endl
            << "    Replace:                                   " << parameters_.replace << std::endl;
        if (m > 1) {
            info.os()
                << "    Inter-shift" << std::endl
                << "        Block maximum length:                  " << parameters_.inter_shift_block_maximum_length << std::endl
                << "    Inter-swap" << std::endl
                << "        Block maximum length:                  " << parameters_.inter_swap_block_maximum_length << std::endl
                << "    Inter-two-opt:                             " << parameters_.swap_tails << std::endl
                << "    Inter-two-opt-reverse:                     " << parameters_.split << std::endl
                << "    Inter-shift-reverse" << std::endl
                << "        Block maximum length:                  " << parameters_.inter_shift_reverse_block_maximum_length << std::endl
                << "    Inter-swap-star:                           " << parameters_.inter_swap_star << std::endl;
        }
        if (maximum_number_of_modes_ > 1) {
            info.os()
                << "    Shfit-change-mode:                         " << parameters_.shift_change_mode << std::endl
                << "    Mode-swap:                                 " << parameters_.mode_swap << std::endl
                << "    Swap-with-modes:                           " << parameters_.swap_with_modes << std::endl;
        }
        info.os()
            << "Perturbations" << std::endl
            << "    Double-bridge" << std::endl
            << "        Number of perturbations:               " << parameters_.double_bridge_number_of_perturbations << std::endl
            << "    Ruin-and-recreate" << std::endl
            << "        Number of perturbations:               " << parameters_.ruin_and_recreate_number_of_perturbations << std::endl
            << "        Number of elements removed:            " << parameters_.ruin_number_of_elements_removed << std::endl
            << "        Ruin random weight:                    " << parameters_.ruin_random_weight << std::endl
            << "        Ruin nearest weight:                   " << parameters_.ruin_nearest_weight << std::endl
            << "        Ruin adjacent string removal weight:   " << parameters_.ruin_adjacent_string_removal_weight << std::endl
            << "            Maximum string cardinality:        " << parameters_.ruin_adjacent_string_removal_maximum_string_cardinality << std::endl
            << "            Split rate:                        " << parameters_.ruin_adjacent_string_removal_split_rate << std::endl
            << "            Beta:                              " << parameters_.ruin_adjacent_string_removal_beta << std::endl
            << "        Recreate random:                       " << parameters_.recreate_random_weight << std::endl
            << "        Recreate best:                         " << parameters_.recreate_best_weight << std::endl
            << "    Force-add:                                 " << parameters_.force_add << std::endl;
        info.os()
            << "Crossovers" << std::endl
            << "    OX:                                        " << parameters_.crossover_ox_weight << std::endl
            << "    SBOX:                                      " << parameters_.crossover_sbox_weight << std::endl
            << "    SJOX:                                      " << parameters_.crossover_sjox_weight << std::endl
            << "    SREX1:                                     " << parameters_.crossover_srex1_weight << std::endl
            << "    SREX2:                                     " << parameters_.crossover_srex2_weight << std::endl;
    }

    void print_statistics(
            optimizationtools::Info& info) const
    {
        info.os() << "General" << std::endl;
        info.os()
            << "    "
            << std::left << std::setw(28) << "Initial solution:"
            << number_of_initial_solution_calls_
            << " / " << initial_solution_time_ << "s"
            << " / " << initial_solution_time_ / number_of_initial_solution_calls_ << "s"
            << std::endl;
        info.os()
            << "    "
            << std::left << std::setw(28) << "Crossover time:"
            << number_of_crossover_calls_
            << " / " << crossover_time_ << "s"
            << " / " << crossover_time_ / number_of_crossover_calls_ << "s"
            << std::endl;
        info.os()
            << "    "
            << std::left << std::setw(28) << "Local search time:"
            << number_of_local_search_calls_
            << " / " << local_search_time_ << "s"
            << " / " << local_search_time_ / number_of_local_search_calls_ << "s"
            << std::endl;
        info.os()
            << "    "
            << std::left << std::setw(28) << "Local search iterations:"
            << number_of_local_search_iterations_
            << " / " << (double)number_of_local_search_iterations_ / number_of_local_search_calls_ 
            << std::endl;
        info.os()
            << "    "
            << std::left << std::setw(28) << "Temporary structures time:"
            << compute_temporary_structures_time_ << "s"
            << std::endl;
        info.os() << "Neighborhoods" << std::endl;
        for (int a = 0; a < (int)neighborhoods_.size(); ++a) {
            for (int k1 = 0; k1 < (int)neighborhoods_[a].size(); ++k1) {
                for (int k2 = 0; k2 < (int)neighborhoods_[a][k1].size(); ++k2) {
                    std::string s = neighborhood2string((Neighborhoods)a, k1, k2);
                    const Neighborhood& neighborhood = neighborhoods_[a][k1][k2];
                    if (neighborhood.number_of_explorations == 0)
                        continue;
                    info.os()
                        << "    "
                        << std::left << std::setw(28) << s
                        << neighborhood.number_of_explorations
                        << " / " << neighborhood.number_of_successes
                        << " / " << (double)neighborhood.number_of_successes / neighborhood.number_of_explorations * 100 << "%"
                        << " / " << neighborhood.time << "s"
                        << std::endl;
                    info.add_to_json(
                            "Algorithm",
                            s + "NumberOfExplorations",
                            neighborhood.number_of_explorations);
                    info.add_to_json(
                            "Algorithm",
                            s + "NumberOfSuccesses",
                            neighborhood.number_of_successes);

                }
            }
        }
    }

private:

    /*
     * SequencePos number_of_sequences()
     */

    template<typename, typename T>
    struct HasNumberOfSequencesMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasNumberOfSequencesMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().number_of_sequences(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    SequencePos number_of_sequences(
            std::false_type) const
    {
        return 1;
    }

    SequencePos number_of_sequences(
            std::true_type) const
    {
        return sequencing_scheme_.number_of_sequences();;
    }

    SequencePos number_of_sequences() const
    {
        return number_of_sequences(
                std::integral_constant<
                    bool,
                    HasNumberOfSequencesMethod<SequencingScheme, SequencePos()>::value>());
    }

    /*
     * Mode number_of_modes(ElementId)
     */

    template<typename, typename T>
    struct HasNumberOfModesMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasNumberOfModesMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().number_of_modes(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    SequencePos number_of_modes(
            ElementId,
            std::false_type) const
    {
        return 1;
    }

    SequencePos number_of_modes(
            ElementId j,
            std::true_type) const
    {
        return sequencing_scheme_.number_of_modes(j);;
    }

    Mode number_of_modes(ElementId j) const
    {
        return number_of_modes(
                j,
                std::integral_constant<
                    bool,
                    HasNumberOfModesMethod<SequencingScheme, Mode(ElementId)>::value>());
    }

    /*
     * void append(SequenceData&, ElementId)
     * void append(SequenceData&, ElementId, ModeId)
     */

    template<typename, typename T>
    struct HasAppendMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasAppendMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().append(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    void append(
            SequenceData& sequence_data,
            SequenceElement se,
            std::true_type) const
    {
        sequencing_scheme_.append(sequence_data, se.element_id, se.mode);;
    }

    void append(
            SequenceData& sequence_data,
            SequenceElement se,
            std::false_type) const
    {
        sequencing_scheme_.append(sequence_data, se.element_id);;
    }

    void append(
            SequenceData& sequence_data,
            SequenceElement se) const
    {
        return append(
                sequence_data,
                se,
                std::integral_constant<
                    bool,
                    HasAppendMethod<SequencingScheme, void(SequenceData&, ElementId, Mode)>::value>());
    }

    template<typename, typename T>
    struct HasSequenceDataInitMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasSequenceDataInitMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().sequence_data_init(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    SequenceData sequence_data_init(
            SequenceElement se,
            std::true_type) const
    {
        return sequencing_scheme_.sequence_data_init(se.element_id, se.mode);;
    }

    SequenceData sequence_data_init(
            SequenceElement se,
            std::false_type) const
    {
        return sequencing_scheme_.sequence_data_init(se.element_id);;
    }

    SequenceData sequence_data_init(
            SequenceElement se) const
    {
        return sequence_data_init(
                se,
                std::integral_constant<
                    bool,
                    HasSequenceDataInitMethod<SequencingScheme, SequenceData(ElementId, Mode)>::value>());
    }

    /*
     * SequenceData empty_sequence_data(SequenceId)
     */

    template<typename, typename T>
    struct HasEmptySequenceDataMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasEmptySequenceDataMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().empty_sequence_data(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    SequenceData empty_sequence_data(
            SequenceId,
            std::false_type) const
    {
        return SequenceData();
    }

    SequenceData empty_sequence_data(
            SequenceId i,
            std::true_type) const
    {
        return sequencing_scheme_.empty_sequence_data(i);;
    }

    SequenceData empty_sequence_data(SequenceId i) const
    {
        return empty_sequence_data(
                i,
                std::integral_constant<
                    bool,
                    HasEmptySequenceDataMethod<SequencingScheme, SequenceData(SequenceId)>::value>());
    }

    /*
     * diff(
     *         const GlobalCost&,
     *         const GlobalCost&)
     */

    template<typename, typename T>
    struct HasDiffMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasDiffMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().diff(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    GlobalCost diff(
            const GlobalCost& gc1,
            const GlobalCost& gc2,
            std::false_type) const
    {
        return gc1 - gc2;
    }

    GlobalCost diff(
            const GlobalCost& gc1,
            const GlobalCost& gc2,
            std::true_type) const
    {
        return sequencing_scheme_.diff(gc1, gc2);
    }

    GlobalCost diff(
            const GlobalCost& gc1,
            const GlobalCost& gc2) const
    {
        return diff(
                gc1,
                gc2,
                std::integral_constant<
                    bool,
                    HasDiffMethod<SequencingScheme, GlobalCost(
                        const GlobalCost&,
                        const GlobalCost&)>::value>());
    }

    /*
     * merge(
     *         const GlobalCost&,
     *         const GlobalCost&)
     */

    template<typename, typename T>
    struct HasMergeMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasMergeMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().merge(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    GlobalCost merge(
            const GlobalCost& gc1,
            const GlobalCost& gc2,
            std::false_type) const
    {
        return gc1 + gc2;
    }

    GlobalCost merge(
            const GlobalCost& gc1,
            const GlobalCost& gc2,
            std::true_type) const
    {
        return sequencing_scheme_.merge(gc1, gc2);
    }

    /**
     * Return the merged global cost of two global costs.
     *
     * If SequencingScheme does not implement a method
     * 'GlobalCost global_cost(const GlobalCost&, const GlobalCost&)',
     * then it returns the sum of the two global costs, i.e., the global cost
     * of two sequences is the sum of the global cost of each sequence. If it
     * is not the case, for example when the objective is the makespan,
     * SequencingScheme needs to implement a custom method.
     */
    GlobalCost merge(
            const GlobalCost& gc1,
            const GlobalCost& gc2) const
    {
        return merge(
                gc1,
                gc2,
                std::integral_constant<
                    bool,
                    HasMergeMethod<SequencingScheme, GlobalCost(
                        const GlobalCost&,
                        const GlobalCost&)>::value>());
    }

    /*
     * append(SequenceData&, bool, SequenceElement, ...)
     */

    bool append(
            SequenceData& sequence_data,
            bool,
            SequenceElement se,
            const GlobalCost& cutoff,
            const GlobalCost* gci,
            std::false_type) const
    {
        append(sequence_data, se);;
        // Check early termination.
        if (gci == nullptr) {
            if (!strictly_better(bound(sequence_data), cutoff))
                return false;
        } else {
            if (!strictly_better(merge(*gci, bound(sequence_data)), cutoff))
                return false;
        }
        return true;
    }

    bool append(
            SequenceData& sequence_data,
            bool empty,
            SequenceElement se,
            const GlobalCost& cutoff,
            const GlobalCost* gci,
            std::true_type) const
    {
        if (empty) {
            sequence_data = sequence_data_init(se);
        } else {
            sequencing_scheme_.concatenate(
                    sequence_data,
                    sequence_data_init(se));
        }
        // Check early termination.
        if (gci == nullptr) {
            if (!strictly_better(bound(sequence_data), cutoff))
                return false;
        } else {
            if (!strictly_better(merge(*gci, bound(sequence_data)), cutoff))
                return false;
        }
        return true;
    }

    bool append(
            SequenceData& sequence_data,
            bool empty,
            SequenceElement se,
            const GlobalCost& cutoff,
            const GlobalCost* gci) const
    {
        return append(
                sequence_data,
                empty,
                se,
                cutoff,
                gci,
                std::integral_constant<
                    bool,
                    HasConcatenateMethod<SequencingScheme, bool(SequenceData&, const SequenceData&)>::value>());
    }

    void append(
            SequenceData& sequence_data,
            bool,
            SequenceElement se,
            std::false_type) const
    {
        append(sequence_data, se);;
    }

    void append(
            SequenceData& sequence_data,
            bool empty,
            SequenceElement se,
            std::true_type) const
    {
        if (empty) {
            sequence_data = sequence_data_init(se);
        } else {
            sequencing_scheme_.concatenate(
                    sequence_data,
                    sequence_data_init(se));
        }
    }

    void append(
            SequenceData& sequence_data,
            bool empty,
            SequenceElement se) const
    {
        return append(
                sequence_data,
                empty,
                se,
                std::integral_constant<
                    bool,
                    HasConcatenateMethod<SequencingScheme, bool(SequenceData&, const SequenceData&)>::value>());
    }

    void append(
            Sequence& sequence,
            const SequenceElement& se)
    {
        append(sequence.data, sequence.elements.empty(), se);
        sequence.elements.push_back(se);
    }

    /*
     * concatenate(
     *         SequenceData&,
     *         bool,
     *         const Sequence&,
     *         ElementPos,
     *         ElementPos,
     *         bool,
     *         ...)
     */

    template<typename, typename T>
    struct HasConcatenateMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasConcatenateMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().concatenate(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    inline const SequenceData& get_sequence_data(
            const Sequence& sequence,
            ElementPos pos_1,
            ElementPos pos_2,
            bool reverse) const
    {
        return (!reverse)?
            sequence_datas_cur_2_[sequence.sequence_id][pos_1][pos_2]:
            sequence_datas_cur_2_[sequence.sequence_id][pos_2][pos_1];
    }

    inline void concatenate(
            SequenceData& sequence_data,
            bool,
            const Sequence& sequence,
            ElementPos pos_1,
            ElementPos pos_2,
            bool reverse,
            std::false_type) const
    {
        if (!reverse) {
            auto it = sequence.elements.begin() + pos_1;
            const auto it_end = sequence.elements.begin() + pos_2;
            for (;; ++it) {
                // Add next element to sequence_data.
                append(sequence_data, *it);
                // End condition.
                if (it == it_end)
                    break;
            }
        } else {
            auto it = sequence.elements.begin() + pos_2;
            const auto it_end = sequence.elements.begin() + pos_1;
            for (;; --it) {
                // Add next element to sequence_data.
                append(sequence_data, *it);
                // End condition.
                if (it == it_end)
                    break;
            }
        }
    }

    inline bool concatenate(
            SequenceData& sequence_data,
            bool,
            const Sequence& sequence,
            ElementPos pos_1,
            ElementPos pos_2,
            bool reverse,
            const GlobalCost& cutoff,
            const GlobalCost* gci,
            std::false_type) const
    {
        if (!reverse) {
            auto it = sequence.elements.begin() + pos_1;
            const auto it_end = sequence.elements.begin() + pos_2;
            for (;; ++it) {
                // Add next element to sequence_data.
                append(sequence_data, *it);
                //std::cout << "j " << it->j
                //    << " gc " << to_string(global_cost(sequence_data))
                //    << " bnd " << to_string(bound(sequence_data))
                //    << " cutoff " << to_string(cutoff)
                //    << std::endl;
                // Check early termination.
                if (gci == nullptr) {
                    if (!strictly_better(bound(sequence_data), cutoff))
                        return false;
                } else {
                    if (!strictly_better(merge(*gci, bound(sequence_data)), cutoff))
                        return false;
                }
                // End condition.
                if (it == it_end)
                    break;
            }
        } else {
            auto it = sequence.elements.begin() + pos_2;
            const auto it_end = sequence.elements.begin() + pos_1;
            for (;; --it) {
                // Add next element to sequence_data.
                append(sequence_data, *it);
                //std::cout << "j " << it->j
                //    << " gc " << to_string(global_cost(sequence_data))
                //    << " bnd " << to_string(bound(sequence_data))
                //    << std::endl;
                // Check early termination.
                if (gci == nullptr) {
                    if (!strictly_better(bound(sequence_data), cutoff))
                        return false;
                } else {
                    if (!strictly_better(merge(*gci, bound(sequence_data)), cutoff))
                        return false;
                }
                // End condition.
                if (it == it_end)
                    break;
            }
        }
        return true;
    }

    inline void concatenate(
            SequenceData& sequence_data,
            bool empty,
            const Sequence& sequence,
            ElementPos pos_1,
            ElementPos pos_2,
            bool reverse,
            std::true_type) const
    {
        if (empty) {
            sequence_data = get_sequence_data(sequence, pos_1, pos_2, reverse);
        } else {
            sequencing_scheme_.concatenate(
                    sequence_data,
                    get_sequence_data(sequence, pos_1, pos_2, reverse));
        }
    }

    inline bool concatenate(
            SequenceData& sequence_data,
            bool empty,
            const Sequence& sequence,
            ElementPos pos_1,
            ElementPos pos_2,
            bool reverse,
            const GlobalCost& cutoff,
            const GlobalCost* gci,
            std::true_type) const
    {
        if (empty) {
            sequence_data = get_sequence_data(sequence, pos_1, pos_2, reverse);
        } else {
            sequencing_scheme_.concatenate(
                    sequence_data,
                    get_sequence_data(sequence, pos_1, pos_2, reverse));
        }
        // Check early termination.
        if (gci == nullptr) {
            if (!strictly_better(bound(sequence_data), cutoff))
                return false;
        } else {
            if (!strictly_better(merge(*gci, bound(sequence_data)), cutoff))
                return false;
        }
        return true;
    }

    inline void concatenate(
            SequenceData& sequence_data,
            bool empty,
            const Sequence& sequence,
            ElementPos pos_1,
            ElementPos pos_2,
            bool reverse) const
    {
        if (pos_1 > pos_2)
            return;
        concatenate(
                sequence_data,
                empty,
                sequence,
                pos_1,
                pos_2,
                reverse,
                std::integral_constant<
                    bool,
                    HasConcatenateMethod<SequencingScheme, bool(SequenceData&, const SequenceData&)>::value>());
    }

    inline bool concatenate(
            SequenceData& sequence_data,
            bool empty,
            const Sequence& sequence,
            ElementPos pos_1,
            ElementPos pos_2,
            bool reverse,
            const GlobalCost& cutoff,
            const GlobalCost* gci) const
    {
        if (pos_1 > pos_2)
            return true;
        return concatenate(
                sequence_data,
                empty,
                sequence,
                pos_1,
                pos_2,
                reverse,
                cutoff,
                gci,
                std::integral_constant<
                    bool,
                    HasConcatenateMethod<SequencingScheme, bool(SequenceData&, const SequenceData&)>::value>());
    }

    /*
     * void compute_sorted_neighbors()
     */

    template<typename, typename T>
    struct HasDistanceMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasDistanceMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().distance(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    void compute_sorted_neighbors(
            std::false_type)
    {
        ElementPos n = sequencing_scheme_.number_of_elements();
        neighbors_ = std::vector<ElementId>(n);
        std::iota(neighbors_.begin(), neighbors_.end(), 0);
    }

    void compute_sorted_neighbors(
            std::true_type)
    {
        ElementPos n = sequencing_scheme_.number_of_elements();
        sorted_predecessors_ = std::vector<std::vector<ElementId>>(n);
        sorted_successors_ = std::vector<std::vector<ElementId>>(n);
        std::vector<ElementId> neighbors(n);
        std::iota(neighbors.begin(), neighbors.end(), 0);
        for (ElementId element_id = 0; element_id < n; ++element_id) {
            // Successors.
            sorted_successors_[element_id] = neighbors;
            std::sort(
                    sorted_successors_[element_id].begin(),
                    sorted_successors_[element_id].end(),
                    [this, element_id](ElementId j1, ElementId j2) -> bool
                    {
                        return sequencing_scheme_.distance(element_id, j1)
                            < sequencing_scheme_.distance(element_id, j2);
                    });
            // Predecessors.
            sorted_predecessors_[element_id] = neighbors;
            std::sort(
                    sorted_predecessors_[element_id].begin(),
                    sorted_predecessors_[element_id].end(),
                    [this, element_id](ElementId j1, ElementId j2) -> bool
                    {
                        return sequencing_scheme_.distance(j1, element_id)
                            < sequencing_scheme_.distance(j2, element_id);
                    });

            //std::cout << "Closest neighbors of " << j << ":";
            //for (ElementId j2: sorted_neighbors_[j])
            //    std::cout << " " << j2 << "," << sequencing_scheme_.distance(j, j2);
            //std::cout << std::endl;
            //std::cout << "Other neighbors of " << j << ":";
            //for (auto it = neighbors.begin() + 100; it != neighbors.end(); ++it)
            //    std::cout << " " << *it << "," << sequencing_scheme_.distance(j, *it);
        }
    }

    void compute_sorted_neighbors()
    {
        ElementPos n = sequencing_scheme_.number_of_elements();
        if (n <= 25) {
            compute_sorted_neighbors(std::false_type{});
        } else {
            compute_sorted_neighbors(
                    std::integral_constant<
                    bool,
                    HasDistanceMethod<SequencingScheme, double(ElementId, ElementId)>::value>());
        }
    }

    inline void compute_sequence_datas_cur(
            const Solution&,
            SequenceId,
            SequenceId,
            std::false_type)
    {
    }

    inline void compute_sequence_datas_cur(
            const Solution& solution,
            SequenceId i1,
            SequenceId i2,
            std::true_type)
    {
        SequencePos m = number_of_sequences();
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (i1 != -1 && sequence_id != i1 && sequence_id != i2)
                continue;
            const auto& elements = solution.sequences[sequence_id].elements;
            ElementPos seq_size = (ElementPos)elements.size();
            while ((ElementPos)sequence_datas_cur_2_[sequence_id].size()  < seq_size)
                sequence_datas_cur_2_[sequence_id].push_back({});
            for (ElementPos pos_1 = 0; pos_1 < seq_size; ++pos_1) {
                while ((ElementPos)sequence_datas_cur_2_[sequence_id][pos_1].size() < seq_size)
                    sequence_datas_cur_2_[sequence_id][pos_1].push_back({});
                sequence_datas_cur_2_[sequence_id][pos_1][pos_1] = sequence_data_init(elements[pos_1]);
            }
            for (ElementPos size = 2; size <= seq_size; ++size) {
                for (ElementPos pos_1 = 0; pos_1 + size - 1 < seq_size; ++pos_1) {
                    ElementPos pos_2 = pos_1 + size - 1;
                    sequence_datas_cur_2_[sequence_id][pos_1][pos_2]
                        = sequence_datas_cur_2_[sequence_id][pos_1][pos_2 - 1];
                    sequencing_scheme_.concatenate(
                            sequence_datas_cur_2_[sequence_id][pos_1][pos_2],
                            sequence_data_init(elements[pos_2]));
                    sequence_datas_cur_2_[sequence_id][pos_2][pos_1]
                        = sequence_datas_cur_2_[sequence_id][pos_2][pos_1 + 1];
                    sequencing_scheme_.concatenate(
                            sequence_datas_cur_2_[sequence_id][pos_2][pos_1],
                            sequence_data_init(elements[pos_1]));
                }
            }
        }
    }

    /**
     * Compute temporary structures.
     *
     * If 'i1 == -1', the temporary structures for all sequences are computed.
     * If 'i1 != -1', only the temporary structures of 'i1' and 'i2' are
     * updated.
     */
    inline void compute_temporary_structures(
            const Solution& solution,
            SequenceId i1 = -1,
            SequenceId i2 = -1)
    {
        auto begin = std::chrono::steady_clock::now();
        SequencePos m = number_of_sequences();
        // Update global_costs_cur_.
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
            global_costs_cur_[sequence_id] = sequencing_scheme_.global_cost(solution.sequences[sequence_id].data);
        // Update elements_cur_.
        std::fill(elements_cur_.begin(), elements_cur_.end(), SolutionElement());
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            const auto& elements = solution.sequences[sequence_id].elements;
            ElementPos seq_size = elements.size();
            for (ElementPos pos = 0; pos < seq_size; ++pos) {
                ElementId j = elements[pos].element_id;
                elements_cur_[j].sequence_id = sequence_id;
                elements_cur_[j].pos = pos;
                elements_cur_[j].mode = elements[pos].mode;
            }
        }
        // Update sequence_datas_cur_1_.
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            const auto& elements = solution.sequences[sequence_id].elements;
            ElementPos seq_size = elements.size();
            while ((ElementPos)sequence_datas_cur_1_[sequence_id].size() <= seq_size)
                sequence_datas_cur_1_[sequence_id].push_back({});
            sequence_datas_cur_1_[sequence_id][0] = empty_sequence_data(sequence_id);
            for (ElementPos pos = 0; pos < seq_size; ++pos) {
                sequence_datas_cur_1_[sequence_id][pos + 1] = sequence_datas_cur_1_[sequence_id][pos];
                append(
                        sequence_datas_cur_1_[sequence_id][pos + 1], (pos == 0),
                        elements[pos]);
            }
        }
        // Update sequence_datas_cur_2_.
        compute_sequence_datas_cur(
                solution,
                i1,
                i2,
                std::integral_constant<
                    bool,
                    HasConcatenateMethod<SequencingScheme, bool(SequenceData&, const SequenceData&)>::value>());

        if (parameters_.linking_constraints) {
            // Update partial_global_costs_1_.
            if (m > 1) {
                std::vector<GlobalCost> front(m);
                front[0] = global_costs_cur_[0];
                for (SequenceId i = 1; i < m; ++i)
                    front[i] = merge(front[i - 1], global_costs_cur_[i]);

                std::vector<GlobalCost> back(m);
                back[0] = global_costs_cur_[m - 1];
                for (SequenceId i = m - 2; i >= 0; --i)
                    back[i] = merge(back[i + 1], global_costs_cur_[i]);

                partial_global_costs_cur_1_[0] = back[1];
                partial_global_costs_cur_1_[m - 1] = front[m - 2];
                for (SequenceId i = 1; i < m - 1; ++i) {
                    partial_global_costs_cur_1_[i]
                        = merge(front[i - 1], back[i + 1]);
                }

                // Update partial_global_costs_2_.
                if (m > 2) {

                    partial_global_costs_cur_2_[0][2] = global_costs_cur_[1];
                    for (SequenceId i2 = 3; i2 < m; ++i2) {
                        partial_global_costs_cur_2_[0][i2] = merge(
                                partial_global_costs_cur_2_[0][i2 - 1],
                                global_costs_cur_[i2 - 1]);
                    }
                    for (SequenceId i1 = 1; i1 < m - 1; ++i1) {
                        partial_global_costs_cur_2_[i1][i1 + 1] = front[i1 - 1];
                        for (SequenceId i2 = i1 + 2; i2 < m; ++i2) {
                            partial_global_costs_cur_2_[i1][i2] = merge(
                                    partial_global_costs_cur_2_[i1][i2 - 1],
                                    global_costs_cur_[i2 - 1]);
                        }
                    }

                    // partial_global_costs_cur_2_[i1][m - 1] already have the
                    // right values.
                    for (SequenceId i1 = 0; i1 < m; ++i1) {
                        for (SequenceId i2 = i1 + 1; i2 < m - 1; ++i2) {
                            if (i1 == 0 && i2 == 1) {
                                partial_global_costs_cur_2_[i1][i2] = back[2];
                            } else {
                                partial_global_costs_cur_2_[i1][i2] = merge(
                                        partial_global_costs_cur_2_[i1][i2],
                                        back[i2 + 1]);
                                partial_global_costs_cur_2_[i2][i1]
                                    = partial_global_costs_cur_2_[i1][i2];
                            }
                        }
                    }
                }
            }
        }

        auto end = std::chrono::steady_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
        compute_temporary_structures_time_ += time_span.count();
    }

    /**
     * Compute and return the global cost of a solution.
     */
    inline void compute_global_cost(
            Solution& solution) const
    {
        solution.global_cost = sequencing_scheme_.global_cost(solution.sequences[0].data);
        for (SequenceId i = 1; i < number_of_sequences(); ++i) {
            solution.global_cost = merge(
                    solution.global_cost,
                    sequencing_scheme_.global_cost(solution.sequences[i].data));
        }
    }

    /*
     * Explore neighborhoods.
     */

    inline void explore_shift(
            const Solution& solution,
            ElementPos block_size,
            bool reverse = false)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = (!reverse)?
            neighborhoods_[int(Neighborhoods::Shift)][block_size][0]:
            neighborhoods_[int(Neighborhoods::ShiftReverse)][block_size][0];

        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!neighborhood.modified_sequences[sequence_id])
                continue;
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            if (!(parameters_.linking_constraints && m > 1))
                gc = global_costs_cur_[sequence_id];

            for (ElementPos block_pos = 0; block_pos <= seq_size - block_size; ++block_pos) {
                // Loop through all new positions.
                ElementPos pos_min = std::max(
                        (ElementPos)0,
                        block_pos - parameters_.shift_maximum_distance);
                ElementPos pos_max = std::min(
                        seq_size - block_size,
                        block_pos + parameters_.shift_maximum_distance);
                for (ElementPos pos_new = pos_min; pos_new <= pos_max; ++pos_new) {
                    if (pos_new == block_pos)
                        continue;

                    const GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                        &partial_global_costs_cur_1_[sequence_id]: nullptr;
                    GlobalCost gc_tmp;

                    if (block_pos > pos_new) {
                        SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos_new];
                        bool ok = concatenate(
                                sequence_data, (pos_new == 0),
                                sequence, block_pos, block_pos + block_size - 1, reverse,
                                gc, gcm);
                        if (!ok)
                            continue;
                        ok = concatenate(
                                sequence_data, false,
                                sequence, pos_new, block_pos - 1, false,
                                gc, gcm);
                        if (!ok)
                            continue;
                        ok = concatenate(
                                sequence_data, false,
                                sequence, block_pos + block_size, seq_size - 1, false,
                                gc, gcm);
                        if (!ok)
                            continue;
                        gc_tmp = (parameters_.linking_constraints && m > 1)?
                            merge(
                                    sequencing_scheme_.global_cost(sequence_data),
                                    partial_global_costs_cur_1_[sequence_id]):
                            sequencing_scheme_.global_cost(sequence_data);
                    } else {
                        SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][block_pos];
                        bool ok = concatenate(
                                sequence_data, (block_pos == 0),
                                sequence, block_pos + block_size, pos_new + block_size - 1, false,
                                gc, gcm);
                        if (!ok)
                            continue;
                        ok = concatenate(
                                sequence_data, false,
                                sequence, block_pos, block_pos + block_size - 1, reverse,
                                gc, gcm);
                        if (!ok)
                            continue;
                        ok = concatenate(
                                sequence_data, false,
                                sequence, pos_new + block_size, seq_size - 1, false,
                                gc, gcm);
                        if (!ok)
                            continue;
                        gc_tmp = (parameters_.linking_constraints && m > 1)?
                            merge(
                                    sequencing_scheme_.global_cost(sequence_data),
                                    partial_global_costs_cur_1_[sequence_id]):
                            sequencing_scheme_.global_cost(sequence_data);
                    }
                    if (strictly_better(gc_tmp, gc)) {
                        Move move;
                        move.type = (!reverse)?
                            Neighborhoods::Shift:
                            Neighborhoods::ShiftReverse;
                        move.k1 = block_size;
                        move.sequence_id_1 = sequence_id;
                        move.pos_1 = block_pos;
                        move.pos_2 = pos_new;
                        move.global_cost = diff(gc_tmp, gc);
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_swap(
            const Solution& solution,
            Counter block_size_1,
            Counter block_size_2)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Swap)][block_size_1][block_size_2];
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!neighborhood.modified_sequences[sequence_id])
                continue;
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            if (!(parameters_.linking_constraints && m > 1))
                gc = global_costs_cur_[sequence_id];

            // Loop through all pairs.
            for (ElementPos pos = 0; pos < seq_size; ++pos) {

                // block 1 is at [pos, pos + block_size_1[.
                ElementPos pos_1 = pos;
                ElementPos pos_2_max = std::min(
                        seq_size - block_size_2,
                        pos + block_size_1 + parameters_.swap_maximum_distance);
                for (ElementPos pos_2 = pos + block_size_1; pos_2 <= pos_2_max; ++pos_2) {

                    GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                        &partial_global_costs_cur_1_[sequence_id]: nullptr;
                    SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos];
                    bool ok = concatenate(
                            sequence_data, (pos == 0),
                            sequence, pos_2, pos_2 + block_size_2 - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data, false,
                            sequence, pos_1 + block_size_1, pos_2 - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data, false,
                            sequence, pos_1, pos_1 + block_size_1 - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data, false,
                            sequence, pos_2 + block_size_2, seq_size - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp = (parameters_.linking_constraints && m > 1)?
                        merge(
                                sequencing_scheme_.global_cost(sequence_data),
                                partial_global_costs_cur_1_[sequence_id]):
                        sequencing_scheme_.global_cost(sequence_data);

                    if (strictly_better(gc_tmp, gc)) {
                        //std::cout << to_string(gc_tmp)
                        //    << " " << to_string(gc)
                        //    << std::endl;
                        Move move;
                        move.type = Neighborhoods::Swap;
                        move.k1 = block_size_1;
                        move.k2 = block_size_2;
                        move.sequence_id_1 = sequence_id;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = diff(gc_tmp, gc);
                        neighborhood.improving_moves.push_back(move);
                    }
                }

                if (block_size_1 == block_size_2)
                    continue;

                // block 2 is at [pos, pos + block_size_2[.
                ElementPos pos_2 = pos;
                ElementPos pos_1_max = std::min(
                        seq_size - block_size_1,
                        pos + block_size_2 + parameters_.swap_maximum_distance);
                for (ElementPos pos_1 = pos + block_size_2; pos_1 <= pos_1_max; ++pos_1) {

                    GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                        &partial_global_costs_cur_1_[sequence_id]: nullptr;
                    SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos];
                    bool ok = concatenate(
                            sequence_data, (pos == 0),
                            sequence, pos_1, pos_1 + block_size_1 - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data, false,
                            sequence, pos_2 + block_size_2, pos_1 - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data, false,
                            sequence, pos_2, pos_2 + block_size_2 - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data, false,
                            sequence, pos_1 + block_size_1, seq_size - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp = (parameters_.linking_constraints && m > 1)?
                        merge(
                                sequencing_scheme_.global_cost(sequence_data),
                                partial_global_costs_cur_1_[sequence_id]):
                        sequencing_scheme_.global_cost(sequence_data);

                    if (strictly_better(gc_tmp, gc)) {
                        Move move;
                        move.type = Neighborhoods::Swap;
                        move.k1 = block_size_1;
                        move.k2 = block_size_2;
                        move.sequence_id_1 = sequence_id;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = diff(gc_tmp, gc);
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_reverse(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Reverse)][0][0];
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!neighborhood.modified_sequences[sequence_id])
                continue;
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            if (!(parameters_.linking_constraints && m > 1))
                gc = global_costs_cur_[sequence_id];

            // Loop through all pairs.
            for (ElementPos pos_1 = 0; pos_1 < seq_size; ++pos_1) {
                ElementPos pos_max = std::min(
                        seq_size,
                        pos_1 + parameters_.reverse_maximum_length);
                for (ElementPos pos_2 = pos_1 + 2; pos_2 < pos_max; ++pos_2) {

                    GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                        &partial_global_costs_cur_1_[sequence_id]: nullptr;
                    SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos_1];
                    bool ok = concatenate(
                            sequence_data, (pos_1 == 0),
                            sequence, pos_1, pos_2, true,
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data, false,
                            sequence, pos_2 + 1, seq_size - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp = (parameters_.linking_constraints && m > 1)?
                        merge(
                                sequencing_scheme_.global_cost(sequence_data),
                                partial_global_costs_cur_1_[sequence_id]):
                        sequencing_scheme_.global_cost(sequence_data);

                    if (strictly_better(gc_tmp, gc)) {
                        Move move;
                        move.type = Neighborhoods::Reverse;
                        move.sequence_id_1 = sequence_id;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = diff(gc_tmp, gc);
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_add(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Add)][0][0];

        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!neighborhood.modified_sequences[sequence_id])
                continue;
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            if (!(parameters_.linking_constraints && m > 1))
                gc = global_costs_cur_[sequence_id];

            // Loop through all new positions.
            for (ElementPos pos = 0; pos <= seq_size; ++pos) {

                for (ElementId element_id = 0; element_id < n; ++element_id) {
                    if (elements_cur_[element_id].mode != -1)
                        continue;

                    for (Mode mode = 0; mode < number_of_modes(element_id); ++mode) {

                        GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                            &partial_global_costs_cur_1_[sequence_id]: nullptr;
                        SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos];
                        bool ok = append(
                                sequence_data, (pos == 0),
                                {element_id, mode},
                                gc, gcm);
                        if (!ok)
                            continue;
                        ok = concatenate(
                                sequence_data, false,
                                sequence, pos, seq_size - 1, false,
                                gc, gcm);
                        if (!ok)
                            continue;
                        GlobalCost gc_tmp = (parameters_.linking_constraints && m > 1)?
                            merge(
                                    sequencing_scheme_.global_cost(sequence_data),
                                    partial_global_costs_cur_1_[sequence_id]):
                            sequencing_scheme_.global_cost(sequence_data);

                        if (strictly_better(gc_tmp, gc)) {
                            Move move;
                            move.type = Neighborhoods::Add;
                            move.sequence_id_1 = sequence_id;
                            move.element_id = element_id;
                            move.mode = mode;
                            move.pos_1 = pos;
                            move.global_cost = diff(gc_tmp, gc);
                            neighborhood.improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    inline void explore_remove(
            const Solution& solution,
            const Perturbation& perturbation)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Remove)][0][0];
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!neighborhood.modified_sequences[sequence_id])
                continue;
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                &partial_global_costs_cur_1_[sequence_id]: nullptr;

            // Loop through all new positions.
            for (ElementPos pos = 0; pos < seq_size; ++pos) {

                ElementId j = sequence.elements[pos].element_id;
                if (perturbation.type == Perturbations::ForceAdd
                        && j == perturbation.force_add_j)
                    continue;

                SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos];
                bool ok = concatenate(
                        sequence_data, (pos == 0),
                        sequence, pos + 1, seq_size - 1, false,
                        gc, gcm);
                if (!ok)
                    continue;
                GlobalCost gc_tmp = (parameters_.linking_constraints && m > 1)?
                    merge(
                            sequencing_scheme_.global_cost(sequence_data),
                            partial_global_costs_cur_1_[sequence_id]):
                    sequencing_scheme_.global_cost(sequence_data);

                if (strictly_better(gc_tmp, gc)) {
                    Move move;
                    move.type = Neighborhoods::Remove;
                    move.sequence_id_1 = sequence_id;
                    move.pos_1 = pos;
                    move.global_cost = diff(gc_tmp, gc);
                    neighborhood.improving_moves.push_back(move);
                }
            }
        }
    }

    inline void explore_replace(
            const Solution& solution,
            const Perturbation& perturbation)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Replace)][0][0];

        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!neighborhood.modified_sequences[sequence_id])
                continue;
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            if (!(parameters_.linking_constraints && m > 1))
                gc = global_costs_cur_[sequence_id];

            // Loop through all new positions.
            for (ElementPos pos = 0; pos < seq_size; ++pos) {

                if (perturbation.type == Perturbations::ForceAdd
                        && sequence.elements[pos].element_id == perturbation.force_add_j)
                    continue;

                for (ElementId element_id = 0; element_id < n; ++element_id) {
                    if (elements_cur_[element_id].mode != -1)
                        continue;

                    for (Mode mode = 0; mode < number_of_modes(element_id); ++mode) {

                        GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                            &partial_global_costs_cur_1_[sequence_id]: nullptr;
                        SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos];
                        bool ok = append(
                                sequence_data, (pos == 0),
                                {element_id, mode},
                                gc, gcm);
                        if (!ok)
                            continue;
                        ok = concatenate(
                                sequence_data, false,
                                sequence, pos + 1, seq_size - 1, false,
                                gc, gcm);
                        if (!ok)
                            continue;
                        GlobalCost gc_tmp = (parameters_.linking_constraints && m > 1)?
                            merge(
                                    sequencing_scheme_.global_cost(sequence_data),
                                    partial_global_costs_cur_1_[sequence_id]):
                            sequencing_scheme_.global_cost(sequence_data);

                        if (strictly_better(gc_tmp, gc)) {
                            Move move;
                            move.type = Neighborhoods::Replace;
                            move.sequence_id_1 = sequence_id;
                            move.element_id = element_id;
                            move.mode = mode;
                            move.pos_1 = pos;
                            move.global_cost = diff(gc_tmp, gc);
                            neighborhood.improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    inline void explore_inter_shift(
            const Solution& solution,
            ElementPos block_size,
            bool reverse = false)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = (!reverse)?
            neighborhoods_[int(Neighborhoods::InterShift)][block_size][0]:
            neighborhoods_[int(Neighborhoods::InterShiftReverse)][block_size][0];
        ElementPos granularity = 25 + number_of_local_search_calls_ / 256;

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (ElementPos pos_1 = 0; pos_1 < seq_1_size; ++pos_1) {
                if (pos_1 + block_size - 1 >= seq_1_size)
                    continue;
                ElementId j1 = sequence_1.elements[pos_1].element_id;
                const auto it_begin = (!sorted_predecessors_.empty() && granularity < n)?
                    sorted_predecessors_[j1].begin(): neighbors_.begin();
                const auto it_end = (!sorted_predecessors_.empty() && granularity < n)?
                    sorted_predecessors_[j1].begin() + granularity: neighbors_.end();

                SequenceData sequence_data_1 = sequence_datas_cur_1_[i1][pos_1];
                concatenate(
                        sequence_data_1, (pos_1 == 0),
                        sequence_1, pos_1 + block_size, seq_1_size - 1, false);
                GlobalCost gc_tmp_1 = sequencing_scheme_.global_cost(sequence_data_1);

                for (auto it = it_begin; it != it_end; ++it) {
                    ElementId j2 = *it;
                    SequenceId i2 = elements_cur_[j2].sequence_id;
                    ElementPos pos_2 = elements_cur_[j2].pos;
                    if (j2 == j1)
                        continue;
                    if (i2 == -1)
                        continue;
                    if (!neighborhood.modified_sequences[i1]
                            && !neighborhood.modified_sequences[i2])
                        continue;
                    if (i1 == i2)
                        continue;

                    const auto& sequence_2 = solution.sequences[i2];
                    SequencePos seq_2_size = sequence_2.elements.size();

                    if (!(parameters_.linking_constraints && m > 1))
                        gc = merge(
                                global_costs_cur_[i1],
                                global_costs_cur_[i2]);
                    GlobalCost gcm = (parameters_.linking_constraints && m > 2)?
                        merge(
                            gc_tmp_1,
                            partial_global_costs_cur_2_[i1][i2]):
                        gc_tmp_1;

                    SequenceData sequence_data_2 = sequence_datas_cur_1_[i2][pos_2 + 1];
                    bool ok = concatenate(
                            sequence_data_2, false,
                            sequence_1, pos_1, pos_1 + block_size - 1, reverse,
                            gc, &gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data_2, false,
                            sequence_2, pos_2 + 1, seq_2_size - 1, false,
                            gc, &gcm);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp = merge(
                            sequencing_scheme_.global_cost(sequence_data_2),
                            gcm);

                    if (strictly_better(gc_tmp, gc)) {
                        Move move;
                        move.type = (!reverse)?
                            Neighborhoods::InterShift:
                            Neighborhoods::InterShiftReverse;
                        move.k1 = block_size;
                        move.sequence_id_1 = i1;
                        move.sequence_id_2 = i2;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2 + 1;
                        move.global_cost = diff(gc_tmp, gc);
                        neighborhood.improving_moves.push_back(move);
                    }

                    // If j2 is the first element of its sequence, we also try
                    // to insert the block before j2.
                    if (pos_2 == 0) {
                        SequenceData sequence_data_2 = sequence_datas_cur_1_[i2][pos_2];
                        bool ok = concatenate(
                                sequence_data_2, true,
                                sequence_1, pos_1, pos_1 + block_size - 1, reverse,
                                gc, &gcm);
                        if (!ok)
                            continue;
                        ok = concatenate(
                                sequence_data_2, false,
                                sequence_2, pos_2, seq_2_size - 1, false,
                                gc, &gcm);
                        if (!ok)
                            continue;
                        GlobalCost gc_tmp = merge(
                                sequencing_scheme_.global_cost(sequence_data_2),
                                gcm);

                        if (strictly_better(gc_tmp, gc)) {
                            Move move;
                            move.type = (!reverse)?
                                Neighborhoods::InterShift:
                                Neighborhoods::InterShiftReverse;
                            move.k1 = block_size;
                            move.sequence_id_1 = i1;
                            move.sequence_id_2 = i2;
                            move.pos_1 = pos_1;
                            move.pos_2 = pos_2;
                            move.global_cost = diff(gc_tmp, gc);
                            neighborhood.improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    inline void explore_inter_shift_1(
            const Solution& solution)
    {
        if (parameters_.linking_constraints) {
            throw std::logic_error(
                    "Neighborhood 'Inter-shift-1'"
                    " requires no linking constraints.");
        }

        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::InterShift)][1][0];

        // Clean obsolete memory.
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!neighborhood.modified_sequences[sequence_id])
                continue;
            std::fill(
                    inter_shift_1_best_positions_[sequence_id].begin(),
                    inter_shift_1_best_positions_[sequence_id].end(),
                    std::make_pair(-2, GlobalCost()));
        }

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (ElementPos pos_1 = 0; pos_1 < seq_1_size; ++pos_1) {
                ElementId j = sequence_1.elements[pos_1].element_id;

                SequenceData sequence_data_1 = sequence_datas_cur_1_[i1][pos_1];
                concatenate(
                        sequence_data_1, (pos_1 == 0),
                        sequence_1, pos_1 + 1, seq_1_size - 1, false);
                GlobalCost gc_tmp_1 = sequencing_scheme_.global_cost(sequence_data_1);

                for (SequenceId i2 = 0; i2 < m; ++i2) {
                    if (i1 == i2)
                        continue;
                    if (!neighborhood.modified_sequences[i1]
                            && !neighborhood.modified_sequences[i2])
                        continue;

                    const auto& sequence_2 = solution.sequences[i2];
                    SequencePos seq_2_size = sequence_2.elements.size();

                    GlobalCost gc = merge(
                            global_costs_cur_[i1],
                            global_costs_cur_[i2]);

                    if (inter_shift_1_best_positions_[i2][j].first == -2) {
                        ElementPos pos_2_best = -1;
                        GlobalCost gc_best;
                        for (ElementPos pos_2 = 0; pos_2 <= seq_2_size; ++pos_2) {
                            SequenceData sequence_data_2 = sequence_datas_cur_1_[i2][pos_2];
                            append(
                                    sequence_data_2, (pos_2 == 0),
                                    sequence_1.elements[pos_1]);
                            concatenate(
                                    sequence_data_2, false,
                                    sequence_2, pos_2, seq_2_size - 1, false);
                            GlobalCost gc_tmp = sequencing_scheme_.global_cost(sequence_data_2);
                            if (pos_2_best == -1 || !(gc_tmp >= gc_best)) {
                                pos_2_best = pos_2;
                                gc_best = gc_tmp;
                            }
                        }
                        //std::cout << "pos_2_best " << pos_2_best
                        //    << " gc_best " << to_string(gc_best)
                        //    << std::endl;
                        inter_shift_1_best_positions_[i2][j] = {pos_2_best, gc_best};
                    }

                    GlobalCost gc_tmp = merge(
                            inter_shift_1_best_positions_[i2][j].second,
                            gc_tmp_1);
                    if (strictly_better(gc_tmp, gc)) {
                        Move move;
                        move.type = Neighborhoods::InterShift;
                        move.k1 = 1;
                        move.sequence_id_1 = i1;
                        move.sequence_id_2 = i2;
                        move.pos_1 = pos_1;
                        move.pos_2 = inter_shift_1_best_positions_[i2][j].first;
                        move.global_cost = diff(gc_tmp, gc);
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_swap_tails(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        GlobalCost gc = global_cost(solution);
        GlobalCost gc_tmp;
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::SwapTails)][0][0];
        ElementPos granularity = 25 + number_of_local_search_calls_ / 256;

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (ElementPos pos_1 = 0; pos_1 < seq_1_size; ++pos_1) {
                ElementId j1 = sequence_1.elements[pos_1].element_id;
                const auto it_begin = (!sorted_successors_.empty() && granularity < n)?
                    sorted_successors_[j1].begin(): neighbors_.begin();
                const auto it_end = (!sorted_successors_.empty() && granularity < n)?
                    sorted_successors_[j1].begin() + granularity: neighbors_.end();

                for (auto it = it_begin; it != it_end; ++it) {
                    ElementId j2 = *it;
                    SequenceId i2 = elements_cur_[j2].sequence_id;
                    ElementPos pos_2 = elements_cur_[j2].pos;
                    if (j2 == j1)
                        continue;
                    if (i2 == -1)
                        continue;
                    if (!neighborhood.modified_sequences[i1]
                            && !neighborhood.modified_sequences[i2])
                        continue;
                    if (i1 == i2)
                        continue;
                    //if (j1 > j2)
                    //    continue;

                    const auto& sequence_2 = solution.sequences[i2];
                    SequencePos seq_2_size = sequence_2.elements.size();

                    if (!(parameters_.linking_constraints && m > 2))
                        gc = merge(
                                global_costs_cur_[i1],
                                global_costs_cur_[i2]);
                    GlobalCost* gcm = (parameters_.linking_constraints && m > 2)?
                        &partial_global_costs_cur_2_[i1][i2]: nullptr;

                    SequenceData sequence_data_1 = sequence_datas_cur_1_[i1][pos_1 + 1];
                    bool ok = concatenate(
                            sequence_data_1, false,
                            sequence_2, pos_2 + 1, seq_2_size - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp_1 = (parameters_.linking_constraints && m > 2)?
                        merge(
                                sequencing_scheme_.global_cost(sequence_data_1),
                                partial_global_costs_cur_2_[i1][i2]):
                        sequencing_scheme_.global_cost(sequence_data_1);

                    SequenceData sequence_data_2 = sequence_datas_cur_1_[i2][pos_2 + 1];
                    ok = concatenate(
                            sequence_data_2, false,
                            sequence_1, pos_1 + 1, seq_1_size - 1, false,
                            gc, &gc_tmp_1);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp = merge(
                            sequencing_scheme_.global_cost(sequence_data_2),
                            gc_tmp_1);

                    if (strictly_better(gc_tmp, gc)) {
                        Move move;
                        move.type = Neighborhoods::SwapTails;
                        move.sequence_id_1 = i1;
                        move.sequence_id_2 = i2;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = diff(gc_tmp, gc);
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_split(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Split)][0][0];
        ElementPos granularity = 25 + number_of_local_search_calls_ / 256;

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (ElementPos pos_1 = 0; pos_1 < seq_1_size; ++pos_1) {
                ElementId j1 = sequence_1.elements[pos_1].element_id;
                const auto it_begin = (!sorted_successors_.empty() && granularity < n)?
                    sorted_successors_[j1].begin(): neighbors_.begin();
                const auto it_end = (!sorted_successors_.empty() && granularity < n)?
                    sorted_successors_[j1].begin() + granularity: neighbors_.end();

                for (auto it = it_begin; it != it_end; ++it) {
                    ElementId j2 = *it;
                    SequenceId i2 = elements_cur_[j2].sequence_id;
                    ElementPos pos_2 = elements_cur_[j2].pos;
                    if (j2 == j1)
                        continue;
                    if (i2 == -1)
                        continue;
                    if (!neighborhood.modified_sequences[i1]
                            && !neighborhood.modified_sequences[i2])
                        continue;
                    if (i1 == i2)
                        continue;
                    //if (j1 > j2)
                    //    continue;

                    const auto& sequence_2 = solution.sequences[i2];
                    SequencePos seq_2_size = sequence_2.elements.size();

                    if (!(parameters_.linking_constraints && m > 1))
                        gc = merge(
                                global_costs_cur_[i1],
                                global_costs_cur_[i2]);
                    GlobalCost* gcm = (parameters_.linking_constraints && m > 2)?
                        &partial_global_costs_cur_2_[i1][i2]: nullptr;

                    SequenceData sequence_data_1 = sequence_datas_cur_1_[i1][pos_1 + 1];
                    bool ok = concatenate(
                            sequence_data_1, false,
                            sequence_2, 0, pos_2, true,
                            gc, gcm);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp_1 = (parameters_.linking_constraints && m > 2)?
                        merge(
                                sequencing_scheme_.global_cost(sequence_data_1),
                                partial_global_costs_cur_2_[i1][i2]):
                        sequencing_scheme_.global_cost(sequence_data_1);

                    SequenceData sequence_data_2 = sequence_datas_cur_1_[i2][0];
                    ok = concatenate(
                            sequence_data_2, true,
                            sequence_1, pos_1 + 1, seq_1_size - 1, true,
                            gc, &gc_tmp_1);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data_2, (pos_1 + 1 > seq_1_size - 1),
                            sequence_2, pos_2 + 1, seq_2_size - 1, false,
                            gc, &gc_tmp_1);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp = merge(
                            sequencing_scheme_.global_cost(sequence_data_2),
                            gc_tmp_1);

                    if (strictly_better(gc_tmp, gc)) {
                        Move move;
                        move.type = Neighborhoods::Split;
                        move.sequence_id_1 = i1;
                        move.sequence_id_2 = i2;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = diff(gc_tmp, gc);
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_inter_swap(
            const Solution& solution,
            ElementPos block_size_1,
            ElementPos block_size_2)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::InterSwap)][block_size_1][block_size_2];
        ElementPos granularity = 25 + number_of_local_search_calls_ / 256;

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (ElementPos pos_1 = 0; pos_1 < seq_1_size; ++pos_1) {
                ElementId j1 = sequence_1.elements[pos_1].element_id;
                if (pos_1 + block_size_1 - 1 >= seq_1_size)
                    continue;
                const auto it_begin = (!sorted_predecessors_.empty() && granularity < n)?
                    sorted_predecessors_[j1].begin(): neighbors_.begin();
                const auto it_end = (!sorted_predecessors_.empty() && granularity < n)?
                    sorted_predecessors_[j1].begin() + granularity: neighbors_.end();

                for (auto it = it_begin; it != it_end; ++it) {
                    ElementId j2 = *it;
                    SequenceId i2 = elements_cur_[j2].sequence_id;
                    ElementPos pos_2 = elements_cur_[j2].pos;
                    if (j2 == j1)
                        continue;
                    if (i2 == -1)
                        continue;
                    if (!neighborhood.modified_sequences[i1]
                            && !neighborhood.modified_sequences[i2])
                        continue;
                    if (i1 == i2)
                        continue;
                    //if (block_size_1 == block_size_2 && j1 > j2)
                    //    continue;

                    const auto& sequence_2 = solution.sequences[i2];
                    SequencePos seq_2_size = sequence_2.elements.size();
                    if (pos_2 + block_size_2 >= seq_2_size)
                        continue;

                    if (!(parameters_.linking_constraints && m > 1))
                        gc = merge(
                                global_costs_cur_[i1],
                                global_costs_cur_[i2]);

                    SequenceData sequence_data_1 = sequence_datas_cur_1_[i1][pos_1];
                    SequenceData sequence_data_2 = sequence_datas_cur_1_[i2][pos_2 + 1];

                    GlobalCost gcm = (parameters_.linking_constraints && m > 2)?
                        merge(
                            bound(sequence_data_2),
                            partial_global_costs_cur_2_[i1][i2]):
                        bound(sequence_data_2);

                    bool ok = concatenate(
                            sequence_data_1, (pos_1 == 0),
                            sequence_2, pos_2 + 1, pos_2 + block_size_2, false,
                            gc, &gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data_1, false,
                            sequence_1, pos_1 + block_size_1, seq_1_size - 1, false,
                            gc, &gcm);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp_1 = (parameters_.linking_constraints && m > 2)?
                        merge(
                                sequencing_scheme_.global_cost(sequence_data_1),
                                partial_global_costs_cur_2_[i1][i2]):
                        sequencing_scheme_.global_cost(sequence_data_1);

                    ok = concatenate(
                            sequence_data_2, false,
                            sequence_1, pos_1, pos_1 + block_size_1 - 1, false,
                            gc, &gc_tmp_1);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data_2, false,
                            sequence_2, pos_2 + block_size_2 + 1, seq_2_size - 1, false,
                            gc, &gc_tmp_1);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp = merge(
                            sequencing_scheme_.global_cost(sequence_data_2),
                            gc_tmp_1);

                    if (ok && strictly_better(gc_tmp, gc)) {
                        Move move;
                        move.type = Neighborhoods::InterSwap;
                        move.k1 = block_size_1;
                        move.k2 = block_size_2;
                        move.sequence_id_1 = i1;
                        move.sequence_id_2 = i2;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2 + 1;
                        move.global_cost = diff(gc_tmp, gc);
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_inter_swap_star(
            const Solution& solution)
    {
        if (parameters_.linking_constraints) {
            throw std::logic_error(
                    "Neighborhood 'Inter-swap-star'"
                    " requires no linking constraints.");
        }

        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::InterSwapStar)][0][0];

        // Clean obsolete memory.
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!neighborhood.modified_sequences[sequence_id])
                continue;
            for (ElementId element_id = 0; element_id < n; ++element_id) {
                std::fill(
                        inter_swap_star_best_positions_[sequence_id][element_id].begin(),
                        inter_swap_star_best_positions_[sequence_id][element_id].end(),
                        -2);
            }
        }

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (SequenceId i2 = i1 + 1; i2 < m; ++i2) {
                if (!neighborhood.modified_sequences[i1]
                        && !neighborhood.modified_sequences[i2])
                    continue;

                GlobalCost gc = merge(
                        global_costs_cur_[i1],
                        global_costs_cur_[i2]);

                const auto& sequence_2 = solution.sequences[i2];
                SequencePos seq_2_size = sequence_2.elements.size();

                // Compute the 3 best position of the elements from the first
                // sequence into the second sequence.
                for (ElementPos pos_1 = 0; pos_1 <= seq_1_size - 1; ++pos_1) {
                    ElementId j = sequence_1.elements[pos_1].element_id;
                    if (inter_swap_star_best_positions_[i2][j][0] == -2) {
                        ElementPos pos_best_1 = -1;
                        ElementPos pos_best_2 = -1;
                        ElementPos pos_best_3 = -1;
                        GlobalCost gc_best_1;
                        GlobalCost gc_best_2;
                        GlobalCost gc_best_3;
                        for (ElementPos pos_2 = 0; pos_2 <= seq_2_size; ++pos_2) {
                            SequenceData sequence_data = sequence_datas_cur_1_[i2][pos_2];
                            bool ok = append(
                                    sequence_data, (pos_2 == 0),
                                    sequence_1.elements[pos_1],
                                    gc, nullptr);
                            if (!ok)
                                continue;
                            ok = concatenate(
                                    sequence_data, false,
                                    sequence_2, pos_2, seq_2_size - 1, false,
                                    gc, nullptr);
                            if (!ok)
                                continue;
                            GlobalCost gci2_tmp = sequencing_scheme_.global_cost(sequence_data);
                            if (pos_best_1 == -1 || strictly_better(gci2_tmp, gc_best_1)) {
                                pos_best_3 = pos_best_2;
                                pos_best_2 = pos_best_1;
                                pos_best_1 = pos_2;
                                gc_best_3 = gc_best_2;
                                gc_best_2 = gc_best_1;
                                gc_best_1 = gci2_tmp;
                            } else if (pos_best_2 == -1 || strictly_better(gci2_tmp, gc_best_2)) {
                                pos_best_3 = pos_best_2;
                                pos_best_2 = pos_2;
                                gc_best_3 = gc_best_2;
                                gc_best_2 = gci2_tmp;
                            } else if (pos_best_3 == -1 || strictly_better(gci2_tmp, gc_best_3)) {
                                pos_best_3 = pos_2;
                                gc_best_3 = gci2_tmp;
                            }
                        }
                        //std::cout << "pos_1 " << pos_1
                        //    << " pos_best_1 " << pos_best_1
                        //    << " gc_best_1 " << to_string(gc_best_1)
                        //    << " pos_best_2 " << pos_best_2
                        //    << " gc_best_2 " << to_string(gc_best_2)
                        //    << " pos_best_3 " << pos_best_3
                        //    << " gc_best_3 " << to_string(gc_best_3)
                        //    << std::endl;
                        inter_swap_star_best_positions_[i2][j]
                            = {pos_best_1, pos_best_2, pos_best_3};
                    }
                }

                // Calcul the 3 best position of the elements from the second
                // sequence into the first sequence.
                for (ElementPos pos_2 = 0; pos_2 <= seq_2_size - 1; ++pos_2) {
                    ElementId j = sequence_2.elements[pos_2].element_id;
                    if (inter_swap_star_best_positions_[i1][j][0] == -2) {
                        ElementPos pos_best_1 = -1;
                        ElementPos pos_best_2 = -1;
                        ElementPos pos_best_3 = -1;
                        GlobalCost gc_best_1;
                        GlobalCost gc_best_2;
                        GlobalCost gc_best_3;
                        for (ElementPos pos_1 = 0; pos_1 <= seq_1_size; ++pos_1) {
                            SequenceData sequence_data = sequence_datas_cur_1_[i1][pos_1];
                            bool ok = append(
                                    sequence_data, (pos_1 == 0),
                                    sequence_2.elements[pos_2],
                                    gc, nullptr);
                            if (!ok)
                                continue;
                            ok = concatenate(
                                    sequence_data, false,
                                    sequence_1, pos_1, seq_1_size - 1, false,
                                    gc, nullptr);
                            if (!ok)
                                continue;
                            GlobalCost gci1_tmp = sequencing_scheme_.global_cost(sequence_data);
                            if (pos_best_1 == -1 || strictly_better(gci1_tmp, gc_best_1)) {
                                pos_best_3 = pos_best_2;
                                pos_best_2 = pos_best_1;
                                pos_best_1 = pos_1;
                                gc_best_3 = gc_best_2;
                                gc_best_2 = gc_best_1;
                                gc_best_1 = gci1_tmp;
                            } else if (pos_best_2 == -1 || strictly_better(gci1_tmp, gc_best_2)) {
                                pos_best_3 = pos_best_2;
                                pos_best_2 = pos_1;
                                gc_best_3 = gc_best_2;
                                gc_best_2 = gci1_tmp;
                            } else if (pos_best_3 == -1 || strictly_better(gci1_tmp, gc_best_3)) {
                                pos_best_3 = pos_1;
                                gc_best_3 = gci1_tmp;
                            }
                        }
                        //std::cout << "pos_2 " << pos_2
                        //    << " pos_best_1 " << pos_best_1
                        //    << " gc_best_1 " << to_string(gc_best_1)
                        //    << " pos_best_2 " << pos_best_2
                        //    << " gc_best_2 " << to_string(gc_best_2)
                        //    << " pos_best_3 " << pos_best_3
                        //    << " gc_best_3 " << to_string(gc_best_3)
                        //    << std::endl;
                        inter_swap_star_best_positions_[i1][j]
                            = {pos_best_1, pos_best_2, pos_best_3};
                    }
                }

                for (ElementPos pos_1 = 0; pos_1 <= seq_1_size - 1; ++pos_1) {
                    ElementId j1 = sequence_1.elements[pos_1].element_id;

                    for (ElementPos pos_2 = 0; pos_2 <= seq_2_size - 1; ++pos_2) {
                        ElementId j2 = sequence_2.elements[pos_2].element_id;

                        //std::cout
                        //    << "i1 " << i1
                        //    << " i2 " << i2
                        //    << " seq_1_size " << seq_1_size
                        //    << " seq_2_size " << seq_2_size
                        //    << " pos_1 " << pos_1
                        //    << " pos_2 " << pos_2
                        //    << std::endl;

                        // Find pos_1_new.
                        ElementPos pos_1_new = -1;
                        GlobalCost gci2_tmp;
                        for (int a = 0; a < 4; ++a) {
                            ElementPos p = (a < 3)? inter_swap_star_best_positions_[i2][j1][a]: pos_2;
                            if (p != -1) {
                                // pos_2: removed positino
                                // p: added position
                                GlobalCost gci2;
                                if (pos_2 < p) {
                                    SequenceData sequence_data = sequence_datas_cur_1_[i2][pos_2];
                                    bool ok = concatenate(
                                            sequence_data, (pos_2 == 0),
                                            sequence_2, pos_2 + 1, p - 1, false,
                                            gc, nullptr);
                                    if (!ok)
                                        continue;
                                    ok = append(
                                            sequence_data, (pos_2 == 0 && pos_2 + 1 > p - 1),
                                            sequence_1.elements[pos_1],
                                            gc, nullptr);
                                    if (!ok)
                                        continue;
                                    ok = concatenate(
                                            sequence_data, false,
                                            sequence_2, p, seq_2_size - 1, false,
                                            gc, nullptr);
                                    if (!ok)
                                        continue;
                                    gci2 = sequencing_scheme_.global_cost(sequence_data);
                                } else {
                                    SequenceData sequence_data = sequence_datas_cur_1_[i2][p];
                                    bool ok = append(
                                            sequence_data, (p == 0),
                                            sequence_1.elements[pos_1],
                                            gc, nullptr);
                                    if (!ok)
                                        continue;
                                    ok = concatenate(
                                            sequence_data, false,
                                            sequence_2, p, pos_2 - 1, false,
                                            gc, nullptr);
                                    if (!ok)
                                        continue;
                                    ok = concatenate(
                                            sequence_data, false,
                                            sequence_2, pos_2 + 1, seq_2_size - 1, false,
                                            gc, nullptr);
                                    if (!ok)
                                        continue;
                                    gci2 = sequencing_scheme_.global_cost(sequence_data);
                                }
                                //std::cout << "p " << p
                                //    << " gc " << to_string(gci2)
                                //    << std::endl;
                                if ((pos_1_new == -1 || strictly_better(gci2, gci2_tmp))) {
                                    gci2_tmp = gci2;
                                    pos_1_new = p;
                                }
                            }
                        }
                        if (pos_1_new == -1)
                            continue;

                        // Find pos_2_new.
                        ElementPos pos_2_new = -1;
                        GlobalCost gci1_tmp;
                        for (int a = 0; a < 4; ++a) {
                            ElementPos p = (a < 3)? inter_swap_star_best_positions_[i1][j2][a]: pos_1;
                            if (p != -1) {
                                // pos_1: removed positino
                                // p: added position
                                GlobalCost gci1;
                                if (pos_1 < p) {
                                    SequenceData sequence_data = sequence_datas_cur_1_[i1][pos_1];
                                    bool ok = concatenate(
                                            sequence_data, (pos_1 == 0),
                                            sequence_1, pos_1 + 1, p - 1, false,
                                            gc, nullptr);
                                    if (!ok)
                                        continue;
                                    ok = append(
                                            sequence_data, (pos_1 == 0 && pos_1 + 1 > p - 1),
                                            sequence_2.elements[pos_2],
                                            gc, nullptr);
                                    if (!ok)
                                        continue;
                                    ok = concatenate(
                                            sequence_data, false,
                                            sequence_1, p, seq_1_size - 1, false,
                                            gc, nullptr);
                                    if (!ok)
                                        continue;
                                    gci1 = sequencing_scheme_.global_cost(sequence_data);
                                } else {
                                    SequenceData sequence_data = sequence_datas_cur_1_[i1][p];
                                    bool ok = append(
                                            sequence_data, (p == 0),
                                            sequence_2.elements[pos_2],
                                            gc, &gci2_tmp);
                                    if (!ok)
                                        continue;
                                    ok = concatenate(
                                            sequence_data, false,
                                            sequence_1, p, pos_1 - 1, false,
                                            gc, &gci2_tmp);
                                    if (!ok)
                                        continue;
                                    ok = concatenate(
                                            sequence_data, false,
                                            sequence_1, pos_1 + 1, seq_1_size - 1, false,
                                            gc, &gci2_tmp);
                                    if (!ok)
                                        continue;
                                    gci1 = sequencing_scheme_.global_cost(sequence_data);
                                }
                                //std::cout << "p " << p
                                //    << " gc " << to_string(gci1)
                                //    << std::endl;
                                if ((pos_2_new == -1 || strictly_better(gci1, gci1_tmp))) {
                                    gci1_tmp = gci1;
                                    pos_2_new = p;
                                }
                            }
                        }
                        if (pos_2_new == -1)
                            continue;
                        GlobalCost gc_tmp = merge(
                                gci1_tmp,
                                gci2_tmp);

                        //std::cout << "pos_1 " << pos_2
                        //    << " pos_2 " << pos_2
                        //    << " pos_1_new " << pos_1_new
                        //    << " pos_2_new " << pos_2_new
                        //    << std::endl;

                        if (strictly_better(gc_tmp, gc)) {
                            Move move;
                            move.type = Neighborhoods::InterSwapStar;
                            move.sequence_id_1 = i1;
                            move.sequence_id_2 = i2;
                            move.pos_1 = pos_1;
                            move.pos_2 = pos_1_new;
                            move.pos_3 = pos_2;
                            move.pos_4 = pos_2_new;
                            move.global_cost = diff(gc_tmp, gc);
                            neighborhood.improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    inline void explore_shift_change_mode(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::ShiftChangeMode)][0][0];
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!neighborhood.modified_sequences[sequence_id])
                continue;
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            if (!(parameters_.linking_constraints && m > 1))
                gc = global_costs_cur_[sequence_id];

            for (ElementPos block_pos = 0; block_pos <= seq_size - 1; ++block_pos) {
                ElementId j = sequence.elements[block_pos].element_id;
                // Loop through all new positions.
                ElementPos pos_min = std::max(
                        (ElementPos)0,
                        block_pos - parameters_.shift_maximum_distance);
                ElementPos pos_max = std::min(
                        seq_size - 1,
                        block_pos + parameters_.shift_maximum_distance);
                for (ElementPos pos_new = pos_min; pos_new <= pos_max; ++pos_new) {
                    for (Mode mode = 0; mode < number_of_modes(j); ++mode) {

                        GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                            &partial_global_costs_cur_1_[sequence_id]: nullptr;
                        GlobalCost gc_tmp;
                        if (block_pos > pos_new) {
                            SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos_new];
                            bool ok = append(
                                    sequence_data, (pos_new == 0),
                                    {j, mode},
                                    gc, gcm);
                            if (!ok)
                                continue;
                            ok = concatenate(
                                    sequence_data, false,
                                    sequence, pos_new, block_pos - 1, false,
                                    gc, gcm);
                            if (!ok)
                                continue;
                            ok = concatenate(
                                    sequence_data, false,
                                    sequence, block_pos + 1, seq_size - 1, false,
                                    gc, gcm);
                            if (!ok)
                                continue;
                            gc_tmp = (parameters_.linking_constraints && m > 1)?
                                merge(
                                        sequencing_scheme_.global_cost(sequence_data),
                                        partial_global_costs_cur_1_[sequence_id]):
                                sequencing_scheme_.global_cost(sequence_data);
                        } else {
                            SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][block_pos];
                            bool ok = concatenate(
                                    sequence_data, (block_pos == 0),
                                    sequence, block_pos + 1, pos_new + 1 - 1, false,
                                    gc, gcm);
                            if (!ok)
                                continue;
                            ok = append(
                                    sequence_data, false,
                                    {j, mode},
                                    gc, gcm);
                            if (!ok)
                                continue;
                            ok = concatenate(
                                    sequence_data, false,
                                    sequence, pos_new + 1, seq_size - 1, false,
                                    gc, gcm);
                            if (!ok)
                                continue;
                            gc_tmp = (parameters_.linking_constraints && m > 1)?
                                merge(
                                        sequencing_scheme_.global_cost(sequence_data),
                                        partial_global_costs_cur_1_[sequence_id]):
                                sequencing_scheme_.global_cost(sequence_data);
                        }
                        if (strictly_better(gc_tmp, gc)) {
                            Move move;
                            move.type = Neighborhoods::ShiftChangeMode;
                            move.sequence_id_1 = sequence_id;
                            move.pos_1 = block_pos;
                            move.pos_2 = pos_new;
                            move.mode = mode;
                            move.global_cost = diff(gc_tmp, gc);
                            neighborhood.improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    inline void explore_mode_swap(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::ModeSwap)][0][0];
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!neighborhood.modified_sequences[sequence_id])
                continue;
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            if (!(parameters_.linking_constraints && m > 1))
                gc = global_costs_cur_[sequence_id];

            // Loop through all pairs.
            Counter pos_max = seq_size - 1;
            for (ElementPos pos_1 = 0; pos_1 <= pos_max; ++pos_1) {
                ElementId j1 = sequence.elements[pos_1].element_id;
                ElementId mode_1 = sequence.elements[pos_1].mode;

                ElementPos pos_2_max = std::min(
                        pos_max,
                        pos_1 + parameters_.swap_maximum_distance);
                for (ElementPos pos_2 = pos_1 + 1; pos_2 < pos_2_max; ++pos_2) {
                    ElementId j2 = sequence.elements[pos_2].element_id;
                    ElementId mode_2 = sequence.elements[pos_2].mode;
                    if (mode_1 == mode_2)
                        continue;

                    GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                        &partial_global_costs_cur_1_[sequence_id]: nullptr;
                    SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos_1];
                    bool ok = append(
                            sequence_data, (pos_1 == 0),
                            {j1, mode_2},
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data, false,
                            sequence, pos_1 + 1, pos_2 - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = append(
                            sequence_data, false,
                            {j2, mode_1},
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data, false,
                            sequence, pos_2 + 1, seq_size - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp = (parameters_.linking_constraints && m > 1)?
                        merge(
                                sequencing_scheme_.global_cost(sequence_data),
                                partial_global_costs_cur_1_[sequence_id]):
                        sequencing_scheme_.global_cost(sequence_data);

                    if (strictly_better(gc_tmp, gc)) {
                        Move move;
                        move.type = Neighborhoods::ModeSwap;
                        move.sequence_id_1 = sequence_id;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = diff(gc_tmp, gc);
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_swap_with_modes(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::SwapWithModes)][0][0];
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            if (!neighborhood.modified_sequences[sequence_id])
                continue;
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            if (!(parameters_.linking_constraints && m > 1))
                gc = global_costs_cur_[sequence_id];

            // Loop through all pairs.
            Counter pos_max = seq_size - 1;
            for (ElementPos pos_1 = 0; pos_1 <= pos_max; ++pos_1) {
                ElementId j1 = sequence.elements[pos_1].element_id;
                ElementId mode_1 = sequence.elements[pos_1].mode;

                ElementPos pos_2_max = std::min(
                        pos_max,
                        pos_1 + parameters_.swap_maximum_distance);
                for (ElementPos pos_2 = pos_1 + 1; pos_2 < pos_2_max; ++pos_2) {
                    ElementId j2 = sequence.elements[pos_2].element_id;
                    ElementId mode_2 = sequence.elements[pos_2].mode;
                    if (mode_1 == mode_2)
                        continue;

                    GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                        &partial_global_costs_cur_1_[sequence_id]: nullptr;
                    SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos_1];
                    bool ok = append(
                            sequence_data, (pos_1 == 0),
                            {j2, mode_1},
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data, false,
                            sequence, pos_1 + 1, pos_2 - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = append(
                            sequence_data, false,
                            {j1, mode_2},
                            gc, gcm);
                    if (!ok)
                        continue;
                    ok = concatenate(
                            sequence_data, false,
                            sequence, pos_2 + 1, seq_size - 1, false,
                            gc, gcm);
                    if (!ok)
                        continue;
                    GlobalCost gc_tmp = (parameters_.linking_constraints && m > 1)?
                        merge(
                                sequencing_scheme_.global_cost(sequence_data),
                                partial_global_costs_cur_1_[sequence_id]):
                        sequencing_scheme_.global_cost(sequence_data);

                    if (ok && strictly_better(gc_tmp, gc)) {
                        Move move;
                        move.type = Neighborhoods::SwapWithModes;
                        move.sequence_id_1 = sequence_id;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = diff(gc_tmp, gc);
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline std::vector<Move> explore_add_2(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        GlobalCost gc = global_cost(solution);
        std::vector<Move> improving_moves;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            if (!(parameters_.linking_constraints && m > 1))
                gc = global_costs_cur_[sequence_id];

            // Loop through all new positions.
            for (ElementPos pos = 0; pos <= seq_size; ++pos) {

                for (ElementId element_id = 0; element_id < n; ++element_id) {
                    if (elements_cur_[element_id].mode != -1)
                        continue;

                    for (Mode mode = 0; mode < number_of_modes(element_id); ++mode) {
                        GlobalCost gc_tmp;
                        bool accept = false;
                        if (!parameters_.add_remove) {
                            SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos];
                            append(
                                    sequence_data, (pos == 0),
                                    {element_id, mode});
                            concatenate(
                                    sequence_data, false,
                                    sequence, pos, seq_size - 1, false);
                            gc_tmp = (parameters_.linking_constraints && m > 1)?
                                merge(
                                        sequencing_scheme_.global_cost(sequence_data),
                                        partial_global_costs_cur_1_[sequence_id]):
                                sequencing_scheme_.global_cost(sequence_data);
                            accept = true;
                        } else {
                            GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                                &partial_global_costs_cur_1_[sequence_id]: nullptr;
                            SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos];
                            bool ok = append(
                                    sequence_data, (pos == 0),
                                    {element_id, mode},
                                    gc, gcm);
                            if (!ok)
                                continue;
                            ok = concatenate(
                                    sequence_data, false,
                                    sequence, pos, seq_size - 1, false,
                                    gc, gcm);
                            if (!ok)
                                continue;
                            gc_tmp = (parameters_.linking_constraints && m > 1)?
                                merge(
                                        sequencing_scheme_.global_cost(sequence_data),
                                        partial_global_costs_cur_1_[sequence_id]):
                                sequencing_scheme_.global_cost(sequence_data);
                            accept = (strictly_better(gc_tmp, gc));
                        }

                        if (accept) {
                            Move move;
                            move.type = Neighborhoods::Add;
                            move.sequence_id_1 = sequence_id;
                            move.element_id = element_id;
                            move.mode = mode;
                            move.pos_1 = pos;
                            move.global_cost = diff(gc_tmp, gc);
                            improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
        return improving_moves;
    }

    inline std::vector<Move> explore_add(
            const Solution& solution,
            ElementId j,
            SequenceId i_old = -1,
            ElementId j_prec_old = -2,
            Mode mode_old = -1)
    {
        //std::cout << "explore_add j " << j << std::endl;
        std::vector<Move> improving_moves;
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            if (!(parameters_.linking_constraints && m > 1))
                gc = global_costs_cur_[sequence_id];

            // Loop through all new positions.
            for (ElementPos pos = 0; pos <= seq_size; ++pos) {

                for (Mode mode = 0; mode < number_of_modes(j); ++mode) {

                    if (j_prec_old == -1
                            && sequence_id == i_old
                            && (seq_size > 0 || m > 1)
                            && pos == 0
                            && mode == mode_old)
                        continue;
                    if (pos > 0
                            && j_prec_old == sequence.elements[pos - 1].element_id
                            && mode == mode_old)
                        continue;

                    GlobalCost gc_tmp;
                    bool accept = false;
                    if (!parameters_.add_remove) {
                        SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos];
                        append(
                                sequence_data, (pos == 0),
                                {j, mode});
                        concatenate(
                                sequence_data, false,
                                sequence, pos, seq_size - 1, false);
                        gc_tmp = (parameters_.linking_constraints && m > 1)?
                            merge(
                                    sequencing_scheme_.global_cost(sequence_data),
                                    partial_global_costs_cur_1_[sequence_id]):
                            sequencing_scheme_.global_cost(sequence_data);
                        accept = true;
                    } else {
                        GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                            &partial_global_costs_cur_1_[sequence_id]: nullptr;
                        SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][pos];
                        bool ok = append(
                                sequence_data, (pos == 0),
                                {j, mode},
                                gc, gcm);
                        if (!ok)
                            continue;
                        ok = concatenate(
                                sequence_data, false,
                                sequence, pos, seq_size - 1, false,
                                gc, gcm);
                        if (!ok)
                            continue;
                        gc_tmp = (parameters_.linking_constraints && m > 1)?
                            merge(
                                    sequencing_scheme_.global_cost(sequence_data),
                                    partial_global_costs_cur_1_[sequence_id]):
                            sequencing_scheme_.global_cost(sequence_data);
                        accept = (strictly_better(gc_tmp, gc));
                    }

                    //std::cout << to_string(gc_tmp) << std::endl;
                    if (accept) {
                        Move move;
                        move.type = Neighborhoods::Add;
                        move.sequence_id_1 = sequence_id;
                        move.element_id = j;
                        move.mode = mode;
                        move.pos_1 = pos;
                        move.global_cost = diff(gc_tmp, gc);
                        improving_moves.push_back(move);
                    }
                }
            }
        }
        return improving_moves;
    }

    inline std::vector<Move> explore_add_end(
            const Solution& solution)
    {
        std::vector<Move> improving_moves;
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        GlobalCost gc = global_cost(solution);
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id) {
            const auto& sequence = solution.sequences[sequence_id];
            SequencePos seq_size = sequence.elements.size();
            if (!(parameters_.linking_constraints && m > 1))
                gc = global_costs_cur_[sequence_id];

            // Loop through all new positions.
            for (ElementId element_id = 0; element_id < n; ++element_id) {
                if (elements_cur_[element_id].mode != -1)
                    continue;

                for (Mode mode = 0; mode < number_of_modes(element_id); ++mode) {
                    GlobalCost gc_tmp;
                    bool accept = false;
                    if (!parameters_.add_remove) {
                        SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][seq_size];
                        append(
                                sequence_data, (seq_size == 0),
                                {element_id, mode});
                        gc_tmp = (parameters_.linking_constraints && m > 1)?
                            merge(
                                    sequencing_scheme_.global_cost(sequence_data),
                                    partial_global_costs_cur_1_[sequence_id]):
                            sequencing_scheme_.global_cost(sequence_data);
                        accept = true;
                    } else {
                        GlobalCost* gcm = (parameters_.linking_constraints && m > 1)?
                            &partial_global_costs_cur_1_[sequence_id]: nullptr;
                        SequenceData sequence_data = sequence_datas_cur_1_[sequence_id][seq_size];
                        bool ok = append(
                                sequence_data, (seq_size == 0),
                                {element_id, mode},
                                gc, gcm);
                        if (!ok)
                            continue;
                        gc_tmp = (parameters_.linking_constraints && m > 1)?
                            merge(
                                    sequencing_scheme_.global_cost(sequence_data),
                                    partial_global_costs_cur_1_[sequence_id]):
                            sequencing_scheme_.global_cost(sequence_data);
                        accept = (strictly_better(gc_tmp, gc));
                    }

                    if (accept) {
                        Move move;
                        move.type = Neighborhoods::Add;
                        move.sequence_id_1 = sequence_id;
                        move.element_id = element_id;
                        move.mode = mode;
                        move.pos_1 = seq_size;
                        move.global_cost = diff(gc_tmp, gc);
                        improving_moves.push_back(move);
                    }
                }
            }
        }
        return improving_moves;
    }

    inline void apply_move(
            Solution& solution,
            const Move& move)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();

        ElementPos n_old = 0;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
            n_old += solution.sequences[sequence_id].elements.size();
        if (n_old > n) {
            throw std::logic_error(
                    "Too many element before applying move "
                    + neighborhood2string(move.type, move.k1, move.k2)
                    + ": " + std::to_string(n_old)
                    + " / " + std::to_string(n) + ".");
        }

        switch (move.type) {
        case Neighborhoods::Shift: {
            ElementPos block_size = move.k1;
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements = solution.sequences[move.sequence_id_1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp = empty_sequence(move.sequence_id_1);
            if (move.pos_1 > move.pos_2) {
                for (ElementPos p = 0; p < move.pos_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1; p < move.pos_1 + block_size; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            } else {
                for (ElementPos p = 0; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size; p < move.pos_2 + block_size; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1; p < move.pos_1 + block_size; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2 + block_size; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::Swap: {
            ElementPos block_size_1 = move.k1;
            ElementPos block_size_2 = move.k2;
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements = solution.sequences[move.sequence_id_1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp = empty_sequence(move.sequence_id_1);
            if (move.pos_1 < move.pos_2) {
                for (ElementPos p = 0; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2; p < move.pos_2 + block_size_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size_1; p < move.pos_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1; p < move.pos_1 + block_size_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2 + block_size_2; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            } else {
                for (ElementPos p = 0; p < move.pos_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1; p < move.pos_1 + block_size_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2 + block_size_2; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2; p < move.pos_2 + block_size_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size_1; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::Reverse: {
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements = solution.sequences[move.sequence_id_1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp = empty_sequence(move.sequence_id_1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp, elements[p]);
            for (ElementPos p = move.pos_2; p >= move.pos_1; --p)
                append(sequence_tmp, elements[p]);
            for (ElementPos p = move.pos_2 + 1; p < seq_size; ++p)
                append(sequence_tmp, elements[p]);
            compute_global_cost(solution_tmp_);
            solution = solution_tmp_;
            break;
        } case Neighborhoods::ShiftReverse: {
            ElementPos block_size = move.k1;
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements = solution.sequences[move.sequence_id_1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp = empty_sequence(move.sequence_id_1);
            if (move.pos_1 > move.pos_2) {
                for (ElementPos p = 0; p < move.pos_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size - 1; p >= move.pos_1; --p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            } else {
                for (ElementPos p = 0; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size; p < move.pos_2 + block_size; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size - 1; p >= move.pos_1; --p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2 + block_size; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::Add: {
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements = solution.sequences[move.sequence_id_1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp = empty_sequence(move.sequence_id_1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp, elements[p]);
            append(sequence_tmp, {move.element_id, move.mode});
            for (ElementPos p = move.pos_1; p < seq_size; ++p)
                append(sequence_tmp, elements[p]);
            break;
        } case Neighborhoods::Remove: {
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements = solution.sequences[move.sequence_id_1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp = empty_sequence(move.sequence_id_1);
            for (ElementPos p = 0; p < seq_size; ++p) {
                if (p == move.pos_1)
                    continue;
                append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::Replace: {
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements = solution.sequences[move.sequence_id_1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp = empty_sequence(move.sequence_id_1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp, elements[p]);
            append(sequence_tmp, {move.element_id, move.mode});
            for (ElementPos p = move.pos_1 + 1; p < seq_size; ++p)
                append(sequence_tmp, elements[p]);
            break;
        } case Neighborhoods::SwapTails: {
            //std::cout << "Apply SwapTails"
            //    << " i1 " << move.sequence_id_1
            //    << " i2 " << move.sequence_id_2
            //    << " pos_1 " << move.pos_1
            //    << " pos_2 " << move.pos_2
            //    << " gc " << to_string(move.global_cost)
            //    << std::endl;
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1 && sequence_id != move.sequence_id_2)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements_1 = solution.sequences[move.sequence_id_1].elements;
            const auto& elements_2 = solution.sequences[move.sequence_id_2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp_1 = empty_sequence(move.sequence_id_1);
            for (ElementPos p = 0; p <= move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_2 + 1; p <= seq_2_size - 1; ++p)
                append(sequence_tmp_1, elements_2[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.sequence_id_2];
            sequence_tmp_2 = empty_sequence(move.sequence_id_2);
            for (ElementPos p = 0; p <= move.pos_2; ++p)
                append(sequence_tmp_2, elements_2[p]);
            for (ElementPos p = move.pos_1 + 1; p <= seq_1_size - 1; ++p)
                append(sequence_tmp_2, elements_1[p]);
            break;
        } case Neighborhoods::Split: {
            //std::cout << "Apply Split" << std::endl;
            //std::cout << "i1 " << move.sequence_id_1
            //    << " i2 " << move.sequence_id_2
            //    << " pos_1 " << move.pos_1
            //    << " pos_2 " << move.pos_2
            //    << " gc " << to_string(move.global_cost)
            //    << std::endl;
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1 && sequence_id != move.sequence_id_2)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements_1 = solution.sequences[move.sequence_id_1].elements;
            const auto& elements_2 = solution.sequences[move.sequence_id_2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp_1 = empty_sequence(move.sequence_id_1);
            for (ElementPos p = 0; p <= move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_2; p >= 0; --p)
                append(sequence_tmp_1, elements_2[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.sequence_id_2];
            sequence_tmp_2 = empty_sequence(move.sequence_id_2);
            for (ElementPos p = seq_1_size - 1; p >= move.pos_1 + 1; --p)
                append(sequence_tmp_2, elements_1[p]);
            for (ElementPos p = move.pos_2 + 1; p <= seq_2_size - 1; ++p)
                append(sequence_tmp_2, elements_2[p]);
            //std::cout << "End apply Split" << std::endl;
            break;
        } case Neighborhoods::InterShift: {
            ElementPos block_size = move.k1;
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1 && sequence_id != move.sequence_id_2)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements_1 = solution.sequences[move.sequence_id_1].elements;
            const auto& elements_2 = solution.sequences[move.sequence_id_2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp_1 = empty_sequence(move.sequence_id_1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_1 + block_size; p < seq_1_size; ++p)
                append(sequence_tmp_1, elements_1[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.sequence_id_2];
            sequence_tmp_2 = empty_sequence(move.sequence_id_2);
            for (ElementPos p = 0; p < move.pos_2; ++p)
                append(sequence_tmp_2, elements_2[p]);
            for (ElementPos p = move.pos_1; p < move.pos_1 + block_size; ++p)
                append(sequence_tmp_2, elements_1[p]);
            for (ElementPos p = move.pos_2; p < seq_2_size; ++p)
                append(sequence_tmp_2, elements_2[p]);
            break;
        } case Neighborhoods::InterSwap: {
            //std::cout << "Apply InterSwap"
            //    << " k1 " << move.k1
            //    << " k2 " << move.k2
            //    << " i1 " << move.sequence_id_1
            //    << " i2 " << move.sequence_id_2
            //    << " pos_1 " << move.pos_1
            //    << " pos_2 " << move.pos_2
            //    << " gc " << to_string(move.global_cost)
            //    << std::endl;
            ElementPos block_size_1 = move.k1;
            ElementPos block_size_2 = move.k2;
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1 && sequence_id != move.sequence_id_2)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements_1 = solution.sequences[move.sequence_id_1].elements;
            const auto& elements_2 = solution.sequences[move.sequence_id_2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp_1 = empty_sequence(move.sequence_id_1);
            for (ElementPos p = 0; p <= move.pos_1 - 1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_2; p <= move.pos_2 + block_size_2 - 1; ++p)
                append(sequence_tmp_1, elements_2[p]);
            for (ElementPos p = move.pos_1 + block_size_1; p <= seq_1_size - 1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.sequence_id_2];
            sequence_tmp_2 = empty_sequence(move.sequence_id_2);
            for (ElementPos p = 0; p <= move.pos_2 - 1; ++p)
                append(sequence_tmp_2, elements_2[p]);
            for (ElementPos p = move.pos_1; p <= move.pos_1 + block_size_1 - 1; ++p)
                append(sequence_tmp_2, elements_1[p]);
            for (ElementPos p = move.pos_2 + block_size_2; p <= seq_2_size - 1; ++p)
                append(sequence_tmp_2, elements_2[p]);
            break;
        } case Neighborhoods::InterShiftReverse: {
            //std::cout << "Apply InterShiftReverse"
            //    << " k1 " << move.k1
            //    << " i1 " << move.sequence_id_1
            //    << " i2 " << move.sequence_id_2
            //    << " pos_1 " << move.pos_1
            //    << " pos_2 " << move.pos_2
            //    << " gc " << to_string(move.global_cost)
            //    << std::endl;
            ElementPos block_size = move.k1;
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1 && sequence_id != move.sequence_id_2)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements_1 = solution.sequences[move.sequence_id_1].elements;
            const auto& elements_2 = solution.sequences[move.sequence_id_2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp_1 = empty_sequence(move.sequence_id_1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_1 + block_size; p < seq_1_size; ++p)
                append(sequence_tmp_1, elements_1[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.sequence_id_2];
            sequence_tmp_2 = empty_sequence(move.sequence_id_2);
            for (ElementPos p = 0; p < move.pos_2; ++p)
                append(sequence_tmp_2, elements_2[p]);
            for (ElementPos p = move.pos_1 + block_size - 1; p >= move.pos_1; --p)
                append(sequence_tmp_2, elements_1[p]);
            for (ElementPos p = move.pos_2; p < seq_2_size; ++p)
                append(sequence_tmp_2, elements_2[p]);
            break;
        } case Neighborhoods::ShiftChangeMode: {
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements = solution.sequences[move.sequence_id_1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp = empty_sequence(move.sequence_id_1);
            if (move.pos_1 > move.pos_2) {
                for (ElementPos p = 0; p < move.pos_2; ++p)
                    append(sequence_tmp, elements[p]);
                append(sequence_tmp, {elements[move.pos_1].element_id, move.mode});
                for (ElementPos p = move.pos_2; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + 1; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            } else {
                for (ElementPos p = 0; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + 1; p < move.pos_2 + 1; ++p)
                    append(sequence_tmp, elements[p]);
                append(sequence_tmp, {elements[move.pos_1].element_id, move.mode});
                for (ElementPos p = move.pos_2 + 1; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::ModeSwap: {
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements = solution.sequences[move.sequence_id_1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp = empty_sequence(move.sequence_id_1);
            for (ElementPos p = 0; p < seq_size; ++p) {
                SequenceElement se = elements[p];
                if (p == move.pos_1) {
                    se.mode = elements[move.pos_2].mode;
                } else if (p == move.pos_2) {
                    se.mode = elements[move.pos_1].mode;
                }
                append(sequence_tmp, se);
            }
            break;
        } case Neighborhoods::SwapWithModes: {
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements = solution.sequences[move.sequence_id_1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.sequence_id_1];
            sequence_tmp = empty_sequence(move.sequence_id_1);
            for (ElementPos p = 0; p < seq_size; ++p) {
                SequenceElement se = elements[p];
                if (p == move.pos_1) {
                    se.element_id = elements[move.pos_2].element_id;
                } else if (p == move.pos_2) {
                    se.element_id = elements[move.pos_1].element_id;
                }
                append(sequence_tmp, se);
            }
            break;
        } case Neighborhoods::InterSwapStar: {
            //std::cout << "i1 " << move.sequence_id_1
            //    << " i2 " << move.sequence_id_2
            //    << " pos_1_old " << move.pos_1
            //    << " pos_2_old " << move.pos_3
            //    << " pos_1_new " << move.pos_2
            //    << " pos_2_new " << move.pos_4
            //    << " gc " << to_string(move.global_cost)
            //    << std::endl;
            for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
                if (sequence_id != move.sequence_id_1 && sequence_id != move.sequence_id_2)
                    solution_tmp_.sequences[sequence_id] = solution.sequences[sequence_id];
            const auto& elements_1 = solution.sequences[move.sequence_id_1].elements;
            const auto& elements_2 = solution.sequences[move.sequence_id_2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.sequence_id_1];
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.sequence_id_2];
            sequence_tmp_1 = empty_sequence(move.sequence_id_1);
            ElementPos pos_1_old = move.pos_1;
            ElementPos pos_1_new = move.pos_2;
            ElementPos pos_2_old = move.pos_3;
            ElementPos pos_2_new = move.pos_4;
            if (pos_1_old < pos_2_new) {
                for (ElementPos p = 0; p < pos_1_old; ++p)
                    append(sequence_tmp_1, elements_1[p]);
                for (ElementPos p = pos_1_old + 1; p < pos_2_new; ++p)
                    append(sequence_tmp_1, elements_1[p]);
                append(sequence_tmp_1, elements_2[pos_2_old]);
                for (ElementPos p = pos_2_new; p < seq_1_size; ++p)
                    append(sequence_tmp_1, elements_1[p]);
            } else {
                for (ElementPos p = 0; p < pos_2_new; ++p)
                    append(sequence_tmp_1, elements_1[p]);
                append(sequence_tmp_1, elements_2[pos_2_old]);
                for (ElementPos p = pos_2_new; p < pos_1_old; ++p)
                    append(sequence_tmp_1, elements_1[p]);
                for (ElementPos p = pos_1_old + 1; p < seq_1_size; ++p)
                    append(sequence_tmp_1, elements_1[p]);
            }
            sequence_tmp_2 = empty_sequence(move.sequence_id_2);
            if (pos_2_old < pos_1_new) {
                for (ElementPos p = 0; p < pos_2_old; ++p)
                    append(sequence_tmp_2, elements_2[p]);
                for (ElementPos p = pos_2_old + 1; p < pos_1_new; ++p)
                    append(sequence_tmp_2, elements_2[p]);
                append(sequence_tmp_2, elements_1[pos_1_old]);
                for (ElementPos p = pos_1_new; p < seq_2_size; ++p)
                    append(sequence_tmp_2, elements_2[p]);
            } else {
                for (ElementPos p = 0; p < pos_1_new; ++p)
                    append(sequence_tmp_2, elements_2[p]);
                append(sequence_tmp_2, elements_1[pos_1_old]);
                for (ElementPos p = pos_1_new; p < pos_2_old; ++p)
                    append(sequence_tmp_2, elements_2[p]);
                for (ElementPos p = pos_2_old + 1; p < seq_2_size; ++p)
                    append(sequence_tmp_2, elements_2[p]);
            }
            break;
        }
        }

        ElementPos n_new = 0;
        for (SequenceId sequence_id = 0; sequence_id < m; ++sequence_id)
            n_new += solution_tmp_.sequences[sequence_id].elements.size();
        if (n_new > n) {
            throw std::logic_error(
                    "Too many elements after applying move "
                    + neighborhood2string(move.type, move.k1, move.k2)
                    + ": " + std::to_string(n_new)
                    + " / " + std::to_string(n) + ".");
        }
        if (move.type != Neighborhoods::Add
                && move.type != Neighborhoods::Remove
                && n_new != n_old) {
            throw std::logic_error(
                    "Wrong number of elements after applying move "
                    + neighborhood2string(move.type, move.k1, move.k2)
                    + "; old: " + std::to_string(n_old)
                    + "; new " + std::to_string(n_new) + ".");
        }
        if (move.type == Neighborhoods::Add && n_new != n_old + 1) {
            throw std::logic_error(
                    "Wrong number of elements after applying move "
                    + neighborhood2string(move.type, move.k1, move.k2)
                    + "; old: " + std::to_string(n_old)
                    + "; new " + std::to_string(n_new) + ".");
        }
        if (move.type == Neighborhoods::Remove && n_new != n_old - 1) {
            throw std::logic_error(
                    "Wrong number of elements after applying move "
                    + neighborhood2string(move.type, move.k1, move.k2)
                    + "; old: " + std::to_string(n_old)
                    + "; new " + std::to_string(n_new) + ".");
        }

        solution_tmp_.modified_sequences = solution.modified_sequences;
        compute_global_cost(solution_tmp_);
        solution = solution_tmp_;
    }

    /*
     * Private attributes.
     */

    /** Sequencing scheme. */
    SequencingScheme& sequencing_scheme_;

    /** Parameters. */
    Parameters parameters_;

    /**
     * Maximum number of modes.
     *
     * Comupted in the constructor.
     */
    Mode maximum_number_of_modes_ = 1;

    /** Structure storing neighborhood related information. */
    std::vector<std::vector<std::vector<Neighborhood>>> neighborhoods_
        = std::vector<std::vector<std::vector<Neighborhood>>>(16);

    /** inter_shift_1_best_positions_[i][j] = {pos, GlobalCost}. */
    std::vector<std::vector<std::pair<ElementPos, GlobalCost>>> inter_shift_1_best_positions_;

    /** inter_swap_star_best_positions_[i][j] = {pos_1, pos_2, pos_3}. */
    std::vector<std::vector<std::vector<ElementPos>>> inter_swap_star_best_positions_;

    /*
     * Sorted predecessors/successors.
     */

    /** Structure storing the sorted predecessors for each element. */
    std::vector<std::vector<ElementId>> sorted_predecessors_;

    /** Structure storing the sorted successors for each element. */
    std::vector<std::vector<ElementId>> sorted_successors_;

    /**
     * Vector which contains all element ids.
     *
     * It is used to loop through neighbors when sorted_neighbors_ is not
     * used.
     */
    std::vector<ElementId> neighbors_;

    /*
     * Statistics.
     */

    /** Number of calls to a crossover operator. */
    Counter number_of_crossover_calls_ = 0;
    /** Time spent in crossover operators. */
    double crossover_time_ = 0.0;
    /** Number of calls to the local search method. */
    Counter number_of_local_search_calls_ = 0;
    /** Time spent in the local search method. */
    double local_search_time_ = 0.0;
    /** Number of local search iterations. */
    Counter number_of_local_search_iterations_ = 0;
    /** Number of calls to an initial solution generator. */
    Counter number_of_initial_solution_calls_ = 0;
    /** Time spent generating initial solutions. */
    double initial_solution_time_ = 0.0;
    /** Time spent in method 'compute_temporary_structures'. */
    double compute_temporary_structures_time_ = 0.0;

    /*
     * Temporary structures.
     */

    std::vector<std::vector<SequenceData>> sequence_datas_cur_1_;
    std::vector<std::vector<std::vector<SequenceData>>> sequence_datas_cur_2_;
    std::vector<GlobalCost> global_costs_cur_;
    std::vector<GlobalCost> partial_global_costs_cur_1_;
    std::vector<std::vector<GlobalCost>> partial_global_costs_cur_2_;
    std::vector<SolutionElement> elements_cur_;
    Solution solution_cur_;
    Solution solution_tmp_;

};

}

}

