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

    InterTwoOpt,
    InterTwoOptReverse,
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
    case Neighborhoods::InterTwoOpt:
        return "Inter-two-opt";
    case Neighborhoods::InterTwoOptReverse:
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

    bool inter_two_opt = false;
    bool inter_two_opt_reverse = false;
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

    Counter double_bridge_number_of_perturbations = 0;
    Counter ruin_and_recreate_number_of_perturbations = 0;
    ElementPos ruin_and_recreate_number_of_elements_removed = 4;
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

    using SequenceData = typename SequencingScheme::SequenceData;

    struct SolutionElement
    {
        SequenceId i = -1;
        ElementPos pos = -1;
        Mode mode = -1;
    };

    struct SequenceElement
    {
        ElementId j;
        Mode mode;

        bool operator==(const SequenceElement& rhs) const
        {
            return j == rhs.j && mode == rhs.mode;
        }
    };

    struct Sequence
    {
        SequenceId i = -1;
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

    struct Move0
    {
        Neighborhoods type;
        ElementPos k1 = -1;
        ElementPos k2 = -1;

        SequenceId i1 = -1;
        SequenceId i2 = -1;
        ElementId j = -1;
        ElementPos pos_1 = -1;
        ElementPos pos_2 = -1;
        ElementPos pos_3 = -1;
        ElementPos pos_4 = -1;
        Mode mode = -1;

        GlobalCost global_cost = worst<GlobalCost>();
    };

    struct Neighborhood
    {
        std::vector<Move0> improving_moves = {};
        // modified_sequences[i] == true iff sequence i has changed since last
        // neighborhood exploration.
        std::vector<bool> modified_sequences;
        Counter number_of_explorations = 0;
        Counter number_of_successes = 0;
        double time = 0.0;
    };


    void append(
            Sequence& sequence,
            SequenceElement se) const
    {
        append(sequence.data, se);
        sequence.elements.push_back(se);
    }


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
                    optimizationtools::hash_combine(hash_tmp, hasher_j(se.j));
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
        sequence_datas_cur_2_ = std::vector<std::vector<std::vector<SequenceData>>>(m);
        global_costs_cur_ = std::vector<GlobalCost>(m);
        if (m > 1)
            partial_global_costs_cur_1_ = std::vector<GlobalCost>(m);
        if (m > 2)
            partial_global_costs_cur_2_ = std::vector<std::vector<GlobalCost>>(
                    m, std::vector<GlobalCost>(m));
        elements_cur_ = std::vector<SolutionElement>(n);
        solution_cur_ = empty_solution();
        solution_tmp_ = empty_solution();

        compute_closest_neighbors();

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

        neighborhoods_[int(Neighborhoods::InterTwoOpt)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::InterTwoOptReverse)] = {{{Neighborhood()}}};
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
    }

    LocalScheme(const LocalScheme& sequencing_scheme):
        LocalScheme(sequencing_scheme.sequencing_scheme_, sequencing_scheme.parameters_) { }

    virtual ~LocalScheme() { }

    inline Sequence empty_sequence(SequenceId i) const
    {
        Sequence sequence;
        sequence.i = i;
        sequence.data = empty_sequence_data(i);
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
        //    std::cout << " " << elements[pos].j;
        //std::cout << std::endl;

        std::uniform_int_distribution<SequenceId> d_i(0, m - 1);
        SequenceId i = d_i(generator);

        // Compute edge costs.
        std::vector<std::vector<GlobalCost>> edges(n + 1);
        for (ElementPos pos = 0; pos <= n; ++pos)
            edges[pos] = std::vector<GlobalCost>(pos);
        for (ElementPos pos_1 = 0; pos_1 < seq_size; ++pos_1) {
            SequenceData sequence_data = empty_sequence_data(i);
            for (ElementPos pos_2 = pos_1 + 1; pos_2 <= seq_size; ++pos_2) {
                append(sequence_data, elements[pos_2 - 1]);
                edges[pos_2][pos_1] = sequencing_scheme_.global_cost(sequence_data);
            }
        }

        // Run Bellman algorithm.
        std::vector<std::vector<ElementPos>> prev(
                m + 1, std::vector<ElementPos>(n + 1, -1));
        std::vector<std::vector<GlobalCost>> distance(
                m + 1, std::vector<GlobalCost>(n + 1, worst<GlobalCost>()));
        for (SequenceId i = 0; i < m; ++i) {
            for (ElementPos pos_1 = 0; pos_1 < seq_size; ++pos_1) {
                for (ElementPos pos_2 = pos_1 + 1; pos_2 <= seq_size; ++pos_2) {
                    if (distance[i + 1][pos_2] > distance[i][pos_1] + edges[pos_2][pos_1]) {
                        distance[i + 1][pos_2] = distance[i][pos_1] + edges[pos_2][pos_1];
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
            for (ElementPos p = pos; p < pos_prev; ++p)
                append(solution.sequences[i_cur], elements[p]);
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
                    Move0 move_best;
                    move_best.global_cost = worst<GlobalCost>();
                    for (const Move0 move: improving_moves)
                        if (dominates(move.global_cost, move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    compute_temporary_structures(solution, move_best.i1, move_best.i2);
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
                    Move0 move_best;
                    move_best.global_cost = worst<GlobalCost>();
                    for (const Move0 move: improving_moves)
                        if (dominates(move.global_cost, move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    compute_temporary_structures(solution, move_best.i1, move_best.i2);
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
                    Move0 move_best;
                    move_best.global_cost = worst<GlobalCost>();
                    for (const Move0 move: improving_moves)
                        if (dominates(move.global_cost, move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    compute_temporary_structures(solution, move_best.i1, move_best.i2);
                } else {
                    break;
                }
            }
            return solution;
        } case 4: {
            // Random sequence splitted.
            std::vector<SequenceElement> elements(n);
            for (ElementId j = 0; j < n; ++j) {
                std::uniform_int_distribution<SequenceId> d_mode(0, number_of_modes(j) - 1);
                Mode mode = d_mode(generator);
                elements[j] = {j, mode};
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
        initial_solution_time += time_span.count();
        number_of_initial_solution_calls++;
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
        Solution solution;
        switch (x) {
        case 0: {
            solution = crossover_ox(solution_parent_1, solution_parent_2, generator);
            break;
        } case 1: {
            solution = crossover_sjox(solution_parent_1, solution_parent_2, generator);
            break;
        } case 2: {
            solution = crossover_sbox(solution_parent_1, solution_parent_2, generator);
            break;
        } case 3: {
            solution = crossover_srex1(solution_parent_1, solution_parent_2, generator);
            break;
        } case 4: {
            solution = crossover_srex2(solution_parent_1, solution_parent_2, generator);
            break;
        } default: {
            solution = crossover_ox(solution_parent_1, solution_parent_2, generator);
            break;
        }
        }

        SequenceId m = number_of_sequences();
        solution.modified_sequences = std::vector<bool>(m, true);
        for (SequenceId i1 = 0; i1 < m; ++i1) {
            for (SequenceId i2 = 0; i2 < m; ++i2) {
                if (solution.sequences[i1].elements
                        == solution_parent_1.sequences[i2].elements
                        && sequencing_scheme_.global_cost(solution.sequences[i1].data)
                        == sequencing_scheme_.global_cost(solution_parent_1.sequences[i2].data)) {
                    solution.modified_sequences[i1] = false;
                    break;
                }
                if (solution.sequences[i1].elements
                        == solution_parent_2.sequences[i2].elements
                        && sequencing_scheme_.global_cost(solution.sequences[i1].data)
                        == sequencing_scheme_.global_cost(solution_parent_2.sequences[i2].data)) {
                    solution.modified_sequences[i1] = false;
                    break;
                }
            }
        }

        auto end = std::chrono::steady_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
        crossover_time += time_span.count();
        number_of_crossover_calls++;
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
        for (SequenceId i = 0; i < m; ++i) {
            elements_parent_1.insert(
                    elements_parent_1.begin(),
                    solution_parent_1.sequences[i].elements.begin(),
                    solution_parent_1.sequences[i].elements.end());
        }
        std::vector<SequenceElement> elements_parent_2;
        for (SequenceId i = 0; i < m; ++i) {
            elements_parent_2.insert(
                    elements_parent_2.begin(),
                    solution_parent_2.sequences[i].elements.begin(),
                    solution_parent_2.sequences[i].elements.end());
        }

        ElementPos seq_1_size = elements_parent_1.size();

        std::vector<ElementPos> edges = optimizationtools::bob_floyd<ElementPos>(
                2, seq_1_size + 1, generator);
        std::sort(edges.begin(), edges.end());

        ElementPos pos_1 = edges[0];
        ElementPos pos_2 = edges[1];

        std::vector<uint8_t> in_substring(seq_1_size, false);
        for (ElementPos pos = pos_1; pos < pos_2; ++pos) {
            ElementId j = elements_parent_1[pos].j;
            in_substring[j] = true;
        }

        std::vector<SequenceElement> elements;

        for (ElementPos pos = 0; pos < seq_1_size; ++pos) {
            if ((ElementPos)elements.size() == pos_1)
                for (ElementPos p = pos_1; p < pos_2; ++p)
                    elements.push_back(elements_parent_1[p]);
            ElementId j = elements_parent_2[pos].j;
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
        for (SequenceId i = 0; i < m; ++i) {
            elements_parent_1.insert(
                    elements_parent_1.begin(),
                    solution_parent_1.sequences[i].elements.begin(),
                    solution_parent_1.sequences[i].elements.end());
        }
        std::vector<SequenceElement> elements_parent_2;
        for (SequenceId i = 0; i < m; ++i) {
            elements_parent_2.insert(
                    elements_parent_2.begin(),
                    solution_parent_2.sequences[i].elements.begin(),
                    solution_parent_2.sequences[i].elements.end());
        }

        ElementPos seq_1_size = elements_parent_1.size();
        std::vector<ElementPos> positions(seq_1_size, -1);

        std::vector<SequenceElement> elements;

        // Add elements from parent_1 up to a given cut point.
        std::uniform_int_distribution<ElementPos> d_point(1, seq_1_size);
        ElementPos pos_0 = d_point(generator);
        for (ElementPos pos = 0; pos < pos_0; ++pos) {
            ElementId j = elements_parent_1[pos].j;
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
                ElementId j1 = elements_parent_1[p].j;
                ElementId j2 = elements_parent_2[p].j;
                if (j1 == j2) {
                    positions[j1] = p;
                    elements.push_back(elements_parent_1[p]);
                    continue;
                }
                break;
            }
            ElementId j = elements_parent_2[pos].j;
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
        for (SequenceId i = 0; i < m; ++i) {
            elements_parent_1.insert(
                    elements_parent_1.begin(),
                    solution_parent_1.sequences[i].elements.begin(),
                    solution_parent_1.sequences[i].elements.end());
        }
        std::vector<SequenceElement> elements_parent_2;
        for (SequenceId i = 0; i < m; ++i) {
            elements_parent_2.insert(
                    elements_parent_2.begin(),
                    solution_parent_2.sequences[i].elements.begin(),
                    solution_parent_2.sequences[i].elements.end());
        }

        ElementPos seq_1_size = elements_parent_1.size();
        std::vector<ElementPos> positions(seq_1_size, -1);

        std::vector<SequenceElement> elements;

        // Add elements from parent_1 up to a given cut point.
        std::uniform_int_distribution<ElementPos> d_point(1, seq_1_size);
        ElementPos pos_0 = d_point(generator);
        for (ElementPos pos = 0; pos < pos_0; ++pos) {
            ElementId j = elements_parent_1[pos].j;
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
                ElementId j1 = elements_parent_1[p].j;
                ElementId j2 = elements_parent_2[p].j;
                if (p <= seq_1_size - 1) {
                    ElementId j1_next = elements_parent_1[p + 1].j;
                    ElementId j2_next = elements_parent_2[p + 1].j;
                    if (j1 == j2 && j1_next == j2_next) {
                        positions[j1] = p;
                        elements.push_back(elements_parent_1[p + 1]);
                        continue;
                    }
                }
                if (p >= 1) {
                    ElementId j1_prev = elements_parent_1[p - 1].j;
                    ElementId j2_prev = elements_parent_2[p - 1].j;
                    if (j1 == j2 && j1_prev == j2_prev) {
                        positions[j1] = p;
                        elements.push_back(elements_parent_1[p - 1]);
                        continue;
                    }
                }
                break;
            }
            ElementId j = elements_parent_2[pos].j;
            if (positions[j] != -1)
                continue;
            positions[j] = elements.size();
            elements.push_back(elements_parent_2[pos]);
        }

        return split(elements, generator);
    }

    inline Solution crossover_srex1(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        //std::cout << "crossover_srex1 start" << std::endl;
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        std::uniform_real_distribution<double> d01(0, 1);

        Solution solution = empty_solution();
        solution.modified_sequences = std::vector<bool>(m, false);

        // Draw sequences to inherite from the first parent.
        std::vector<bool> sequences(m, false);
        std::vector<bool> added(n, false);
        ElementId n_cur = 0;

        // Add sequences from the first parent.
        for (SequenceId i = 0; i < m; ++i) {
            if (d01(generator) > 0.5)
                continue;
            solution.sequences[i] = solution_parent_1.sequences[i];
            sequences[i] = true;
            for (const auto& se: solution.sequences[i].elements) {
                added[se.j] = true;
                n_cur++;
            }
            //std::cout << "i " << to_string(global_cost(solution.sequences[i])) << std::endl;
        }
        //std::cout << "n_cur_1 " << n_cur << " / " << n << std::endl;

        // Add sequences from the second parent.
        for (SequenceId i = 0; i < m; ++i) {
            if (sequences[i])
                continue;
            // Find the sequence from the second parent with the highest number
            // of unscheduled elements.
            SequenceId i2_best = -1;
            ElementPos seq_size_best = -1;
            for (SequenceId i2 = 0; i2 < m; ++i2) {
                ElementPos seq_size = 0;
                for (const auto& se: solution_parent_2.sequences[i2].elements)
                    if (!added[se.j])
                        seq_size++;
                if (seq_size <= (SequencePos)solution_parent_2.sequences[i2].elements.size() / 2)
                    continue;
                if (i2_best == 0 || seq_size_best < seq_size) {
                    i2_best = i2;
                    seq_size_best = seq_size;
                }
            }
            if (i2_best == -1)
                break;
            //std::cout << "i2_best " << i2_best << " seq_size_best " << seq_size_best << std::endl;

            solution.sequences[i] = empty_sequence(i);
            for (const auto& se: solution_parent_2.sequences[i2_best].elements) {
                if (added[se.j]) {
                    solution.modified_sequences[i] = true;
                    continue;
                }
                append(solution.sequences[i], se);
                added[se.j] = true;
                n_cur++;
            }
            //std::cout << "i " << i
            //    << " gc " << to_string(global_cost(solution.sequences[i]))
            //    << std::endl;
        }
        //std::cout << "n_cur_2 " << n_cur << " / " << n << std::endl;

        // Add remaining elements at random positions.
        compute_temporary_structures(solution);
        for (ElementId j = 0; j < n; ++j) {
            if (added[j])
                continue;
            auto improving_moves = explore_add(solution, j);
            if (!improving_moves.empty()) {
                std::shuffle(
                        improving_moves.begin(),
                        improving_moves.end(),
                        generator);
                Move0 move_best;
                move_best.global_cost = worst<GlobalCost>();
                for (const Move0 move: improving_moves)
                    if (dominates(move.global_cost, move_best.global_cost))
                        move_best = move;
                apply_move(solution, move_best);
                solution.modified_sequences[move_best.i1] = true;
                compute_temporary_structures(solution, move_best.i1, move_best.i2);
            }
        }

        compute_global_cost(solution);
        //std::cout << "crossover_srex1 end" << std::endl;
        return solution;
    }

    inline Solution crossover_srex2(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        std::uniform_real_distribution<double> d01(0, 1);

        Solution solution = empty_solution();
        solution.modified_sequences = std::vector<bool>(m, false);

        // Draw sequences to inherite from the first parent.
        std::vector<bool> sequences(m, false);
        std::vector<bool> added(n, false);

        // Add sequences from the first parent.
        for (SequenceId i = 0; i < m; ++i) {
            if (d01(generator) > 0.5)
                continue;
            solution.sequences[i] = solution_parent_1.sequences[i];
            sequences[i] = true;
            for (const auto& se: solution.sequences[i].elements)
                added[se.j] = true;
        }

        // Add sequences from the second parent.
        for (SequenceId i = 0; i < m; ++i) {
            if (sequences[i])
                continue;
            solution.sequences[i] = empty_sequence(i);
            for (const auto& se: solution_parent_2.sequences[i].elements) {
                if (added[se.j]) {
                    solution.modified_sequences[i] = true;
                    continue;
                }
                append(solution.sequences[i], se);
                added[se.j] = true;
            }
        }

        // Add remaining elements at random positions.
        for (ElementId j = 0; j < n; ++j) {
            if (added[j])
                continue;
            auto improving_moves = explore_add(solution, j);
            if (!improving_moves.empty()) {
                std::shuffle(
                        improving_moves.begin(),
                        improving_moves.end(),
                        generator);
                Move0 move_best;
                move_best.global_cost = worst<GlobalCost>();
                for (const Move0 move: improving_moves)
                    if (dominates(move.global_cost, move_best.global_cost))
                        move_best = move;
                apply_move(solution, move_best);
                solution.modified_sequences[move_best.i1] = true;
                compute_temporary_structures(solution, move_best.i1, move_best.i2);
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
        for (SequenceId i = 0; i < m; ++i) {
            const auto& elements = solution_1.sequences[i].elements;
            ElementPos seq_size = elements.size();
            for (ElementPos pos = 0; pos < seq_size - 1; ++pos) {
                ElementId j = elements[pos].j;
                ElementId j_next = elements[pos + 1].j;
                next_1[j] = j_next;
            }
        }

        // Compare with second solution.
        ElementPos d = 0;
        for (SequenceId i = 0; i < m; ++i) {
            const auto& elements = solution_2.sequences[i].elements;
            ElementPos seq_size = elements.size();
            for (ElementPos pos = 0; pos < seq_size - 1; ++pos) {
                ElementId j = elements[pos].j;
                ElementId j_next = elements[pos + 1].j;
                if (j_next != next_1[j])
                    d++;
            }
        }
        return d;
    }

    /*
     * Local search.
     */

    struct Move
    {
        Move(): type(Perturbations::None), global_cost(worst<GlobalCost>()) { }

        /** Type of perturbation. */
        Perturbations type;
        /** ForceAdd: element to add. */
        ElementId force_add_j = -1;
        /** Global cost of the move. */
        GlobalCost global_cost;
    };

    struct MoveHasher
    {
        std::hash<ElementPos> hasher;

        inline bool hashable(const Move& move) const
        {
            if (move.type == Perturbations::ForceAdd)
                return true;
            return false;
        }

        inline bool operator()(
                const Move& move_1,
                const Move& move_2) const
        {
            if (move_1.type == Perturbations::ForceAdd
                    && move_2.type == Perturbations::ForceAdd
                    && move_1.force_add_j == move_2.force_add_j)
                return true;
            return false;
        }

        inline std::size_t operator()(
                const Move& move) const
        {
            return hasher(move.force_add_j);
        }
    };

    inline MoveHasher move_hasher() const { return MoveHasher(); }

    /*
     * Perturbations.
     */

    inline std::vector<Move> perturbations(
            const Solution& solution,
            std::mt19937_64& generator)
    {
        std::vector<Move> moves;

        // Double-bridge.
        for (Counter perturbation = 0;
                perturbation < parameters_.double_bridge_number_of_perturbations;
                ++perturbation) {
            Move move;
            move.type = Perturbations::DoubleBridge;
            move.global_cost = global_cost(solution);
            moves.push_back(move);
        }

        // Ruin-and-recreate.
        for (Counter perturbation = 0;
                perturbation < parameters_.ruin_and_recreate_number_of_perturbations;
                ++perturbation) {
            Move move;
            move.type = Perturbations::RuinAndRecreate;
            move.global_cost = global_cost(solution);
            moves.push_back(move);
        }

        // Force-add.
        if (parameters_.force_add) {
            std::vector<uint8_t> contains(sequencing_scheme_.number_of_elements(), 0);
            for (SequenceId i = 0; i < number_of_sequences(); ++i)
                for (SequenceElement se: solution.sequences[i].elements)
                    contains[se.j] = 1;
            for (ElementId j = 0; j < sequencing_scheme_.number_of_elements(); ++j) {
                if (contains[j])
                    continue;
                Move move;
                move.type = Perturbations::ForceAdd;
                move.force_add_j = j;
                move.global_cost = global_cost(solution);
                moves.push_back(move);
            }
        }

        std::shuffle(moves.begin(), moves.end(), generator);
        return moves;
    }

    inline void apply_move(
            Solution& solution,
            Move& move,
            std::mt19937_64& generator)
    {
        SequenceId m = number_of_sequences();

        switch (move.type) {
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

        } case Perturbations::RuinAndRecreate: {
            //std::cout << "RuinAndRecreate" << std::endl;
            //std::cout << to_string(solution.global_cost) << std::endl;
            ElementId n = sequencing_scheme_.number_of_elements();
            std::vector<SequenceId> sequences(n, -1);
            std::vector<Mode> modes(n, -1);
            std::vector<ElementPos> positions(n, -1);
            optimizationtools::IndexedSet elts(n);
            for (SequenceId i = 0; i < m; ++i) {
                const auto& elements = solution.sequences[i].elements;
                ElementPos seq_size = (ElementPos)elements.size();
                for (ElementPos pos = 0; pos < seq_size; ++pos) {
                    ElementId j = elements[pos].j;
                    elts.add(j);
                    sequences[j] = i;
                    modes[j] = elements[pos].mode;
                    positions[j] = pos;
                }
            }

            Solution solution_cur_ = empty_solution();
            solution_cur_.modified_sequences = std::vector<bool>(m, false);

            ElementPos number_of_elements_removed = std::min(
                    parameters_.ruin_and_recreate_number_of_elements_removed,
                    elts.size());
            elts.shuffle_in(generator);
            while (elts.size() > number_of_elements_removed) {
                ElementId j = *(elts.begin() + (elts.size() - 1));
                elts.remove(j);
                solution_cur_.modified_sequences[sequences[j]] = true;
            }
            for (SequenceId i = 0; i < m; ++i) {
                const auto& elements = solution.sequences[i].elements;
                SequencePos seq_size = elements.size();
                Sequence& sequence_cur = solution_cur_.sequences[i];
                for (ElementPos pos = 0; pos < seq_size; ++pos) {
                    ElementId j = elements[pos].j;
                    if (!elts.contains(j))
                        append(sequence_cur, elements[pos]);
                }
            }
            // Add back removed elements.
            compute_global_cost(solution_cur_);
            compute_temporary_structures(solution_cur_);
            for (ElementId j: elts) {
                //std::cout << to_string(solution_cur_.global_cost) << std::endl;
                ElementId j_prec_old = -1;
                if (positions[j] > 0) {
                    const auto& elements = solution.sequences[sequences[j]].elements;
                    j_prec_old = elements[positions[j] - 1].j;
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
                    Move0 move_best;
                    move_best.global_cost = worst<GlobalCost>();
                    for (const Move0 move: improving_moves)
                        if (dominates(move.global_cost, move_best.global_cost))
                            move_best = move;
                    apply_move(solution_cur_, move_best);
                    solution_cur_.modified_sequences[move_best.i1] = true;
                    compute_temporary_structures(solution_cur_, move_best.i1, move_best.i2);
                }
            }
            solution = solution_cur_;
            //std::cout << to_string(solution.global_cost) << std::endl;
            break;

        } case Perturbations::ForceAdd: {
            ElementId j = move.force_add_j;
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
            Move0 move;
            move.type = Neighborhoods::Add;
            move.i1 = i;
            move.j = j;
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
            const Move& perturbation = Move())
    {
        //std::cout << "local_search" << std::endl;
        //if (tabu.j != -1)
        //    std::cout << "j " << tabu.j << " j_prev " << tabu.j_prev
        //        << std::endl;
        //print(std::cout, solution);
        //std::cout << to_string(global_cost(solution)) << std::endl;

        auto local_search_begin = std::chrono::steady_clock::now();
        SequencePos m = number_of_sequences();
        //ElementPos n = sequencing_scheme_.number_of_elements();

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
            if (parameters_.inter_two_opt)
                neighborhoods.push_back({Neighborhoods::InterTwoOpt, 0, 0});
            if (parameters_.inter_two_opt_reverse)
                neighborhoods.push_back({Neighborhoods::InterTwoOptReverse, 0, 0});
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
                        for (SequenceId i = 0; i < m; ++i)
                            neighborhood.modified_sequences[i] = solution.modified_sequences[i];
                    }
                }
            }
        }

        compute_temporary_structures(solution);

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
                //for (SequenceId i = 0; i < m; ++i)
                //    std::cout << " " << neighborhood.modified_sequences[i];
                //std::cout << std::endl;

                // Remove moves which have changed from improving_moves.
                for (auto it = neighborhood.improving_moves.begin();
                        it != neighborhood.improving_moves.end();) {
                    if (neighborhood.modified_sequences[it->i1]
                            || (it->i2 != -1 && neighborhood.modified_sequences[it->i2])) {
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
                case Neighborhoods::InterTwoOpt:
                    explore_inter_two_opt(solution);
                    break;
                case Neighborhoods::InterTwoOptReverse:
                    explore_inter_two_opt_reverse(solution);
                    break;
                case Neighborhoods::InterShift:
                    explore_inter_shift(solution, k1);
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
                    Move0 move_best;
                    move_best.global_cost = worst<GlobalCost>();
                    for (const Move0 move: neighborhood.improving_moves)
                        if (dominates(move.global_cost, move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    compute_temporary_structures(solution, move_best.i1, move_best.i2);
                    // Check new current solution cost.
                    if (global_cost(solution) - gc_old != move_best.global_cost) {
                        //std::cout
                        //    << "k1 " << move_best.k1
                        //    << " k2 " << move_best.k2
                        //    << " i1 " << move_best.i1
                        //    << " i2 " << move_best.i2
                        //    << " pos_1 " << move_best.pos_1
                        //    << " pos_2 " << move_best.pos_2
                        //    << std::endl;
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

                    // Update modified sequences.
                    for (int a = 0; a < (int)neighborhoods_.size(); ++a) {
                        for (int k1 = 0; k1 < (int)neighborhoods_[a].size(); ++k1) {
                            for (int k2 = 0; k2 < (int)neighborhoods_[a][k1].size(); ++k2) {
                                Neighborhood& neighborhood_2 = neighborhoods_[a][k1][k2];
                                neighborhood_2.modified_sequences[move_best.i1] = true;
                                if (move_best.i2 != -1)
                                    neighborhood_2.modified_sequences[move_best.i2] = true;
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
        local_search_time += time_span.count();
        number_of_local_search_calls++;
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
                    os << " " << se.j;
                } else {
                    os << " " << se.j << " " << se.mode;
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
                    file << se.j + offset << " ";
                } else {
                    file << se.j + offset << " " << se.mode << " ";
                }
            }
        } else {
            for (SequenceId i = 0; i < m; ++i) {
                file << solution.sequences[i].elements.size() << std::endl;
                for (SequenceElement se: solution.sequences[i].elements) {
                    if (maximum_number_of_modes_ == 1) {
                        file << se.j + offset << " ";
                    } else {
                        file << se.j + offset << " " << se.mode << " ";
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
                << "    Inter-two-opt:                             " << parameters_.inter_two_opt << std::endl
                << "    Inter-two-opt-reverse:                     " << parameters_.inter_two_opt_reverse << std::endl
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
            << "        Number of elements removed:            " << parameters_.ruin_and_recreate_number_of_elements_removed << std::endl
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
            << number_of_initial_solution_calls
            << " / " << initial_solution_time << "s"
            << " / " << initial_solution_time / number_of_initial_solution_calls << "s"
            << std::endl;
        info.os()
            << "    "
            << std::left << std::setw(28) << "Crossover time:"
            << number_of_crossover_calls
            << " / " << crossover_time << "s"
            << " / " << crossover_time / number_of_crossover_calls << "s"
            << std::endl;
        info.os()
            << "    "
            << std::left << std::setw(28) << "Local search time:"
            << number_of_local_search_calls
            << " / " << local_search_time << "s"
            << " / " << local_search_time / number_of_local_search_calls << "s"
            << std::endl;
        info.os()
            << "    "
            << std::left << std::setw(28) << "Apply search time:"
            << number_of_apply_move_calls
            << " / " << apply_move_time << "s"
            << " / " << apply_move_time / number_of_apply_move_calls << "s"
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

    void append(
            SequenceData& sequence_data,
            SequenceElement se,
            std::true_type,
            std::false_type,
            std::false_type,
            std::false_type) const
    {
        sequencing_scheme_.append(sequence_data, se.j, se.mode);;
    }

    void append(
            SequenceData& sequence_data,
            SequenceElement se,
            std::false_type,
            std::true_type,
            std::false_type,
            std::false_type) const
    {
        sequencing_scheme_.append(sequence_data, se.j);;
    }

    void append(
            SequenceData& sequence_data,
            SequenceElement se,
            std::false_type,
            std::false_type,
            std::true_type,
            std::false_type) const
    {
        if (sequence_data.number_of_locations == 0) {
            sequence_data = sequencing_scheme_.sequence_data_init(se.j, se.mode);
        } else {
            sequencing_scheme_.concatenate(sequence_data, sequencing_scheme_.sequence_data_init(se.j, se.mode));
        }
    }

    void append(
            SequenceData& sequence_data,
            SequenceElement se,
            std::false_type,
            std::false_type,
            std::false_type,
            std::true_type) const
    {
        if (sequence_data.number_of_locations == 0) {
            sequence_data = sequencing_scheme_.sequence_data_init(se.j);
        } else {
            sequencing_scheme_.concatenate(
                    sequence_data,
                    sequencing_scheme_.sequence_data_init(se.j));
        }
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
                    HasAppendMethod<SequencingScheme, void(SequenceData&, ElementId, Mode)>::value>(),
                std::integral_constant<
                    bool,
                    HasAppendMethod<SequencingScheme, void(SequenceData&, ElementId)>::value>(),
                std::integral_constant<
                    bool,
                    HasSequenceDataInitMethod<SequencingScheme, SequenceData(ElementId, Mode)>::value>(),
                std::integral_constant<
                    bool,
                    HasSequenceDataInitMethod<SequencingScheme, SequenceData(ElementId)>::value>());
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
     * global_cost_merge(
     *         const GlobalCost&,
     *         const GlobalCost&)
     */

    template<typename, typename T>
    struct HasGlobalCostMergeMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasGlobalCostMergeMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().global_cost_merge(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    GlobalCost global_cost_merge(
            const GlobalCost& gc1,
            const GlobalCost& gc2,
            std::false_type) const
    {
        return gc1 + gc2;
    }

    GlobalCost global_cost_merge(
            const GlobalCost& gc1,
            const GlobalCost& gc2,
            std::true_type) const
    {
        return sequencing_scheme_.global_cost_merge(gc1, gc2);
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
    GlobalCost global_cost_merge(
            const GlobalCost& gc1,
            const GlobalCost& gc2) const
    {
        return global_cost_merge(
                gc1,
                gc2,
                std::integral_constant<
                    bool,
                    HasGlobalCostMergeMethod<SequencingScheme, GlobalCost(
                        const GlobalCost&,
                        const GlobalCost&)>::value>());
    }

    /*
     * bool concatenate(
     *         SequenceData&,
     *         const SequenceData&)
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

    inline const SequenceData& get_sequence_data(const SubSequence& sub_sequence) const
    {
        return (!sub_sequence.reverse)?
            sequence_datas_cur_2_[sub_sequence.sequence->i][sub_sequence.pos_1][sub_sequence.pos_2]:
            sequence_datas_cur_2_[sub_sequence.sequence->i][sub_sequence.pos_2][sub_sequence.pos_1];
    }

    inline GlobalCost global_cost_concatenate_1(
            const SubSequence* sub_sequences,
            ElementPos size,
            const GlobalCost* gci,
            std::false_type) const
    {
        SequenceData sequence_data;
        ElementPos p = 0;
        for (;;) {
            if (sub_sequences[p].sequence != nullptr) {
                if (sub_sequences[p].pos_1 <= sub_sequences[p].pos_2) {
                    sequence_data = (!sub_sequences[p].reverse)?
                        sequence_datas_cur_2_[sub_sequences[p].sequence->i][sub_sequences[p].pos_1][sub_sequences[p].pos_2]:
                        sequence_datas_cur_2_[sub_sequences[p].sequence->i][sub_sequences[p].pos_2][sub_sequences[p].pos_1];
                    break;
                }
            } else {
                sequence_data = empty_sequence_data(sub_sequences[0].sequence->i);
                append(sequence_data, {sub_sequences[p].j, sub_sequences[p].mode});
                break;
            }
            p++;
            if (p == size) {
                sequence_data = empty_sequence_data(sub_sequences[0].sequence->i);
                return (gci == nullptr)?
                    sequencing_scheme_.global_cost(sequence_data):
                    global_cost_merge(sequencing_scheme_.global_cost(sequence_data), *gci);
            }
        }
        p++;
        for (; p < size; ++p) {
            if (sub_sequences[p].sequence != nullptr) {
                if (sub_sequences[p].pos_1 > sub_sequences[p].pos_2)
                    continue;
                if (!sub_sequences[p].reverse) {
                    auto it = sub_sequences[p].sequence->elements.begin() + sub_sequences[p].pos_1;
                    const auto it_end = sub_sequences[p].sequence->elements.begin() + sub_sequences[p].pos_2;
                    for (;; ++it) {
                        // Add next element to sequence_data.
                        append(sequence_data, *it);
                        // End condition.
                        if (it == it_end)
                            break;
                    }
                } else {
                    auto it = sub_sequences[p].sequence->elements.begin() + sub_sequences[p].pos_2;
                    const auto it_end = sub_sequences[p].sequence->elements.begin() + sub_sequences[p].pos_1;
                    for (;; --it) {
                        // Add next element to sequence_data.
                        append(sequence_data, *it);
                        // End condition.
                        if (it == it_end)
                            break;
                    }
                }
            } else {
                append(sequence_data, {sub_sequences[p].j, sub_sequences[p].mode});
            }
        }
        return (gci == nullptr)?
            sequencing_scheme_.global_cost(sequence_data):
            global_cost_merge(sequencing_scheme_.global_cost(sequence_data), *gci);
    }

    inline std::pair<bool, GlobalCost> global_cost_concatenate_1(
            const SubSequence* sub_sequences,
            ElementPos size,
            GlobalCost cutoff,
            const GlobalCost* gci,
            std::false_type) const
    {
        //std::cout << "global_cost_concatenate_1"
        //    << " size " << size
        //    << " gci " << gci
        //    << std::endl;
        SequenceData sequence_data;
        ElementPos p = 0;
        for (;;) {
            if (sub_sequences[p].sequence != nullptr) {
                if (sub_sequences[p].pos_1 <= sub_sequences[p].pos_2) {
                    sequence_data = (!sub_sequences[p].reverse)?
                        sequence_datas_cur_2_[sub_sequences[p].sequence->i][sub_sequences[p].pos_1][sub_sequences[p].pos_2]:
                        sequence_datas_cur_2_[sub_sequences[p].sequence->i][sub_sequences[p].pos_2][sub_sequences[p].pos_1];
                    break;
                }
            } else {
                sequence_data = empty_sequence_data(sub_sequences[0].sequence->i);
                append(sequence_data, {sub_sequences[p].j, sub_sequences[p].mode});
                break;
            }
            p++;
            if (p == size) {
                sequence_data = empty_sequence_data(sub_sequences[0].sequence->i);
                return {true, (gci == nullptr)?
                    sequencing_scheme_.global_cost(sequence_data):
                    global_cost_merge(sequencing_scheme_.global_cost(sequence_data), *gci)};
            }
        }
        p++;
        for (; p < size; ++p) {
            if (sub_sequences[p].sequence != nullptr) {
                //std::cout << "p " << p
                //    << " seq " << sub_sequences[p].sequence
                //    << " size " << sub_sequences[p].sequence->elements.size()
                //    << " i " << sub_sequences[p].sequence->i
                //    << " pos_1 " << sub_sequences[p].pos_1
                //    << " pos_2 " << sub_sequences[p].pos_2
                //    << " rev " << sub_sequences[p].reverse
                //    << std::endl;
                if (sub_sequences[p].pos_1 > sub_sequences[p].pos_2)
                    continue;
                if (!sub_sequences[p].reverse) {
                    auto it = sub_sequences[p].sequence->elements.begin() + sub_sequences[p].pos_1;
                    const auto it_end = sub_sequences[p].sequence->elements.begin() + sub_sequences[p].pos_2;
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
                            if (bound(sequence_data) >= cutoff)
                                return {false, GlobalCost()};
                        } else {
                            if (global_cost_merge(*gci, bound(sequence_data)) >= cutoff)
                                return {false, GlobalCost()};
                        }
                        // End condition.
                        if (it == it_end)
                            break;
                    }
                } else {
                    auto it = sub_sequences[p].sequence->elements.begin() + sub_sequences[p].pos_2;
                    const auto it_end = sub_sequences[p].sequence->elements.begin() + sub_sequences[p].pos_1;
                    for (;; --it) {
                        // Add next element to sequence_data.
                        append(sequence_data, *it);
                        //std::cout << "j " << it->j
                        //    << " gc " << to_string(global_cost(sequence_data))
                        //    << " bnd " << to_string(bound(sequence_data))
                        //    << std::endl;
                        // Check early termination.
                        if (gci == nullptr) {
                            if (bound(sequence_data) >= cutoff)
                                return {false, GlobalCost()};
                        } else {
                            if (global_cost_merge(*gci, bound(sequence_data)) >= cutoff)
                                return {false, GlobalCost()};
                        }
                        // End condition.
                        if (it == it_end)
                            break;
                    }
                }
            } else {
                append(sequence_data, {sub_sequences[p].j, sub_sequences[p].mode});
                //std::cout << "j " << sub_sequences[p].j
                //    << " gc " << to_string(global_cost(sequence_data))
                //    << " bnd " << to_string(bound(sequence_data))
                //    << std::endl;
                // Check early termination.
                if (gci == nullptr) {
                    if (bound(sequence_data) >= cutoff)
                        return {false, GlobalCost()};
                } else {
                    if (global_cost_merge(*gci, bound(sequence_data)) >= cutoff)
                        return {false, GlobalCost()};
                }
            }
        }
        return {true, (gci == nullptr)?
            sequencing_scheme_.global_cost(sequence_data):
            global_cost_merge(sequencing_scheme_.global_cost(sequence_data), *gci)};
    }

    inline GlobalCost global_cost_concatenate_1(
            const SubSequence* sub_sequences,
            ElementPos size,
            const GlobalCost* gci,
            std::true_type) const
    {
        SequenceData sequence_data;
        ElementPos p = 0;
        for (;;) {
            if (sub_sequences[p].sequence != nullptr) {
                if (sub_sequences[p].pos_1 <= sub_sequences[p].pos_2) {
                    sequence_data = (!sub_sequences[p].reverse)?
                        sequence_datas_cur_2_[sub_sequences[p].sequence->i][sub_sequences[p].pos_1][sub_sequences[p].pos_2]:
                        sequence_datas_cur_2_[sub_sequences[p].sequence->i][sub_sequences[p].pos_2][sub_sequences[p].pos_1];
                    break;
                }
            } else {
                sequence_data = empty_sequence_data(sub_sequences[0].sequence->i);
                append(sequence_data, {sub_sequences[p].j, sub_sequences[p].mode});
                break;
            }
            p++;
            if (p == size) {
                sequence_data = empty_sequence_data(sub_sequences[0].sequence->i);
                return (gci == nullptr)?
                    sequencing_scheme_.global_cost(sequence_data):
                    global_cost_merge(sequencing_scheme_.global_cost(sequence_data), *gci);
            }
        }
        p++;
        for (; p < size; ++p) {
            if (sub_sequences[p].sequence != nullptr) {
                if (sub_sequences[p].pos_1 <= sub_sequences[p].pos_2)
                    sequencing_scheme_.concatenate(sequence_data, get_sequence_data(sub_sequences[p]));
            } else {
                append(sequence_data, {sub_sequences[p].j, sub_sequences[p].mode});
            }
        }
        return (gci == nullptr)?
            sequencing_scheme_.global_cost(sequence_data):
            global_cost_merge(sequencing_scheme_.global_cost(sequence_data), *gci);
    }

    inline std::pair<bool, GlobalCost> global_cost_concatenate_1(
            const SubSequence* sub_sequences,
            ElementPos size,
            GlobalCost,
            const GlobalCost* gci,
            std::true_type) const
    {
        //std::cout << "global_cost_concatenate_1"
        //    << " size " << size
        //    << " gci " << gci
        //    << std::endl;
        SequenceData sequence_data;
        ElementPos p = 0;
        for (;;) {
            if (sub_sequences[p].sequence != nullptr) {
                if (sub_sequences[p].pos_1 <= sub_sequences[p].pos_2) {
                    sequence_data = (!sub_sequences[p].reverse)?
                        sequence_datas_cur_2_[sub_sequences[p].sequence->i][sub_sequences[p].pos_1][sub_sequences[p].pos_2]:
                        sequence_datas_cur_2_[sub_sequences[p].sequence->i][sub_sequences[p].pos_2][sub_sequences[p].pos_1];
                    break;
                }
            } else {
                sequence_data = empty_sequence_data(sub_sequences[0].sequence->i);
                append(sequence_data, {sub_sequences[p].j, sub_sequences[p].mode});
                break;
            }
            p++;
            if (p == size) {
                sequence_data = empty_sequence_data(sub_sequences[0].sequence->i);
                return {true, (gci == nullptr)?
                    sequencing_scheme_.global_cost(sequence_data):
                    global_cost_merge(sequencing_scheme_.global_cost(sequence_data), *gci)};
            }
        }
        p++;
        for (; p < size; ++p) {
            if (sub_sequences[p].sequence != nullptr) {
                if (sub_sequences[p].pos_1 <= sub_sequences[p].pos_2)
                    sequencing_scheme_.concatenate(sequence_data, get_sequence_data(sub_sequences[p]));
            } else {
                append(sequence_data, {sub_sequences[p].j, sub_sequences[p].mode});
            }
        }
        return {true, (gci == nullptr)?
            sequencing_scheme_.global_cost(sequence_data):
            global_cost_merge(sequencing_scheme_.global_cost(sequence_data), *gci)};
    }

    inline GlobalCost global_cost_concatenate_1(
            const SubSequence* sub_sequences,
            ElementPos size,
            const GlobalCost* gci) const
    {
        return global_cost_concatenate_1(
                sub_sequences,
                size,
                gci,
                std::integral_constant<
                    bool,
                    HasConcatenateMethod<SequencingScheme, bool(SequenceData&, const SequenceData&)>::value>());
    }

    inline std::pair<bool, GlobalCost> global_cost_concatenate_1(
            const SubSequence* sub_sequences,
            ElementPos size,
            GlobalCost cutoff,
            GlobalCost* gci) const
    {
        return global_cost_concatenate_1(
                sub_sequences,
                size,
                cutoff,
                gci,
                std::integral_constant<
                    bool,
                    HasConcatenateMethod<SequencingScheme, bool(SequenceData&, const SequenceData&)>::value>());
    }

    inline std::pair<bool, GlobalCost> global_cost_concatenate_2(
            const SubSequence* sub_sequences_1,
            ElementPos size_1,
            const SubSequence* sub_sequences_2,
            ElementPos size_2,
            GlobalCost cutoff,
            GlobalCost* gci) const
    {
        auto bgc = global_cost_concatenate_1(
                sub_sequences_1, size_1, cutoff, gci);
        if (!bgc.first)
            return bgc;
        return global_cost_concatenate_1(
                sub_sequences_2, size_2, cutoff, &bgc.second);
    }

    /*
     * void compute_closest_neighbors()
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

    void compute_closest_neighbors(
            std::false_type)
    {
        ElementPos n = sequencing_scheme_.number_of_elements();
        neighbors_ = std::vector<ElementId>(n);
        std::iota(neighbors_.begin(), neighbors_.end(), 0);
    }

    void compute_closest_neighbors(
            std::true_type)
    {
        ElementPos n = sequencing_scheme_.number_of_elements();
        closest_neighbors_ = std::vector<std::vector<ElementId>>(n);
        std::vector<ElementId> neighbors(n);
        std::iota(neighbors.begin(), neighbors.end(), 0);
        for (ElementId j = 0; j < n; ++j) {
            std::nth_element(
                    neighbors.begin(),
                    neighbors.begin() + neighborhood_size_,
                    neighbors.end(),
                    [this, j](ElementId j1, ElementId j2) -> bool
                    {
                        return sequencing_scheme_.distance(j, j1)
                            < sequencing_scheme_.distance(j, j2);
                    });
            closest_neighbors_[j].insert(
                    closest_neighbors_[j].begin(),
                    neighbors.begin(),
                    neighbors.begin() + neighborhood_size_);

            //std::cout << "Closest neighbors of " << j << ":";
            //for (ElementId j2: closest_neighbors_[j])
            //    std::cout << " " << j2 << "," << sequencing_scheme_.distance(j, j2);
            //std::cout << std::endl;
            //std::cout << "Other neighbors of " << j << ":";
            //for (auto it = neighbors.begin() + 100; it != neighbors.end(); ++it)
            //    std::cout << " " << *it << "," << sequencing_scheme_.distance(j, *it);
        }
    }

    void compute_closest_neighbors()
    {
        ElementPos n = sequencing_scheme_.number_of_elements();
        if (n <= neighborhood_size_) {
            compute_closest_neighbors(std::false_type{});
        } else {
            compute_closest_neighbors(
                    std::integral_constant<
                    bool,
                    HasDistanceMethod<SequencingScheme, double(ElementId, ElementId)>::value>());
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
            SequenceId i2 = -2)
    {
        SequencePos m = number_of_sequences();
        // Update global_costs_cur_.
        for (SequenceId i = 0; i < m; ++i)
            global_costs_cur_[i] = sequencing_scheme_.global_cost(solution.sequences[i].data);
        // Update elements_cur_.
        std::fill(elements_cur_.begin(), elements_cur_.end(), SolutionElement());
        for (SequenceId i = 0; i < m; ++i) {
            const auto& elements = solution.sequences[i].elements;
            ElementPos seq_size = elements.size();
            for (ElementPos pos = 0; pos < seq_size; ++pos) {
                ElementId j = elements[pos].j;
                elements_cur_[j].i = i;
                elements_cur_[j].pos = pos;
                elements_cur_[j].mode = elements[pos].mode;
            }
        }
        // Update sequence_datas_cur_2_.
        for (SequenceId i = 0; i < m; ++i) {
            if (i1 != -1 && i != i1 && i != i2)
                continue;
            const auto& elements = solution.sequences[i].elements;
            ElementPos seq_size = (ElementPos)elements.size();
            while ((ElementPos)sequence_datas_cur_2_[i].size()  < seq_size)
                sequence_datas_cur_2_[i].push_back({});
            for (ElementPos pos_1 = 0; pos_1 < seq_size; ++pos_1) {
                while ((ElementPos)sequence_datas_cur_2_[i][pos_1].size() < seq_size)
                    sequence_datas_cur_2_[i][pos_1].push_back({});
                sequence_datas_cur_2_[i][pos_1][pos_1] = empty_sequence_data(i);
                append(sequence_datas_cur_2_[i][pos_1][pos_1], elements[pos_1]);
            }
            for (ElementPos size = 2; size <= seq_size; ++size) {
                for (ElementPos pos_1 = 0; pos_1 + size - 1 < seq_size; ++pos_1) {
                    ElementPos pos_2 = pos_1 + size - 1;
                    sequence_datas_cur_2_[i][pos_1][pos_2] = sequence_datas_cur_2_[i][pos_1][pos_2 - 1];
                    append(sequence_datas_cur_2_[i][pos_1][pos_2], elements[pos_2]);
                    sequence_datas_cur_2_[i][pos_2][pos_1] = sequence_datas_cur_2_[i][pos_2][pos_1 + 1];
                    append(sequence_datas_cur_2_[i][pos_2][pos_1], elements[pos_1]);
                }
            }
        }
        // Update partial_global_costs_1_.
        if (m > 1) {
            for (SequenceId i = 0; i < m; ++i) {
                partial_global_costs_cur_1_[i]
                    = compute_partial_global_cost_1(i);
            }
        }
        // Update partial_global_costs_2_.
        // TODO optimize.
        if (m > 2) {
            for (SequenceId i1 = 0; i1 < m; ++i1) {
                for (SequenceId i2 = i1 + 1; i2 < m; ++i2) {
                    partial_global_costs_cur_2_[i1][i2]
                        = compute_partial_global_cost_2(i1, i2);
                    partial_global_costs_cur_2_[i2][i1]
                        = partial_global_costs_cur_2_[i1][i2];
                }
            }
        }
    }

    /**
     * Compute and return the global cost of a solution.
     */
    inline void compute_global_cost(
            Solution& solution) const
    {
        solution.global_cost = sequencing_scheme_.global_cost(solution.sequences[0].data);
        for (SequenceId i = 1; i < number_of_sequences(); ++i) {
            solution.global_cost = global_cost_merge(
                    solution.global_cost,
                    sequencing_scheme_.global_cost(solution.sequences[i].data));
        }
    }

    /**
     * Compute and return the global cost of a solution without sequence
     * 'i_out'.
     */
    inline GlobalCost compute_partial_global_cost_1(
            SequenceId i_out) const
    {
        if (number_of_sequences() == 1)
            return GlobalCost();
        if (i_out != 0) {
            GlobalCost gc = global_costs_cur_[0];
            for (SequenceId i = 1; i < number_of_sequences(); ++i)
                if (i != i_out)
                    gc = global_cost_merge(gc, global_costs_cur_[i]);
            return gc;
        } else {
            GlobalCost gc = global_costs_cur_[1];
            for (SequenceId i = 2; i < number_of_sequences(); ++i)
                gc = global_cost_merge(gc, global_costs_cur_[i]);
            return gc;
        }
    }

    /**
     * Compute and return the global cost of a solution without sequences
     * 'i_out_1' and 'i_out_2'.
     */
    inline GlobalCost compute_partial_global_cost_2(
            SequenceId i_out_1,
            SequenceId i_out_2) const
    {
        if (number_of_sequences() == 1) {
            throw std::logic_error(
                    "'compute_global_cost(solution, i_out_1, i_out_2)'"
                    " requires at least 2 sequences.");
        }
        if (number_of_sequences() == 2)
            return GlobalCost();
        if (i_out_1 != 0 && i_out_2 != 0) {
            GlobalCost gc = global_costs_cur_[0];
            for (SequenceId i = 1; i < number_of_sequences(); ++i)
                if (i != i_out_1 && i != i_out_2)
                    gc = global_cost_merge(gc, global_costs_cur_[i]);
            return gc;
        } else if (i_out_1 != 1 && i_out_2 != 1) {
            GlobalCost gc = global_costs_cur_[1];
            for (SequenceId i = 2; i < number_of_sequences(); ++i)
                if (i != i_out_1 && i != i_out_2)
                    gc = global_cost_merge(gc, global_costs_cur_[i]);
            return gc;
        } else {
            GlobalCost gc = global_costs_cur_[2];
            for (SequenceId i = 3; i < number_of_sequences(); ++i)
                gc = global_cost_merge(gc, global_costs_cur_[i]);
            return gc;
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

        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();

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
                    std::pair<bool, GlobalCost> bgc_tmp;
                    if (block_pos > pos_new) {
                        SubSequence sub_sequences[] = {
                            SubSequence(sequence, 0, pos_new - 1),
                            SubSequence(sequence, block_pos, block_pos + block_size - 1, reverse),
                            SubSequence(sequence, pos_new, block_pos - 1),
                            SubSequence(sequence, block_pos + block_size, seq_size - 1)};
                        bgc_tmp = global_cost_concatenate_1(
                                sub_sequences, 4, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                    } else {
                        SubSequence sub_sequences[] = {
                            SubSequence(sequence, 0, block_pos - 1),
                            SubSequence(sequence, block_pos + block_size, pos_new + block_size - 1),
                            SubSequence(sequence, block_pos, block_pos + block_size - 1, reverse),
                            SubSequence(sequence, pos_new + block_size, seq_size - 1)};
                        bgc_tmp = global_cost_concatenate_1(
                                sub_sequences, 4, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                    }
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        Move0 move;
                        move.type = (!reverse)?
                            Neighborhoods::Shift:
                            Neighborhoods::ShiftReverse;
                        move.k1 = block_size;
                        move.i1 = i;
                        move.pos_1 = block_pos;
                        move.pos_2 = pos_new;
                        move.global_cost = bgc_tmp.second - gc;
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
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();

            // Loop through all pairs.
            for (ElementPos pos = 0; pos < seq_size; ++pos) {

                // block 1 is at [pos, pos + block_size_1[.
                ElementPos pos_1 = pos;
                ElementPos pos_2_max = std::min(
                        seq_size - block_size_2,
                        pos + block_size_1 + parameters_.swap_maximum_distance);
                for (ElementPos pos_2 = pos + block_size_1; pos_2 <= pos_2_max; ++pos_2) {
                    SubSequence sub_sequences[] = {
                        SubSequence(sequence, 0, pos - 1),
                        SubSequence(sequence, pos_2, pos_2 + block_size_2 - 1),
                        SubSequence(sequence, pos_1 + block_size_1, pos_2 - 1),
                        SubSequence(sequence, pos_1, pos_1 + block_size_1 - 1),
                        SubSequence(sequence, pos_2 + block_size_2, seq_size - 1)};
                    auto bgc_tmp = global_cost_concatenate_1(
                            sub_sequences, 5, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        //std::cout << to_string(gc_tmp)
                        //    << " " << to_string(gc)
                        //    << std::endl;
                        Move0 move;
                        move.type = Neighborhoods::Swap;
                        move.k1 = block_size_1;
                        move.k2 = block_size_2;
                        move.i1 = i;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = bgc_tmp.second - gc;
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
                    SubSequence sub_sequences[] = {
                        SubSequence(sequence, 0, pos - 1),
                        SubSequence(sequence, pos_1, pos_1 + block_size_1 - 1),
                        SubSequence(sequence, pos_2 + block_size_2, pos_1 - 1),
                        SubSequence(sequence, pos_2, pos_2 + block_size_2 - 1),
                        SubSequence(sequence, pos_1 + block_size_1, seq_size - 1)};
                    auto bgc_tmp = global_cost_concatenate_1(
                            sub_sequences, 5, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        Move0 move;
                        move.type = Neighborhoods::Swap;
                        move.k1 = block_size_1;
                        move.k2 = block_size_2;
                        move.i1 = i;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = bgc_tmp.second - gc;
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
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();

            // Loop through all pairs.
            for (ElementPos pos_1 = 0; pos_1 < seq_size; ++pos_1) {
                ElementPos pos_max = std::min(
                        seq_size,
                        pos_1 + parameters_.reverse_maximum_length);
                for (ElementPos pos_2 = pos_1 + 2; pos_2 < pos_max; ++pos_2) {
                    SubSequence sub_sequences[] = {
                        SubSequence(sequence, 0, pos_1 - 1),
                        SubSequence(sequence, pos_1, pos_2, true),
                        SubSequence(sequence, pos_2 + 1, seq_size - 1)};
                    auto bgc_tmp = global_cost_concatenate_1(
                            sub_sequences, 3, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        Move0 move;
                        move.type = Neighborhoods::Reverse;
                        move.i1 = i;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = bgc_tmp.second - gc;
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

        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();

            // Loop through all new positions.
            for (ElementPos pos = 0; pos <= seq_size; ++pos) {

                for (ElementId j = 0; j < n; ++j) {
                    if (elements_cur_[j].mode != -1)
                        continue;

                    for (Mode mode = 0; mode < number_of_modes(j); ++mode) {
                        SubSequence sub_sequences[] = {
                            SubSequence(sequence, 0, pos - 1),
                            SubSequence(j, mode),
                            SubSequence(sequence, pos, seq_size - 1)};
                        auto bgc_tmp = global_cost_concatenate_1(
                                sub_sequences, 3, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                        if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                            Move0 move;
                            move.type = Neighborhoods::Add;
                            move.i1 = i;
                            move.j = j;
                            move.mode = mode;
                            move.pos_1 = pos;
                            move.global_cost = bgc_tmp.second - gc;
                            neighborhood.improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    inline void explore_remove(
            const Solution& solution,
            const Move& perturbation)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Remove)][0][0];
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();

            // Loop through all new positions.
            for (ElementPos pos = 0; pos < seq_size; ++pos) {

                ElementId j = sequence.elements[pos].j;
                if (perturbation.type == Perturbations::ForceAdd
                        && j == perturbation.force_add_j)
                    continue;

                SubSequence sub_sequences[] = {
                    SubSequence(sequence, 0, pos - 1),
                    SubSequence(sequence, pos + 1, seq_size - 1)};
                auto bgc_tmp = global_cost_concatenate_1(
                        sub_sequences, 2, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                    Move0 move;
                    move.type = Neighborhoods::Remove;
                    move.i1 = i;
                    move.pos_1 = pos;
                    move.global_cost = bgc_tmp.second - gc;
                    neighborhood.improving_moves.push_back(move);
                }
            }
        }
    }

    inline void explore_replace(
            const Solution& solution,
            const Move& perturbation)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Add)][0][0];

        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();

            // Loop through all new positions.
            for (ElementPos pos = 0; pos < seq_size; ++pos) {

                if (perturbation.type == Perturbations::ForceAdd
                        && sequence.elements[pos].j == perturbation.force_add_j)
                    continue;

                for (ElementId j = 0; j < n; ++j) {
                    if (elements_cur_[j].mode != -1)
                        continue;

                    for (Mode mode = 0; mode < number_of_modes(j); ++mode) {
                        SubSequence sub_sequences[] = {
                            SubSequence(sequence, 0, pos - 1),
                            SubSequence(j, mode),
                            SubSequence(sequence, pos + 1, seq_size - 1)};
                        auto bgc_tmp = global_cost_concatenate_1(
                                sub_sequences, 3, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                        if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                            Move0 move;
                            move.type = Neighborhoods::Replace;
                            move.i1 = i;
                            move.j = j;
                            move.mode = mode;
                            move.pos_1 = pos;
                            move.global_cost = bgc_tmp.second - gc;
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
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = (!reverse)?
            neighborhoods_[int(Neighborhoods::InterShift)][block_size][0]:
            neighborhoods_[int(Neighborhoods::InterShiftReverse)][block_size][0];

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (ElementPos pos_1 = 0; pos_1 < seq_1_size; ++pos_1) {
                if (pos_1 + block_size - 1 >= seq_1_size)
                    continue;
                ElementId j1 = sequence_1.elements[pos_1].j;
                const auto& neighbors = (closest_neighbors_.empty())? neighbors_: closest_neighbors_[j1];

                SubSequence sub_sequences[] = {
                    SubSequence(sequence_1, 0, pos_1 - 1),
                    SubSequence(sequence_1, pos_1 + block_size, seq_1_size - 1)};
                auto bgc1_tmp = global_cost_concatenate_1(
                        sub_sequences, 2, gc, nullptr);
                if (!bgc1_tmp.first)
                    continue;

                for (ElementId j2: neighbors) {
                    SequenceId i2 = elements_cur_[j2].i;
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

                    GlobalCost gci = (m == 2)? bgc1_tmp.second:
                        global_cost_merge(
                            bgc1_tmp.second,
                            partial_global_costs_cur_2_[i1][i2]);
                    SubSequence sub_sequences[] = {
                        SubSequence(sequence_2, 0, pos_2 - 1),
                        SubSequence(sequence_1, pos_1, pos_1 + block_size - 1, reverse),
                        SubSequence(sequence_2, pos_2, seq_2_size - 1)};
                    auto bgc_tmp = global_cost_concatenate_1(
                            sub_sequences, 3, gc, &gci);
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        Move0 move;
                        move.type = (!reverse)?
                            Neighborhoods::InterShift:
                            Neighborhoods::InterShiftReverse;
                        move.k1 = block_size;
                        move.i1 = i1;
                        move.i2 = i2;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = bgc_tmp.second - gc;
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_inter_two_opt(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::InterTwoOpt)][0][0];

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (ElementPos pos_1 = 0; pos_1 < seq_1_size; ++pos_1) {
                ElementId j1 = sequence_1.elements[pos_1].j;
                const auto& neighbors = (closest_neighbors_.empty())? neighbors_: closest_neighbors_[j1];

                for (ElementId j2: neighbors) {
                    SequenceId i2 = elements_cur_[j2].i;
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
                    if (j1 > j2)
                        continue;

                    const auto& sequence_2 = solution.sequences[i2];
                    SequencePos seq_2_size = sequence_2.elements.size();

                    SubSequence sub_sequences_1[] = {
                        SubSequence(sequence_1, 0, pos_1),
                        SubSequence(sequence_2, pos_2 + 1, seq_2_size - 1)};
                    SubSequence sub_sequences_2[] = {
                        SubSequence(sequence_2, 0, pos_2),
                        SubSequence(sequence_1, pos_1 + 1, seq_1_size - 1)};
                    auto bgc_tmp = global_cost_concatenate_2(
                            sub_sequences_1, 2, sub_sequences_2, 2,
                            gc, ((m > 1)? &partial_global_costs_cur_2_[i1][i2]: nullptr));
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        Move0 move;
                        move.type = Neighborhoods::InterTwoOpt;
                        move.i1 = i1;
                        move.i2 = i2;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = bgc_tmp.second - gc;
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_inter_two_opt_reverse(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::InterTwoOptReverse)][0][0];

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (ElementPos pos_1 = 0; pos_1 < seq_1_size; ++pos_1) {
                ElementId j1 = sequence_1.elements[pos_1].j;
                const auto& neighbors = (closest_neighbors_.empty())? neighbors_: closest_neighbors_[j1];

                for (ElementId j2: neighbors) {
                    SequenceId i2 = elements_cur_[j2].i;
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
                    if (j1 > j2)
                        continue;

                    const auto& sequence_2 = solution.sequences[i2];
                    SequencePos seq_2_size = sequence_2.elements.size();

                    SubSequence sub_sequences_1[] = {
                        SubSequence(sequence_1, 0, pos_1),
                        SubSequence(sequence_2, 0, pos_2, true)};
                    SubSequence sub_sequences_2[] = {
                        SubSequence(sequence_1, pos_1 + 1, seq_1_size - 1, true),
                        SubSequence(sequence_2, pos_2 + 1, seq_2_size - 1)};
                    auto bgc_tmp = global_cost_concatenate_2(
                            sub_sequences_1, 2, sub_sequences_2, 2,
                            gc, ((m > 1)? &partial_global_costs_cur_2_[i1][i2]: nullptr));
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        Move0 move;
                        move.type = Neighborhoods::InterTwoOptReverse;
                        move.i1 = i1;
                        move.i2 = i2;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = bgc_tmp.second - gc;
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
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::InterSwap)][block_size_1][block_size_2];

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (ElementPos pos_1 = 0; pos_1 < seq_1_size; ++pos_1) {
                ElementId j1 = sequence_1.elements[pos_1].j;
                const auto& neighbors = (closest_neighbors_.empty())? neighbors_: closest_neighbors_[j1];
                if (pos_1 + block_size_1 - 1 >= seq_1_size)
                    continue;

                for (ElementId j2: neighbors) {
                    SequenceId i2 = elements_cur_[j2].i;
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
                    if (block_size_1 == block_size_2 && j1 > j2)
                        continue;

                    const auto& sequence_2 = solution.sequences[i2];
                    SequencePos seq_2_size = sequence_2.elements.size();
                    if (pos_2 + block_size_2 - 1 >= seq_2_size)
                        continue;

                    SubSequence sub_sequences_1[] = {
                        SubSequence(sequence_1, 0, pos_1 - 1),
                        SubSequence(sequence_2, pos_2, pos_2 + block_size_2 - 1),
                        SubSequence(sequence_1, pos_1 + block_size_1, seq_1_size - 1)};
                    SubSequence sub_sequences_2[] = {
                        SubSequence(sequence_2, 0, pos_2 - 1),
                        SubSequence(sequence_1, pos_1, pos_1 + block_size_1 - 1),
                        SubSequence(sequence_2, pos_2 + block_size_2, seq_2_size - 1)};
                    auto bgc_tmp = global_cost_concatenate_2(
                            sub_sequences_1, 3, sub_sequences_2, 3,
                            gc, ((m > 1)? &partial_global_costs_cur_2_[i1][i2]: nullptr));
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        Move0 move;
                        move.type = Neighborhoods::InterSwap;
                        move.k1 = block_size_1;
                        move.k2 = block_size_2;
                        move.i1 = i1;
                        move.i2 = i2;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = bgc_tmp.second - gc;
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_inter_swap_star(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::InterSwapStar)][0][0];

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (SequenceId i2 = i1 + 1; i2 < m; ++i2) {
                if (!neighborhood.modified_sequences[i1]
                        && !neighborhood.modified_sequences[i2])
                    continue;

                const auto& sequence_2 = solution.sequences[i2];
                SequencePos seq_2_size = sequence_2.elements.size();

                // Calcul the 3 best position of the elements from the first
                // sequence into the second sequence.
                std::vector<std::vector<ElementPos>> best_positions_1(seq_1_size);
                for (ElementPos pos_1 = 0; pos_1 <= seq_1_size - 1; ++pos_1) {
                    ElementPos pos_best_1 = -1;
                    ElementPos pos_best_2 = -1;
                    ElementPos pos_best_3 = -1;
                    GlobalCost gc_best_1 = worst<GlobalCost>();
                    GlobalCost gc_best_2 = worst<GlobalCost>();
                    GlobalCost gc_best_3 = worst<GlobalCost>();
                    for (ElementPos pos_2 = 0; pos_2 <= seq_2_size; ++pos_2) {
                        SubSequence sub_sequences[] = {
                            SubSequence(sequence_2, 0, pos_2 - 1),
                            SubSequence(sequence_1, pos_1, pos_1),
                            SubSequence(sequence_2, pos_2, seq_2_size - 1)};
                        auto gci2_tmp = global_cost_concatenate_1(
                                sub_sequences, 3, nullptr);
                        if (pos_best_1 == -1 || !(gci2_tmp >= gc_best_1)) {
                            pos_best_3 = pos_best_2;
                            pos_best_2 = pos_best_1;
                            pos_best_1 = pos_2;
                            gc_best_3 = gc_best_2;
                            gc_best_2 = gc_best_1;
                            gc_best_1 = gci2_tmp;
                        } else if (pos_best_2 == -1 || !(gci2_tmp >= gc_best_2)) {
                            pos_best_3 = pos_best_2;
                            pos_best_2 = pos_2;
                            gc_best_3 = gc_best_2;
                            gc_best_2 = gci2_tmp;
                        } else if (pos_best_3 == -1 || !(gci2_tmp >= gc_best_3)) {
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
                    best_positions_1[pos_1] = {pos_best_1, pos_best_2, pos_best_3};
                }

                // Calcul the 3 best position of the elements from the second
                // sequence into the first sequence.
                std::vector<std::vector<ElementPos>> best_positions_2(seq_2_size);
                for (ElementPos pos_2 = 0; pos_2 <= seq_2_size - 1; ++pos_2) {
                    ElementPos pos_best_1 = -1;
                    ElementPos pos_best_2 = -1;
                    ElementPos pos_best_3 = -1;
                    GlobalCost gc_best_1 = worst<GlobalCost>();
                    GlobalCost gc_best_2 = worst<GlobalCost>();
                    GlobalCost gc_best_3 = worst<GlobalCost>();
                    for (ElementPos pos_1 = 0; pos_1 <= seq_1_size; ++pos_1) {
                        SubSequence sub_sequences[] = {
                            SubSequence(sequence_1, 0, pos_1 - 1),
                            SubSequence(sequence_2, pos_2, pos_2),
                            SubSequence(sequence_1, pos_1, seq_1_size - 1)};
                        auto gci1_tmp = global_cost_concatenate_1(
                                sub_sequences, 3, nullptr);
                        if (pos_best_1 == -1 || !(gci1_tmp >= gc_best_1)) {
                            pos_best_3 = pos_best_2;
                            pos_best_2 = pos_best_1;
                            pos_best_1 = pos_1;
                            gc_best_3 = gc_best_2;
                            gc_best_2 = gc_best_1;
                            gc_best_1 = gci1_tmp;
                        } else if (pos_best_2 == -1 || !(gci1_tmp >= gc_best_2)) {
                            pos_best_3 = pos_best_2;
                            pos_best_2 = pos_1;
                            gc_best_3 = gc_best_2;
                            gc_best_2 = gci1_tmp;
                        } else if (pos_best_3 == -1 || !(gci1_tmp >= gc_best_3)) {
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
                    best_positions_2[pos_2] = {pos_best_1, pos_best_2, pos_best_3};
                }

                for (ElementPos pos_1 = 0; pos_1 <= seq_1_size - 1; ++pos_1) {

                    for (ElementPos pos_2 = 0; pos_2 <= seq_2_size - 1; ++pos_2) {

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
                        GlobalCost gci2_tmp = worst<GlobalCost>();
                        for (int a = 0; a < 4; ++a) {
                            ElementPos p = (a < 3)? best_positions_1[pos_1][a]: pos_2;
                            if (p != -1) {
                                SubSequence sub_sequences[] = {
                                    SubSequence(sequence_2, 0, pos_2 - 1),
                                    SubSequence(sequence_1, p, p),
                                    SubSequence(sequence_2, pos_2, seq_2_size - 1)};
                                auto bgci2 = global_cost_concatenate_1(
                                        sub_sequences, 3, gc, &partial_global_costs_cur_2_[i1][i2]);
                                //std::cout << "p " << p
                                //    << " gc " << to_string(gci2)
                                //    << std::endl;
                                if (bgci2.first && !(bgci2.second >= gci2_tmp)) {
                                    gci2_tmp = bgci2.second;
                                    pos_1_new = p;
                                }
                            }
                        }
                        if (pos_1_new == -1)
                            continue;

                        // Find pos_2_new.
                        ElementPos pos_2_new = -1;
                        GlobalCost gc_tmp = worst<GlobalCost>();
                        for (int a = 0; a < 4; ++a) {
                            ElementPos p = (a < 3)? best_positions_2[pos_2][a]: pos_1;
                            if (p != -1) {
                                SubSequence sub_sequences[] = {
                                    SubSequence(sequence_1, 0, pos_1 - 1),
                                    SubSequence(sequence_2, p, p),
                                    SubSequence(sequence_1, pos_1, seq_1_size - 1)};
                                auto bgc = global_cost_concatenate_1(
                                        sub_sequences, 3, gc, &gci2_tmp);
                                //std::cout << "p " << p
                                //    << " gc " << to_string(gci1)
                                //    << std::endl;
                                if (bgc.first && !(bgc.second >= gc_tmp)) {
                                    gc_tmp = bgc.second;
                                    pos_2_new = p;
                                }
                            }
                        }
                        if (pos_2_new == -1)
                            continue;

                        //std::cout << "pos_1 " << pos_2
                        //    << " pos_2 " << pos_2
                        //    << " pos_1_new " << pos_1_new
                        //    << " pos_2_new " << pos_2_new
                        //    << std::endl;

                        if (!(gc_tmp >= gc)) {
                            Move0 move;
                            move.type = Neighborhoods::InterSwapStar;
                            move.i1 = i1;
                            move.i2 = i2;
                            move.pos_1 = pos_1;
                            move.pos_2 = pos_1_new;
                            move.pos_3 = pos_2;
                            move.pos_4 = pos_2_new;
                            move.global_cost = gc_tmp - gc;
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
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();

            for (ElementPos block_pos = 0; block_pos <= seq_size - 1; ++block_pos) {
                ElementId j = sequence.elements[block_pos].j;
                // Loop through all new positions.
                ElementPos pos_min = std::max(
                        (ElementPos)0,
                        block_pos - parameters_.shift_maximum_distance);
                ElementPos pos_max = std::min(
                        seq_size - 1,
                        block_pos + parameters_.shift_maximum_distance);
                for (ElementPos pos_new = pos_min; pos_new <= pos_max; ++pos_new) {
                    for (Mode mode = 0; mode < number_of_modes(j); ++mode) {
                        std::pair<bool, GlobalCost> bgc_tmp;
                        if (block_pos > pos_new) {
                            SubSequence sub_sequences[] = {
                                SubSequence(sequence, 0, pos_new - 1),
                                SubSequence(j, mode),
                                SubSequence(sequence, pos_new, block_pos - 1),
                                SubSequence(sequence, block_pos + 1, seq_size - 1)};
                            bgc_tmp = global_cost_concatenate_1(
                                    sub_sequences, 4, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                        } else {
                            SubSequence sub_sequences[] = {
                                SubSequence(sequence, 0, block_pos - 1),
                                SubSequence(sequence, block_pos + 1, pos_new),
                                SubSequence(j, mode),
                                SubSequence(sequence, pos_new + 1, seq_size - 1)};
                            bgc_tmp = global_cost_concatenate_1(
                                    sub_sequences, 4, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                        }
                        if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                            Move0 move;
                            move.type = Neighborhoods::ShiftChangeMode;
                            move.i1 = i;
                            move.pos_1 = block_pos;
                            move.pos_2 = pos_new;
                            move.mode = mode;
                            move.global_cost = bgc_tmp.second - gc;
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
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();

            // Loop through all pairs.
            Counter pos_max = seq_size - 1;
            for (ElementPos pos_1 = 0; pos_1 <= pos_max; ++pos_1) {
                ElementId j1 = sequence.elements[pos_1].j;
                ElementId mode_1 = sequence.elements[pos_1].mode;

                ElementPos pos_2_max = std::min(
                        pos_max,
                        pos_1 + parameters_.swap_maximum_distance);
                for (ElementPos pos_2 = pos_1 + 1; pos_2 < pos_2_max; ++pos_2) {
                    ElementId j2 = sequence.elements[pos_2].j;
                    ElementId mode_2 = sequence.elements[pos_2].mode;
                    if (mode_1 == mode_2)
                        continue;

                    SubSequence sub_sequences[] = {
                        SubSequence(sequence, 0, pos_1 - 1),
                        SubSequence(j1, mode_2),
                        SubSequence(sequence, pos_1 + 1, pos_2 - 1),
                        SubSequence(j2, mode_1),
                        SubSequence(sequence, pos_2 + 1, seq_size - 1)};
                    auto bgc_tmp = global_cost_concatenate_1(
                            sub_sequences, 5, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        Move0 move;
                        move.type = Neighborhoods::ModeSwap;
                        move.i1 = i;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = bgc_tmp.second - gc;
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
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();

            // Loop through all pairs.
            Counter pos_max = seq_size - 1;
            for (ElementPos pos_1 = 0; pos_1 <= pos_max; ++pos_1) {
                ElementId j1 = sequence.elements[pos_1].j;
                ElementId mode_1 = sequence.elements[pos_1].mode;

                ElementPos pos_2_max = std::min(
                        pos_max,
                        pos_1 + parameters_.swap_maximum_distance);
                for (ElementPos pos_2 = pos_1 + 1; pos_2 < pos_2_max; ++pos_2) {
                    ElementId j2 = sequence.elements[pos_2].j;
                    ElementId mode_2 = sequence.elements[pos_2].mode;
                    if (mode_1 == mode_2)
                        continue;

                    SubSequence sub_sequences[] = {
                        SubSequence(sequence, 0, pos_1 - 1),
                        SubSequence(j2, mode_1),
                        SubSequence(sequence, pos_1 + 1, pos_2 - 1),
                        SubSequence(j1, mode_2),
                        SubSequence(sequence, pos_2 + 1, seq_size - 1)};
                    auto bgc_tmp = global_cost_concatenate_1(
                            sub_sequences, 5, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        Move0 move;
                        move.type = Neighborhoods::SwapWithModes;
                        move.i1 = i;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = bgc_tmp.second - gc;
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline std::vector<Move0> explore_add_2(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        GlobalCost gc = global_cost(solution);
        std::vector<Move0> improving_moves;
        for (SequenceId i = 0; i < m; ++i) {
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();

            // Loop through all new positions.
            for (ElementPos pos = 0; pos <= seq_size; ++pos) {

                for (ElementId j = 0; j < n; ++j) {
                    if (elements_cur_[j].mode != -1)
                        continue;

                    for (Mode mode = 0; mode < number_of_modes(j); ++mode) {
                        SubSequence sub_sequences[] = {
                            SubSequence(sequence, 0, pos - 1),
                            SubSequence(j, mode),
                            SubSequence(sequence, pos, seq_size - 1)};
                        auto bgc_tmp = global_cost_concatenate_1(
                                sub_sequences, 3, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                        if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                            Move0 move;
                            move.type = Neighborhoods::Add;
                            move.i1 = i;
                            move.j = j;
                            move.mode = mode;
                            move.pos_1 = pos;
                            move.global_cost = bgc_tmp.second - gc;
                            improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
        return improving_moves;
    }

    inline std::vector<Move0> explore_add(
            const Solution& solution,
            ElementId j,
            SequenceId i_old = -1,
            ElementId j_prec_old = -2,
            Mode mode_old = -1)
    {
        //std::cout << "explore_add j " << j << std::endl;
        std::vector<Move0> improving_moves;
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        for (SequenceId i = 0; i < m; ++i) {
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();

            // Loop through all new positions.
            for (ElementPos pos = 0; pos <= seq_size; ++pos) {

                for (Mode mode = 0; mode < number_of_modes(j); ++mode) {

                    if (j_prec_old == -1
                            && i == i_old
                            && (seq_size > 0 || m > 1)
                            && pos == 0
                            && mode == mode_old)
                        continue;
                    if (pos > 0
                            && j_prec_old == sequence.elements[pos - 1].j
                            && mode == mode_old)
                        continue;

                    SubSequence sub_sequences[] = {
                        SubSequence(sequence, 0, pos - 1),
                        SubSequence(j, mode),
                        SubSequence(sequence, pos, seq_size - 1)};
                    auto bgc_tmp = global_cost_concatenate_1(
                            sub_sequences, 3, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                    //std::cout << to_string(gc_tmp) << std::endl;
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        Move0 move;
                        move.type = Neighborhoods::Add;
                        move.i1 = i;
                        move.j = j;
                        move.mode = mode;
                        move.pos_1 = pos;
                        move.global_cost = bgc_tmp.second - gc;
                        improving_moves.push_back(move);
                    }
                }
            }
        }
        return improving_moves;
    }

    inline std::vector<Move0> explore_add_end(
            const Solution& solution)
    {
        std::vector<Move0> improving_moves;
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        GlobalCost gc = global_cost(solution);
        for (SequenceId i = 0; i < m; ++i) {
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            ElementPos pos = seq_size;

            // Loop through all new positions.
            for (ElementId j = 0; j < n; ++j) {
                if (elements_cur_[j].mode != -1)
                    continue;

                for (Mode mode = 0; mode < number_of_modes(j); ++mode) {
                    SubSequence sub_sequences[] = {
                        SubSequence(sequence, 0, pos - 1),
                        SubSequence(j, mode),
                        SubSequence(sequence, pos, seq_size - 1)};
                    auto bgc_tmp = global_cost_concatenate_1(
                            sub_sequences, 3, gc, ((m > 1)? &partial_global_costs_cur_1_[i]: nullptr));
                    if (bgc_tmp.first && !(bgc_tmp.second >= gc)) {
                        Move0 move;
                        move.type = Neighborhoods::Add;
                        move.i1 = i;
                        move.j = j;
                        move.mode = mode;
                        move.pos_1 = pos;
                        move.global_cost = bgc_tmp.second - gc;
                        improving_moves.push_back(move);
                    }
                }
            }
        }
        return improving_moves;
    }

    inline void apply_move(
            Solution& solution,
            const Move0& move)
    {
        auto begin = std::chrono::steady_clock::now();
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();

        ElementPos n_old = 0;
        for (SequenceId i = 0; i < m; ++i)
            n_old += solution.sequences[i].elements.size();
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
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
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
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
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
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
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
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
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
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp, elements[p]);
            append(sequence_tmp, {move.j, move.mode});
            for (ElementPos p = move.pos_1; p < seq_size; ++p)
                append(sequence_tmp, elements[p]);
            break;
        } case Neighborhoods::Remove: {
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            for (ElementPos p = 0; p < seq_size; ++p) {
                if (p == move.pos_1)
                    continue;
                append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::Replace: {
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp, elements[p]);
            append(sequence_tmp, {move.j, move.mode});
            for (ElementPos p = move.pos_1 + 1; p < seq_size; ++p)
                append(sequence_tmp, elements[p]);
            break;
        } case Neighborhoods::InterTwoOpt: {
            //std::cout << "Apply InterTwoOpt"
            //    << " i1 " << move.i1
            //    << " i2 " << move.i2
            //    << " pos_1 " << move.pos_1
            //    << " pos_2 " << move.pos_2
            //    << " gc " << to_string(move.global_cost)
            //    << std::endl;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1 && i != move.i2)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements_1 = solution.sequences[move.i1].elements;
            const auto& elements_2 = solution.sequences[move.i2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.i1];
            sequence_tmp_1 = empty_sequence(move.i1);
            for (ElementPos p = 0; p <= move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_2 + 1; p <= seq_2_size - 1; ++p)
                append(sequence_tmp_1, elements_2[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.i2];
            sequence_tmp_2 = empty_sequence(move.i2);
            for (ElementPos p = 0; p <= move.pos_2; ++p)
                append(sequence_tmp_2, elements_2[p]);
            for (ElementPos p = move.pos_1 + 1; p <= seq_1_size - 1; ++p)
                append(sequence_tmp_2, elements_1[p]);
            break;
        } case Neighborhoods::InterTwoOptReverse: {
            //std::cout << "Apply InterTwoOptReverse" << std::endl;
            //std::cout << "i1 " << move.i1
            //    << " i2 " << move.i2
            //    << " pos_1 " << move.pos_1
            //    << " pos_2 " << move.pos_2
            //    << " gc " << to_string(move.global_cost)
            //    << std::endl;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1 && i != move.i2)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements_1 = solution.sequences[move.i1].elements;
            const auto& elements_2 = solution.sequences[move.i2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.i1];
            sequence_tmp_1 = empty_sequence(move.i1);
            for (ElementPos p = 0; p <= move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_2; p >= 0; --p)
                append(sequence_tmp_1, elements_2[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.i2];
            sequence_tmp_2 = empty_sequence(move.i2);
            for (ElementPos p = seq_1_size - 1; p >= move.pos_1 + 1; --p)
                append(sequence_tmp_2, elements_1[p]);
            for (ElementPos p = move.pos_2 + 1; p <= seq_2_size - 1; ++p)
                append(sequence_tmp_2, elements_2[p]);
            //std::cout << "End apply InterTwoOptReverse" << std::endl;
            break;
        } case Neighborhoods::InterShift: {
            ElementPos block_size = move.k1;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1 && i != move.i2)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements_1 = solution.sequences[move.i1].elements;
            const auto& elements_2 = solution.sequences[move.i2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.i1];
            sequence_tmp_1 = empty_sequence(move.i1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_1 + block_size; p < seq_1_size; ++p)
                append(sequence_tmp_1, elements_1[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.i2];
            sequence_tmp_2 = empty_sequence(move.i2);
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
            //    << " i1 " << move.i1
            //    << " i2 " << move.i2
            //    << " pos_1 " << move.pos_1
            //    << " pos_2 " << move.pos_2
            //    << " gc " << to_string(move.global_cost)
            //    << std::endl;
            ElementPos block_size_1 = move.k1;
            ElementPos block_size_2 = move.k2;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1 && i != move.i2)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements_1 = solution.sequences[move.i1].elements;
            const auto& elements_2 = solution.sequences[move.i2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.i1];
            sequence_tmp_1 = empty_sequence(move.i1);
            for (ElementPos p = 0; p <= move.pos_1 - 1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_2; p <= move.pos_2 + block_size_2 - 1; ++p)
                append(sequence_tmp_1, elements_2[p]);
            for (ElementPos p = move.pos_1 + block_size_1; p <= seq_1_size - 1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.i2];
            sequence_tmp_2 = empty_sequence(move.i2);
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
            //    << " i1 " << move.i1
            //    << " i2 " << move.i2
            //    << " pos_1 " << move.pos_1
            //    << " pos_2 " << move.pos_2
            //    << " gc " << to_string(move.global_cost)
            //    << std::endl;
            ElementPos block_size = move.k1;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1 && i != move.i2)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements_1 = solution.sequences[move.i1].elements;
            const auto& elements_2 = solution.sequences[move.i2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.i1];
            sequence_tmp_1 = empty_sequence(move.i1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_1 + block_size; p < seq_1_size; ++p)
                append(sequence_tmp_1, elements_1[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.i2];
            sequence_tmp_2 = empty_sequence(move.i2);
            for (ElementPos p = 0; p < move.pos_2; ++p)
                append(sequence_tmp_2, elements_2[p]);
            for (ElementPos p = move.pos_1 + block_size - 1; p >= move.pos_1; --p)
                append(sequence_tmp_2, elements_1[p]);
            for (ElementPos p = move.pos_2; p < seq_2_size; ++p)
                append(sequence_tmp_2, elements_2[p]);
            break;
        } case Neighborhoods::ShiftChangeMode: {
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            if (move.pos_1 > move.pos_2) {
                for (ElementPos p = 0; p < move.pos_2; ++p)
                    append(sequence_tmp, elements[p]);
                append(sequence_tmp, {elements[move.pos_1].j, move.mode});
                for (ElementPos p = move.pos_2; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + 1; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            } else {
                for (ElementPos p = 0; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + 1; p < move.pos_2 + 1; ++p)
                    append(sequence_tmp, elements[p]);
                append(sequence_tmp, {elements[move.pos_1].j, move.mode});
                for (ElementPos p = move.pos_2 + 1; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::ModeSwap: {
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
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
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            for (ElementPos p = 0; p < seq_size; ++p) {
                SequenceElement se = elements[p];
                if (p == move.pos_1) {
                    se.j = elements[move.pos_2].j;
                } else if (p == move.pos_2) {
                    se.j = elements[move.pos_1].j;
                }
                append(sequence_tmp, se);
            }
            break;
        } case Neighborhoods::InterSwapStar: {
            //std::cout << "i1 " << move.i1
            //    << " i2 " << move.i2
            //    << " pos_1_old " << move.pos_1
            //    << " pos_2_old " << move.pos_3
            //    << " pos_1_new " << move.pos_2
            //    << " pos_2_new " << move.pos_4
            //    << " gc " << to_string(move.global_cost)
            //    << std::endl;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1 && i != move.i2)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements_1 = solution.sequences[move.i1].elements;
            const auto& elements_2 = solution.sequences[move.i2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.i1];
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.i2];
            sequence_tmp_1 = empty_sequence(move.i1);
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
            sequence_tmp_2 = empty_sequence(move.i2);
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
        for (SequenceId i = 0; i < m; ++i)
            n_new += solution_tmp_.sequences[i].elements.size();
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
        auto end = std::chrono::steady_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
        apply_move_time += time_span.count();
        number_of_apply_move_calls++;
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

    /**
     * Structure storing the closest neighbors for each element.
     *
     * Then, in the local search, only moves between close elements are
     * considered. In routing problems, this significantly reduces the
     * computational effort spent in the local search without degrading the
     * final solution too much.
     *
     * This is only used if the number of elements in greater than 100 and if
     * the SequencingScheme implements a 'double distance(ElementId,
     * ElementId)' method.
     */
    std::vector<std::vector<ElementId>> closest_neighbors_;

    /**
     * Size of the neighborhood.
     */
    ElementPos neighborhood_size_ = 100;

    /**
     * Vector which contains all element ids.
     *
     * It is used to loop through neighbors when closest_neighbors_ is not
     * used.
     */
    std::vector<ElementId> neighbors_;

    /** Number of calls to a crossover operator. */
    Counter number_of_crossover_calls = 0;
    /** Time spent in crossover operators. */
    double crossover_time = 0.0;
    /** Number of calls to the local search method. */
    Counter number_of_local_search_calls = 0;
    /** Time spent in the local search method. */
    double local_search_time = 0.0;
    /** Number of calls to the apply_move method. */
    Counter number_of_apply_move_calls = 0;
    /** Time spent in the apply_move method. */
    double apply_move_time = 0.0;
    /** Number of calls to an initial solution generator. */
    Counter number_of_initial_solution_calls = 0;
    /** Time spent generating initial solutions. */
    double initial_solution_time = 0.0;

    /*
     * Temporary structures.
     */

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

