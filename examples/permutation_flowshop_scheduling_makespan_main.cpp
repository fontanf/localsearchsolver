/**
 * Permutation flow shop scheduling problem, makespan
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/permutation_flowshop_scheduling_makespan.hpp
 *
 * TODO
 */

#include "read_args.hpp"

#include "localsearchsolver/common.hpp"

#include "orproblems/scheduling/permutation_flowshop_scheduling_makespan.hpp"

#include "optimizationtools/utils/utils.hpp"

using namespace localsearchsolver;
using namespace orproblems::permutation_flowshop_scheduling_makespan;

class LocalScheme
{

public:

    struct Parameters
    {
        JobPos block_size_max = 8;
        bool shuffle_neighborhood_order = true;
        Counter number_of_perturbations = 10;
    };

    /*
     * Constructors and destructor
     */

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters),
        positions1_(instance.number_of_jobs()),
        positions2_(instance.number_of_jobs()),
        times_(instance_.number_of_machines(), 0),
        heads_(instance.number_of_jobs() + 1),
        tails_(instance.number_of_jobs() + 1),
        completion_times_(instance.number_of_jobs() + 1)
    {
        std::iota(positions1_.begin(), positions1_.end(), 0);
        std::iota(positions2_.begin(), positions2_.end(), 0);
        for (JobId job_id = 0;
                job_id < instance_.number_of_jobs() + 1;
                ++job_id) {
            heads_[job_id] = std::vector<Time>(instance_.number_of_machines(), 0);
            tails_[job_id] = std::vector<Time>(instance_.number_of_machines(), 0);
            completion_times_[job_id] = std::vector<Time>(instance_.number_of_machines(), 0);
        }
    }

    /*
     * Global cost
     */

    /** Global cost: <Makespan>. */
    using GlobalCost = std::tuple<Time>;

    inline Time&       makespan(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time  makespan(const GlobalCost& global_cost) { return std::get<0>(global_cost); }

    /*
     * Solutions
     */

    struct Solution
    {
        /** Jobs. */
        std::vector<JobId> jobs;

        /** Makespan. */
        Time makespan = 0;
    };

    inline Solution empty_solution() const
    {
        Solution solution;
        return solution;
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator)
    {
        Solution solution = empty_solution();
        std::vector<JobId> jobs(instance_.number_of_jobs());
        std::iota(jobs.begin(), jobs.end(), 0);
        std::shuffle(jobs.begin(), jobs.end(), generator);
        return compact2solution(jobs);
    }

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            solution.makespan,
        };
    }

    /*
     * Local search
     */

    struct Perturbation;

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Perturbation& = Perturbation())
    {
        std::vector<Counter> neighborhoods;
        for (JobPos block_size = 1;
                block_size <= parameters_.block_size_max;
                ++block_size) {
            neighborhoods.push_back(block_size);
        }

        Counter it = 0;
        (void)it;
        for (;; ++it) {
            //std::cout << "it " << it
            //    << " c " << to_string(global_cost(solution))
            //    << std::endl;
            //print(std::cout, solution);

            if (parameters_.shuffle_neighborhood_order)
                std::shuffle(neighborhoods.begin(), neighborhoods.end(), generator);
            bool improved = false;
            // Loop through neighborhoods.
            for (Counter block_size: neighborhoods) {
                std::shuffle(positions1_.begin(), positions1_.end(), generator);
                std::shuffle(positions2_.begin(), positions2_.end(), generator);
                JobPos pos_best = -1;
                JobPos pos_new_best = -1;
                GlobalCost c_best = global_cost(solution);
                for (JobPos pos: positions1_) {
                    if (pos > (JobPos)solution.jobs.size() - block_size)
                        continue;
                    compute_structures(solution, pos, block_size);
                    for (JobPos pos_new: positions2_) {
                        if (pos == pos_new || pos_new > (JobPos)solution.jobs.size() - block_size)
                            continue;
                        Time makespan = 0;
                        for (MachineId machine_id = 0;
                                machine_id < instance_.number_of_machines();
                                ++machine_id) {
                            makespan = std::max(makespan,
                                    completion_times_[pos_new][machine_id]
                                             + tails_[pos_new][machine_id]);
                        }
                        GlobalCost c = {makespan};
                        if (c >= c_best)
                            continue;
                        if (pos_best != -1 && !dominates(c, c_best))
                            continue;
                        pos_best = pos;
                        pos_new_best = pos_new;
                        c_best = c;
                    }
                }
                if (pos_best != -1) {
                    improved = true;
                    assert(makespan(c_best) < solution.makespan);
                    // Apply best perturbation.
                    //std::cout << "pos_best " << pos_best
                    //    << " pos_new_best " << pos_new_best
                    //    << " size " << block_size
                    //    << std::endl;
                    std::vector<JobId> jobs;
                    if (pos_best > pos_new_best) {
                        for (JobPos p = 0; p < pos_new_best; ++p)
                            jobs.push_back(solution.jobs[p]);
                        for (JobPos p = pos_best; p < pos_best + block_size; ++p)
                            jobs.push_back(solution.jobs[p]);
                        for (JobPos p = pos_new_best; p < pos_best; ++p)
                            jobs.push_back(solution.jobs[p]);
                        for (JobPos p = pos_best + block_size; p < (JobPos)solution.jobs.size(); ++p)
                            jobs.push_back(solution.jobs[p]);
                    } else {
                        for (JobPos p = 0; p < pos_best; ++p)
                            jobs.push_back(solution.jobs[p]);
                        for (JobPos p = pos_best + block_size; p < pos_new_best + block_size; ++p)
                            jobs.push_back(solution.jobs[p]);
                        for (JobPos p = pos_best; p < pos_best + block_size; ++p)
                            jobs.push_back(solution.jobs[p]);
                        for (JobPos p = pos_new_best + block_size; p < (JobPos)solution.jobs.size(); ++p)
                            jobs.push_back(solution.jobs[p]);
                    }
                    assert((JobPos)jobs.size() <= instance_.number_of_jobs());
                    compute(solution, jobs);
                    if (solution.makespan != makespan(c_best)) {
                        std::cout << "pos_best " << pos_best
                            << " pos_new_best " << pos_new_best
                            << " size " << block_size
                            << std::endl;
                        std::cout << makespan(c_best) << std::endl;
                        std::cout << solution.makespan << std::endl;
                        for (MachineId machine_id = 0;
                                machine_id < instance_.number_of_machines();
                                ++machine_id) {
                            std::cout << "machine_id " << machine_id
                                << " " << heads_[((pos_new_best <= pos_best)? pos_new_best: pos_new_best - block_size)][machine_id]
                                << " " << completion_times_[((pos_new_best <= pos_best)? pos_new_best: pos_new_best - block_size)][machine_id]
                                << " " << tails_[((pos_new_best <= pos_best)? pos_new_best: pos_new_best - block_size)][machine_id]
                                << std::endl;
                        }
                        solution_format(std::cout, solution, 1);
                    }
                    assert(solution.makespan == makespan(c_best));
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
     * Iterated local search
     */

    struct Perturbation
    {
        Perturbation(): pos_1(-1), global_cost(worst<GlobalCost>()) { }

        JobPos pos_1;
        JobPos pos_2;
        JobPos pos_3;
        JobPos pos_4;
        GlobalCost global_cost;
    };

    inline std::vector<Perturbation> perturbations(
            const Solution& solution,
            std::mt19937_64& generator)
    {
        std::vector<Perturbation> perturbations;
        for (Counter perturbation_id = 0;
                perturbation_id < parameters_.number_of_perturbations;
                ++perturbation_id) {
            std::vector<JobPos> edges = optimizationtools::bob_floyd<JobPos>(
                    4,
                    solution.jobs.size() + 1,
                    generator);
            std::sort(edges.begin(), edges.end());
            Perturbation perturbation;
            perturbation.pos_1 = edges[0];
            perturbation.pos_2 = edges[1];
            perturbation.pos_3 = edges[2];
            perturbation.pos_4 = edges[3];
            assert(perturbation.pos_1 >= 0);
            assert(perturbation.pos_4 <= (JobPos)solution.jobs.size());
            perturbation.global_cost = global_cost(solution);
            perturbations.push_back(perturbation);
        }
        return perturbations;
    }

    inline void apply_perturbation(
            Solution& solution,
            const Perturbation& perturbation,
            std::mt19937_64&)
    {
        std::vector<JobId> jobs;
        for (JobPos pos = 0; pos < perturbation.pos_1; ++pos)
            jobs.push_back(solution.jobs[pos]);
        for (JobPos pos = perturbation.pos_3; pos < perturbation.pos_4; ++pos)
            jobs.push_back(solution.jobs[pos]);
        for (JobPos pos = perturbation.pos_2; pos < perturbation.pos_3; ++pos)
            jobs.push_back(solution.jobs[pos]);
        for (JobPos pos = perturbation.pos_1; pos < perturbation.pos_2; ++pos)
            jobs.push_back(solution.jobs[pos]);
        for (JobPos pos = perturbation.pos_4; pos < (JobPos)solution.jobs.size(); ++pos)
            jobs.push_back(solution.jobs[pos]);
        assert((JobPos)jobs.size() <= instance_.number_of_jobs());
        compute(solution, jobs);
    }

    /*
     * Best first local search
     */

    /**
     * solution[n] is the first job scheduled.
     * solution[j] is the index of the job scheduled after job j, n if job
     * j is the last job, or -1 if job j is not in the solution.
     */
    using CompactSolution = std::vector<JobId>;

    struct CompactSolutionHasher
    {
        std::hash<JobId> hasher;

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
            for (JobId job_id: *compact_solution)
                optimizationtools::hash_combine(hash, hasher(job_id));
            return hash;
        }
    };

    inline CompactSolutionHasher compact_solution_hasher() const { return CompactSolutionHasher(); }

    CompactSolution solution2compact(const Solution& solution)
    {
        return solution.jobs;
    }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        auto solution = empty_solution();
        compute(solution, compact_solution);
        return solution;
    }

    struct PerturbationHasher
    {
        inline bool hashable(const Perturbation&) const { return false; }
        inline bool operator()(const Perturbation&, const Perturbation&) const { return false; }
        inline std::size_t operator()(const Perturbation&) const { return 0; }
    };

    inline PerturbationHasher perturbation_hasher() const { return PerturbationHasher(); }

    /*
     * Outputs
     */

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Permutation flow shop scheduling problem, makespan" << std::endl;
        instance_.format(os, verbosity_level);
    }

    void solution_format(
            std::ostream& os,
            const Solution& solution,
            int verbosity_level)
    {
        (void)verbosity_level;
        os << "jobs:";
        for (JobId job_id: solution.jobs)
            os << " " << job_id;
        os << std::endl;
        os << "makespan: " << solution.makespan << std::endl;
    }

    void solution_write(
            const Solution& solution,
            const std::string& certificate_path) const
    {
        if (certificate_path.empty())
            return;
        std::ofstream file(certificate_path);
        if (!file.good()) {
            throw std::runtime_error(
                    "Unable to open file \"" + certificate_path + "\".");
        }

        for (JobId job_id: solution.jobs)
            file << job_id << " ";
    }

private:

    /*
     * Manipulate solutions
     */

    inline void compute(
            Solution& solution,
            const std::vector<JobId>& jobs)
    {
        solution.jobs = jobs;
        std::fill(times_.begin(), times_.end(), 0);
        for (JobId job_id: solution.jobs) {
            times_[0] = times_[0] + instance_.processing_time(job_id, 0);
            for (MachineId machine_id = 1;
                    machine_id < instance_.number_of_machines();
                    ++machine_id) {
                if (times_[machine_id - 1] > times_[machine_id]) {
                    times_[machine_id] = times_[machine_id - 1]
                        + instance_.processing_time(job_id, machine_id);
                } else {
                    times_[machine_id] = times_[machine_id]
                        + instance_.processing_time(job_id, machine_id);
                }
            }
        }
        solution.makespan = times_[instance_.number_of_machines() - 1];
    }

    /*
     * Evaluate moves
     */

    inline void compute_structures(
            const Solution& solution,
            JobPos pos,
            JobPos size)
    {
        // Compute heads_.
        for (JobPos pos_new = 0; pos_new < (JobPos)solution.jobs.size() - size; ++pos_new) {
            JobId job_id = solution.jobs[((pos_new < pos)? pos_new: pos_new + size)];
            heads_[pos_new + 1][0] = heads_[pos_new][0]
                + instance_.processing_time(job_id, 0);
            for (MachineId machine_id = 1;
                    machine_id < instance_.number_of_machines();
                    ++machine_id) {
                if (heads_[pos_new + 1][machine_id - 1] > heads_[pos_new][machine_id]) {
                    heads_[pos_new + 1][machine_id] = heads_[pos_new + 1][machine_id - 1]
                        + instance_.processing_time(job_id, machine_id);
                } else {
                    heads_[pos_new + 1][machine_id] = heads_[pos_new][machine_id]
                        + instance_.processing_time(job_id, machine_id);
                }
            }
        }

        // Compute completion_times_.
        for (JobPos pos_new = 0; pos_new <= (JobPos)solution.jobs.size() - size; ++pos_new) {
            for (MachineId machine_id = 0;
                    machine_id < instance_.number_of_machines();
                    ++machine_id) {
                completion_times_[pos_new][machine_id] = heads_[pos_new][machine_id];
            }
            for (JobPos pos_0 = pos; pos_0 < pos + size; ++pos_0) {
                JobId job_id_0 = solution.jobs[pos_0];
                completion_times_[pos_new][0] = completion_times_[pos_new][0]
                    + instance_.processing_time(job_id_0, 0);
                for (MachineId machine_id = 1;
                        machine_id < instance_.number_of_machines();
                        ++machine_id) {
                    if (completion_times_[pos_new][machine_id] > completion_times_[pos_new][machine_id - 1]) {
                        completion_times_[pos_new][machine_id] = completion_times_[pos_new][machine_id]
                            + instance_.processing_time(job_id_0, machine_id);
                    } else {
                        completion_times_[pos_new][machine_id] = completion_times_[pos_new][machine_id - 1]
                            + instance_.processing_time(job_id_0, machine_id);
                    }
                }
            }
        }

        // Update tails_.
        for (MachineId machine_id = instance_.number_of_machines() - 1;
                machine_id >= 0;
                --machine_id) {
            tails_[solution.jobs.size() - size][machine_id] = 0;
        }
        for (JobPos pos_new = solution.jobs.size() - size - 1; pos_new >= 0; --pos_new) {
            JobId job_id = solution.jobs[((pos_new < pos)? pos_new: pos_new + size)];
            assert(job_id >= 0);
            assert(job_id < instance_.number_of_jobs());
            tails_[pos_new][instance_.number_of_machines() - 1]
                = tails_[pos_new + 1][instance_.number_of_machines() - 1]
                + instance_.processing_time(job_id, instance_.number_of_machines() - 1);
            for (MachineId machine_id = instance_.number_of_machines() - 2;
                    machine_id >= 0;
                    --machine_id) {
                if (tails_[pos_new][machine_id + 1] > tails_[pos_new + 1][machine_id]) {
                    tails_[pos_new][machine_id] = tails_[pos_new][machine_id + 1]
                        + instance_.processing_time(job_id, machine_id);
                } else {
                    tails_[pos_new][machine_id] = tails_[pos_new + 1][machine_id]
                        + instance_.processing_time(job_id, machine_id);
                }
            }
        }
    }

    /*
     * Private attributes
     */

    /** Instance. */
    const Instance& instance_;

    /** Parmaeters. */
    Parameters parameters_;

    std::vector<JobPos> positions1_;

    std::vector<JobPos> positions2_;

    std::vector<Time> times_;

    std::vector<std::vector<Time>> heads_;

    std::vector<std::vector<Time>> tails_;

    std::vector<std::vector<Time>> completion_times_;

};

int main(int argc, char *argv[])
{
    // Create command line options.
    boost::program_options::options_description desc = setup_args();
    desc.add_options()
        ("block-size-max", boost::program_options::value<JobPos>(), "")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;;
        throw "";
    }
    try {
        boost::program_options::notify(vm);
    } catch (const boost::program_options::required_option& e) {
        std::cout << desc << std::endl;;
        throw "";
    }

    // Create instance.
    InstanceBuilder instance_builder;
    instance_builder.read(
            vm["input"].as<std::string>(),
            vm["format"].as<std::string>());
    const Instance instance = instance_builder.build();

    // Create local scheme.
    LocalScheme::Parameters parameters;
    if (vm.count("block-size-max"))
        parameters.block_size_max = vm["block-size-max"].as<Counter>();
    LocalScheme local_scheme(instance, parameters);

    // Run algorithm.
    std::string algorithm = vm["algorithm"].as<std::string>();
    auto output =
        (algorithm == "multi-start-local-search")?
        run_multi_start_local_search(local_scheme, vm):
        (algorithm == "iterated-local-search")?
        run_iterated_local_search(local_scheme, vm):
        run_best_first_local_search(local_scheme, vm);

    // Run checker.
    if (vm["print-checker"].as<int>() > 0
            && vm["certificate"].as<std::string>() != "") {
        std::cout << std::endl
            << "Checker" << std::endl
            << "-------" << std::endl;
        instance.check(
                vm["certificate"].as<std::string>(),
                std::cout,
                vm["print-checker"].as<int>());
    }

    return 0;
}
