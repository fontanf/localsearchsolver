#include "examples/read_args.hpp"
#include "localsearchsolver/read_args.hpp"

#include <boost/program_options.hpp>

using namespace localsearchsolver;

namespace po = boost::program_options;

template <typename LocalScheme>
SolutionPool<LocalScheme> run_astarlocalsearch(
        const std::vector<std::string>& algorithm_args,
        LocalScheme& local_scheme,
        const optimizationtools::Info& info)
{
    std::vector<char*> algorithm_argv;
    for (Counter i = 0; i < (Counter)algorithm_args.size(); ++i)
        algorithm_argv.push_back(const_cast<char*>(algorithm_args[i].c_str()));

    auto parameters = read_astar_args<LocalScheme>(algorithm_argv);
    parameters.info = info;
    if (parameters.initial_solution_ids.empty())
        parameters.initial_solution_ids = std::vector<Counter>(parameters.number_of_threads_2, 0);
    return a_star_local_search(local_scheme, parameters).solution_pool;
}

template <typename LocalScheme>
SolutionPool<LocalScheme> run_geneticlocalsearch(
        const std::vector<std::string>& algorithm_args,
        LocalScheme& local_scheme,
        const optimizationtools::Info& info)
{
    std::vector<char*> algorithm_argv;
    for (Counter i = 0; i < (Counter)algorithm_args.size(); ++i)
        algorithm_argv.push_back(const_cast<char*>(algorithm_args[i].c_str()));

    auto parameters = read_genetic_args<LocalScheme>(algorithm_argv);
    parameters.info = info;
    return genetic_local_search(local_scheme, parameters).solution_pool;
}

int main(int argc, char *argv[])
{

    // Parse program options

    std::string problem = "knapsackwithconflicts";
    std::string instance_path = "";
    std::string output_path = "";
    std::string certificate_path = "";
    std::string format = "";
    std::string algorithm = "astarlocalsearch";
    std::string local_scheme_parameters = "";
    double time_limit = std::numeric_limits<double>::infinity();
    Seed seed = 0;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("problem,p", po::value<std::string>(&problem)->required(), "set problem (required)")
        ("input,i", po::value<std::string>(&instance_path)->required(), "set input path (required)")
        ("format,f", po::value<std::string>(&format), "set input file format")
        ("output,o", po::value<std::string>(&output_path), "set JSON output path")
        ("certificate,c", po::value<std::string>(&certificate_path), "set certificate path")
        ("algorithm_args,a", po::value<std::string>(&algorithm), "set algorithm")
        ("localscheme,l", po::value<std::string>(&local_scheme_parameters), "set local scheme parameters")
        ("time-limit,t", po::value<double>(&time_limit), "set time limit in seconds")
        ("seed,s", po::value<Seed>(&seed), "set seed")
        ("only-write-at-the-end,e", "only write output and certificate files at the end")
        ("verbose,v", "")
        ("print-instance", "")
        ("print-solution", "")
        ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;;
        return 1;
    }
    try {
        po::notify(vm);
    } catch (const po::required_option& e) {
        std::cout << desc << std::endl;;
        return 1;
    }

    std::vector<std::string> algorithm_args
        = boost::program_options::split_unix(algorithm);

    optimizationtools::Info info = optimizationtools::Info()
        .set_verbose(vm.count("verbose"))
        .set_timelimit(time_limit)
        .set_outputfile(output_path)
        .set_onlywriteattheend(vm.count("only-write-at-the-end"))
        ;

    // Run algorithm

    std::vector<std::string> local_scheme_args
        = boost::program_options::split_unix(local_scheme_parameters);
    std::vector<char*> local_scheme_argv;
    std::string dummy = "dummy";
    local_scheme_argv.push_back(const_cast<char*>(dummy.c_str()));
    for (Counter i = 0; i < (Counter)local_scheme_args.size(); ++i)
        local_scheme_argv.push_back(const_cast<char*>(local_scheme_args[i].c_str()));

    if (problem == "knapsackwithconflicts") {
        knapsackwithconflicts::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        auto parameters_local_scheme = read_knapsackwithconflicts_args(local_scheme_argv);
        knapsackwithconflicts::LocalScheme local_scheme(instance, parameters_local_scheme);
        auto solution_pool = (algorithm_args[0] == "astarlocalsearch")?
            run_astarlocalsearch(algorithm_args, local_scheme, info):
            run_geneticlocalsearch(algorithm_args, local_scheme, info);
        local_scheme.write(solution_pool.best(), certificate_path);
        if (vm.count("print-solution"))
            local_scheme.print(std::cout, solution_pool.best());

    } else if (problem == "multidimensionalmultiplechoiceknapsack") {
        multidimensionalmultiplechoiceknapsack::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        auto parameters_local_scheme = read_multidimensionalmultiplechoiceknapsack_args(local_scheme_argv);
        multidimensionalmultiplechoiceknapsack::LocalScheme local_scheme(instance, parameters_local_scheme);
        auto solution_pool = run_astarlocalsearch(algorithm_args, local_scheme, info);
        local_scheme.write(solution_pool.best(), certificate_path);
        if (vm.count("print-solution"))
            local_scheme.print(std::cout, solution_pool.best());

    } else if (problem == "quadraticassignment") {
        quadraticassignment::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        auto parameters_local_scheme = read_quadraticassignment_args(local_scheme_argv);
        quadraticassignment::LocalScheme local_scheme(instance, parameters_local_scheme);
        auto solution_pool = (algorithm_args[0] == "astarlocalsearch")?
            run_astarlocalsearch(algorithm_args, local_scheme, info):
            run_geneticlocalsearch(algorithm_args, local_scheme, info);
        local_scheme.write(solution_pool.best(), certificate_path);
        if (vm.count("print-solution"))
            local_scheme.print(std::cout, solution_pool.best());

    } else if (problem == "travellingsalesman") {
        travellingsalesman::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        auto parameters_local_scheme = read_travellingsalesman_args(local_scheme_argv);
        travellingsalesman::LocalScheme local_scheme(instance, parameters_local_scheme);
        auto solution_pool = run_astarlocalsearch(algorithm_args, local_scheme, info);
        local_scheme.write(solution_pool.best(), certificate_path);
        if (vm.count("print-solution"))
            local_scheme.print(std::cout, solution_pool.best());

    } else if (problem == "schedulingwithsdsttwt") {
        schedulingwithsdsttwt::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        auto parameters_local_scheme = read_schedulingwithsdsttwt_args(local_scheme_argv);
        schedulingwithsdsttwt::LocalScheme local_scheme(instance, parameters_local_scheme);
        auto solution_pool = run_astarlocalsearch(algorithm_args, local_scheme, info);
        local_scheme.write(solution_pool.best(), certificate_path);
        if (vm.count("print-solution"))
            local_scheme.print(std::cout, solution_pool.best());

    } else if (problem == "permutationflowshopschedulingmakespan") {
        permutationflowshopschedulingmakespan::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        auto parameters_local_scheme = read_permutationflowshopschedulingmakespan_args(local_scheme_argv);
        permutationflowshopschedulingmakespan::LocalScheme local_scheme(instance, parameters_local_scheme);
        auto solution_pool = run_astarlocalsearch(algorithm_args, local_scheme, info);
        local_scheme.write(solution_pool.best(), certificate_path);
        if (vm.count("print-solution"))
            local_scheme.print(std::cout, solution_pool.best());

    } else if (problem == "permutationflowshopschedulingtt") {
        permutationflowshopschedulingtt::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        auto parameters_local_scheme = read_permutationflowshopschedulingtt_args(local_scheme_argv);
        permutationflowshopschedulingtt::LocalScheme local_scheme(instance, parameters_local_scheme);
        auto solution_pool = run_astarlocalsearch(algorithm_args, local_scheme, info);
        local_scheme.write(solution_pool.best(), certificate_path);
        if (vm.count("print-solution"))
            local_scheme.print(std::cout, solution_pool.best());

    } else {
        std::cerr << "\033[31m" << "ERROR, unknown problem: '" << problem << "'.\033[0m" << std::endl;
        return 1;
    }

    return 0;
}

