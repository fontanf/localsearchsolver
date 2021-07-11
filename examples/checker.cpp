#include "examples/knapsackwithconflicts.hpp"
#include "examples/multidimensionalmultiplechoiceknapsack.hpp"
#include "examples/quadraticassignment.hpp"
#include "examples/travellingsalesman.hpp"
#include "examples/schedulingwithsdsttwt.hpp"
#include "examples/permutationflowshopschedulingmakespan.hpp"
#include "examples/permutationflowshopschedulingtt.hpp"

#include "optimizationtools/indexed_set.hpp"

#include <boost/program_options.hpp>

using namespace localsearchsolver;

namespace po = boost::program_options;

int main(int argc, char *argv[])
{

    // Parse program options

    std::string problem = "";
    std::string instance_path = "";
    std::string format = "";
    std::string certificate_path = "";

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("problem,p", po::value<std::string>(&problem)->required(), "set problem (required)")
        ("input,i", po::value<std::string>(&instance_path)->required(), "set input path (required)")
        ("certificate,c", po::value<std::string>(&certificate_path), "set certificate path")
        ("format,f", po::value<std::string>(&format), "set input file format (default: orlibrary)")
        ("print-instance", "")
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

    if (problem == "knapsackwithconflicts") {
        knapsackwithconflicts::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        instance.check(certificate_path);

    } else if (problem == "multidimensionalmultiplechoiceknapsack") {
        multidimensionalmultiplechoiceknapsack::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        instance.check(certificate_path);

    } else if (problem == "quadraticassignment") {
        quadraticassignment::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        instance.check(certificate_path);

    } else if (problem == "travellingsalesman") {
        travellingsalesman::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        instance.check(certificate_path);

    } else if (problem == "schedulingwithsdsttwt") {
        schedulingwithsdsttwt::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        instance.check(certificate_path);

    } else if (problem == "permutationflowshopschedulingmakespan") {
        permutationflowshopschedulingmakespan::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        instance.check(certificate_path);

    } else if (problem == "permutationflowshopschedulingtt") {
        permutationflowshopschedulingtt::Instance instance(instance_path, format);
        if (vm.count("print-instance"))
            std::cout << instance << std::endl;
        instance.check(certificate_path);

    } else {
        std::cerr << "\033[31m" << "ERROR, unknown problem: '" << problem << "'.\033[0m" << std::endl;
        return 1;
    }

    return 0;
}

