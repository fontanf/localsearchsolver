config_setting(name = "cplex_build",       values = {"define": "cplex=true"},       visibility = ["//visibility:public"])
config_setting(name = "gurobi_build",      values = {"define": "gurobi=true"},      visibility = ["//visibility:public"])

load("//examples:variables.bzl",
         "CPLEX_COPTS", "CPLEX_DEP",
         "GUROBI_COPTS", "GUROBI_DEP")

cc_library(
        name = "knapsack_with_conflicts",
        hdrs = ["knapsack_with_conflicts.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:knapsack_with_conflicts",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "knapsack_with_conflicts_main",
        srcs = ["knapsack_with_conflicts_main.cpp"],
        deps = [
                ":knapsack_with_conflicts",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "multidimensional_multiple_choice_knapsack",
        hdrs = ["multidimensional_multiple_choice_knapsack.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:multidimensional_multiple_choice_knapsack",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "multidimensional_multiple_choice_knapsack_main",
        srcs = ["multidimensional_multiple_choice_knapsack_main.cpp"],
        deps = [
                ":multidimensional_multiple_choice_knapsack",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "quadratic_assignment",
        hdrs = ["quadratic_assignment.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:quadratic_assignment",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "quadratic_assignment_main",
        srcs = ["quadratic_assignment_main.cpp"],
        deps = [
                ":quadratic_assignment",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "sequential_ordering",
        hdrs = ["sequential_ordering.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:sequential_ordering",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "sequential_ordering_main",
        srcs = ["sequential_ordering_main.cpp"],
        deps = [
                ":sequential_ordering",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "scheduling_with_sdst_twt",
        hdrs = ["scheduling_with_sdst_twt.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:scheduling_with_sdst_twt",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "scheduling_with_sdst_twt_main",
        srcs = ["scheduling_with_sdst_twt_main.cpp"],
        deps = [
                ":scheduling_with_sdst_twt",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "order_acceptance_and_scheduling",
        hdrs = ["order_acceptance_and_scheduling.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:order_acceptance_and_scheduling",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "order_acceptance_and_scheduling_main",
        srcs = ["order_acceptance_and_scheduling_main.cpp"],
        deps = [
                ":order_acceptance_and_scheduling",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "job_sequencing_and_tool_switching",
        hdrs = ["job_sequencing_and_tool_switching.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:job_sequencing_and_tool_switching",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "job_sequencing_and_tool_switching_main",
        srcs = ["job_sequencing_and_tool_switching_main.cpp"],
        deps = [
                ":job_sequencing_and_tool_switching",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "time_dependent_orienteering",
        hdrs = ["time_dependent_orienteering.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:time_dependent_orienteering",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "time_dependent_orienteering_main",
        srcs = ["time_dependent_orienteering_main.cpp"],
        deps = [
                ":time_dependent_orienteering",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "permutation_flowshop_scheduling_makespan",
        hdrs = ["permutation_flowshop_scheduling_makespan.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:permutation_flowshop_scheduling_makespan",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "permutation_flowshop_scheduling_makespan_main",
        srcs = ["permutation_flowshop_scheduling_makespan_main.cpp"],
        deps = [
                ":permutation_flowshop_scheduling_makespan",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "permutation_flowshop_scheduling_tct",
        hdrs = ["permutation_flowshop_scheduling_tct.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:permutation_flowshop_scheduling_tct",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "permutation_flowshop_scheduling_tct_main",
        srcs = ["permutation_flowshop_scheduling_tct_main.cpp"],
        deps = [
                ":permutation_flowshop_scheduling_tct",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "permutation_flowshop_scheduling_tt",
        hdrs = ["permutation_flowshop_scheduling_tt.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:permutation_flowshop_scheduling_tt",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "permutation_flowshop_scheduling_tt_main",
        srcs = ["permutation_flowshop_scheduling_tt_main.cpp"],
        deps = [
                ":permutation_flowshop_scheduling_tt",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "distributed_pfss_tct",
        hdrs = ["distributed_pfss_tct.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:distributed_pfss_tct",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "distributed_pfss_tct_main",
        srcs = ["distributed_pfss_tct_main.cpp"],
        deps = [
                ":distributed_pfss_tct",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "traveling_repairman",
        hdrs = ["traveling_repairman.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:traveling_repairman",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "traveling_repairman_main",
        srcs = ["traveling_repairman_main.cpp"],
        deps = [
                ":traveling_repairman",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "capacitated_vehicle_routing",
        hdrs = ["capacitated_vehicle_routing.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:capacitated_vehicle_routing",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "capacitated_vehicle_routing_main",
        srcs = ["capacitated_vehicle_routing_main.cpp"],
        deps = [
                ":capacitated_vehicle_routing",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)

cc_binary(
        name = "capacitated_vehicle_routing_dimacs_main",
        srcs = ["capacitated_vehicle_routing_dimacs_main.cpp"],
        deps = [
                ":capacitated_vehicle_routing",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "vehicle_routing_with_time_windows",
        hdrs = ["vehicle_routing_with_time_windows.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:vehicle_routing_with_time_windows",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "vehicle_routing_with_time_windows_main",
        srcs = ["vehicle_routing_with_time_windows_main.cpp"],
        deps = [
                ":vehicle_routing_with_time_windows",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)

cc_binary(
        name = "vehicle_routing_with_time_windows_dimacs_main",
        srcs = ["vehicle_routing_with_time_windows_dimacs_main.cpp"],
        deps = [
                ":vehicle_routing_with_time_windows",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "traveling_salesman_with_release_dates",
        hdrs = ["traveling_salesman_with_release_dates.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:traveling_salesman_with_release_dates",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "traveling_salesman_with_release_dates_main",
        srcs = ["traveling_salesman_with_release_dates_main.cpp"],
        deps = [
                ":traveling_salesman_with_release_dates",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "batch_scheduling_total_weighted_tardiness",
        hdrs = ["batch_scheduling_total_weighted_tardiness.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:batch_scheduling_total_weighted_tardiness",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "batch_scheduling_total_weighted_tardiness_main",
        srcs = ["batch_scheduling_total_weighted_tardiness_main.cpp"],
        deps = [
                ":batch_scheduling_total_weighted_tardiness",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)
