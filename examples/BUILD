config_setting(name = "cplex_build",       values = {"define": "cplex=true"},       visibility = ["//visibility:public"])
config_setting(name = "gurobi_build",      values = {"define": "gurobi=true"},      visibility = ["//visibility:public"])

load("//examples:variables.bzl",
         "CPLEX_COPTS", "CPLEX_DEP",
         "GUROBI_COPTS", "GUROBI_DEP")

cc_library(
        name = "knapsackwithconflicts",
        hdrs = ["knapsackwithconflicts.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:knapsackwithconflicts",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "knapsackwithconflicts_main",
        srcs = ["knapsackwithconflicts_main.cpp"],
        deps = [
                ":knapsackwithconflicts",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "multidimensionalmultiplechoiceknapsack",
        hdrs = ["multidimensionalmultiplechoiceknapsack.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:multidimensionalmultiplechoiceknapsack",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "multidimensionalmultiplechoiceknapsack_main",
        srcs = ["multidimensionalmultiplechoiceknapsack_main.cpp"],
        deps = [
                ":multidimensionalmultiplechoiceknapsack",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "quadraticassignment",
        hdrs = ["quadraticassignment.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:quadraticassignment",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "quadraticassignment_main",
        srcs = ["quadraticassignment_main.cpp"],
        deps = [
                ":quadraticassignment",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "schedulingwithsdsttwt",
        hdrs = ["schedulingwithsdsttwt.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:schedulingwithsdsttwt",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "schedulingwithsdsttwt_main",
        srcs = ["schedulingwithsdsttwt_main.cpp"],
        deps = [
                ":schedulingwithsdsttwt",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "orderacceptanceandscheduling",
        hdrs = ["orderacceptanceandscheduling.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:orderacceptanceandscheduling",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "orderacceptanceandscheduling_main",
        srcs = ["orderacceptanceandscheduling_main.cpp"],
        deps = [
                ":orderacceptanceandscheduling",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "permutationflowshopschedulingmakespan",
        hdrs = ["permutationflowshopschedulingmakespan.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:permutationflowshopschedulingmakespan",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "permutationflowshopschedulingmakespan_main",
        srcs = ["permutationflowshopschedulingmakespan_main.cpp"],
        deps = [
                ":permutationflowshopschedulingmakespan",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "permutationflowshopschedulingtt",
        hdrs = ["permutationflowshopschedulingtt.hpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:permutationflowshopschedulingtt",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "permutationflowshopschedulingtt_main",
        srcs = ["permutationflowshopschedulingtt_main.cpp"],
        deps = [
                ":permutationflowshopschedulingtt",
                "//localsearchsolver:read_args",
        ],
        linkopts = select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)


cc_library(
        name = "roadef2020",
        hdrs = ["roadef2020.hpp"],
        srcs = ["roadef2020.cpp"],
        deps = [
                "//localsearchsolver:localsearchsolver",
                "@orproblems//orproblems:roadef2020",
                "@optimizationtools//optimizationtools:indexed_collections",
        ] + CPLEX_DEP + GUROBI_DEP,
        copts = [
                #"-fprofile-generate=/home/florian/Dev/localsearchsolver/profile_dir",
                #"-fprofile-use=/home/florian/Dev/localsearchsolver/profile_dir",
                #"-fprofile-correction",
        ] + CPLEX_COPTS + GUROBI_COPTS,
        linkopts = [
                #"-fprofile-generate=/home/florian/Dev/localsearchsolver/profile_dir",
                #"-fprofile-use=/home/florian/Dev/localsearchsolver/profile_dir",
                #"-fprofile-correction",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "roadef2020_main",
        srcs = ["roadef2020_main.cpp"],
        deps = [
                ":roadef2020",
                "@boost//:program_options",
        ],
        copts = [
                #"-fprofile-generate=/home/florian/Dev/localsearchsolver/profile_dir",
                #"-fprofile-use=/home/florian/Dev/localsearchsolver/profile_dir",
                #"-fprofile-correction",
        ],
        linkopts = [
                #"-fprofile-generate=/home/florian/Dev/localsearchsolver/profile_dir",
                #"-fprofile-use=/home/florian/Dev/localsearchsolver/profile_dir",
                #"-fprofile-correction",
        ] + select({
                "@bazel_tools//src/conditions:windows": [],
                "//conditions:default":                 ["-lpthread"],
        }),
)

