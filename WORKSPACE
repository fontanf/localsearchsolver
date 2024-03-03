load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

git_repository(
    name = "googletest",
    remote = "https://github.com/google/googletest.git",
    commit = "58d77fa8070e8cec2dc1ed015d66b454c8d78850",
    shallow_since = "1656350095 -0400",
)

git_repository(
    name = "com_github_nelhage_rules_boost",
    remote = "https://github.com/nelhage/rules_boost",
    commit = "e83dfef18d91a3e35c8eac9b9aeb1444473c0efd",
    shallow_since = "1671181466 +0000",
)
load("@com_github_nelhage_rules_boost//:boost/boost.bzl", "boost_deps")
boost_deps()

http_archive(
    name = "json",
    build_file_content = """
cc_library(
        name = "json",
        hdrs = ["single_include/nlohmann/json.hpp"],
        visibility = ["//visibility:public"],
        strip_include_prefix = "single_include/"
)
""",
    urls = ["https://github.com/nlohmann/json/releases/download/v3.7.3/include.zip"],
    sha256 = "87b5884741427220d3a33df1363ae0e8b898099fbc59f1c451113f6732891014",
)

git_repository(
    name = "optimizationtools",
    remote = "https://github.com/fontanf/optimizationtools.git",
    commit = "0bd9be63666a4ddf66c2ca89ca984bda59824e01",
)

local_repository(
    name = "optimizationtools_",
    path = "../optimizationtools/",
)

git_repository(
    name = "orproblems",
    remote = "https://github.com/fontanf/orproblems.git",
    commit = "fe9dcd0f6be65154a53f205cb30485adfb07095f",
)

local_repository(
    name = "orproblems_",
    path = "../orproblems/",
)

git_repository(
    name = "travelingsalesmansolver",
    remote = "https://github.com/fontanf/travelingsalesmansolver.git",
    commit = "97a12111e2d8dcdb2072570b09667ec9073f6000",
)

local_repository(
    name = "travelingsalesmansolver_",
    path = "../travelingsalesmansolver/",
)

http_archive(
    name = "simdjson",
    build_file_content = """
cc_library(
        name = "simdjson",
        hdrs = ["simdjson-0.9.2/singleheader/simdjson.h"],
        srcs = ["simdjson-0.9.2/singleheader/simdjson.cpp"],
        visibility = ["//visibility:public"],
        strip_include_prefix = "simdjson-0.9.2/singleheader/"
)
""",
    urls = ["https://github.com/simdjson/simdjson/archive/refs/tags/v0.9.2.zip"],
    sha256 = "4d4e4a9d64b55839ec6c80f0ba54261139543b7ffcf84bceede681ccb95fe94b",
)

new_local_repository(
    name = "gurobi",
    path = "/home/florian/Programmes/gurobi811/linux64/",
    build_file_content = """
cc_library(
    name = "gurobi",
    hdrs = [
            "include/gurobi_c.h",
            "include/gurobi_c++.h",
    ],
    strip_include_prefix = "include/",
    srcs = [
            "lib/libgurobi_c++.a",
            "lib/libgurobi81.so",
    ],
    visibility = ["//visibility:public"],
)
""",
)

new_local_repository(
    name = "cplex",
    path = "/opt/ibm/ILOG/CPLEX_Studio129/",
    build_file_content = """
cc_library(
    name = "concert",
    hdrs = glob(["concert/include/ilconcert/**/*.h"], exclude_directories = 0),
    strip_include_prefix = "concert/include/",
    srcs = ["concert/lib/x86-64_linux/static_pic/libconcert.a"],
    linkopts = [
            "-lm",
            "-lpthread",
            "-ldl",
    ],
    visibility = ["//visibility:public"],
)
cc_library(
    name = "cplex",
    hdrs = glob(["cplex/include/ilcplex/*.h"]),
    strip_include_prefix = "cplex/include/",
    srcs = [
            "cplex/lib/x86-64_linux/static_pic/libilocplex.a",
            "cplex/lib/x86-64_linux/static_pic/libcplex.a",
    ],
    deps = [":concert"],
    visibility = ["//visibility:public"],
)
cc_library(
    name = "cpoptimizer",
    hdrs = glob(["cpoptimizer/include/ilcp/*.h"]),
    strip_include_prefix = "cpoptimizer/include/",
    srcs = ["cpoptimizer/lib/x86-64_linux/static_pic/libcp.a"],
    deps = [":cplex"],
    visibility = ["//visibility:public"],
)
""",
)
