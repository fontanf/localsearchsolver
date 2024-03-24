import sys
import os
sys.path.insert(1, os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..",
    "bazel-localsearchsolver",
    "external",
    "optimizationtools",
    "scripts"))
import process_tests

process_tests.process_tests(
        os.path.join("test_results_ref", "scheduling_with_sdst_twt"),
        os.path.join("test_results", "scheduling_with_sdst_twt"),
        ["Value"],
        ["Time"])

process_tests.process_tests(
        os.path.join("test_results_ref", "team_orienteering"),
        os.path.join("test_results", "team_orienteering"),
        ["Value"],
        ["Time"])
