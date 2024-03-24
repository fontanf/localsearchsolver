import run_tests
import os


commands = run_tests.generate_commands(
        os.path.join(
            "bazel-bin",
            "examples",
            "scheduling_with_sdst_twt_main")
        + " --algorithm best-first-local-search"
        + " --maximum-number-of-nodes 1024",
        os.path.join(
            "..",
            "ordata",
            "scheduling",
            "scheduling_with_sdst_twt",
            "data.csv"),
        os.path.join("test_results", "scheduling_with_sdst_twt"))


if __name__ == "__main__":
    for command in commands:
        run_tests.run(command)
