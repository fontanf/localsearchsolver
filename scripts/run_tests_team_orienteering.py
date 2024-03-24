import run_tests
import os


commands = run_tests.generate_commands(
        os.path.join(
            "bazel-bin",
            "examples",
            "team_orienteering_main")
        + " --algorithm best-first-local-search"
        + " --maximum-number-of-nodes 1024",
        os.path.join(
            "..",
            "ordata",
            "routing",
            "team_orienteering",
            "data.csv"),
        os.path.join("test_results", "team_orienteering"),
        "row['Dataset'] in"
        " ['chao1996_4', 'chao1996_5', 'chao1996_6', 'chao1996_7']")


if __name__ == "__main__":
    for command in commands:
        run_tests.run(command)
