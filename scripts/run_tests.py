import argparse
import sys
import os

parser = argparse.ArgumentParser(description='')
parser.add_argument('directory')
parser.add_argument(
        "-t", "--tests",
        type=str,
        nargs='*',
        help='')

args = parser.parse_args()


if args.tests is None or "capacitated-vehicle-routing" in args.tests:
    print("Capacitated vehicle routing problem")
    print("-----------------------------------")
    print()

    capacitated_vehicle_routing_data = [
            (os.path.join("uchoa2014", "X", "X-n101-k25.vrp"), "cvrplib"),
            (os.path.join("uchoa2014", "X", "X-n106-k14.vrp"), "cvrplib"),
            (os.path.join("uchoa2014", "X", "X-n110-k13.vrp"), "cvrplib"),
            (os.path.join("uchoa2014", "X", "X-n115-k10.vrp"), "cvrplib"),
            (os.path.join("uchoa2014", "X", "X-n120-k6.vrp"), "cvrplib"),
            (os.path.join("uchoa2014", "X", "X-n125-k30.vrp"), "cvrplib"),
            (os.path.join("uchoa2014", "X", "X-n129-k18.vrp"), "cvrplib"),
            (os.path.join("uchoa2014", "X", "X-n134-k13.vrp"), "cvrplib"),
            (os.path.join("uchoa2014", "X", "X-n139-k10.vrp"), "cvrplib"),
            (os.path.join("uchoa2014", "X", "X-n143-k7.vrp"), "cvrplib")]
    capacitated_vehicle_routing_main = os.path.join(
            "bazel-bin",
            "examples",
            "capacitated_vehicle_routing_main")
    for instance, instance_format in capacitated_vehicle_routing_data:
        instance_path = os.path.join(
                "data",
                "capacitated_vehicle_routing",
                instance)
        json_output_path = os.path.join(
                args.directory,
                "capacitated_vehicle_routing",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                capacitated_vehicle_routing_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + " --format \"" + instance_format + "\""
                + "  --algorithm genetic-local-search"
                + " --maximum-number-of-iterations 64"
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()


if args.tests is None or "knapsack-with-conflicts" in args.tests:
    print("Knapsack problem with conflicts")
    print("-------------------------------")
    print()

    knapsack_with_conflicts_data = [
            (os.path.join("hifi2006", "I1 - I10", "1I1"), "hifi2006"),
            (os.path.join("hifi2006", "I1 - I10", "2I2"), "hifi2006"),
            (os.path.join("hifi2006", "I1 - I10", "3I3"), "hifi2006"),
            (os.path.join("hifi2006", "I1 - I10", "4I4"), "hifi2006"),
            (os.path.join("hifi2006", "I1 - I10", "5I5"), "hifi2006"),
            (os.path.join("hifi2006", "I1 - I10", "6I1"), "hifi2006"),
            (os.path.join("hifi2006", "I1 - I10", "7I2"), "hifi2006"),
            (os.path.join("hifi2006", "I1 - I10", "8I3"), "hifi2006"),
            (os.path.join("hifi2006", "I1 - I10", "9I4"), "hifi2006"),
            (os.path.join("hifi2006", "I1 - I10", "10I5"), "hifi2006")]
    knapsack_with_conflicts_main = os.path.join(
            "bazel-bin",
            "examples",
            "knapsack_with_conflicts_main")
    for instance, instance_format in knapsack_with_conflicts_data:
        instance_path = os.path.join(
                "data",
                "knapsack_with_conflicts",
                instance)
        json_output_path = os.path.join(
                args.directory,
                "knapsack_with_conflicts",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                knapsack_with_conflicts_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + " --format \"" + instance_format + "\""
                + "  --algorithm genetic-local-search"
                + " --maximum-number-of-iterations 64"
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()


if args.tests is None or "job-sequencing-and-tool-switching" in args.tests:
    print("Job sequencing and tool switching")
    print("---------------------------------")
    print()

    job_sequencing_and_tool_switching_data = [
            (os.path.join("crama1994", "Tabela1", "s1n001.txt"), ""),
            (os.path.join("crama1994", "Tabela1", "s2n001.txt"), ""),
            (os.path.join("crama1994", "Tabela1", "s3n001.txt"), ""),
            (os.path.join("crama1994", "Tabela2", "s1n001.txt"), ""),
            (os.path.join("crama1994", "Tabela2", "s2n001.txt"), ""),
            (os.path.join("crama1994", "Tabela2", "s3n001.txt"), ""),
            (os.path.join("crama1994", "Tabela3", "s1n001.txt"), ""),
            (os.path.join("crama1994", "Tabela3", "s2n001.txt"), ""),
            (os.path.join("crama1994", "Tabela3", "s3n001.txt"), ""),
            (os.path.join("crama1994", "Tabela4", "s1n001.txt"), ""),
            (os.path.join("crama1994", "Tabela4", "s2n001.txt"), ""),
            (os.path.join("crama1994", "Tabela4", "s3n001.txt"), "")]
    job_sequencing_and_tool_switching_main = os.path.join(
            "bazel-bin",
            "examples",
            "job_sequencing_and_tool_switching_main")
    for instance, instance_format in job_sequencing_and_tool_switching_data:
        instance_path = os.path.join(
                "data",
                "job_sequencing_and_tool_switching",
                instance)
        json_output_path = os.path.join(
                args.directory,
                "job_sequencing_and_tool_switching",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                job_sequencing_and_tool_switching_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + " --format \"" + instance_format + "\""
                + "  --algorithm genetic-local-search"
                + " --maximum-number-of-iterations 64"
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()


if args.tests is None or "permutation-flowshop-scheduling-tct" in args.tests:
    print("Permutation flowshop scheduling problem, total completion time")
    print("--------------------------------------------------------------")
    print()

    permutation_flowshop_scheduling_tct_data = [
            (os.path.join("taillard1993", "tai20_5_0.txt"), "default"),
            (os.path.join("taillard1993", "tai20_5_1.txt"), "default"),
            (os.path.join("taillard1993", "tai20_5_2.txt"), "default"),
            (os.path.join("taillard1993", "tai20_5_3.txt"), "default"),
            (os.path.join("taillard1993", "tai20_5_4.txt"), "default"),
            (os.path.join("taillard1993", "tai20_5_5.txt"), "default"),
            (os.path.join("taillard1993", "tai20_5_6.txt"), "default"),
            (os.path.join("taillard1993", "tai20_5_7.txt"), "default"),
            (os.path.join("taillard1993", "tai20_5_8.txt"), "default"),
            (os.path.join("taillard1993", "tai20_5_9.txt"), "default")]
    permutation_flowshop_scheduling_tct_main = os.path.join(
            "bazel-bin",
            "examples",
            "permutation_flowshop_scheduling_tct_main")
    for instance, instance_format in permutation_flowshop_scheduling_tct_data:
        instance_path = os.path.join(
                "data",
                "flowshop_scheduling",
                instance)
        json_output_path = os.path.join(
                args.directory,
                "permutation_flowshop_scheduling_tct",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                permutation_flowshop_scheduling_tct_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + " --format \"" + instance_format + "\""
                + "  --algorithm genetic-local-search"
                + " --maximum-number-of-iterations 64"
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()


if args.tests is None or "permutation-flowshop-scheduling-makespan" in args.tests:
    print("Permutation flowshop scheduling problem, makespan")
    print("-------------------------------------------------")
    print()

    permutation_flowshop_scheduling_makespan_data = [
            (os.path.join("vallada2015", "Small", "VFR10_5_1_Gap.txt"), "default"),
            (os.path.join("vallada2015", "Small", "VFR10_5_2_Gap.txt"), "default"),
            (os.path.join("vallada2015", "Small", "VFR10_5_3_Gap.txt"), "default"),
            (os.path.join("vallada2015", "Small", "VFR10_5_4_Gap.txt"), "default"),
            (os.path.join("vallada2015", "Small", "VFR10_5_5_Gap.txt"), "default"),
            (os.path.join("vallada2015", "Small", "VFR10_5_6_Gap.txt"), "default"),
            (os.path.join("vallada2015", "Small", "VFR10_5_7_Gap.txt"), "default"),
            (os.path.join("vallada2015", "Small", "VFR10_5_8_Gap.txt"), "default"),
            (os.path.join("vallada2015", "Small", "VFR10_5_9_Gap.txt"), "default"),
            (os.path.join("vallada2015", "Small", "VFR10_5_10_Gap.txt"), "default")]
    permutation_flowshop_scheduling_makespan_main = os.path.join(
            "bazel-bin",
            "examples",
            "permutation_flowshop_scheduling_makespan_main")
    for instance, instance_format in permutation_flowshop_scheduling_makespan_data:
        instance_path = os.path.join(
                "data",
                "flowshop_scheduling",
                instance)
        json_output_path = os.path.join(
                args.directory,
                "permutation_flowshop_scheduling_makespan",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                permutation_flowshop_scheduling_makespan_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + " --format \"" + instance_format + "\""
                + "  --algorithm best-first-local-search"
                + " --maximum-number-of-nodes 64"
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()


if args.tests is None or "scheduling-with-sdst-twt" in args.tests:
    print("Single machine scheduling problem with sequence-dependent setup times, total weighted tardiness")
    print("-----------------------------------------------------------------------------------------------")
    print()

    scheduling_with_sdst_twt_data = [
            (os.path.join("cicirello2005", "wt_sds_10.instance"), ""),
            (os.path.join("cicirello2005", "wt_sds_20.instance"), ""),
            (os.path.join("cicirello2005", "wt_sds_30.instance"), ""),
            (os.path.join("cicirello2005", "wt_sds_40.instance"), ""),
            (os.path.join("cicirello2005", "wt_sds_50.instance"), ""),
            (os.path.join("cicirello2005", "wt_sds_60.instance"), ""),
            (os.path.join("cicirello2005", "wt_sds_70.instance"), ""),
            (os.path.join("cicirello2005", "wt_sds_80.instance"), ""),
            (os.path.join("cicirello2005", "wt_sds_90.instance"), ""),
            (os.path.join("cicirello2005", "wt_sds_100.instance"), ""),
            (os.path.join("cicirello2005", "wt_sds_110.instance"), ""),
            (os.path.join("cicirello2005", "wt_sds_120.instance"), "")]
    scheduling_with_sdst_twt_main = os.path.join(
            "bazel-bin",
            "examples",
            "scheduling_with_sdst_twt_main")
    for instance, instance_format in scheduling_with_sdst_twt_data:
        instance_path = os.path.join(
                "data",
                "scheduling_with_sdst_twt",
                instance)
        json_output_path = os.path.join(
                args.directory,
                "scheduling_with_sdst_twt",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                scheduling_with_sdst_twt_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + " --format \"" + instance_format + "\""
                + "  --algorithm genetic-local-search"
                + " --maximum-number-of-iterations 128"
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()


if args.tests is None or "sequential-ordering" in args.tests:
    print("Sequential ordering problem")
    print("---------------------------")
    print()

    sequential_ordering_data = [
            (os.path.join("soplib", "R.200.100.1.sop"), "soplib"),
            (os.path.join("soplib", "R.200.100.15.sop"), "soplib"),
            (os.path.join("soplib", "R.200.100.30.sop"), "soplib"),
            (os.path.join("soplib", "R.200.100.60.sop"), "soplib"),
            (os.path.join("soplib", "R.200.1000.1.sop"), "soplib"),
            (os.path.join("soplib", "R.200.1000.15.sop"), "soplib"),
            (os.path.join("soplib", "R.200.1000.30.sop"), "soplib"),
            (os.path.join("soplib", "R.200.1000.60.sop"), "soplib")]
    sequential_ordering_main = os.path.join(
            "bazel-bin",
            "examples",
            "sequential_ordering_main")
    for instance, instance_format in sequential_ordering_data:
        instance_path = os.path.join(
                "data",
                "sequential_ordering",
                instance)
        json_output_path = os.path.join(
                args.directory,
                "sequential_ordering",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                sequential_ordering_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + " --format \"" + instance_format + "\""
                + "  --algorithm best-first-local-search"
                + " --maximum-number-of-nodes 16"
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()


if args.tests is None or "team-orienteering" in args.tests:
    print("Team orienteering problem")
    print("-------------------------")
    print()

    team_orienteering_data = [
            (os.path.join("chao1996", "Set_102_234", "p7.2.a.txt"), ""),
            (os.path.join("chao1996", "Set_102_234", "p7.2.b.txt"), ""),
            (os.path.join("chao1996", "Set_102_234", "p7.2.c.txt"), ""),
            (os.path.join("chao1996", "Set_102_234", "p7.3.d.txt"), ""),
            (os.path.join("chao1996", "Set_102_234", "p7.3.e.txt"), ""),
            (os.path.join("chao1996", "Set_102_234", "p7.3.f.txt"), ""),
            (os.path.join("chao1996", "Set_102_234", "p7.4.i.txt"), ""),
            (os.path.join("chao1996", "Set_102_234", "p7.4.j.txt"), ""),
            (os.path.join("chao1996", "Set_102_234", "p7.4.k.txt"), "")]
    team_orienteering_main = os.path.join(
            "bazel-bin",
            "examples",
            "team_orienteering_main")
    for instance, instance_format in team_orienteering_data:
        instance_path = os.path.join(
                "data",
                "team_orienteering",
                instance)
        json_output_path = os.path.join(
                args.directory,
                "team_orienteering",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                team_orienteering_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + " --format \"" + instance_format + "\""
                + "  --algorithm best-first-local-search"
                + " --maximum-number-of-nodes 128"
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()


if args.tests is None or "time-dependent-orienteering" in args.tests:
    print("Time-dependent orienteering problem")
    print("-----------------------------------")
    print()

    time_dependent_orienteering_data = [
            (os.path.join("verbeeck2014", "dataset 1", "OP_instances", "p1.1.a.txt"), ""),
            (os.path.join("verbeeck2014", "dataset 2", "OP_instances", "p2.1.a.txt"), ""),
            (os.path.join("verbeeck2014", "dataset 3", "OP_instances", "p3.1.a.txt"), ""),
            (os.path.join("verbeeck2014", "dataset 4", "OP_instances", "p4.1.a.txt"), ""),
            (os.path.join("verbeeck2014", "dataset 5", "OP_instances", "p5.1.a.txt"), ""),
            (os.path.join("verbeeck2014", "dataset 6", "OP_instances", "p6.1.a.txt"), ""),
            (os.path.join("verbeeck2014", "dataset 7", "OP_instances", "p7.1.a.txt"), "")]
    time_dependent_orienteering_main = os.path.join(
            "bazel-bin",
            "examples",
            "time_dependent_orienteering_main")
    for instance, instance_format in time_dependent_orienteering_data:
        instance_path = os.path.join(
                "data",
                "time_dependent_orienteering",
                instance)
        json_output_path = os.path.join(
                args.directory,
                "time_dependent_orienteering",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                time_dependent_orienteering_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + " --format \"" + instance_format + "\""
                + "  --algorithm best-first-local-search"
                + " --maximum-number-of-nodes 128"
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()


if args.tests is None or "vehicle-routing-with-time-windows" in args.tests:
    print("Vehicle routing problem with time-windows")
    print("-----------------------------------------")
    print()

    vehicle_routing_with_time_windows_data = [
            (os.path.join("solomon1987", "C101.txt"), ""),
            (os.path.join("solomon1987", "C102.txt"), ""),
            (os.path.join("solomon1987", "C103.txt"), ""),
            (os.path.join("solomon1987", "R101.txt"), ""),
            (os.path.join("solomon1987", "R102.txt"), ""),
            (os.path.join("solomon1987", "R103.txt"), ""),
            (os.path.join("solomon1987", "RC101.txt"), ""),
            (os.path.join("solomon1987", "RC102.txt"), ""),
            (os.path.join("solomon1987", "RC103.txt"), "")]
    vehicle_routing_with_time_windows_main = os.path.join(
            "bazel-bin",
            "examples",
            "vehicle_routing_with_time_windows_main")
    for instance, instance_format in vehicle_routing_with_time_windows_data:
        instance_path = os.path.join(
                "data",
                "vehicle_routing_with_time_windows",
                instance)
        json_output_path = os.path.join(
                args.directory,
                "vehicle_routing_with_time_windows",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                vehicle_routing_with_time_windows_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + " --format \"" + instance_format + "\""
                + "  --algorithm genetic-local-search"
                + " --maximum-number-of-iterations 64"
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()
