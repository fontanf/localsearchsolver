import os
import csv


def generate_commands(
        main,
        datacsv_path,
        output_directory,
        instance_filter=True):

    commands = []

    reader = csv.DictReader(open(datacsv_path))  # noqa: F841
    rows_filtered = eval("filter(lambda row: %s, reader)" % (instance_filter))

    # Loop through instances.
    for row in rows_filtered:

        instance_path = os.path.join(
                os.path.dirname(datacsv_path),
                row["Path"])

        instance_format = row["Format"]

        json_output_path = os.path.join(
                output_directory,
                "outputs",
                row["Dataset"],
                row["Path"] + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))

        certificate_path = os.path.join(
                output_directory,
                "certificates",
                row["Dataset"],
                row["Path"] + "_solution.txt")
        if not os.path.exists(os.path.dirname(certificate_path)):
            os.makedirs(os.path.dirname(certificate_path))

        command = (
                main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + "  --format \"" + instance_format + "\""
                + "  --output \"" + json_output_path + "\""
                + "  --certificate \"" + certificate_path + "\""
                )

        commands.append(command)

    return commands


def run(command):
    print()
    print(command)
    os.system(command)
