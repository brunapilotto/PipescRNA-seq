import argparse
import os
import shutil
import subprocess
import time


def menu():
    parser = argparse.ArgumentParser(
        description="Find long non-coddings in scRNA-seq data"
        )
    parser.add_argument(
        "--fastq",
        help="path to fastq file",
        required=True,
        type=str
        )
    parser.add_argument(
        "--output",
        help="output folder",
        required=True,
        type=str
        )
    args = parser.parse_args()
    return args


def set_work_env(output_folder):
    exec_path = os.path.expanduser(output_folder)
    os.makedirs(exec_path, exist_ok=True)
    os.chdir(exec_path)


if __name__ == "__main__":
    start = time.time()
    options = menu()
    # Set work environment directory
    set_work_env(options.output)
