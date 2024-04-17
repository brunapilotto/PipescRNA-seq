import argparse
import os
import shutil
import subprocess
import time

from termcolor import colored
from rich.console import Console

from utils.check_input import CheckInput


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


def main():
    start = time.time()
    options = menu()
    console = Console()
    with console.status("[bold green]Searching fastqs..."):
        check_input = CheckInput(options.fastq)
        fastqs = check_input.check_fastq()

    # Set work environment directory
    set_work_env(options.output)


if __name__ == "__main__":
    main()
