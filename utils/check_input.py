import os

from pathlib import Path
from termcolor import colored


class CheckInput:
    def __init__(self, input_path):
        self.input = Path(os.path.abspath(input_path))

    def check_fastq(self):
        fastqs = sorted(
            str(fastq_path)
            for fastq_path in self.input.glob("*.fastq")
            if fastq_path.is_file()
        )
        for fastq in fastqs:
            print(f"Found fastq file: {colored(fastq, "green")}")

        if not fastqs:
            print(colored(f"No fastq files found in {self.input}", "red"))
            exit()
        return fastqs
