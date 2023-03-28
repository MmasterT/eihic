#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create the run_config.yaml
"""

# authorship
__author__ = "Mariano Olivera Fedi"
__maintainer__ = "Mariano Olivera Fedi"
__email__ = "mariano.olivera.fedi@hotmail.com"

# import libraries
import argparse
from pathlib import Path
import sys
import os
import yaml
import csv

# import pandas as pd
from eihic import (
    DEFAULT_CONFIG_FILE,
    DEFAULT_HPC_CONFIG_FILE,
    DEFAULT_SAMPLE_CSV_FILE,
)

# get script name
script = Path(sys.argv[0]).name

cwd = os.getcwd()


class HI_CCONFIGURE:
    @staticmethod
    def check_exits(file_path):
        if not Path(file_path).exists():
            raise FileNotFoundError(f"File not found: '{file_path}'")
        return str(Path(file_path).resolve())

    def __init__(self, args):
        self.args = args
        self.args.samples_csv = HI_CCONFIGURE.check_exits(self.args.samples_csv)
        self.run_config = dict()
        self.args.output = str(Path(self.args.output).resolve())
        self.args.logs = str(Path(self.args.output).resolve() / "logs")
        self.run_config_file = str()
        self.args.curation = bool(self.args.curation)
        self.args.bwa_mem2 = bool(self.args.bwa_mem2)

    def process_run_config(self):
        with open(DEFAULT_CONFIG_FILE, "r") as fh:
            try:
                self.run_config = yaml.safe_load(fh)
            except yaml.YAMLError as err:
                print(err)
        if not self.run_config:
            raise ValueError(
                f"No information processed from run_config file - '{self.args.run_config}'"
            )

    def process_samples_csv(self):
        data = {}
        with open(self.args.samples_csv, newline='') as f:
            reader = csv.reader(f)
            file = list(reader)


        #  Forward reads and sample reads should be =
        if len(file[0]) != len(file[1]):
            print("Error: the number pair-end samples are not the same for R1 and R2")
            print(data)
            exit() 

        # Save the data to a list: R1 forward reads, R2 reverse reads, and sample name/ organism name
       
        data["R1"] = file[0]
        data["R2"] = file[1]
        data["reference"] = file[2][0]
        data["organism"] = file[3][0]
        
        if self.args.curation:
            data["long_reads"] = file[4]
            if len(file) != 5:
                print("Error: Your data sample_csv file does not have the right amount of fields. Please check the documentation and your file.")
                print(data)
                exit() 
        

        # add details
        for sample in data["R1"]:

            if not Path(sample).is_file():
                raise ValueError(
                    f"Sample '{sample}' R1 path does not exist {sample} ."
            )
        for sample in data["R2"]:

            if not Path(sample).is_file():
                raise ValueError(
                    f"Sample '{sample}' R1 path does not exist {sample} ."
            )

        # link the reads
        read_dir = Path(self.args.output).joinpath("reads/short_reads")
        Path(read_dir).mkdir(parents=True, exist_ok=True)

        # r1_cmd = f"cd {read_dir} && ln -s {R1[sample]} x[SID].R1.{ext}"
        for sample in range(len(data["R1"])):
            r1_name = os.path.basename(os.path.abspath(data["R1"][sample]))
            r2_name = os.path.basename(os.path.abspath(data["R2"][sample]))
            r1_path = Path(read_dir).joinpath(r1_name)
            r2_path = Path(read_dir).joinpath(r2_name)

            if self.args.force_reconfiguration:
                if Path(r1_path).is_symlink():
                    print(f"Unlinking '{r1_path}' with --force")
                    r1_path.unlink()
                if Path(r2_path).is_symlink():
                    print(f"Unlinking '{r2_path}' with --force")
                    r2_path.unlink()

            Path(r1_path).symlink_to(os.path.abspath(data["R1"][sample]))
            Path(r2_path).symlink_to(os.path.abspath(data["R2"][sample]))
        
        read_dir = Path(self.args.output).joinpath("reads/long_reads")
        Path(read_dir).mkdir(parents=True, exist_ok=True)
        for sample in range(len(data["long_reads"])):
            sample_name = os.path.basename(os.path.abspath(data["long_reads"][sample]))
            sample_path = Path(read_dir).joinpath(sample_name)
            

            if self.args.force_reconfiguration:
                if Path(sample_path).is_symlink():
                    print(f"Unlinking '{sample_path}' with --force")
                    sample_path.unlink()

            Path(sample_path).symlink_to(os.path.abspath(data["long_reads"][sample]))
        
        self.run_config["input_samples"] = data

        #Create reference dir
        reference_dir = Path(self.args.output).joinpath("reference/genome")
        Path(reference_dir).mkdir(parents=True, exist_ok=True)
	
        reference_path = Path(reference_dir).joinpath(data["organism"] + ".fasta")
	
        if self.args.force_reconfiguration:
            if reference_path.is_symlink():
                print(f"Unlinking '{reference_path}' with --force")
                reference_path.unlink()
        Path(reference_dir).joinpath(data["organism"] + ".fasta").symlink_to(os.path.abspath(data["reference"]))

        #Creates workflow and tmp directories
        tmp_dir = Path(self.args.output).joinpath("tmp")
        workflow_dir = Path(self.args.output).joinpath("workflow")

        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
        Path(workflow_dir).mkdir(parents=True, exist_ok=True)

    def write_run_config(self):
        # output directory
        self.run_config["samples_csv"] = self.args.samples_csv
        self.run_config["output"] = self.args.output
        self.run_config["logs"] = self.args.logs

        Path(self.args.logs).mkdir(parents=True, exist_ok=True)

        # process samples csv to yaml
        self.process_samples_csv()

        # modify any additional defaults
        self.run_config["hpc_config"] = DEFAULT_HPC_CONFIG_FILE
        self.run_config["jira"]["jira_id"] = self.args.jira
        
        # write bwa-mem2 option
        if self.args.bwa_mem2:
            self.run_config["bwa-mem2"] = 'True'
        else:
            self.run_config["bwa-mem2"] = 'False'
        
        # write the new run config file
        with open(self.run_config_file, "w") as fh:
            yaml.dump(self.run_config, fh, sort_keys=False)

    def run(self):
        self.process_run_config()
        self.run_config_file = os.path.join(self.args.output, "run_config.yaml")
        self.write_run_config()
        print(f"\nGreat! Created run_config file: '{self.run_config_file}'\n")


def main():

    parser = argparse.ArgumentParser(description="Script to create the run_config.yaml")
    parser.add_argument(
        "-s",
	"--samples_csv",
        help=f"Provide sample information in csv format. Please refer to the sample file is here: {DEFAULT_CONFIG_FILE}, for more information above the tsv format. A template is provided here {DEFAULT_SAMPLE_CSV_FILE}.  (default: %(default)s)",
    )
    parser.add_argument(
        "-c",
        "--curation",
        action="store_true",
        help="Use this flag if you have your final scaffold and want to run the curation pipeline. Bare in mind that the sample_csv file requires an extra field with the long reads file paths as a fifth line. (default: %(default)s)"
    )
    parser.add_argument(
        "--jira",
        help="Provide JIRA id for posting job summary. E.g., PPBFX-611 (default: %(default)s)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=os.path.join(cwd, "output"),
        help="Provide output directory (default: %(default)s)",
    )
    parser.add_argument(
        "-f",
        "--force-reconfiguration",
        action="store_true",
        help="Force reconfiguration (default: %(default)s)",
    )
    parser.add_argument(
        "-bm2",
        "--bwa_mem2",
        action="store_true",
        help="Use bwa-mem2 insted of bwa mem. This option use a lot more RAM, use with precaution. (default: %(default)s)",
    )
    args = parser.parse_args()
    HI_CCONFIGURE(args).run()


if __name__ == "__main__":
    main()