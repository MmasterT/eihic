# -*- coding: utf-8 -*-
"""
EI Hi-C Pipeline
"""

# authorship
__author__ = "Marino Olivera Fedi"
__maintainer__ = "Mariano Olivera Fedi"
__email__ = "mariano.olivera.fedi@hotmail.com"

# import libraries
import argparse
import os
import sys
import subprocess
import yaml
import pkg_resources
from pathlib import Path
import datetime
import time

from eihic import __version__

from eihic.scripts.jiracomms import JiraInfo
from eihic.scripts.hi_c_configure import HI_CCONFIGURE

from eihic import (
    DEFAULT_CONFIG_FILE,
    DEFAULT_HPC_CONFIG_FILE,
    DEFAULT_SAMPLE_CSV_FILE,
)

# Request min version of snakemake
from snakemake.utils import min_version

min_version("7.10")

NOW = datetime.datetime.fromtimestamp(time.time()).strftime("%Y%m%d_%H%M%S")

# get script name
script = os.path.basename(sys.argv[0])
script_dir = pkg_resources.resource_filename("eihic", "workflow")
cwd = os.getcwd()


def command_configure(args):
    print("Running configure..")
    # check if we have a config file
    run_config = False
    run_config_file = os.path.join(args.output, "run_config.yaml")
    try:
        run_config = os.path.exists(run_config_file)
    except:
        pass
    if run_config is False or args.force_reconfiguration:
        HI_CCONFIGURE(args).run()
    elif run_config:
        print(
            f"\nWARNING: Configuration file '{run_config_file}' already present. Please set --force-reconfiguration/-f to override this.\n"
        )


def command_run(args):
    HI_C(args).run()


class HI_C:
    @staticmethod
    def check_exits(file_path):
        if not Path(file_path).exists():
            raise FileNotFoundError(f"File not found: '{file_path}'")
        return str(Path(file_path).resolve())

    def __init__(self, args):
        print("Initialising pipeline")
        self.args = args
        self.run_config = HI_C.check_exits(args.run_config)
        self.hpc_config = HI_C.check_exits(args.hpc_config)
        self.jobs = args.jobs
        self.latency_wait = args.latency_wait
        self.no_posting = args.no_posting_off
        self.verbose = args.verbose
        self.dry_run = args.dry_run
        self.library = args.library
        self.curation = args.curation
        self.loaded_run_config = yaml.load(
            open(self.run_config), Loader=yaml.SafeLoader
        )
        # replace with new HPC config
        if self.hpc_config and self.hpc_config != DEFAULT_HPC_CONFIG_FILE:
            self.loaded_run_config["hpc_config"] = self.hpc_config
        self.run_config = self.run_config.replace(".yaml", ".{}.yaml".format(NOW))
        # write the new run config file
        with open(self.run_config, "w") as fh:
            yaml.dump(self.loaded_run_config, fh, sort_keys=False)

        self.jira_id = self.loaded_run_config["jira"]["jira_id"]
        self.output = self.loaded_run_config["output"]
        self.samples_csv = self.loaded_run_config["samples_csv"]
        self.logs = self.loaded_run_config["logs"]
        # Load the config file
        self.config = yaml.load(open(DEFAULT_CONFIG_FILE), Loader=yaml.SafeLoader)

        # Gets JIRA ticket from server (or makes one from args provided if args are set appropriately)
        if self.jira_id:
            JiraInfo(self.jira_id).initialise(config=self.config)

    def run(self):
        print("Running the pipeline..")
        cmd = None
        if self.dry_run:
            print("Enabling dry run..")
            if self.curation:
                self.library = "curation_" + self.library
            cmd = (
                f"snakemake --snakefile {script_dir}/{self.library}.smk"
                f" --configfile {self.run_config} --latency-wait {self.latency_wait} --jobs {self.jobs} --cluster-config {self.hpc_config}"
                f" --config notify={self.no_posting} verbose={self.verbose}"
                f" --drmaa ' -p {{cluster.partition}} -c {{cluster.cores}} --mem={{cluster.memory}} -J {{cluster.name}} ' -np --reason "
            )
            print(cmd)

        elif self.library == "arima":
            cmd = (
                f"snakemake --snakefile {script_dir}/arima.smk"
                f" --configfile {self.run_config} --latency-wait {self.latency_wait} --jobs {self.jobs} --nolock --cluster-config {self.hpc_config}"
                f" --config notify={self.no_posting} verbose={self.verbose}"
                f" --drmaa ' -C AVX-512 -p {{cluster.partition}} -c {{cluster.cores}} --mem={{cluster.memory}} -J {{cluster.name}} -o {self.logs}/{{rule}}.%N.%j.cluster.log' --printshellcmds --reason "
            )
        elif self.library == "omni-c" and self.curation:
            cmd = (
                f"snakemake --snakefile {script_dir}/curation_omni-c.smk"
                f" --configfile {self.run_config} --latency-wait {self.latency_wait} --jobs {self.jobs} --nolock --cluster-config {self.hpc_config}"
                f" --config notify={self.no_posting} verbose={self.verbose}"
                f" --drmaa ' -C AVX-512 -p {{cluster.partition}} -c {{cluster.cores}} --mem={{cluster.memory}} -J {{cluster.name}} -o {self.logs}/{{rule}}.%N.%j.cluster.log' --printshellcmds --reason "
            )
        elif self.library == "omni-c" and self.curation == False:
            cmd = (
                f"snakemake --snakefile {script_dir}/omni-c.smk"
                f" --configfile {self.run_config} --latency-wait {self.latency_wait} --jobs {self.jobs} --nolock --cluster-config {self.hpc_config}"
                f" --config notify={self.no_posting} verbose={self.verbose}"
                f" --drmaa ' -C AVX-512 -p {{cluster.partition}} -c {{cluster.cores}} --mem={{cluster.memory}} -J {{cluster.name}} -o {self.logs}/{{rule}}.%N.%j.cluster.log' --printshellcmds --reason "
            )


        # for universal_newlines - https://stackoverflow.com/a/4417735
        p = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True
        )
        # https://www.programcreek.com/python/example/50/subprocess.Popen
        (result, error) = p.communicate()
        exit_code = p.returncode
        print(f"\nRESULTS:\n{result}\n\nERRORS:\n{error}\n\nEXIT_CODE:\n{exit_code}\n")
        if exit_code:
            raise subprocess.CalledProcessError(exit_code, cmd)

        if self.dry_run:
            print("Dry run completed successfully!\n")

def main():
    parser = argparse.ArgumentParser(
        prog="EI HI-C", description="EI HI-C Pipeline", add_help=True,
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + __version__
    )
    subparsers = parser.add_subparsers()

    # configure
    parser_configure = subparsers.add_parser("configure", help="see `configure -h`")
    parser_configure.add_argument(
        "-s",
    "--samples_csv",
        help=f"Provide sample information in tab-separated format. Please refer to the sample file: {DEFAULT_CONFIG_FILE} for more information above the csv format. A template is provided here {DEFAULT_SAMPLE_CSV_FILE}.  (default: %(default)s)",
    )
    parser_configure.add_argument(
        "-c",
        "--curation",
        action="store_true",
        help="Use this flag if you have your final scaffold and want to run the curation pipeline. Bare in mind that the sample_csv file requires an extra field with the long reads file paths as a fifth line. (default: %(default)s)"
    )
    parser_configure.add_argument(
        "--jira",
        help="Provide JIRA id for posting job summary. E.g., PPBFX-611 (default: %(default)s)",
    )
    parser_configure.add_argument(
        "-o",
        "--output",
        default=os.path.join(cwd, "output"),
        help="Provide output directory (default: %(default)s)",
    )
    parser_configure.add_argument(
        "-f",
        "--force-reconfiguration",
        action="store_true",
        help="Force reconfiguration (default: %(default)s)",
    )
    parser_configure.add_argument(
        "-bm2",
        "--bwa_mem2",
        action="store_true",
        help="Use bwa-mem2 insted of bwa mem. This option use a lot more RAM, use with precaution. (default: %(default)s)",
    )
    parser_configure.set_defaults(handler=command_configure)

    # run
    parser_run = subparsers.add_parser("run", help="see `run -h`")
    parser_run.add_argument(
        "run_config",
        help=f"Provide run configuration YAML. Run 'eihic configure -h' to generate the run configuration YAML file. (Description template file is here: {DEFAULT_CONFIG_FILE})",
    )
    parser_run.add_argument(
        "--library",
        default="omni-c",
        help="This pipeline supports the following library protocols for hi-c data: arima (2 enzymes protocol), and  omni-c. Provide the name (arima or omni-c) as the second positional argument after the sample file (default: %(default)s)",
    )
    parser_run.add_argument(
        "-c",
        "--curation",
        action="store_true",
        help="This flag run the steps generate all the required Hi-C contact matrices and tracks for its manual curation step. (default: %(default)s)",
    )
    parser_run.add_argument(
        "--hpc_config",
        default=DEFAULT_HPC_CONFIG_FILE,
        help="Provide HPC configuration YAML (default: %(default)s)",
    )
    parser_run.add_argument(
        "--jobs",
        "-j",
        default=100,
        help="Use at most N CPU cluster/cloud jobs in parallel (default: %(default)s)",
    )
    parser_run.add_argument(
        "--latency_wait",
        default=120,
        help="Wait given seconds if an output file of a job is not present after the job finished (default: %(default)s)",
    )
    parser_run.add_argument(
        "--no_posting_off",
        action="store_false",
        help="Use this flag if you want to post comments to JIRA tickets (default: %(default)s)",
    )
    parser_run.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose mode for debugging (default: %(default)s)",
    )
    parser_run.add_argument(
        "-np", "--dry_run", action="store_true", help="Dry run (default: %(default)s)"
    )
    parser_run.set_defaults(handler=command_run)

    args = parser.parse_args()
    if hasattr(args, "handler"):
        args.handler(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()