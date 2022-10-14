import pkg_resources

__title__ = "eihic"
__author__ = "Mariano Olivera Fedi"
__email__ = "mariano.Olivera-Fedi@earlham.ac.uk"
__copyright__ = "Copyright 2022 Earlham Institute"
__version__ = pkg_resources.require("eihic")[0].version

DEFAULT_CONFIG_FILE = pkg_resources.resource_filename("eihic.etc", "run_config.yaml")
DEFAULT_HPC_CONFIG_FILE = pkg_resources.resource_filename(
    "eihic.etc", "hpc_config.json"
)
DEFAULT_SAMPLE_CSV_FILE = pkg_resources.resource_filename("eihic.etc", "samples.csv")