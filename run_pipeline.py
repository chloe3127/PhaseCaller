"""
Main script that runs everything, needs a config file to run.

Use with:

python run_pipeline.py --gene <gene_name> --input <input_file> --output <output_folder>
"""

import logging
import shutil
import importlib.util
from collections import Counter
from pathlib import Path
import argparse



def load_config(config_path):
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    spec = importlib.util.spec_from_file_location(config_path.stem, config_path)

    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load config file from {config_path}")

    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)

    return config

def get_configs_to_run(gene_args, config_dir = "configs"):
    #lets you run multiple at once by creating a list of all config files matching gene arguments - saves time by using same blast db

    config_dir = Path(config_dir)

    if "all" in gene_args:
        config_paths = sorted(p for p in config_dir.glob("config_*.py") if p.name != "config_template.py" and p.name != "config_example.py")

    else: 
        config_paths = []
        for gene in gene_args:
            config_path = config_dir / f"config_{gene}.py"

            if not config_path.exists():
                raise FileNotFoundError(f"Config file not found: {config_path}")
            config_paths.append(config_path)

    return config_paths

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--gene", nargs = "+", required = True, help = "Gene name(s) to run the pipeline for. Use 'all' to run for all available genes"
    )

    parser.add_argument(
        "--input", required = True, help = "path to input file/folder/zip with genome sequences"
    )

    parser.add_argument(
        "--output", required = True, help = "path to output folder for results"
    )

    parser.add_argument(
        "--keep_blast", action = "store_true", help = "Keep intermediate BLAST files after running"
    )

    return parser.parse_args()

def main():
    args = parse_args()

    config_paths = get_configs_to_run(args.gene)
    configs = [load_config(p) for p in config_paths]
    base_output_dir = Path(args.output)
    base_output_dir.mkdir(parents=True, exist_ok=True)
    