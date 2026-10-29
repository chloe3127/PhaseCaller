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

from pipeline import pipeline_prep, report_maker, plots
from pipeline.blast import blast_parser, blast_runner, blast_prep


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
    pipeline_prep.setup_logging(base_output_dir, "PhaseCaller")
    logger = logging.getLogger("PhaseCaller")

    logger.info("Starting PhaseCaller pipeline for %s", ", ".join(args.gene))

    output_dirs_by_config = []
    for config in configs:
        t_base_output_dir = str(base_output_dir)
        output_dirs_by_config.append((config, pipeline_prep.create_output_dirs(config, t_base_output_dir)))

    blast_temp = base_output_dir / "blast_temp" 
    combined_fasta, contig_id_map, isolate_ids = blast_prep.create_combined_genome_fasta(args.input, blast_temp)

    db_prefix = blast_prep.make_blast_db(combined_fasta)
    contig_counts = Counter(internal_id.rsplit("_ctg", 1)[0] for internal_id in contig_id_map)

    pipeline_success = False
    try:
        for config, output_dirs in output_dirs_by_config:
            logger.info("Running pipeline for gene: %s", config.gene_name)

            gene_ref_fasta = blast_prep.create_gene_ref_fasta([config], output_dirs["blast"])
            gene_start_ref_fasta = blast_prep.create_gene_start_ref_fasta([config], output_dirs["blast"])
            gene_end_ref_fasta = blast_prep.create_gene_end_ref_fasta(config, output_dirs["blast"])
            primer_ref_fasta = blast_prep.create_primer_ref_fasta([config], output_dirs["blast"])

            blast_outputs = blast_runner.run_all_blast_searches(config, gene_ref_fasta, gene_start_ref_fasta, gene_end_ref_fasta, primer_ref_fasta, db_prefix, output_dirs["blast"])

            blast_results = blast_parser.parse_blast_outputs(blast_outputs, config, contig_id_map)

            blast_results["contig_counts"] = contig_counts

            report_maker.generate_all_reports(config, blast_results, output_dirs["reports"], isolate_ids)
            plots.plot_raw_blast_results(blast_outputs, config, output_dirs["plots"])
            plots.plot_raw_blast_results(blast_outputs, config, output_dirs["plots"])
            plots.plot_primer_mismatch(blast_outputs, config, output_dirs["plots"])





            if len(configs) == 1 and not args.keep_blast:
                shutil.rmtree(blast_temp)
                logger.info("Removed temporary BLAST directory: %s", blast_temp)
            pipeline_success = True
    finally:
        if len(configs) > 1 and pipeline_success and not args.keep_blast:
            shutil.rmtree(blast_temp)
            logger.info("Removed temporary BLAST directory: %s", blast_temp)


if __name__ == "__main__":
    main()



