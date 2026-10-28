"""
Basic pipeline setup and preparation for PhaseCaller pipeline
"""

import os
import sys
import logging
import subprocess
from pathlib import Path

def create_output_dirs(config, base_output_dir):
    gene_name = config.gene_name
    gene_id = config.gene_id
    repeat_location = config.repeat_location

    output_label = gene_id if gene_id else gene_name

    base_output_path = Path(base_output_dir) / output_label

    dirs = {
        "base": base_output_path,
        "blast_temp": Path(base_output_dir) / "blast_temp",
        "blast": base_output_path / "blast",
        "sequences": base_output_path / "sequences",
        "alignments": base_output_path / "alignments",
        "reports": base_output_path / "reports",
        "plots": base_output_path / "plots",
    }

    dirs["frames"] = dirs["sequences"] / "frames"
    dirs["atg_variants"] = dirs["sequences"] / "atg_variants"
    dirs["multicopy"] = dirs["sequences"] / "multicopy"
    dirs["translations"] = dirs["sequences"] / "translations"
    dirs["pcr_fragments"] = dirs["sequences"] / "pcr_fragments"
    dirs["alleles"] = dirs["sequences"] / "alleles"

    if repeat_location.lower() == "promoter": 
        dirs["promoters"] = dirs["sequences"] / "promoters"

    for directory in dirs.values():
        directory.mkdir(parents=True, exist_ok=True)

    return dirs

def setup_logging(output_dir, gene_name):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / f"{gene_name}_pipeline.log"

    logging.basicConfig(
        level = logging.INFO,
        format = "%(asctime)s [%(levelname)s] %(message)s",
        handlers = [logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)]
    )

    return logging.getLogger("pipeline_logger")

logger = logging.getLogger(__name__)


