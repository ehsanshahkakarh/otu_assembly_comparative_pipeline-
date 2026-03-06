"""
Configuration and Path Setup for 16S Census Parser
===================================================

Handles directory paths and logging configuration for the census parser.
"""

import logging
import sys
from pathlib import Path
from dataclasses import dataclass, field


@dataclass
class PathConfig:
    """Configuration for directory paths."""

    # Base directory (16S_censusparse/)
    base_dir: Path = field(default_factory=lambda: Path(__file__).resolve().parents[1])

    @property
    def metadata_dir(self) -> Path:
        """Directory containing input metadata files."""
        return self.base_dir / "metadata"

    @property
    def csv_output_dir(self) -> Path:
        """Directory for CSV output files."""
        return self.base_dir / "output"

    @property
    def log_dir(self) -> Path:
        """Directory for log files."""
        return self.base_dir / "logs"

    @property
    def cache_dir(self) -> Path:
        """Directory for cache files."""
        return self.base_dir / "cache"

    @property
    def input_file(self) -> Path:
        """Path to the input TSV file."""
        return self.metadata_dir / "eukcensus_16S.clusters.97.tsv"

    def ensure_directories(self):
        """Create necessary directories if they don't exist."""
        self.csv_output_dir.mkdir(exist_ok=True)
        self.log_dir.mkdir(exist_ok=True)
        self.cache_dir.mkdir(exist_ok=True)


def setup_directory_paths():
    """
    Set up directory paths for the reorganized 16S_censusparse structure.
    
    Returns:
        Tuple of (script_dir, metadata_dir, csv_output_dir, log_dir)
    """
    config = PathConfig()
    config.ensure_directories()
    
    script_dir = Path(__file__).resolve().parent.parent
    return script_dir, config.metadata_dir, config.csv_output_dir, config.log_dir


def setup_logging(log_dir, output_prefix="eukcensus16S"):
    """
    Set up logging configuration for the enhanced script.
    
    Args:
        log_dir: Directory to store log files
        output_prefix: Prefix for output files to include in log messages
        
    Returns:
        Logger instance
    """
    # Create log file path in the logs directory
    log_file = log_dir / "eukcensus16S_processing.log"
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='a'),  # Append to existing log
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    # Log the start of processing
    logging.info(f"Starting enhanced EukCensus 16S processing: eukcensus_16S.clusters.97.tsv -> {output_prefix}_*")
    logging.info(f"Log file: {log_file}")
    return logging.getLogger(__name__)

