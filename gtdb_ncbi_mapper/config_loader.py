#!/usr/bin/env python3
"""
Configuration Loader for GTDB-NCBI Mapper
==========================================

Loads and validates configuration from config.yaml
"""

import yaml
from pathlib import Path
from typing import Dict, Any


class Config:
    """Configuration manager for GTDB-NCBI Mapper."""
    
    def __init__(self, config_path: Path = None):
        """
        Load configuration from YAML file.
        
        Args:
            config_path: Path to config.yaml. If None, uses default location.
        """
        if config_path is None:
            config_path = Path(__file__).parent / 'config.yaml'
        
        self.config_path = config_path
        self.base_dir = config_path.parent
        
        # Load config
        with open(config_path, 'r') as f:
            self._config = yaml.safe_load(f)
    
    def get(self, *keys, default=None):
        """
        Get configuration value using dot notation.
        
        Example:
            config.get('input', 'gtdb', 'archaea')
        """
        value = self._config
        for key in keys:
            if isinstance(value, dict) and key in value:
                value = value[key]
            else:
                return default
        return value
    
    def get_path(self, *keys, create_dir=False) -> Path:
        """
        Get path from config, resolved relative to config file location.
        
        Args:
            *keys: Configuration keys to traverse
            create_dir: If True, create directory if it doesn't exist
        
        Returns:
            Absolute Path object
        """
        path_str = self.get(*keys)
        if path_str is None:
            raise ValueError(f"Configuration path not found: {'.'.join(keys)}")
        
        # Resolve relative to config file location
        path = self.base_dir / path_str
        
        # Create directory if requested
        if create_dir and not path.exists():
            path.mkdir(parents=True, exist_ok=True)
        
        return path
    
    def get_input_paths(self) -> Dict[str, Path]:
        """Get all input file paths."""
        return {
            'archaea': self.get_path('input', 'gtdb', 'archaea'),
            'bacteria': self.get_path('input', 'gtdb', 'bacteria'),
            'ncbi_assembly': self.get_path('input', 'ncbi', 'assembly_summary')
        }
    
    def get_output_dirs(self, create=True) -> Dict[str, Path]:
        """Get all output directory paths."""
        dirs = {
            'raw': self.get_path('output', 'raw', create_dir=create),
            'mapping_tables': self.get_path('output', 'mapping_tables', create_dir=create),
            'logs': self.get_path('output', 'logs', create_dir=create),
            'archive': self.get_path('output', 'archive', create_dir=create)
        }
        return dirs
    
    def get_output_file(self, category: str, name: str) -> Path:
        """
        Get output file path.
        
        Args:
            category: 'raw', 'mapping', or 'status'
            name: File name key from config
        
        Returns:
            Full path to output file
        """
        filename = self.get('output_files', category, name)
        if filename is None:
            raise ValueError(f"Output file not found in config: {category}.{name}")
        
        # Determine output directory based on category
        if category == 'raw':
            output_dir = self.get_path('output', 'raw', create_dir=True)
        elif category == 'mapping':
            output_dir = self.get_path('output', 'mapping_tables', create_dir=True)
        elif category == 'status':
            output_dir = self.get_path('output', 'archive', create_dir=True)
        else:
            raise ValueError(f"Unknown output category: {category}")
        
        return output_dir / filename
    
    def get_log_file(self, log_name: str) -> Path:
        """Get log file path."""
        log_path = self.get('logging', 'files', log_name)
        if log_path is None:
            raise ValueError(f"Log file not found in config: {log_name}")
        
        # Resolve relative to base directory
        full_path = self.base_dir / log_path
        
        # Create parent directory if needed
        full_path.parent.mkdir(parents=True, exist_ok=True)
        
        return full_path
    
    @property
    def metadata(self) -> Dict[str, Any]:
        """Get pipeline metadata."""
        return self.get('metadata', default={})
    
    @property
    def ncbi_datasets_config(self) -> Dict[str, Any]:
        """Get NCBI datasets configuration."""
        return self.get('ncbi_datasets', default={})
    
    @property
    def processing_config(self) -> Dict[str, Any]:
        """Get processing configuration."""
        return self.get('processing', default={})
    
    def __repr__(self):
        return f"Config(config_path={self.config_path})"


# Convenience function
def load_config(config_path: Path = None) -> Config:
    """Load configuration from YAML file."""
    return Config(config_path)

