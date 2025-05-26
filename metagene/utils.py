#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Common utilities for Metagene

import os
import logging
from pathlib import Path
from typing import Optional, Union, List

def setup_logger(name: str = "metagene", level: int = logging.INFO) -> logging.Logger:
    """
    Set up and configure a logger for the application.
    
    Args:
        name: Name of the logger (default: "metagene")
        level: Logging level (default: INFO)
        
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    
    # Only add handlers if they haven't been added already
    if not logger.handlers:
        logger.setLevel(level)
        
        # Create console handler with formatting
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    logger.propagate = False  # Prevent double logging
    return logger

# Initialize logger
logger = setup_logger()

def get_cache_dir(app_name: str = "metagene") -> Path:
    """
    Get the cache directory path following XDG Base Directory Specification.
    
    Args:
        app_name: Name of the application (default: "metagene")
        
    Returns:
        Path object pointing to the cache directory
    """
    if "XDG_CACHE_HOME" in os.environ:
        cache_home = Path(os.environ["XDG_CACHE_HOME"])
    else:
        cache_home = Path.home() / ".cache"
    
    cache_dir = cache_home / app_name
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    return cache_dir


def ensure_dir(path: Union[str, Path]) -> Path:
    """
    Ensure a directory exists, creating it if necessary.
    
    Args:
        path: Directory path to ensure exists
        
    Returns:
        Path object for the directory
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_file_hash(file_path: Union[str, Path]) -> str:
    """
    Calculate MD5 hash of a file.
    
    Args:
        file_path: Path to the file
        
    Returns:
        MD5 hash as a hexadecimal string
    """
    import hashlib
    
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    
    return hash_md5.hexdigest()


def get_file_size(file_path: Union[str, Path]) -> int:
    """
    Get the size of a file in bytes.
    
    Args:
        file_path: Path to the file
        
    Returns:
        File size in bytes
    """
    return Path(file_path).stat().st_size


def format_file_size(size_bytes: int) -> str:
    """
    Format file size in bytes to human-readable format.
    
    Args:
        size_bytes: Size in bytes
        
    Returns:
        Human-readable file size string
    """
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB" 


def print_analysis_summary(params: dict) -> None:
    """Print a summary of analysis parameters.
    
    Args:
        params: Dictionary containing analysis parameters
    """
    logger.info("\nAnalysis Summary:")
    logger.info("------------------")
    logger.info(f"Input file: {params['input_file']}")
    logger.info(f"Output file: {params['output_file']}")
    
    if params.get('output_score'):
        logger.info(f"Score output file: {params['output_score']}")
    if params.get('output_figure'):
        logger.info(f"Figure output file: {params['output_figure']}")
    if params.get('reference'):
        logger.info(f"Reference: {params['reference']}")
    if params.get('gtf'):
        logger.info(f"GTF file: {params['gtf']}")
    
    logger.info(f"Region: {params.get('region', 'N/A')}")
    logger.info(f"Number of bins: {params.get('bin_number', params.get('bins', 'N/A'))}")
    logger.info(f"Has header: {params.get('with_header', 'N/A')}")
    logger.info(f"Meta columns: {params.get('meta_columns', 'N/A')}")
    logger.info(f"Weight columns: {params.get('weight_columns', 'N/A')}")
    
    if params.get('weight_names'):
        logger.info(f"Weight names: {params['weight_names']}")
    
    logger.info(f"Score transform: {params.get('score_transform', 'N/A')}")
    logger.info(f"Normalize: {params.get('normalize', 'N/A')}")
    
    if params.get('plot'):
        logger.info("\nPlot Settings:")
        if params.get('plot_dir'):
            logger.info(f"Plot directory: {params['plot_dir']}")
        logger.info(f"Plot format: {params.get('plot_format', 'N/A')}")
    
    logger.info("------------------")