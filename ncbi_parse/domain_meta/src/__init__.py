"""
Domain Meta Source Modules
===========================

Tools for splitting species-grouped data by taxonomic domain.
"""

from .domain_splitter import split_by_domain, get_domain_summary, extract_domain

__all__ = [
    'split_by_domain',
    'get_domain_summary',
    'extract_domain'
]

