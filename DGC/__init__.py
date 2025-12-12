"""
DGC - Domain Growth Competition

Competitive domain growth algorithm for protein domain prediction
from AlphaFold structures.
"""

from .domain_growth import (
    NonConvexDomain,
    DomainGrowthPredictor,
    competitive_domain_growth,
    simple_merge,
    boundary_strength_merge,
    size_weighted_merge,
    spatial_merge,
    get_seeds,
    find_plddt_peaks,
    AlphaFoldLoader,
    BenchmarkDomain,
    run_benchmark,
    generate_chimerax_commands,
    save_chimerax_file,
    DOMAIN_COLORS,
    calculate_radius_of_gyration,
    filter_domains_by_quality,
    build_contact_map,
    count_long_range_contacts,
    trim_terminal_extensions,
    calculate_adaptive_threshold,
)

__all__ = [
    "NonConvexDomain",
    "DomainGrowthPredictor",
    "competitive_domain_growth",
    "simple_merge",
    "boundary_strength_merge",
    "size_weighted_merge",
    "spatial_merge",
    "get_seeds",
    "find_plddt_peaks",
    "AlphaFoldLoader",
    "BenchmarkDomain",
    "run_benchmark",
    "generate_chimerax_commands",
    "save_chimerax_file",
    "DOMAIN_COLORS",
    "calculate_radius_of_gyration",
    "filter_domains_by_quality",
    "build_contact_map",
    "count_long_range_contacts",
    "trim_terminal_extensions",
    "calculate_adaptive_threshold",
]
