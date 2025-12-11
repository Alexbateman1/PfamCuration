"""
Interactive Graph Visualization for ABC Domain Predictor

Creates an HTML visualization of the contact graph used by Leiden clustering.
Supports multiple coloring modes to compare ground truth, predictions, and clusters.
"""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class DomainAnnotation:
    """A single domain annotation from ground truth or prediction."""
    domain_id: int
    segments: List[Tuple[int, int]]  # List of (start, end) tuples

    def contains_residue(self, resnum: int) -> bool:
        """Check if this domain contains the given residue number."""
        for start, end in self.segments:
            if start <= resnum <= end:
                return True
        return False


def parse_dev_file(dev_file_path: str) -> Dict[str, List[DomainAnnotation]]:
    """
    Parse a dev file with gold standard domain definitions.

    Format:
        ACCESSION seg1-seg2_seg3-seg4 seg5-seg6 ...

    Example:
        A9WIN4 19-179_273-325 217-253 329-420

    Returns dict mapping accession -> list of DomainAnnotation
    """
    annotations = {}

    with open(dev_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < 1:
                continue

            accession = parts[0]
            domains = []

            for i, domain_str in enumerate(parts[1:], start=1):
                # Parse segments (e.g., "19-179_273-325" -> [(19, 179), (273, 325)])
                segments = []
                for seg_str in domain_str.split('_'):
                    if '-' in seg_str:
                        start, end = seg_str.split('-')
                        segments.append((int(start), int(end)))

                if segments:
                    domains.append(DomainAnnotation(domain_id=i, segments=segments))

            annotations[accession] = domains

    return annotations


def get_residue_domain_map(
    domains: List[DomainAnnotation],
    max_resnum: int
) -> Dict[int, int]:
    """
    Create a mapping from residue number to domain ID.
    Returns dict: resnum -> domain_id (0 = no domain)
    """
    mapping = {}
    for resnum in range(1, max_resnum + 1):
        domain_id = 0
        for domain in domains:
            if domain.contains_residue(resnum):
                domain_id = domain.domain_id
                break
        mapping[resnum] = domain_id
    return mapping


def generate_color_palette(n_colors: int) -> List[str]:
    """Generate n distinct colors for visualization."""
    # Use a colorblind-friendly palette for small n, otherwise generate
    if n_colors <= 10:
        colors = [
            '#1f77b4',  # blue
            '#ff7f0e',  # orange
            '#2ca02c',  # green
            '#d62728',  # red
            '#9467bd',  # purple
            '#8c564b',  # brown
            '#e377c2',  # pink
            '#7f7f7f',  # gray
            '#bcbd22',  # olive
            '#17becf',  # cyan
        ]
        return colors[:n_colors]
    else:
        # Generate colors using HSL
        colors = []
        for i in range(n_colors):
            hue = i / n_colors
            # Convert HSL to hex (simplified)
            r, g, b = _hsl_to_rgb(hue, 0.7, 0.5)
            colors.append(f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}')
        return colors


def _hsl_to_rgb(h: float, s: float, l: float) -> Tuple[float, float, float]:
    """Convert HSL to RGB."""
    if s == 0:
        return l, l, l

    def hue_to_rgb(p, q, t):
        if t < 0: t += 1
        if t > 1: t -= 1
        if t < 1/6: return p + (q - p) * 6 * t
        if t < 1/2: return q
        if t < 2/3: return p + (q - p) * (2/3 - t) * 6
        return p

    q = l * (1 + s) if l < 0.5 else l + s - l * s
    p = 2 * l - q
    r = hue_to_rgb(p, q, h + 1/3)
    g = hue_to_rgb(p, q, h)
    b = hue_to_rgb(p, q, h - 1/3)
    return r, g, b


def create_graph_visualization(
    graph: nx.Graph,
    residues: List,  # List of Residue objects
    cluster_assignments: Dict[int, int],
    ground_truth_domains: Optional[List[DomainAnnotation]] = None,
    predicted_domains: Optional[List] = None,  # List of Domain objects
    output_path: str = "graph_visualization.html",
    title: str = "Contact Graph Visualization",
    uniprot_acc: str = None,
    show_chain_trace: bool = False,
) -> str:
    """
    Create an interactive HTML visualization of the contact graph.

    Parameters:
    -----------
    graph : nx.Graph
        The contact graph
    residues : List[Residue]
        List of residue objects with resnum, plddt, etc.
    cluster_assignments : Dict[int, int]
        Mapping of residue index to cluster ID
    ground_truth_domains : List[DomainAnnotation], optional
        Gold standard domain definitions
    predicted_domains : List[Domain], optional
        ABC predicted domains
    output_path : str
        Path for output HTML file
    title : str
        Title for the visualization
    show_chain_trace : bool
        If True, add backbone edges connecting sequential residues

    Returns:
    --------
    str : Path to the generated HTML file
    """

    # Build residue info
    n_residues = len(residues)
    resnum_to_idx = {res.resnum: i for i, res in enumerate(residues)}

    # Get unique clusters and create color map
    unique_clusters = sorted(set(cluster_assignments.values()))
    n_clusters = len(unique_clusters)
    cluster_colors = generate_color_palette(max(n_clusters, 10))
    cluster_color_map = {c: cluster_colors[i % len(cluster_colors)] for i, c in enumerate(unique_clusters)}

    # Get domain color maps
    if ground_truth_domains:
        n_gt_domains = len(ground_truth_domains)
        gt_colors = generate_color_palette(max(n_gt_domains + 1, 10))
        gt_domain_map = get_residue_domain_map(ground_truth_domains, max(r.resnum for r in residues))
    else:
        gt_domain_map = {}
        gt_colors = ['#cccccc']

    if predicted_domains:
        n_pred_domains = len(predicted_domains)
        pred_colors = generate_color_palette(max(n_pred_domains + 1, 10))
        pred_domain_map = {}
        for domain in predicted_domains:
            for seg_start, seg_end in domain.segments:
                for resnum in range(seg_start, seg_end + 1):
                    pred_domain_map[resnum] = domain.domain_id
    else:
        pred_domain_map = {}
        pred_colors = ['#cccccc']

    # Build nodes data
    nodes_data = []
    for i, res in enumerate(residues):
        cluster_id = cluster_assignments.get(i, -1)
        gt_domain_id = gt_domain_map.get(res.resnum, 0)
        pred_domain_id = pred_domain_map.get(res.resnum, 0)

        # Track if residue is in GT but missed by prediction (false negative)
        is_missed = gt_domain_id > 0 and pred_domain_id == 0

        node = {
            'id': i,
            'resnum': res.resnum,
            'plddt': round(res.plddt, 1),
            'cluster': cluster_id,
            'gt_domain': gt_domain_id,
            'pred_domain': pred_domain_id,
            'is_missed': is_missed,  # In GT but not in prediction
            'cluster_color': cluster_color_map.get(cluster_id, '#cccccc'),
            'gt_color': gt_colors[gt_domain_id % len(gt_colors)] if gt_domain_id > 0 else '#cccccc',
            'pred_color': pred_colors[pred_domain_id % len(pred_colors)] if pred_domain_id > 0 else '#cccccc',
        }
        nodes_data.append(node)

    # Build edges data
    edges_data = []
    for u, v, data in graph.edges(data=True):
        edge = {
            'from': u,
            'to': v,
            'weight': round(data.get('weight', 1.0), 3),
            'distance': round(data.get('distance', 0), 1),
            'hbond': data.get('hbond', False),
            'backbone': False,
        }
        edges_data.append(edge)

    # Add backbone chain trace edges if requested
    if show_chain_trace:
        # Create set of existing edges to avoid duplicates
        existing_edges = {(e['from'], e['to']) for e in edges_data}
        existing_edges.update({(e['to'], e['from']) for e in edges_data})

        # Add backbone edges (connecting sequential residues by index)
        for i in range(len(residues) - 1):
            # Check if residues are actually sequential (no gaps)
            if residues[i + 1].resnum - residues[i].resnum == 1:
                if (i, i + 1) not in existing_edges:
                    edges_data.append({
                        'from': i,
                        'to': i + 1,
                        'weight': 0,
                        'distance': 3.8,  # Typical CA-CA distance
                        'hbond': False,
                        'backbone': True,
                    })

    # Generate HTML
    html = _generate_html(
        nodes_data=nodes_data,
        edges_data=edges_data,
        title=title,
        n_clusters=n_clusters,
        n_gt_domains=len(ground_truth_domains) if ground_truth_domains else 0,
        n_pred_domains=len(predicted_domains) if predicted_domains else 0,
        uniprot_acc=uniprot_acc,
    )

    # Write to file
    output_path = Path(output_path)
    output_path.write_text(html)
    logger.info(f"Graph visualization saved to {output_path}")

    return str(output_path)


def _generate_html(
    nodes_data: List[Dict],
    edges_data: List[Dict],
    title: str,
    n_clusters: int,
    n_gt_domains: int,
    n_pred_domains: int,
    uniprot_acc: str = None,
) -> str:
    """Generate the HTML content with vis.js visualization."""

    nodes_json = json.dumps(nodes_data)
    edges_json = json.dumps(edges_data)

    # Create clickable title with TED link if UniProt accession provided
    if uniprot_acc:
        title_html = f'Contact Graph: <a href="https://ted.cathdb.info/uniprot/{uniprot_acc}" target="_blank" style="color: #0066cc;">{uniprot_acc}</a>'
    else:
        title_html = title

    html = f'''<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        h1 {{
            margin: 0 0 10px 0;
            color: #333;
        }}
        h1 a {{
            text-decoration: none;
        }}
        h1 a:hover {{
            text-decoration: underline;
        }}
        .container {{
            display: flex;
            gap: 20px;
        }}
        #graph {{
            width: 80%;
            height: 800px;
            border: 1px solid #ddd;
            background: white;
            border-radius: 8px;
        }}
        .controls {{
            width: 20%;
            background: white;
            padding: 15px;
            border-radius: 8px;
            border: 1px solid #ddd;
        }}
        .control-group {{
            margin-bottom: 15px;
        }}
        .control-group label {{
            display: block;
            font-weight: bold;
            margin-bottom: 5px;
            color: #555;
        }}
        select, input {{
            width: 100%;
            padding: 8px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }}
        .stats {{
            background: #f9f9f9;
            padding: 10px;
            border-radius: 4px;
            font-size: 13px;
            margin-top: 15px;
        }}
        .stats h3 {{
            margin: 0 0 10px 0;
            font-size: 14px;
        }}
        .stats p {{
            margin: 5px 0;
        }}
        .legend {{
            margin-top: 15px;
            font-size: 12px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            margin: 3px 0;
        }}
        .legend-color {{
            width: 16px;
            height: 16px;
            border-radius: 50%;
            margin-right: 8px;
            border: 1px solid #999;
        }}
        #tooltip {{
            position: absolute;
            background: rgba(0,0,0,0.8);
            color: white;
            padding: 10px;
            border-radius: 4px;
            font-size: 12px;
            pointer-events: none;
            display: none;
            z-index: 1000;
            max-width: 250px;
        }}
        .highlight-controls {{
            margin-top: 15px;
        }}
        .highlight-controls input {{
            width: auto;
            margin-right: 10px;
        }}
    </style>
</head>
<body>
    <h1>{title_html}</h1>
    <div class="container">
        <div id="graph"></div>
        <div class="controls">
            <div class="control-group">
                <label>Color By:</label>
                <select id="colorMode">
                    <option value="cluster">Initial Clusters ({n_clusters})</option>
                    <option value="gt_domain">Ground Truth Domains ({n_gt_domains})</option>
                    <option value="pred_domain">Predicted Domains ({n_pred_domains})</option>
                    <option value="plddt">pLDDT Score</option>
                </select>
            </div>

            <div class="control-group">
                <label>Node Size:</label>
                <select id="sizeMode">
                    <option value="uniform">Uniform</option>
                    <option value="plddt">By pLDDT</option>
                    <option value="degree">By Degree</option>
                </select>
            </div>

            <div class="control-group">
                <label>Edge Filter (min weight):</label>
                <input type="range" id="edgeFilter" min="0" max="1" step="0.05" value="0">
                <span id="edgeFilterValue">0</span>
            </div>

            <div class="control-group">
                <label>Show H-bond edges only:</label>
                <input type="checkbox" id="hbondOnly" style="width: auto;">
            </div>

            <div class="highlight-controls">
                <label>Highlight Residue Range:</label>
                <input type="number" id="highlightStart" placeholder="Start" style="width: 60px;">
                <input type="number" id="highlightEnd" placeholder="End" style="width: 60px;">
                <button onclick="highlightRange()">Highlight</button>
                <button onclick="clearHighlight()">Clear</button>
            </div>

            <div class="stats">
                <h3>Statistics</h3>
                <p>Nodes: {len(nodes_data)}</p>
                <p>Edges: {len(edges_data)}</p>
                <p>Clusters: {n_clusters}</p>
                <p>GT Domains: {n_gt_domains}</p>
                <p>Predicted: {n_pred_domains}</p>
            </div>

            <div id="legend" class="legend"></div>
        </div>
    </div>
    <div id="tooltip"></div>

    <script>
        // Data
        const nodesData = {nodes_json};
        const edgesData = {edges_json};

        // Build vis.js dataset
        let nodes = new vis.DataSet();
        let edges = new vis.DataSet();

        // Initialize nodes
        nodesData.forEach(n => {{
            nodes.add({{
                id: n.id,
                label: String(n.resnum),
                title: `Residue ${{n.resnum}}\\npLDDT: ${{n.plddt}}\\nCluster: ${{n.cluster}}\\nGT Domain: ${{n.gt_domain}}\\nPred Domain: ${{n.pred_domain}}${{n.is_missed ? '\\n⚠️ MISSED (in GT, not predicted)' : ''}}`,
                color: n.cluster_color,
                shape: 'dot',
                resnum: n.resnum,
                plddt: n.plddt,
                cluster: n.cluster,
                gt_domain: n.gt_domain,
                pred_domain: n.pred_domain,
                is_missed: n.is_missed,
                cluster_color: n.cluster_color,
                gt_color: n.gt_color,
                pred_color: n.pred_color,
            }});
        }});

        // Initialize edges
        edgesData.forEach((e, i) => {{
            let edgeWidth, edgeColor, edgeTitle;
            if (e.backbone) {{
                // Backbone chain trace - thicker, dark gray
                edgeWidth = 3;
                edgeColor = '#333333';
                edgeTitle = 'Backbone chain trace';
            }} else {{
                // Contact edge
                edgeWidth = Math.max(0.5, e.weight * 2);
                edgeColor = e.hbond ? '#e74c3c' : '#999999';
                edgeTitle = `Weight: ${{e.weight}}\\nDistance: ${{e.distance}}Å${{e.hbond ? '\\nH-bonded' : ''}}`;
            }}
            edges.add({{
                id: i,
                from: e.from,
                to: e.to,
                weight: e.weight,
                distance: e.distance,
                hbond: e.hbond,
                backbone: e.backbone,
                title: edgeTitle,
                width: edgeWidth,
                color: edgeColor,
            }});
        }});

        // Create network
        const container = document.getElementById('graph');
        const data = {{ nodes: nodes, edges: edges }};
        const options = {{
            nodes: {{
                shape: 'dot',
                size: 8,
                font: {{ size: 10, color: '#333' }},
                borderWidth: 1,
                borderWidthSelected: 3,
            }},
            edges: {{
                smooth: {{ type: 'continuous' }},
                color: {{ inherit: false }},
            }},
            physics: {{
                enabled: true,
                solver: 'forceAtlas2Based',
                forceAtlas2Based: {{
                    gravitationalConstant: -50,
                    centralGravity: 0.01,
                    springLength: 100,
                    springConstant: 0.08,
                }},
                stabilization: {{
                    iterations: 200,
                    updateInterval: 25,
                }},
            }},
            interaction: {{
                hover: true,
                tooltipDelay: 100,
            }},
        }};

        const network = new vis.Network(container, data, options);

        // Color mode change
        document.getElementById('colorMode').addEventListener('change', function() {{
            const mode = this.value;
            nodes.forEach(node => {{
                let color;
                let shape = 'dot';  // Default shape

                if (mode === 'cluster') {{
                    color = node.cluster_color;
                }} else if (mode === 'gt_domain') {{
                    color = node.gt_color;
                }} else if (mode === 'pred_domain') {{
                    color = node.pred_color;
                    // Show missed residues (in GT but not predicted) as grey squares
                    if (node.is_missed) {{
                        shape = 'square';
                        color = '#cccccc';  // Grey - not in any predicted domain
                    }}
                }} else if (mode === 'plddt') {{
                    // AlphaFold standard pLDDT coloring
                    const plddt = node.plddt;
                    if (plddt >= 90) color = '#0053D6';      // Very high - dark blue
                    else if (plddt >= 70) color = '#65CBF3'; // Confident - light blue
                    else if (plddt >= 50) color = '#FFDB13'; // Low - yellow
                    else color = '#FF7D45';                   // Very low - orange
                }}
                nodes.update({{ id: node.id, color: color, shape: shape }});
            }});
            updateLegend(mode);
        }});

        // Size mode change
        document.getElementById('sizeMode').addEventListener('change', function() {{
            const mode = this.value;
            nodes.forEach(node => {{
                let size;
                if (mode === 'uniform') {{
                    size = 8;
                }} else if (mode === 'plddt') {{
                    size = 4 + (node.plddt / 100) * 12;
                }} else if (mode === 'degree') {{
                    const degree = network.getConnectedEdges(node.id).length;
                    size = 4 + Math.min(degree, 20);
                }}
                nodes.update({{ id: node.id, size: size }});
            }});
        }});

        // Edge filter
        document.getElementById('edgeFilter').addEventListener('input', function() {{
            const minWeight = parseFloat(this.value);
            document.getElementById('edgeFilterValue').textContent = minWeight.toFixed(2);
            const hbondOnly = document.getElementById('hbondOnly').checked;

            edges.forEach(edge => {{
                // Backbone edges are always shown (not affected by filters)
                if (edge.backbone) return;
                const hidden = edge.weight < minWeight || (hbondOnly && !edge.hbond);
                edges.update({{ id: edge.id, hidden: hidden }});
            }});
        }});

        // H-bond only filter
        document.getElementById('hbondOnly').addEventListener('change', function() {{
            const hbondOnly = this.checked;
            const minWeight = parseFloat(document.getElementById('edgeFilter').value);

            edges.forEach(edge => {{
                // Backbone edges are always shown (not affected by filters)
                if (edge.backbone) return;
                const hidden = edge.weight < minWeight || (hbondOnly && !edge.hbond);
                edges.update({{ id: edge.id, hidden: hidden }});
            }});
        }});

        // Highlight range
        function highlightRange() {{
            const start = parseInt(document.getElementById('highlightStart').value);
            const end = parseInt(document.getElementById('highlightEnd').value);
            if (isNaN(start) || isNaN(end)) return;

            nodes.forEach(node => {{
                const inRange = node.resnum >= start && node.resnum <= end;
                nodes.update({{
                    id: node.id,
                    borderWidth: inRange ? 3 : 1,
                    borderWidthSelected: inRange ? 5 : 3,
                    font: {{ size: inRange ? 14 : 10 }},
                }});
            }});
        }}

        function clearHighlight() {{
            nodes.forEach(node => {{
                nodes.update({{
                    id: node.id,
                    borderWidth: 1,
                    borderWidthSelected: 3,
                    font: {{ size: 10 }},
                }});
            }});
        }}

        // Legend
        function updateLegend(mode) {{
            const legend = document.getElementById('legend');
            let html = '<h3>Legend</h3>';

            if (mode === 'plddt') {{
                html += `
                    <div class="legend-item"><div class="legend-color" style="background:#0053D6"></div>pLDDT ≥ 90 (very high)</div>
                    <div class="legend-item"><div class="legend-color" style="background:#65CBF3"></div>pLDDT 70-90 (confident)</div>
                    <div class="legend-item"><div class="legend-color" style="background:#FFDB13"></div>pLDDT 50-70 (low)</div>
                    <div class="legend-item"><div class="legend-color" style="background:#FF7D45"></div>pLDDT < 50 (very low)</div>
                `;
            }} else if (mode === 'pred_domain') {{
                html += '<div class="legend-item"><div class="legend-color" style="background:#cccccc"></div>No predicted domain</div>';
                html += '<div class="legend-item"><div class="legend-color" style="background:#cccccc; border-radius:0"></div>■ Missed (in GT, not predicted)</div>';
            }} else {{
                html += '<div class="legend-item"><div class="legend-color" style="background:#cccccc"></div>No domain / Cluster</div>';
            }}

            html += '<br><div class="legend-item"><div class="legend-color" style="background:#333333"></div>Backbone trace</div>';
            html += '<div class="legend-item"><div class="legend-color" style="background:#e74c3c"></div>H-bonded edge</div>';
            html += '<div class="legend-item"><div class="legend-color" style="background:#999999"></div>Regular edge</div>';

            legend.innerHTML = html;
        }}

        // Initial legend
        updateLegend('cluster');

        // Stop physics after stabilization
        network.on('stabilizationIterationsDone', function() {{
            network.setOptions({{ physics: {{ enabled: false }} }});
        }});
    </script>
</body>
</html>'''

    return html


def generate_chimerax_cluster_commands(
    residues: List,
    cluster_assignments: Dict[int, int],
    output_path: str,
    uniprot_acc: str,
    max_clusters: int = 20,
) -> str:
    """
    Generate ChimeraX commands to color structure by Leiden clusters.

    Parameters:
    -----------
    residues : List[Residue]
        List of residue objects
    cluster_assignments : Dict[int, int]
        Mapping of residue index to cluster ID
    output_path : str
        Path for output .cxc file
    uniprot_acc : str
        UniProt accession for fetching structure
    max_clusters : int
        Maximum number of clusters to show with distinct colors (default: 20)

    Returns:
    --------
    str : Path to the generated .cxc file
    """
    from collections import Counter

    # Count cluster sizes
    cluster_counts = Counter(cluster_assignments.values())

    # Sort by size (largest first)
    sorted_clusters = sorted(cluster_counts.items(), key=lambda x: -x[1])

    # Get top N clusters
    top_clusters = sorted_clusters[:max_clusters]
    top_cluster_ids = {c[0] for c in top_clusters}

    # Generate distinct colors for top clusters
    # Using a perceptually distinct color palette
    distinct_colors = [
        (0.12, 0.47, 0.71),   # blue
        (1.00, 0.50, 0.05),   # orange
        (0.17, 0.63, 0.17),   # green
        (0.84, 0.15, 0.16),   # red
        (0.58, 0.40, 0.74),   # purple
        (0.55, 0.34, 0.29),   # brown
        (0.89, 0.47, 0.76),   # pink
        (0.50, 0.50, 0.50),   # gray
        (0.74, 0.74, 0.13),   # olive
        (0.09, 0.75, 0.81),   # cyan
        (0.68, 0.78, 0.91),   # light blue
        (1.00, 0.73, 0.47),   # light orange
        (0.60, 0.87, 0.70),   # light green
        (1.00, 0.60, 0.60),   # light red
        (0.79, 0.70, 0.84),   # light purple
        (0.77, 0.69, 0.63),   # tan
        (0.98, 0.71, 0.85),   # light pink
        (0.78, 0.78, 0.78),   # light gray
        (0.86, 0.86, 0.55),   # light olive
        (0.62, 0.85, 0.90),   # light cyan
    ]

    # Build cluster to color mapping
    cluster_to_color = {}
    for i, (cluster_id, count) in enumerate(top_clusters):
        if i < len(distinct_colors):
            cluster_to_color[cluster_id] = distinct_colors[i]
        else:
            # Generate additional colors if needed
            hue = (i - len(distinct_colors)) / (max_clusters - len(distinct_colors) + 1)
            cluster_to_color[cluster_id] = _hsl_to_rgb(hue, 0.6, 0.5)

    # Build residue to cluster mapping (using residue numbers)
    resnum_to_cluster = {}
    for idx, cluster_id in cluster_assignments.items():
        if idx < len(residues):
            resnum_to_cluster[residues[idx].resnum] = cluster_id

    # Group residues by cluster for efficient commands
    cluster_residues = {}
    for resnum, cluster_id in resnum_to_cluster.items():
        if cluster_id not in cluster_residues:
            cluster_residues[cluster_id] = []
        cluster_residues[cluster_id].append(resnum)

    # Generate ChimeraX commands
    lines = [
        f"# ChimeraX commands to color {uniprot_acc} by Leiden clusters",
        f"# Top {max_clusters} clusters shown with distinct colors",
        f"# Generated by ABC Domain Predictor",
        "",
        f"# Open AlphaFold structure",
        f"alphafold fetch {uniprot_acc}",
        "",
        "# Set background and initial coloring",
        "set bgColor white",
        "color gray target a",
        "",
        f"# Color top {len(top_clusters)} clusters",
    ]

    # Add color commands for each top cluster
    for cluster_id, count in top_clusters:
        if cluster_id in cluster_residues:
            resnums = cluster_residues[cluster_id]
            color = cluster_to_color[cluster_id]
            r, g, b = [int(c * 255) for c in color]

            # Group consecutive residues for cleaner commands
            ranges = _group_consecutive(sorted(resnums))
            range_strs = [f"{s}-{e}" if s != e else str(s) for s, e in ranges]
            residue_spec = ",".join(range_strs)

            lines.append(f"# Cluster {cluster_id}: {count} residues")
            lines.append(f"color :{residue_spec} rgb({r},{g},{b})")

    # Add legend info as comments
    lines.extend([
        "",
        "# Cluster sizes (top 20):",
    ])
    for i, (cluster_id, count) in enumerate(top_clusters):
        lines.append(f"#   {i+1}. Cluster {cluster_id}: {count} residues")

    # Write file
    output_path = Path(output_path)
    output_path.write_text("\n".join(lines))
    logger.info(f"ChimeraX cluster commands saved to {output_path}")

    return str(output_path)


def _group_consecutive(nums: List[int]) -> List[Tuple[int, int]]:
    """Group consecutive numbers into ranges."""
    if not nums:
        return []

    ranges = []
    start = nums[0]
    end = nums[0]

    for n in nums[1:]:
        if n == end + 1:
            end = n
        else:
            ranges.append((start, end))
            start = n
            end = n

    ranges.append((start, end))
    return ranges
