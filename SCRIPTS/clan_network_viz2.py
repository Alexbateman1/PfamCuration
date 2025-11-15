#!/usr/bin/env python3
"""
Visualize clan relationships using multiple comparison methods:
- Foldseek (structural comparison)
- HHblits (sequence profile comparison)
- SCOOP (profile-profile comparison)

This script:
1. Queries the Pfam database to get all families in a specified clan
2. Parses results from selected comparison methods
3. Creates an interactive HTML network visualization showing relationships

Multiple edges between nodes represent different comparison methods.
"""

import mysql.connector
import pandas as pd
import networkx as nx
import re
import sys
import json
from pathlib import Path

class ClanNetworkVisualizer:
    def __init__(self, mysql_config_file):
        """
        Initialize the visualizer.

        Args:
            mysql_config_file: Path to MySQL config file (e.g., ~/.my.cnf)
        """
        self.mysql_config_file = mysql_config_file
        self.connection = None

        # Hardcoded file paths for comparison methods
        self.foldseek_file = '/nfs/production/agb/interpro/users/typhaine/interpro-pfam-curation-tools/pfam_add_clan_search/results/foldseek_all.out'
        self.hhblits_file = '/nfs/production/agb/pfam/users/agb/HHblits/RESULTS/SUMMARY/hhblits_all_vs_all_with_pfam.tsv'
        self.scoop_file = '/nfs/research/agb/research/agb/SCOOP/Latest/sump.E100.output'

        # Map Pfam types to vis.js shapes
        self.type_to_vis_shape = {
            'Domain': 'square',
            'Family': 'dot',
            'Repeat': 'triangle',
            'Coiled-coil': 'ellipse',
            'Motif': 'star',
            'Disordered': 'hexagon'
        }

        # Colorblind-friendly color scheme for each method
        # Using research-backed colorblind-safe palette
        self.method_colors = {
            'foldseek': {
                0: '#0072B2',  # Dark blue (strong)
                1: '#3A8FC7',  # Medium blue
                2: '#56B4E9'   # Light blue (weak)
            },
            'hhblits': {
                0: '#D55E00',  # Dark orange (strong)
                1: '#DD7E33',  # Medium orange
                2: '#E69F00'   # Light orange (weak)
            },
            'scoop': {
                0: '#009E73',  # Dark teal (strong)
                1: '#00B883',  # Medium teal
                2: '#00D9A3'   # Light teal (weak)
            }
        }

    def connect_to_database(self):
        """Connect to the Pfam database using the config file."""
        from configparser import ConfigParser

        config_file = Path(self.mysql_config_file).expanduser()

        if not config_file.exists():
            raise FileNotFoundError(f"MySQL config file not found: {config_file}")

        # Parse MySQL config file using ConfigParser
        config = ConfigParser()
        config.read(str(config_file))

        # Get credentials from [client] section
        if 'client' not in config:
            raise ValueError(f"No [client] section found in {config_file}")

        db_config = {
            'host': config.get('client', 'host', fallback='localhost'),
            'user': config.get('client', 'user'),
            'password': config.get('client', 'password', fallback=''),
            'database': 'pfam_live'
        }

        if config.has_option('client', 'port'):
            db_config['port'] = config.getint('client', 'port')

        try:
            self.connection = mysql.connector.connect(**db_config)
            print(f"Connected to database: {db_config['database']}")
        except mysql.connector.Error as e:
            raise Exception(f"Failed to connect to database: {e}")

    def get_clan_families(self, clan_acc):
        """
        Get all Pfam families in a clan.

        Args:
            clan_acc: Clan accession (e.g., 'CL0004')

        Returns:
            DataFrame with family information
        """
        query = """
            SELECT pfamA.pfamA_acc, pfamA_id, type, model_length, num_full, pfamA.description
            FROM clan
            JOIN clan_membership ON clan.clan_acc = clan_membership.clan_acc
            JOIN pfamA ON clan_membership.pfamA_acc = pfamA.pfamA_acc
            WHERE clan.clan_acc = %s
        """

        cursor = self.connection.cursor()
        cursor.execute(query, (clan_acc,))

        columns = ['pfamA_acc', 'pfamA_id', 'type', 'model_length', 'num_full', 'description']
        families = pd.DataFrame(cursor.fetchall(), columns=columns)
        cursor.close()

        print(f"\nFound {len(families)} families in clan {clan_acc}")
        print(families[['pfamA_acc', 'pfamA_id', 'model_length']])

        return families

    def get_family_info(self, pfam_acc):
        """
        Get information about a Pfam family including clan membership.

        Args:
            pfam_acc: Pfam accession

        Returns:
            Dictionary with family info or None
        """
        query = """
            SELECT pfamA.pfamA_acc, pfamA_id, type, model_length, num_full, description,
                   clan_membership.clan_acc
            FROM pfamA
            LEFT JOIN clan_membership ON pfamA.pfamA_acc = clan_membership.pfamA_acc
            WHERE pfamA.pfamA_acc = %s
        """

        cursor = self.connection.cursor()
        cursor.execute(query, (pfam_acc,))
        result = cursor.fetchone()
        cursor.close()

        if result:
            return {
                'pfamA_acc': result[0],
                'pfamA_id': result[1],
                'type': result[2],
                'model_length': result[3],
                'num_full': result[4],
                'description': result[5],
                'clan_acc': result[6]
            }
        return None

    def categorize_evalue(self, evalue):
        """
        Categorize E-value into significance levels.

        Args:
            evalue: E-value

        Returns:
            Category number (0-2) where 0 is most significant
        """
        if evalue < 1e-10:
            return 0  # Very significant
        elif evalue < 1e-5:
            return 1  # Significant
        else:
            return 2  # Moderate

    def categorize_pvalue(self, pvalue):
        """
        Categorize p-value into significance levels (for HHblits).

        Args:
            pvalue: p-value

        Returns:
            Category number (0-2) where 0 is most significant
        """
        if pvalue < 1e-10:
            return 0  # Very significant
        elif pvalue < 1e-5:
            return 1  # Significant
        else:
            return 2  # Moderate

    def categorize_scoop_score(self, score):
        """
        Categorize SCOOP score into significance levels.

        Args:
            score: SCOOP score

        Returns:
            Category number (0-2) where 0 is most significant
        """
        if score >= 30:
            return 0  # Very significant
        elif score >= 20:
            return 1  # Significant
        else:
            return 2  # Moderate (10-20 range)

    def parse_foldseek_results(self, clan_families, evalue_threshold=1e-3):
        """
        Parse foldseek results and extract matches involving clan families.

        Args:
            clan_families: DataFrame of families in the clan
            evalue_threshold: Maximum E-value to include

        Returns:
            Dictionary mapping (pfam1, pfam2) to {'evalue': value, 'category': cat}
        """
        print(f"\n[FOLDSEEK] Reading results from {self.foldseek_file}...")
        print(f"Filtering for E-value < {evalue_threshold}")

        # Get set of clan family accessions
        clan_pfam_accs = set(clan_families['pfamA_acc'])

        # Define column names for standard foldseek output
        col_names = [
            'query', 'target', 'identity', 'alnlen', 'mismatch', 'gapopen',
            'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits'
        ]

        # Read file in chunks
        chunk_size = 100000
        edges = {}
        total_rows = 0
        kept_rows = 0

        for chunk in pd.read_csv(self.foldseek_file, sep='\t', chunksize=chunk_size,
                                  names=col_names, header=None):
            total_rows += len(chunk)

            # Convert E-value to float and filter
            chunk['evalue'] = pd.to_numeric(chunk['evalue'], errors='coerce')
            chunk = chunk[chunk['evalue'] < evalue_threshold]

            if len(chunk) == 0:
                continue

            # Extract Pfam accessions
            chunk['Pfam_Accession'] = chunk['query'].str.extract(r'^(PF\d+)_', expand=False)
            chunk['Target_Pfam'] = chunk['target'].str.extract(r'^(PF\d+)_', expand=False)

            # Filter for clan families
            mask = (chunk['Pfam_Accession'].isin(clan_pfam_accs)) | (chunk['Target_Pfam'].isin(clan_pfam_accs))
            chunk = chunk[mask]

            if len(chunk) == 0:
                continue

            kept_rows += len(chunk)

            # Build edges dictionary
            for _, row in chunk.iterrows():
                pfam1 = row['Pfam_Accession']
                pfam2 = row['Target_Pfam']

                if pd.isna(pfam1) or pd.isna(pfam2) or pfam1 == pfam2:
                    continue

                # Use sorted tuple as key to avoid duplicates
                edge_key = tuple(sorted([pfam1, pfam2]))
                evalue = row['evalue']
                category = self.categorize_evalue(evalue)

                # Keep best (lowest) E-value if duplicate
                if edge_key not in edges or evalue < edges[edge_key]['evalue']:
                    edges[edge_key] = {'evalue': evalue, 'category': category}

            if total_rows % 1000000 == 0:
                print(f"  Processed {total_rows:,} rows, kept {kept_rows:,}")

        print(f"Total: {total_rows:,} rows -> {len(edges)} unique edges")
        return edges

    def parse_hhblits_results(self, clan_families, pvalue_threshold=0.01):
        """
        Parse HHblits results and extract matches involving clan families.

        Args:
            clan_families: DataFrame of families in the clan
            pvalue_threshold: Maximum p-value to include (column 6)

        Returns:
            Dictionary mapping (pfam1, pfam2) to {'pvalue': value, 'category': cat}
        """
        print(f"\n[HHBLITS] Reading results from {self.hhblits_file}...")
        print(f"Filtering for p-value < {pvalue_threshold}")

        clan_pfam_accs = set(clan_families['pfamA_acc'])

        # Read file in chunks
        chunk_size = 100000
        edges = {}
        total_rows = 0
        kept_rows = 0
        first_chunk = True

        try:
            for chunk in pd.read_csv(self.hhblits_file, sep='\t', chunksize=chunk_size):
                total_rows += len(chunk)

                # Debug: show column names on first chunk
                if first_chunk:
                    print(f"  DEBUG: HHblits file columns: {list(chunk.columns)}")
                    print(f"  DEBUG: First row sample:")
                    if len(chunk) > 0:
                        print(f"    query_family: {chunk.iloc[0]['query_family']}")
                        print(f"    target_pfam_family: {chunk.iloc[0]['target_pfam_family']}")
                        print(f"    p_value: {chunk.iloc[0]['p_value']}")
                    first_chunk = False

                # Convert p-value to float and filter
                chunk['p_value'] = pd.to_numeric(chunk['p_value'], errors='coerce')
                pre_filter_count = len(chunk)
                chunk = chunk[chunk['p_value'] < pvalue_threshold]
                print(f"  DEBUG: p-value filter: {pre_filter_count} -> {len(chunk)} rows")

                if len(chunk) == 0:
                    continue

                # Filter for clan families
                mask = (chunk['query_family'].isin(clan_pfam_accs)) | (chunk['target_pfam_family'].isin(clan_pfam_accs))
                pre_clan_count = len(chunk)
                chunk = chunk[mask]
                print(f"  DEBUG: Clan filter: {pre_clan_count} -> {len(chunk)} rows")

                if len(chunk) == 0:
                    continue

                kept_rows += len(chunk)

                # Build edges dictionary
                for _, row in chunk.iterrows():
                    pfam1 = row['query_family']
                    pfam2 = row['target_pfam_family']

                    if pd.isna(pfam1) or pd.isna(pfam2) or pfam1 == pfam2:
                        continue

                    edge_key = tuple(sorted([pfam1, pfam2]))
                    pvalue = row['p_value']
                    category = self.categorize_pvalue(pvalue)

                    if edge_key not in edges or pvalue < edges[edge_key]['pvalue']:
                        edges[edge_key] = {'pvalue': pvalue, 'category': category}

                if total_rows % 1000000 == 0:
                    print(f"  Processed {total_rows:,} rows, kept {kept_rows:,}")

        except Exception as e:
            print(f"  ERROR parsing HHblits file: {e}")
            import traceback
            traceback.print_exc()
            return {}

        print(f"Total: {total_rows:,} rows -> {len(edges)} unique edges")
        return edges

    def parse_scoop_results(self, clan_families, score_threshold=10.0):
        """
        Parse SCOOP results and extract matches involving clan families.

        Args:
            clan_families: DataFrame of families in the clan
            score_threshold: Minimum SCOOP score to include

        Returns:
            Dictionary mapping (pfam1, pfam2) to {'score': value, 'category': cat}
        """
        print(f"\n[SCOOP] Reading results from {self.scoop_file}...")
        print(f"Filtering for SCOOP score >= {score_threshold}")

        clan_pfam_accs = set(clan_families['pfamA_acc'])

        edges = {}
        total_rows = 0
        kept_rows = 0

        # SCOOP file format: score pfam1 num1 pfam2 num2 direction num3
        with open(self.scoop_file, 'r') as f:
            for line in f:
                total_rows += 1
                parts = line.strip().split()

                if len(parts) < 4:
                    continue

                try:
                    score = float(parts[0])
                    pfam1 = parts[1]
                    pfam2 = parts[3]
                except (ValueError, IndexError):
                    continue

                if score < score_threshold:
                    continue

                # Filter for clan families
                if pfam1 not in clan_pfam_accs and pfam2 not in clan_pfam_accs:
                    continue

                if pfam1 == pfam2:
                    continue

                kept_rows += 1

                edge_key = tuple(sorted([pfam1, pfam2]))
                category = self.categorize_scoop_score(score)

                # Keep best (highest) score if duplicate
                if edge_key not in edges or score > edges[edge_key]['score']:
                    edges[edge_key] = {'score': score, 'category': category}

                if total_rows % 1000000 == 0:
                    print(f"  Processed {total_rows:,} rows, kept {kept_rows:,}")

        print(f"Total: {total_rows:,} rows -> {len(edges)} unique edges")
        return edges

    def combine_method_results(self, foldseek_edges, hhblits_edges, scoop_edges, selected_methods):
        """
        Combine results from different methods into a unified structure.

        Args:
            foldseek_edges: Dictionary of foldseek edges
            hhblits_edges: Dictionary of hhblits edges
            scoop_edges: Dictionary of scoop edges
            selected_methods: List of method names to include

        Returns:
            Dictionary mapping (pfam1, pfam2) to method-specific data
        """
        combined = {}

        if 'foldseek' in selected_methods and foldseek_edges:
            for edge_key, data in foldseek_edges.items():
                if edge_key not in combined:
                    combined[edge_key] = {}
                combined[edge_key]['foldseek'] = data

        if 'hhblits' in selected_methods and hhblits_edges:
            for edge_key, data in hhblits_edges.items():
                if edge_key not in combined:
                    combined[edge_key] = {}
                combined[edge_key]['hhblits'] = data

        if 'scoop' in selected_methods and scoop_edges:
            for edge_key, data in scoop_edges.items():
                if edge_key not in combined:
                    combined[edge_key] = {}
                combined[edge_key]['scoop'] = data

        return combined

    def create_network(self, clan_acc, clan_families, combined_edges):
        """
        Create network graph from combined edges.

        Args:
            clan_acc: Clan accession
            clan_families: DataFrame of clan families
            combined_edges: Dictionary of combined edge data

        Returns:
            NetworkX graph with nodes and combined_edges dict
        """
        G = nx.Graph()

        # Add clan family nodes
        clan_pfam_accs = set(clan_families['pfamA_acc'])

        for _, fam in clan_families.iterrows():
            G.add_node(fam['pfamA_acc'],
                      pfam_id=fam['pfamA_id'],
                      model_length=fam['model_length'],
                      type=fam['type'],
                      in_clan=True,
                      clan_acc=clan_acc,
                      description=fam['description'])

        # Collect all Pfam accessions from edges
        all_pfams = set()
        for edge_key in combined_edges.keys():
            all_pfams.add(edge_key[0])
            all_pfams.add(edge_key[1])

        # Add nodes for external families
        for pfam_acc in all_pfams:
            if pfam_acc not in G:
                pfam_info = self.get_family_info(pfam_acc)
                if pfam_info:
                    G.add_node(pfam_acc,
                              pfam_id=pfam_info['pfamA_id'],
                              model_length=pfam_info['model_length'],
                              type=pfam_info['type'],
                              in_clan=pfam_acc in clan_pfam_accs,
                              clan_acc=pfam_info['clan_acc'],
                              description=pfam_info['description'])
                else:
                    G.add_node(pfam_acc,
                              pfam_id=pfam_acc,
                              model_length=100,
                              type='Domain',
                              in_clan=pfam_acc in clan_pfam_accs,
                              clan_acc=None,
                              description='Unknown')

        # Add placeholder edges (actual edge rendering happens in HTML)
        for edge_key in combined_edges.keys():
            if not G.has_edge(edge_key[0], edge_key[1]):
                G.add_edge(edge_key[0], edge_key[1])

        print(f"\nNetwork statistics:")
        print(f"  Nodes: {G.number_of_nodes()}")
        print(f"  Unique edge pairs: {len(combined_edges)}")
        print(f"  Target clan families: {sum(1 for n, d in G.nodes(data=True) if d.get('in_clan', False))}")

        no_clan_count = sum(1 for n, d in G.nodes(data=True)
                           if not d.get('in_clan', False) and d.get('clan_acc') is None)
        different_clan_count = sum(1 for n, d in G.nodes(data=True)
                                  if not d.get('in_clan', False) and d.get('clan_acc') is not None)

        print(f"  Families with no clan: {no_clan_count}")
        print(f"  Families in different clan: {different_clan_count}")

        return G, combined_edges

    def export_interactive_html(self, G, combined_edges, clan_acc, selected_methods, output_file):
        """
        Export network as interactive HTML using vis.js with multiple edges per node pair.

        Args:
            G: NetworkX graph
            combined_edges: Dictionary of combined edge data
            clan_acc: Clan accession
            selected_methods: List of selected method names
            output_file: Output HTML file path
        """
        # Prepare nodes data
        nodes = []
        node_titles = {}
        for node, data in G.nodes(data=True):
            node_clan = data.get('clan_acc', None)
            in_target_clan = data.get('in_clan', False)

            # Determine color based on clan membership
            if in_target_clan:
                color = '#4CAF50'  # Green
                border_color = '#2E7D32'
            elif node_clan is None:
                color = '#BDBDBD'  # Grey
                border_color = '#757575'
            else:
                color = '#F44336'  # Red
                border_color = '#C62828'

            # Node size based on model length
            model_length = data.get('model_length', 100)
            size = 20 + (model_length / 10)

            # Get shape based on type
            node_type = data.get('type', 'Domain')
            shape = self.type_to_vis_shape.get(node_type, 'dot')

            # Handle ellipse shapes specially
            if shape == 'ellipse':
                font_config = {
                    'size': 18,
                    'color': '#000000',
                    'face': 'Arial',
                    'bold': {'size': 20},
                    'vadjust': 25
                }
                node_config = {
                    'id': node,
                    'label': f"{data.get('pfam_id', node)}\n{node}",
                    'shape': shape,
                    'color': {
                        'background': color,
                        'border': border_color,
                        'highlight': {'background': color, 'border': '#000000'}
                    },
                    'size': size,
                    'font': font_config,
                    'widthConstraint': {'minimum': size * 1.5, 'maximum': size * 1.5},
                    'heightConstraint': {'minimum': size, 'maximum': size}
                }
            else:
                node_config = {
                    'id': node,
                    'label': f"{data.get('pfam_id', node)}\n{node}",
                    'shape': shape,
                    'color': {
                        'background': color,
                        'border': border_color,
                        'highlight': {'background': color, 'border': '#000000'}
                    },
                    'size': size,
                    'font': {'size': 18, 'color': '#000000', 'face': 'Arial', 'bold': {'size': 20}}
                }

            node_titles[node] = (f"<div style='font-family: Arial; padding: 5px;'>"
                        f"<b style='font-size: 16px;'>{data.get('pfam_id', node)}</b><br/>"
                        f"<span style='color: #666; font-size: 14px;'>{node}</span><br/>"
                        f"<b style='font-size: 13px;'>Type:</b> <span style='font-size: 13px;'>{node_type}</span><br/>"
                        f"<b style='font-size: 13px;'>Length:</b> <span style='font-size: 13px;'>{model_length}</span><br/>"
                        f"<b style='font-size: 13px;'>Clan:</b> <span style='font-size: 13px;'>{node_clan if node_clan else 'None'}</span><br/>"
                        f"<span style='font-size: 12px; color: #555;'>{data.get('description', '')}</span>"
                        f"</div>")

            nodes.append(node_config)

        # Prepare edges data - create separate edge for each method
        edges = []
        edge_titles = {}

        # Define offset for each method to separate multiple edges visually
        method_offsets = {
            'foldseek': -0.2,   # Curve slightly to one side
            'hhblits': 0,       # Straight
            'scoop': 0.2        # Curve slightly to other side
        }

        for edge_key, method_data in combined_edges.items():
            pfam1, pfam2 = edge_key

            # Create separate edges for each method present
            for method in selected_methods:
                if method not in method_data:
                    continue

                data = method_data[method]
                category = data['category']
                color = self.method_colors[method][category]

                # Width based on category
                if category == 0:
                    width = 3
                elif category == 1:
                    width = 2
                else:
                    width = 1

                # Create unique edge ID
                edge_id = f"{pfam1}-{pfam2}-{method}"

                # Build tooltip
                if method == 'scoop':
                    score_or_eval = f"SCOOP Score: {data['score']:.2f}"
                elif method == 'hhblits':
                    score_or_eval = f"p-value: {data['pvalue']:.2e}"
                else:
                    score_or_eval = f"E-value: {data['evalue']:.2e}"

                edge_titles[edge_id] = (f"<div style='font-family: Arial; padding: 5px;'>"
                            f"<b style='font-size: 15px;'>{method.upper()}</b><br/>"
                            f"<b style='font-size: 14px;'>{score_or_eval}</b><br/>"
                            f"<b style='font-size: 13px;'>Between:</b><br/>"
                            f"<span style='color: #666; font-size: 13px;'>{pfam1}</span><br/>"
                            f"<span style='color: #666; font-size: 13px;'>{pfam2}</span>"
                            f"</div>")

                # Use different smooth settings to offset edges
                roundness = method_offsets.get(method, 0)
                if roundness == 0:
                    # Straight edge - use False directly for vis.js
                    smooth_config = False
                else:
                    # Slightly curved to separate from other edges
                    smooth_config = {'enabled': True, 'type': 'curvedCW' if roundness > 0 else 'curvedCCW', 'roundness': abs(roundness)}

                edges.append({
                    'id': edge_id,
                    'from': pfam1,
                    'to': pfam2,
                    'value': width,
                    'color': {'color': color, 'highlight': '#000000'},
                    'smooth': smooth_config
                })

        # Debug: count edges per method
        edge_counts = {method: 0 for method in selected_methods}
        for edge in edges:
            for method in selected_methods:
                if edge['id'].endswith(f'-{method}'):
                    edge_counts[method] += 1

        print(f"\nEdge counts by method:")
        for method, count in edge_counts.items():
            print(f"  {method.upper()}: {count} edges")
        print(f"  Total edges in visualization: {len(edges)}")

        # Build legend items for selected methods
        legend_html = self._build_legend_html(selected_methods)

        # Create HTML content
        html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>Clan {clan_acc} Network Visualization</title>
    <script type="text/javascript" src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
    <style type="text/css">
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        #mynetwork {{
            width: 100%;
            height: 800px;
            border: 1px solid #ddd;
            background-color: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            cursor: default;
        }}
        #mynetwork canvas {{
            cursor: pointer !important;
        }}
        h1 {{
            color: #333;
            text-align: center;
            margin-bottom: 10px;
        }}
        .subtitle {{
            text-align: center;
            color: #666;
            margin-bottom: 20px;
        }}
        .legend {{
            background-color: white;
            padding: 15px;
            border: 1px solid #ddd;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .legend-item {{
            display: inline-block;
            margin-right: 20px;
            margin-bottom: 5px;
        }}
        .legend-circle {{
            display: inline-block;
            width: 20px;
            height: 20px;
            border-radius: 50%;
            margin-right: 5px;
            vertical-align: middle;
            border: 2px solid;
        }}
        .legend-line {{
            display: inline-block;
            width: 30px;
            height: 3px;
            margin-right: 5px;
            vertical-align: middle;
        }}
        .method-section {{
            margin-bottom: 10px;
        }}
        .vis-tooltip {{
            position: absolute;
            visibility: hidden;
            padding: 12px;
            white-space: nowrap;
            font-family: Arial;
            font-size: 14px;
            color: #000;
            background-color: rgba(255, 255, 255, 0.98);
            border: 2px solid #333;
            border-radius: 6px;
            box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.4);
            pointer-events: none;
            z-index: 5;
            max-width: 400px;
        }}
    </style>
</head>
<body>
    <h1>Structural Relationships for Clan {clan_acc}</h1>
    <div class="subtitle">
        Interactive Network Visualization ({', '.join([m.upper() for m in selected_methods])})<br/>
        <small>Click nodes to open InterPro | Drag nodes to reposition | Hover for details | Scroll to zoom | Navigation buttons in bottom-left</small>
    </div>

    <div style="text-align: center; margin-bottom: 15px;">
        <button id="editModeBtn" onclick="toggleEditMode()" style="
            padding: 10px 20px;
            font-size: 14px;
            background-color: #4CAF50;
            color: white;
            border: none;
            border-radius: 5px;
            cursor: pointer;
            box-shadow: 0 2px 4px rgba(0,0,0,0.2);
        ">
            🔓 Editing Mode: OFF (Click opens InterPro)
        </button>
        <div id="editModeStatus" style="font-size: 12px; color: #666; margin-top: 5px;">
            Click the button to enable editing mode for easier node repositioning
        </div>
    </div>

    <div style="text-align: center; margin-bottom: 15px; background-color: white; padding: 15px; border: 1px solid #ddd; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
        <label for="sizeSlider" style="font-weight: bold; margin-right: 10px;">Node & Text Scale:</label>
        <input type="range" id="sizeSlider" min="0.3" max="3.0" step="0.1" value="1.0"
               style="width: 300px; vertical-align: middle;">
        <span id="sizeValue" style="margin-left: 10px; font-weight: bold; color: #4CAF50;">100%</span>
        <div style="font-size: 12px; color: #666; margin-top: 5px;">
            Adjust the size of all nodes and text labels
        </div>
    </div>

    <div class="legend">
        <strong>Node Colors:</strong>
        <div class="legend-item">
            <span class="legend-circle" style="background-color: #4CAF50; border-color: #2E7D32;"></span>
            Target clan families
        </div>
        <div class="legend-item">
            <span class="legend-circle" style="background-color: #BDBDBD; border-color: #757575;"></span>
            Families with no clan
        </div>
        <div class="legend-item">
            <span class="legend-circle" style="background-color: #F44336; border-color: #C62828;"></span>
            Families in different clan
        </div>
        <br/><br/>
        {legend_html}
    </div>

    <div id="mynetwork"></div>

    <script type="text/javascript">
        // Editing mode state
        var editingMode = false;

        function toggleEditMode() {{
            editingMode = !editingMode;
            var btn = document.getElementById('editModeBtn');
            var status = document.getElementById('editModeStatus');

            if (editingMode) {{
                btn.textContent = '🔒 Editing Mode: ON (Click moves nodes)';
                btn.style.backgroundColor = '#FF9800';
                status.textContent = 'Editing mode active - click and drag to reposition nodes';
            }} else {{
                btn.textContent = '🔓 Editing Mode: OFF (Click opens InterPro)';
                btn.style.backgroundColor = '#4CAF50';
                status.textContent = 'Click nodes to open their InterPro pages';
            }}
        }}

        // Create network data
        var nodes = new vis.DataSet({json.dumps(nodes, indent=8)});

        var edges = new vis.DataSet({json.dumps(edges, indent=8)});

        // Store titles separately
        var nodeTitles = {json.dumps(node_titles)};
        var edgeTitles = {json.dumps(edge_titles)};

        // Find isolated nodes and position them in a grid
        var connectedNodes = new Set();
        edges.forEach(function(edge) {{
            connectedNodes.add(edge.from);
            connectedNodes.add(edge.to);
        }});

        var isolatedNodes = [];
        nodes.forEach(function(node) {{
            if (!connectedNodes.has(node.id)) {{
                isolatedNodes.push(node.id);
            }}
        }});

        console.log('Found ' + isolatedNodes.length + ' isolated nodes (will be arranged in grid)');

        // Calculate grid dimensions for isolated nodes
        var gridCols = Math.ceil(Math.sqrt(isolatedNodes.length));
        var gridRows = Math.ceil(isolatedNodes.length / gridCols);
        var gridSpacing = 150;
        var gridStartX = -2500;  // Position grid far on the left
        var gridStartY = -400;

        // Position isolated nodes in a grid
        isolatedNodes.forEach(function(nodeId, index) {{
            var col = index % gridCols;
            var row = Math.floor(index / gridCols);
            var x = gridStartX + (col * gridSpacing);
            var y = gridStartY + (row * gridSpacing);

            nodes.update({{
                id: nodeId,
                x: x,
                y: y,
                fixed: {{ x: true, y: true }}
            }});
        }});

        // Create network
        var container = document.getElementById('mynetwork');
        var data = {{
            nodes: nodes,
            edges: edges
        }};

        var options = {{
            nodes: {{
                shape: 'dot',
                borderWidth: 2,
                shadow: true,
                font: {{
                    size: 18,
                    face: 'Arial',
                    color: '#000000'
                }}
            }},
            edges: {{
                width: 1,
                shadow: true
                // Individual edge smooth settings will override
            }},
            physics: {{
                enabled: true,
                stabilization: {{
                    enabled: true,
                    iterations: 200,
                    onlyDynamicEdges: false,
                    fit: true
                }},
                barnesHut: {{
                    gravitationalConstant: -8000,
                    centralGravity: 0.3,
                    springLength: 150,
                    springConstant: 0.04,
                    damping: 0.09,
                    avoidOverlap: 0.5
                }}
            }},
            interaction: {{
                dragNodes: true,
                dragView: true,
                zoomView: true,
                hover: true,
                tooltipDelay: 100,
                navigationButtons: true,
                keyboard: true
            }},
            layout: {{
                improvedLayout: true
            }}
        }};

        var network = new vis.Network(container, data, options);

        // Disable physics after stabilization
        network.on('stabilizationIterationsDone', function() {{
            isolatedNodes.forEach(function(nodeId) {{
                nodes.update({{
                    id: nodeId,
                    fixed: false
                }});
            }});

            network.setOptions({{physics: false}});
            console.log('Physics disabled - nodes will stay in place');
        }});

        setTimeout(function() {{
            isolatedNodes.forEach(function(nodeId) {{
                nodes.update({{
                    id: nodeId,
                    fixed: false
                }});
            }});

            network.setOptions({{physics: false}});
        }}, 3000);

        // Create custom tooltip element
        var tooltip = document.createElement('div');
        tooltip.className = 'vis-tooltip';
        document.body.appendChild(tooltip);

        // Handle node clicks
        network.on('click', function(params) {{
            if (params.nodes.length > 0 && !editingMode) {{
                var nodeId = params.nodes[0];
                var url = 'https://www.ebi.ac.uk/interpro/entry/pfam/' + nodeId + '/';
                window.open(url, '_blank');
            }}
        }});

        // Show custom tooltip on hover
        network.on('hoverNode', function(params) {{
            var nodeId = params.node;
            var title = nodeTitles[nodeId];
            if (title) {{
                tooltip.innerHTML = title;
                tooltip.style.visibility = 'visible';
                tooltip.style.left = params.event.pageX + 10 + 'px';
                tooltip.style.top = params.event.pageY + 10 + 'px';
            }}
        }});

        network.on('hoverEdge', function(params) {{
            var edgeId = params.edge;
            var title = edgeTitles[edgeId];
            if (title) {{
                tooltip.innerHTML = title;
                tooltip.style.visibility = 'visible';
                tooltip.style.left = params.event.pageX + 10 + 'px';
                tooltip.style.top = params.event.pageY + 10 + 'px';
            }}
        }});

        network.on('blurNode', function() {{
            tooltip.style.visibility = 'hidden';
        }});

        network.on('blurEdge', function() {{
            tooltip.style.visibility = 'hidden';
        }});

        network.on('dragStart', function() {{
            tooltip.style.visibility = 'hidden';
        }});

        // Store original node sizes and font sizes for scaling
        var originalNodeData = {{}};
        nodes.forEach(function(node) {{
            originalNodeData[node.id] = {{
                size: node.size,
                fontSize: node.font ? node.font.size : 18,
                boldSize: node.font && node.font.bold ? node.font.bold.size : 20,
                vadjust: node.font && node.font.vadjust ? node.font.vadjust : 0,
                widthMin: node.widthConstraint ? node.widthConstraint.minimum : null,
                widthMax: node.widthConstraint ? node.widthConstraint.maximum : null,
                heightMin: node.heightConstraint ? node.heightConstraint.minimum : null,
                heightMax: node.heightConstraint ? node.heightConstraint.maximum : null
            }};
        }});

        // Handle size slider changes
        var sizeSlider = document.getElementById('sizeSlider');
        var sizeValue = document.getElementById('sizeValue');

        sizeSlider.addEventListener('input', function() {{
            var scale = parseFloat(this.value);
            sizeValue.textContent = Math.round(scale * 100) + '%';

            // Update all nodes with scaled sizes
            var updates = [];
            nodes.forEach(function(node) {{
                var original = originalNodeData[node.id];
                var update = {{
                    id: node.id,
                    size: original.size * scale,
                    font: {{
                        size: Math.round(original.fontSize * scale),
                        color: node.font.color,
                        face: node.font.face,
                        bold: {{ size: Math.round(original.boldSize * scale) }}
                    }}
                }};

                if (original.vadjust) {{
                    update.font.vadjust = Math.round(original.vadjust * scale);
                }}

                if (original.widthMin !== null) {{
                    update.widthConstraint = {{
                        minimum: original.widthMin * scale,
                        maximum: original.widthMax * scale
                    }};
                }}
                if (original.heightMin !== null) {{
                    update.heightConstraint = {{
                        minimum: original.heightMin * scale,
                        maximum: original.heightMax * scale
                    }};
                }}

                updates.push(update);
            }});

            nodes.update(updates);
        }});

        console.log('Network created with ' + nodes.length + ' nodes and ' + edges.length + ' edges');
        console.log('Isolated nodes arranged in grid on left side: ' + isolatedNodes.length);
        console.log('Methods: {', '.join(selected_methods)}');
    </script>
</body>
</html>'''

        # Write to file
        with open(output_file, 'w') as f:
            f.write(html_content)

        print(f"\nSaved interactive HTML to {output_file}")

    def _build_legend_html(self, selected_methods):
        """Build HTML for method-specific legend items."""
        legend_parts = []

        method_names = {
            'foldseek': 'Foldseek (structural)',
            'hhblits': 'HHblits (profile)',
            'scoop': 'SCOOP (profile-profile)'
        }

        for method in selected_methods:
            colors = self.method_colors[method]
            name = method_names[method]

            if method == 'scoop':
                cat0_label = 'Score ≥ 30'
                cat1_label = 'Score 20-30'
                cat2_label = 'Score 10-20'
            elif method == 'hhblits':
                cat0_label = 'p-value < 1e-10'
                cat1_label = 'p-value 1e-10 to 1e-5'
                cat2_label = 'p-value 1e-5 to threshold'
            else:
                cat0_label = 'E-value < 1e-10'
                cat1_label = 'E-value 1e-10 to 1e-5'
                cat2_label = 'E-value 1e-5 to threshold'

            legend_parts.append(f'''
        <div class="method-section">
            <strong>{name}:</strong>
            <div class="legend-item">
                <span class="legend-line" style="background-color: {colors[0]}; height: 3px;"></span>
                {cat0_label}
            </div>
            <div class="legend-item">
                <span class="legend-line" style="background-color: {colors[1]}; height: 2px;"></span>
                {cat1_label}
            </div>
            <div class="legend-item">
                <span class="legend-line" style="background-color: {colors[2]}; height: 1px;"></span>
                {cat2_label}
            </div>
        </div>''')

        return '\n'.join(legend_parts)

    def analyze_clan(self, clan_acc, selected_methods, evalue_threshold=1e-5, hhblits_pvalue_threshold=0.01, scoop_threshold=10.0):
        """
        Complete analysis pipeline for a clan.

        Args:
            clan_acc: Clan accession (e.g., 'CL0004')
            selected_methods: List of methods to use (['foldseek'], ['hhblits'], ['scoop'], or combinations)
            evalue_threshold: Maximum E-value for foldseek
            hhblits_pvalue_threshold: Maximum p-value for hhblits (default 0.01)
            scoop_threshold: Minimum SCOOP score
        """
        print(f"\n{'='*60}")
        print(f"Analyzing clan {clan_acc}")
        print(f"Methods: {', '.join([m.upper() for m in selected_methods])}")
        print(f"{'='*60}")

        # Get clan families
        clan_families = self.get_clan_families(clan_acc)

        if len(clan_families) == 0:
            print(f"No families found for clan {clan_acc}")
            return

        # Parse results from selected methods
        foldseek_edges = {}
        hhblits_edges = {}
        scoop_edges = {}

        if 'foldseek' in selected_methods:
            foldseek_edges = self.parse_foldseek_results(clan_families, evalue_threshold)

        if 'hhblits' in selected_methods:
            hhblits_edges = self.parse_hhblits_results(clan_families, hhblits_pvalue_threshold)

        if 'scoop' in selected_methods:
            scoop_edges = self.parse_scoop_results(clan_families, scoop_threshold)

        # Combine results
        combined_edges = self.combine_method_results(foldseek_edges, hhblits_edges, scoop_edges, selected_methods)

        if len(combined_edges) == 0:
            print(f"No matches found for clan {clan_acc}")
            return

        # Create network
        G, combined_edges = self.create_network(clan_acc, clan_families, combined_edges)

        # Create interactive HTML
        methods_str = '_'.join(selected_methods)
        html_file = f'{clan_acc}_{methods_str}_interactive.html'
        self.export_interactive_html(G, combined_edges, clan_acc, selected_methods, html_file)

        print(f"\nView the interactive HTML: {html_file}")

        # Print statistics
        print(f"\nTop external families by number of connections:")
        external_degrees = [(n, G.degree(n)) for n in G.nodes()
                           if not G.nodes[n].get('in_clan', False)]
        external_degrees.sort(key=lambda x: x[1], reverse=True)

        for pfam_acc, degree in external_degrees[:10]:
            pfam_id = G.nodes[pfam_acc].get('pfam_id', pfam_acc)
            print(f"  {pfam_acc} ({pfam_id}): {degree} connections")

    def close(self):
        """Close database connection."""
        if self.connection:
            self.connection.close()


def main():
    """Main function."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Visualize Pfam clan relationships using multiple comparison methods',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Single method
  %(prog)s CL0004 --foldseek
  %(prog)s CL0004 --hhblits
  %(prog)s CL0004 --scoop

  # Multiple methods
  %(prog)s CL0004 --foldseek --hhblits
  %(prog)s CL0004 --foldseek --scoop

  # All methods
  %(prog)s CL0004 --all

  # Custom thresholds
  %(prog)s CL0004 --all --evalue-threshold 1e-10 --hhblits-pvalue-threshold 0.001 --scoop-threshold 20
        '''
    )

    parser.add_argument('clan', help='Clan accession (e.g., CL0004)')
    parser.add_argument('--foldseek', action='store_true', help='Include foldseek structural comparisons')
    parser.add_argument('--hhblits', action='store_true', help='Include HHblits profile comparisons')
    parser.add_argument('--scoop', action='store_true', help='Include SCOOP profile-profile comparisons')
    parser.add_argument('--all', action='store_true', help='Include all comparison methods')
    parser.add_argument('--evalue-threshold', type=float, default=1e-5,
                       help='E-value threshold for foldseek (default: 1e-5)')
    parser.add_argument('--hhblits-pvalue-threshold', type=float, default=0.01,
                       help='p-value threshold for hhblits (default: 0.01)')
    parser.add_argument('--scoop-threshold', type=float, default=10.0,
                       help='Minimum SCOOP score threshold (default: 10.0)')
    parser.add_argument('--mysql-config', default='~/.my.cnf',
                       help='Path to MySQL config file (default: ~/.my.cnf)')

    args = parser.parse_args()

    # Determine which methods to use
    selected_methods = []
    if args.all:
        selected_methods = ['foldseek', 'hhblits', 'scoop']
    else:
        if args.foldseek:
            selected_methods.append('foldseek')
        if args.hhblits:
            selected_methods.append('hhblits')
        if args.scoop:
            selected_methods.append('scoop')

    if not selected_methods:
        print("Error: Must specify at least one method (--foldseek, --hhblits, --scoop, or --all)")
        sys.exit(1)

    # Create visualizer
    viz = ClanNetworkVisualizer(args.mysql_config)

    try:
        viz.connect_to_database()
        viz.analyze_clan(args.clan, selected_methods, args.evalue_threshold, args.hhblits_pvalue_threshold, args.scoop_threshold)
    finally:
        viz.close()


if __name__ == '__main__':
    main()
