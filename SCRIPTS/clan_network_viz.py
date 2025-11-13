#!/usr/bin/env python3
"""
Visualize clan relationships using foldseek structural comparison data.

This script:
1. Queries the Pfam database to get all families in a specified clan
2. Parses foldseek results to find structural matches
3. Creates a network visualization showing relationships between families

For best layout results, install graphviz support:
  - pygraphviz: pip install pygraphviz (may need system graphviz first)
  - pydot: pip install pydot (easier, recommended)
Script uses neato layout algorithm which works well for small networks.
Falls back to spring layout if graphviz is unavailable.
"""

import mysql.connector
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import re
import sys
import json
from pathlib import Path

class ClanNetworkVisualizer:
    def __init__(self, foldseek_file, mysql_config_file):
        """
        Initialize the visualizer.

        Args:
            foldseek_file: Path to foldseek results file
            mysql_config_file: Path to MySQL config file (e.g., ~/.my.cnf)
        """
        self.foldseek_file = foldseek_file
        self.mysql_config_file = mysql_config_file
        self.connection = None

        # Map Pfam types to shapes for visualization
        self.type_to_shape = {
            'Domain': 's',       # Square (matplotlib) / 'box' (vis.js)
            'Family': 'o',       # Circle (matplotlib) / 'dot' (vis.js)
            'Repeat': '^',       # Triangle (matplotlib) / 'triangle' (vis.js)
            'Coiled-coil': 'D',  # Diamond/Oval (matplotlib) / 'ellipse' (vis.js)
            'Motif': '8',        # Octagon (matplotlib) / 'star' (vis.js) - closest to octagon
            'Disordered': 'h'    # Hexagon (matplotlib) / 'hexagon' (vis.js)
        }

        self.type_to_vis_shape = {
            'Domain': 'square',  # 'square' behaves like 'dot' with label outside
            'Family': 'dot',
            'Repeat': 'triangle',
            'Coiled-coil': 'ellipse',
            'Motif': 'star',
            'Disordered': 'hexagon'
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

        # Debug: show what sections were found
        print(f"DEBUG: Config file: {config_file}")
        print(f"DEBUG: Sections found: {config.sections()}")

        # Get credentials from [client] section
        if 'client' not in config:
            raise ValueError(f"No [client] section found in {config_file}")

        # Debug: show what options are in [client] section
        print(f"DEBUG: Options in [client]: {config.options('client')}")

        db_config = {
            'host': config.get('client', 'host', fallback='localhost'),
            'user': config.get('client', 'user'),
            'password': config.get('client', 'password', fallback=''),
            'database': 'pfam_live'
        }

        # Debug: show config (without exposing password)
        debug_config = {k: ('***' if k == 'password' and v else v) for k, v in db_config.items()}
        print(f"DEBUG: db_config: {debug_config}")
        print(f"DEBUG: Has password: {bool(db_config.get('password'))}")

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
    
    def extract_target_pfam(self, target_full):
        """
        Extract Pfam accession from Target_Full field.
        
        Format: PFxxxxx_UniProtID_domain_resStart-End
        
        Args:
            target_full: Target string from foldseek results
            
        Returns:
            Pfam accession or None
        """
        match = re.match(r'(PF\d+)_', target_full)
        if match:
            return match.group(1)
        return None
    
    def parse_foldseek_results(self, clan_families, evalue_threshold=1e-3):
        """
        Parse foldseek results and extract matches involving clan families.
        
        Args:
            clan_families: DataFrame of families in the clan
            evalue_threshold: Maximum E-value to include
            
        Returns:
            DataFrame of filtered matches
        """
        print(f"\nReading foldseek results from {self.foldseek_file}...")
        print(f"Filtering for E-value < {evalue_threshold} during read to save memory...")
        
        # Get set of clan family accessions
        clan_pfam_accs = set(clan_families['pfamA_acc'])
        
        # Define column names for standard foldseek output (12 columns, no header)
        col_names = [
            'query', 'target', 'identity', 'alnlen', 'mismatch', 'gapopen',
            'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits'
        ]
        
        # Read file in chunks and filter immediately to save memory
        chunk_size = 100000
        filtered_chunks = []
        total_rows = 0
        kept_rows = 0
        clan_rows = 0
        
        for chunk in pd.read_csv(self.foldseek_file, sep='\t', chunksize=chunk_size, 
                                  names=col_names, header=None):
            total_rows += len(chunk)
            
            # Convert E-value to float and filter immediately
            chunk['evalue'] = pd.to_numeric(chunk['evalue'], errors='coerce')
            chunk = chunk[chunk['evalue'] < evalue_threshold]
            
            if len(chunk) == 0:
                continue
            
            kept_rows += len(chunk)
            
            # Extract Pfam accessions from query and target strings
            chunk['Pfam_Accession'] = chunk['query'].str.extract(r'^(PF\d+)_', expand=False)
            chunk['Target_Pfam'] = chunk['target'].str.extract(r'^(PF\d+)_', expand=False)
            
            # Filter for clan families IMMEDIATELY to save memory
            mask = (chunk['Pfam_Accession'].isin(clan_pfam_accs)) | (chunk['Target_Pfam'].isin(clan_pfam_accs))
            clan_chunk = chunk[mask].copy()
            
            if len(clan_chunk) == 0:
                continue
            
            clan_rows += len(clan_chunk)
            
            # Add remaining columns only for clan-relevant rows
            clan_chunk['Target_Full'] = clan_chunk['target']
            clan_chunk['E-value'] = clan_chunk['evalue']
            
            filtered_chunks.append(clan_chunk)
            
            # Progress indicator showing memory savings
            if total_rows % 1000000 == 0:
                print(f"  Processed {total_rows:,} rows, kept {kept_rows:,} with E-value < {evalue_threshold}, {clan_rows:,} relevant to clan")
        
        if not filtered_chunks:
            print(f"No matches found with E-value < {evalue_threshold}")
            return pd.DataFrame()
        
        # Combine all filtered chunks
        df = pd.concat(filtered_chunks, ignore_index=True)
        
        print(f"\nTotal rows processed: {total_rows:,}")
        print(f"Rows with E-value < {evalue_threshold}: {kept_rows:,}")
        print(f"Rows relevant to clan (kept in memory): {clan_rows:,}")
        print(f"Memory efficiency: {100 * clan_rows / kept_rows:.1f}% of E-value filtered rows kept")
        print(f"Matches involving clan families: {len(df)}")
        
        return df
    
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
                'clan_acc': result[6]  # Will be None if family is not in a clan
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
    
    def create_network(self, clan_acc, clan_families, matches):
        """
        Create network graph from matches.
        
        Args:
            clan_acc: Clan accession
            clan_families: DataFrame of clan families
            matches: DataFrame of foldseek matches
            
        Returns:
            NetworkX graph
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
        
        # Add edges and external family nodes
        for _, match in matches.iterrows():
            query_pfam = match['Pfam_Accession']
            target_pfam = match['Target_Pfam']
            
            if pd.isna(target_pfam) or pd.isna(query_pfam) or query_pfam == target_pfam:
                continue
            
            # Add query node if not already present
            if query_pfam not in G:
                query_info = self.get_family_info(query_pfam)
                if query_info:
                    G.add_node(query_pfam,
                              pfam_id=query_info['pfamA_id'],
                              model_length=query_info['model_length'],
                              type=query_info['type'],
                              in_clan=query_pfam in clan_pfam_accs,
                              clan_acc=query_info['clan_acc'],
                              description=query_info['description'])
                else:
                    # If we can't get info, at least add the node with accession as label
                    G.add_node(query_pfam,
                              pfam_id=query_pfam,
                              model_length=100,
                              type='Domain',  # Default to Domain
                              in_clan=query_pfam in clan_pfam_accs,
                              clan_acc=None,
                              description='Unknown')
            
            # Add target node if not already present
            if target_pfam not in G:
                target_info = self.get_family_info(target_pfam)
                if target_info:
                    G.add_node(target_pfam,
                              pfam_id=target_info['pfamA_id'],
                              model_length=target_info['model_length'],
                              type=target_info['type'],
                              in_clan=target_pfam in clan_pfam_accs,
                              clan_acc=target_info['clan_acc'],
                              description=target_info['description'])
                else:
                    # If we can't get info, at least add the node with accession as label
                    G.add_node(target_pfam,
                              pfam_id=target_pfam,
                              model_length=100,
                              type='Domain',  # Default to Domain
                              in_clan=target_pfam in clan_pfam_accs,
                              clan_acc=None,
                              description='Unknown')
            
            # Add edge
            evalue = match['E-value']
            evalue_cat = self.categorize_evalue(evalue)
            
            # If edge exists, keep the better (lower) E-value
            if G.has_edge(query_pfam, target_pfam):
                existing_evalue = G[query_pfam][target_pfam]['evalue']
                if evalue < existing_evalue:
                    G[query_pfam][target_pfam]['evalue'] = evalue
                    G[query_pfam][target_pfam]['evalue_cat'] = evalue_cat
            else:
                G.add_edge(query_pfam, target_pfam, 
                          evalue=evalue,
                          evalue_cat=evalue_cat)
        
        print(f"\nNetwork statistics:")
        print(f"  Nodes: {G.number_of_nodes()}")
        print(f"  Edges: {G.number_of_edges()}")
        print(f"  Target clan families: {sum(1 for n, d in G.nodes(data=True) if d.get('in_clan', False))}")
        
        # Count families by clan status
        no_clan_count = sum(1 for n, d in G.nodes(data=True) 
                           if not d.get('in_clan', False) and d.get('clan_acc') is None)
        different_clan_count = sum(1 for n, d in G.nodes(data=True) 
                                  if not d.get('in_clan', False) and d.get('clan_acc') is not None)
        
        print(f"  Families with no clan: {no_clan_count}")
        print(f"  Families in different clan: {different_clan_count}")
        
        return G
    
    def visualize_network(self, G, clan_acc, evalue_threshold, output_file=None):
        """
        Visualize the network.
        
        Args:
            G: NetworkX graph
            clan_acc: Clan accession
            evalue_threshold: E-value threshold used for filtering
            output_file: Output file path (optional)
        """
        # Get connected components info
        components = list(nx.connected_components(G))
        print(f"\nNetwork has {len(components)} connected component(s)")
        
        # Weight edges by E-value for layout
        for u, v, data in G.edges(data=True):
            data['weight'] = 1.0 / (data['evalue'] + 1e-100)
        
        # Try to use graphviz neato layout (better for small graphs, no overlap issue)
        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog='neato')
            print("Using graphviz neato layout")
        except:
            try:
                # Fallback to pydot
                pos = nx.nx_pydot.graphviz_layout(G, prog='neato')
                print("Using pydot neato layout")
            except:
                # Final fallback to spring layout
                print("Graphviz not available, using spring layout")
                pos = nx.spring_layout(G, weight='weight', k=0.5, iterations=100, seed=42)
        
        # Create figure with extra space on right for legend
        fig = plt.figure(figsize=(22, 16))  # Wider to accommodate legend
        ax = fig.add_axes([0.05, 0.05, 0.8, 0.87])  # Left, bottom, width, height as fractions
        
        # Categorize nodes by clan membership
        # Green: families in the target clan
        # Grey: families with no clan
        # Red: families with a different clan
        target_clan_nodes = []
        no_clan_nodes = []
        different_clan_nodes = []
        
        for n, d in G.nodes(data=True):
            node_clan = d.get('clan_acc', None)
            if d.get('in_clan', False):
                # This is a family in the target clan
                target_clan_nodes.append(n)
            elif node_clan is None:
                # This family has no clan
                no_clan_nodes.append(n)
            else:
                # This family has a different clan
                different_clan_nodes.append(n)
        
        # Calculate node sizes based on model_length
        # Scale to reasonable display sizes - will be applied to ALL node types
        model_lengths = [G.nodes[n].get('model_length', 100) for n in G.nodes()]
        max_length = max(model_lengths) if model_lengths else 1000
        min_length = min(model_lengths) if model_lengths else 50

        def scale_size(length):
            # Scale between 300 and 3000 for node size
            # This applies to ALL marker shapes (circles, squares, triangles, etc.)
            if max_length == min_length:
                return 1000
            return 300 + (length - min_length) / (max_length - min_length) * 2700

        # Draw edges with different widths based on E-value category
        edge_widths = []
        edge_colors = []
        
        for u, v, data in G.edges(data=True):
            cat = data.get('evalue_cat', 2)
            if cat == 0:
                edge_widths.append(3.0)
                edge_colors.append('darkred')
            elif cat == 1:
                edge_widths.append(2.0)
                edge_colors.append('red')
            else:
                edge_widths.append(1.0)
                edge_colors.append('lightcoral')
        
        nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color=edge_colors, alpha=0.6)

        # Draw nodes by type and clan membership for different shapes
        # Group nodes by type and clan status
        from collections import defaultdict
        nodes_by_type_clan = defaultdict(list)

        for n in G.nodes():
            node_data = G.nodes[n]
            node_type = node_data.get('type', 'Domain')
            node_clan = node_data.get('clan_acc', None)
            in_clan = node_data.get('in_clan', False)

            # Determine clan category
            if in_clan:
                clan_cat = 'target'
            elif node_clan is None:
                clan_cat = 'no_clan'
            else:
                clan_cat = 'different'

            nodes_by_type_clan[(node_type, clan_cat)].append(n)

        # Draw each type/clan combination with appropriate shape and color
        for (node_type, clan_cat), nodes in nodes_by_type_clan.items():
            if not nodes:
                continue

            # Get shape for this type
            marker = self.type_to_shape.get(node_type, 'o')  # Default to circle

            # Get color for clan category
            if clan_cat == 'target':
                color = 'green'
                edge_color = 'darkgreen'
                linewidth = 2
                alpha = 0.8
            elif clan_cat == 'no_clan':
                color = 'lightgrey'
                edge_color = 'grey'
                linewidth = 1
                alpha = 0.6
            else:  # different clan
                color = 'red'
                edge_color = 'darkred'
                linewidth = 2
                alpha = 0.8

            # Get sizes proportional to model_length for this type/clan combination
            # Note: All node types (circles, squares, triangles, etc.) are scaled by model_length
            sizes = [scale_size(G.nodes[n].get('model_length', 100)) for n in nodes]

            # Draw nodes with this type/color combination
            # 's' marker renders as true squares (not rectangles) in matplotlib
            nx.draw_networkx_nodes(G, pos, nodelist=nodes,
                                  node_size=sizes,
                                  node_shape=marker,
                                  node_color=color,
                                  alpha=alpha,
                                  edgecolors=edge_color,
                                  linewidths=linewidth)
        
        # Draw labels with ID and accession on separate lines
        labels = {n: f"{G.nodes[n].get('pfam_id', n)}\n{n}" for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight='bold')
        
        # Add legend with dynamic E-value ranges based on actual threshold
        # Calculate the category boundaries
        cat0_threshold = 1e-10
        cat1_threshold = 1e-5
        cat2_threshold = evalue_threshold
        
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', label='Target clan families',
                      markerfacecolor='green', markersize=12, markeredgecolor='darkgreen', markeredgewidth=2),
            plt.Line2D([0], [0], marker='o', color='w', label='Families with no clan',
                      markerfacecolor='lightgrey', markersize=12, markeredgecolor='grey'),
            plt.Line2D([0], [0], marker='o', color='w', label='Families in different clan',
                      markerfacecolor='red', markersize=12, markeredgecolor='darkred', markeredgewidth=2),
            plt.Line2D([0], [0], color='darkred', linewidth=3, label=f'E-value < 1e-10'),
            plt.Line2D([0], [0], color='red', linewidth=2, label=f'E-value 1e-10 to 1e-5'),
            plt.Line2D([0], [0], color='lightcoral', linewidth=1, label=f'E-value 1e-5 to {evalue_threshold:.0e}'),
            # Add shape legend
            plt.Line2D([0], [0], marker='s', color='w', label='Domain (square)',
                      markerfacecolor='grey', markersize=10, markeredgecolor='black'),
            plt.Line2D([0], [0], marker='o', color='w', label='Family (circle)',
                      markerfacecolor='grey', markersize=10, markeredgecolor='black'),
            plt.Line2D([0], [0], marker='^', color='w', label='Repeat (triangle)',
                      markerfacecolor='grey', markersize=10, markeredgecolor='black'),
            plt.Line2D([0], [0], marker='D', color='w', label='Coiled-coil (diamond)',
                      markerfacecolor='grey', markersize=10, markeredgecolor='black'),
            plt.Line2D([0], [0], marker='8', color='w', label='Motif (octagon)',
                      markerfacecolor='grey', markersize=10, markeredgecolor='black'),
            plt.Line2D([0], [0], marker='h', color='w', label='Disordered (hexagon)',
                      markerfacecolor='grey', markersize=10, markeredgecolor='black'),
        ]
        
        # Place legend in the right margin (using figure coordinates)
        fig.legend(handles=legend_elements, loc='center right', fontsize=10, 
                  bbox_to_anchor=(0.97, 0.5), frameon=True)
        
        # Add title with proper spacing
        fig.suptitle(f'Structural relationships for clan {clan_acc}\n'
                    f'Node size proportional to model length', 
                    fontsize=14, fontweight='bold', y=0.96)
        plt.axis('off')
        
        # Don't use tight_layout as it will override our careful positioning
        # plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300)
            print(f"\nSaved visualization to {output_file}")
        
        plt.show()
    
    def export_interactive_html(self, G, clan_acc, evalue_threshold, output_file):
        """
        Export network as interactive HTML using vis.js.
        
        Args:
            G: NetworkX graph
            clan_acc: Clan accession
            evalue_threshold: E-value threshold used
            output_file: Output HTML file path
        """
        # Prepare nodes data
        nodes = []
        node_titles = {}  # Store titles separately for custom tooltips
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
            size = 20 + (model_length / 10)  # Scale size

            # Get shape based on type
            node_type = data.get('type', 'Domain')
            shape = self.type_to_vis_shape.get(node_type, 'dot')

            # For ellipse shapes, we need special handling to prevent text from being inside
            # All other shapes (dot, square, triangle, hexagon, star) render text outside by default
            if shape == 'ellipse':
                # Position label below the ellipse, not inside it
                font_config = {
                    'size': 18,
                    'color': '#000000',
                    'face': 'Arial',
                    'bold': {'size': 20},
                    'vadjust': 25  # Move label below the ellipse
                }
                # Make ellipse a fixed size (not resizing to fit text)
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
                    'widthConstraint': {'minimum': size * 1.5, 'maximum': size * 1.5},  # Oval width
                    'heightConstraint': {'minimum': size, 'maximum': size}  # Oval height
                }
            else:
                # For other shapes (dot, square, triangle, hexagon, star), use default label positioning
                # These shapes automatically position labels outside and scale properly
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

            # Store title separately (not in node data)
            node_titles[node] = (f"<div style='font-family: Arial; padding: 5px;'>"
                        f"<b style='font-size: 16px;'>{data.get('pfam_id', node)}</b><br/>"
                        f"<span style='color: #666; font-size: 14px;'>{node}</span><br/>"
                        f"<b style='font-size: 13px;'>Type:</b> <span style='font-size: 13px;'>{node_type}</span><br/>"
                        f"<b style='font-size: 13px;'>Length:</b> <span style='font-size: 13px;'>{model_length}</span><br/>"
                        f"<b style='font-size: 13px;'>Clan:</b> <span style='font-size: 13px;'>{node_clan if node_clan else 'None'}</span><br/>"
                        f"<span style='font-size: 12px; color: #555;'>{data.get('description', '')}</span>"
                        f"</div>")

            nodes.append(node_config)
        
        # Prepare edges data
        edges = []
        edge_titles = {}  # Store titles separately for custom tooltips
        for u, v, data in G.edges(data=True):
            evalue = data.get('evalue', 1.0)
            evalue_cat = data.get('evalue_cat', 2)
            
            # Edge color based on E-value category
            if evalue_cat == 0:
                color = '#B71C1C'  # Dark red
                width = 3
            elif evalue_cat == 1:
                color = '#F44336'  # Red
                width = 2
            else:
                color = '#FFCDD2'  # Light coral
                width = 1
            
            # Create unique edge ID
            edge_id = f"{u}-{v}"
            
            # Store title separately (not in edge data)
            edge_titles[edge_id] = (f"<div style='font-family: Arial; padding: 5px;'>"
                        f"<b style='font-size: 15px;'>E-value: {evalue:.2e}</b><br/>"
                        f"<b style='font-size: 13px;'>Between:</b><br/>"
                        f"<span style='color: #666; font-size: 13px;'>{u}</span><br/>"
                        f"<span style='color: #666; font-size: 13px;'>{v}</span>"
                        f"</div>")
            
            edges.append({
                'id': edge_id,
                'from': u,
                'to': v,
                'value': width,
                'color': {'color': color, 'highlight': '#000000'},
                'smooth': {'type': 'continuous'}
            })
        
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
        #editModeBtn:hover {{
            opacity: 0.9;
            transform: scale(1.02);
        }}
        #editModeBtn:active {{
            transform: scale(0.98);
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
        Interactive Network Visualization<br/>
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
        <strong>Legend:</strong>
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
        <div class="legend-item">
            <span class="legend-line" style="background-color: #B71C1C; height: 3px;"></span>
            E-value &lt; 1e-10
        </div>
        <div class="legend-item">
            <span class="legend-line" style="background-color: #F44336; height: 2px;"></span>
            E-value 1e-10 to 1e-5
        </div>
        <div class="legend-item">
            <span class="legend-line" style="background-color: #FFCDD2; height: 1px;"></span>
            E-value 1e-5 to {evalue_threshold:.0e}
        </div>
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
        
        // Store titles separately (not part of vis.js data to prevent built-in tooltips)
        var nodeTitles = {json.dumps(node_titles)};
        var edgeTitles = {json.dumps(edge_titles)};
        
        // Find isolated nodes (nodes with no connections) and position them in a grid
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
        var gridStartX = -2500;  // Position grid far on the left to avoid overlap
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
                fixed: {{ x: true, y: true }}  // Fix isolated nodes during physics
            }});
        }});
        
        // Create a network
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
        
        // Disable physics after initial stabilization to prevent nodes from moving when dragging
        network.on('stabilizationIterationsDone', function() {{
            // Unfix isolated nodes so they can be dragged
            isolatedNodes.forEach(function(nodeId) {{
                nodes.update({{
                    id: nodeId,
                    fixed: false
                }});
            }});
            
            network.setOptions({{physics: false}});
            console.log('Physics disabled - nodes will now stay in place when dragging');
            console.log('Isolated nodes unfixed - can now be repositioned');
        }});
        
        // Also disable after a timeout as backup
        setTimeout(function() {{
            // Unfix isolated nodes
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
        
        // Handle node clicks - open InterPro page (only when not in editing mode)
        network.on('click', function(params) {{
            if (params.nodes.length > 0 && !editingMode) {{
                var nodeId = params.nodes[0];
                var url = 'https://www.ebi.ac.uk/interpro/entry/pfam/' + nodeId + '/';
                window.open(url, '_blank');
                console.log('Opening InterPro page: ' + url);
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

                // Scale vadjust for ellipse shapes
                if (original.vadjust) {{
                    update.font.vadjust = Math.round(original.vadjust * scale);
                }}

                // Scale width/height constraints for ellipse shapes
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

        // Add some info on console
        console.log('Network created with ' + nodes.length + ' nodes and ' + edges.length + ' edges');
        console.log('Isolated nodes arranged in grid on left side: ' + isolatedNodes.length);
        console.log('Toggle editing mode to prevent accidental page opens while repositioning nodes.');
        console.log('Drag nodes to reposition them. Hover over nodes and edges to see details.');
        console.log('Use the size slider to adjust node and text scale.');
    </script>
</body>
</html>'''
        
        # Write to file
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        print(f"\nSaved interactive HTML to {output_file}")
        print(f"Open {output_file} in a web browser to interact with the network")
    
    def analyze_clan(self, clan_acc, evalue_threshold=1e-3, output_file=None):
        """
        Complete analysis pipeline for a clan.
        
        Args:
            clan_acc: Clan accession (e.g., 'CL0004')
            evalue_threshold: Maximum E-value to include
            output_file: Output file path for visualization
        """
        print(f"\n{'='*60}")
        print(f"Analyzing clan {clan_acc}")
        print(f"{'='*60}")
        
        # Get clan families
        clan_families = self.get_clan_families(clan_acc)
        
        if len(clan_families) == 0:
            print(f"No families found for clan {clan_acc}")
            return
        
        # Parse foldseek results
        matches = self.parse_foldseek_results(clan_families, evalue_threshold)
        
        if len(matches) == 0:
            print(f"No matches found for clan {clan_acc}")
            return
        
        # Create network
        G = self.create_network(clan_acc, clan_families, matches)

        # Always create interactive HTML version based on clan accession
        html_file = f'{clan_acc}_interactive.html'
        self.export_interactive_html(G, clan_acc, evalue_threshold, html_file)

        # Optionally create static visualization if output file is specified
        if output_file:
            self.visualize_network(G, clan_acc, evalue_threshold, output_file)
        else:
            # Skip static visualization, show the network briefly
            print("\nSkipping static visualization (no --output specified)")
            print(f"View the interactive HTML: {html_file}")
        
        # Print some interesting statistics
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
        description='Visualize clan relationships using foldseek structural comparison data'
    )
    parser.add_argument('clan_acc', help='Clan accession (e.g., CL0004)')
    parser.add_argument('--foldseek', 
                       default='/nfs/production/agb/interpro/users/typhaine/interpro-pfam-curation-tools/pfam_add_clan_search/results/foldseek_all.out',
                       help='Path to foldseek results file')
    parser.add_argument('--config', 
                       default='~/.my.cnf',
                       help='Path to MySQL config file')
    parser.add_argument('--evalue', type=float, default=1e-3,
                       help='Maximum E-value threshold (default: 1e-3)')
    parser.add_argument('--output', 
                       help='Output file path for visualization')
    
    args = parser.parse_args()
    
    # Create visualizer
    viz = ClanNetworkVisualizer(args.foldseek, args.config)
    
    try:
        # Connect to database
        viz.connect_to_database()
        
        # Analyze clan
        viz.analyze_clan(args.clan_acc, args.evalue, args.output)
        
    finally:
        viz.close()


if __name__ == '__main__':
    main()
