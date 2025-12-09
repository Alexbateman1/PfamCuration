"""
Domain Visualization

Generates:
- Matplotlib-based domain architecture diagrams
- ChimeraX coloring commands
- PyMOL coloring scripts
- HTML reports
"""

from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import numpy as np

if TYPE_CHECKING:
    from .abc_predictor import ABCPrediction, Domain, NDRRegion, Residue

# Domain color palette (10 distinct colors)
DOMAIN_COLORS = [
    "#4A79A7",  # Blue
    "#F28E2C",  # Orange
    "#E15759",  # Red
    "#76B7B2",  # Teal
    "#59A14F",  # Green
    "#EDC949",  # Yellow
    "#AF7AA1",  # Purple
    "#FF9DA7",  # Pink
    "#9C755F",  # Brown
    "#BAB0AB",  # Gray
]

NDR_COLOR = "#CCCCCC"  # Light gray for NDRs


class DomainVisualizer:
    """
    Visualizes domain predictions.

    Generates:
    - Architecture diagrams (matplotlib)
    - ChimeraX commands
    - PyMOL scripts
    """

    def __init__(
        self,
        colors: List[str] = None,
        ndr_color: str = NDR_COLOR,
    ):
        self.colors = colors or DOMAIN_COLORS
        self.ndr_color = ndr_color

    def visualize(
        self,
        prediction: "ABCPrediction",
        output_path: Optional[str] = None,
        residues: Optional[List["Residue"]] = None,
    ) -> str:
        """
        Generate visualization outputs.

        Parameters:
        -----------
        prediction : ABCPrediction
            Prediction result
        output_path : str, optional
            Base path for output files (without extension)
        residues : List[Residue], optional
            Residue data for pLDDT visualization

        Returns:
        --------
        str: ChimeraX commands
        """
        chimerax_commands = self.generate_chimerax_commands(prediction)

        if output_path:
            path = Path(output_path)

            # Save ChimeraX commands
            chimerax_file = path.with_suffix(".cxc")
            with open(chimerax_file, "w") as f:
                f.write(chimerax_commands)

            # Save PyMOL script
            pymol_file = path.with_suffix(".pml")
            pymol_script = self.generate_pymol_script(prediction)
            with open(pymol_file, "w") as f:
                f.write(pymol_script)

            # Generate architecture diagram
            png_file = path.with_suffix(".png")
            self.create_architecture_diagram(prediction, str(png_file), residues)

            # Generate HTML report
            html_file = path.with_suffix(".html")
            self.generate_html_report(prediction, str(html_file))

        return chimerax_commands

    def generate_chimerax_commands(
        self,
        prediction: "ABCPrediction",
    ) -> str:
        """
        Generate ChimeraX coloring commands.

        Returns:
        --------
        str: ChimeraX commands
        """
        lines = [
            f"# ABC Domain Prediction for {prediction.uniprot_acc}",
            f"# {len(prediction.domains)} domains, {len(prediction.ndr_regions)} NDR regions",
            "",
            "# Color by domain",
            "color all white",
        ]

        # Color each domain
        for i, domain in enumerate(prediction.domains):
            color = self.colors[i % len(self.colors)]
            # ChimeraX: color each segment separately for reliability
            lines.append(f"# Domain {domain.domain_id}: {domain.to_chopping_string()}")
            for start, end in domain.segments:
                lines.append(f"color :{start}-{end} {color}")

        # Color NDRs
        lines.append("")
        lines.append("# NDR regions")
        for ndr in prediction.ndr_regions:
            lines.append(f"color :{ndr.start}-{ndr.end} {self.ndr_color}")

        lines.extend([
            "",
            "# Show cartoon representation",
            "hide atoms",
            "show cartoons",
            "",
            "# Optional: show surface",
            "# surface",
            "# transparency 50",
        ])

        return "\n".join(lines)

    def generate_pymol_script(
        self,
        prediction: "ABCPrediction",
    ) -> str:
        """
        Generate PyMOL coloring script.

        Returns:
        --------
        str: PyMOL commands
        """
        lines = [
            f"# ABC Domain Prediction for {prediction.uniprot_acc}",
            f"# {len(prediction.domains)} domains, {len(prediction.ndr_regions)} NDR regions",
            "",
            "hide everything",
            "show cartoon",
            "color white, all",
            "",
        ]

        # Color each domain
        for i, domain in enumerate(prediction.domains):
            color = self.colors[i % len(self.colors)]
            color_name = f"domain_{domain.domain_id}"

            # Convert hex to RGB for PyMOL
            r, g, b = self._hex_to_rgb(color)
            lines.append(f"set_color {color_name}, [{r:.3f}, {g:.3f}, {b:.3f}]")

            sel = self._segments_to_pymol_selection(domain.segments)
            lines.append(f"color {color_name}, {sel}")
            lines.append("")

        # Color NDRs
        r, g, b = self._hex_to_rgb(self.ndr_color)
        lines.append(f"set_color ndr_color, [{r:.3f}, {g:.3f}, {b:.3f}]")
        for ndr in prediction.ndr_regions:
            lines.append(f"color ndr_color, resi {ndr.start}-{ndr.end}")

        lines.extend([
            "",
            "# Show surface (optional)",
            "# show surface",
            "# set transparency, 0.5",
            "",
            "# Orient view",
            "orient",
            "zoom",
        ])

        return "\n".join(lines)

    def create_architecture_diagram(
        self,
        prediction: "ABCPrediction",
        output_path: str,
        residues: Optional[List["Residue"]] = None,
    ):
        """
        Create a matplotlib domain architecture diagram.

        Shows:
        - Protein backbone as a line
        - Domains as colored boxes
        - NDRs as gray boxes
        - pLDDT trace (if residues provided)
        """
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as mpatches
        except ImportError:
            return  # Skip if matplotlib not available

        # Figure setup
        n_rows = 2 if residues else 1
        fig, axes = plt.subplots(n_rows, 1, figsize=(14, 3 * n_rows))

        if n_rows == 1:
            axes = [axes]

        ax_arch = axes[0]

        # Draw protein backbone
        seq_len = prediction.sequence_length
        ax_arch.axhline(y=0.5, color="gray", linewidth=2, zorder=1)

        # Draw domains
        for i, domain in enumerate(prediction.domains):
            color = self.colors[i % len(self.colors)]
            for start, end in domain.segments:
                rect = mpatches.FancyBboxPatch(
                    (start, 0.3),
                    end - start + 1,
                    0.4,
                    boxstyle="round,pad=0.02,rounding_size=0.1",
                    facecolor=color,
                    edgecolor="black",
                    linewidth=1,
                    zorder=2,
                )
                ax_arch.add_patch(rect)

                # Add domain label
                mid = (start + end) / 2
                ax_arch.text(mid, 0.5, f"D{domain.domain_id}", ha="center", va="center",
                           fontsize=9, fontweight="bold", color="white", zorder=3)

        # Draw NDRs
        for ndr in prediction.ndr_regions:
            rect = mpatches.Rectangle(
                (ndr.start, 0.35),
                ndr.end - ndr.start + 1,
                0.3,
                facecolor=self.ndr_color,
                edgecolor="gray",
                linewidth=0.5,
                alpha=0.7,
                zorder=1,
            )
            ax_arch.add_patch(rect)

        ax_arch.set_xlim(0, seq_len + 1)
        ax_arch.set_ylim(0, 1)
        ax_arch.set_xlabel("Residue Position")
        ax_arch.set_title(f"ABC Domain Prediction: {prediction.uniprot_acc}")
        ax_arch.set_yticks([])

        # Draw pLDDT trace if available
        if residues and len(axes) > 1:
            ax_plddt = axes[1]
            positions = [r.resnum for r in residues]
            plddt_values = [r.plddt for r in residues]

            ax_plddt.fill_between(positions, plddt_values, alpha=0.3, color="blue")
            ax_plddt.plot(positions, plddt_values, color="blue", linewidth=1)

            ax_plddt.axhline(y=70, color="orange", linestyle="--", alpha=0.7, label="pLDDT=70")
            ax_plddt.set_xlim(0, seq_len + 1)
            ax_plddt.set_ylim(0, 100)
            ax_plddt.set_xlabel("Residue Position")
            ax_plddt.set_ylabel("pLDDT")
            ax_plddt.legend(loc="lower right")

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close()

    def generate_html_report(
        self,
        prediction: "ABCPrediction",
        output_path: str,
    ):
        """
        Generate an HTML report with all prediction details.
        """
        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>ABC Domain Prediction: {prediction.uniprot_acc}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #333; }}
        .domain-box {{
            display: inline-block;
            padding: 5px 10px;
            margin: 2px;
            border-radius: 4px;
            color: white;
            font-weight: bold;
        }}
        .ndr-box {{
            display: inline-block;
            padding: 5px 10px;
            margin: 2px;
            border-radius: 4px;
            background-color: #CCCCCC;
            color: #333;
        }}
        table {{ border-collapse: collapse; margin: 10px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #4A79A7; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .params {{ background-color: #f9f9f9; padding: 10px; border-radius: 4px; }}
    </style>
</head>
<body>
    <h1>ABC Domain Prediction: {prediction.uniprot_acc}</h1>

    <h2>Summary</h2>
    <p>Sequence length: {prediction.sequence_length} residues</p>
    <p>Domains found: {len(prediction.domains)}</p>
    <p>NDR regions: {len(prediction.ndr_regions)}</p>

    <h2>Domain Architecture</h2>
    <div style="margin: 10px 0;">
"""
        # Add domain boxes
        for i, domain in enumerate(prediction.domains):
            color = self.colors[i % len(self.colors)]
            html += f'        <span class="domain-box" style="background-color:{color}">D{domain.domain_id}: {domain.to_chopping_string()}</span>\n'

        for ndr in prediction.ndr_regions:
            html += f'        <span class="ndr-box">NDR: {ndr.start}-{ndr.end}</span>\n'

        html += """    </div>

    <h2>Domain Details</h2>
    <table>
        <tr>
            <th>Domain</th>
            <th>Segments</th>
            <th>Size</th>
            <th>Discontinuous</th>
            <th>Avg pLDDT</th>
            <th>Rg (Å)</th>
            <th>Contact Ratio</th>
            <th>Quality Score</th>
        </tr>
"""
        for domain in prediction.domains:
            m = domain.quality_metrics
            html += f"""        <tr>
            <td>D{domain.domain_id}</td>
            <td>{domain.to_chopping_string()}</td>
            <td>{domain.size}</td>
            <td>{'Yes' if domain.is_discontinuous else 'No'}</td>
            <td>{m.get('avg_plddt', 'N/A'):.1f}</td>
            <td>{m.get('radius_of_gyration', 'N/A'):.1f}</td>
            <td>{m.get('contact_density_ratio', 'N/A'):.2f}</td>
            <td>{m.get('quality_score', 'N/A'):.1f}</td>
        </tr>
"""
        html += """    </table>

    <h2>NDR Regions</h2>
"""
        if prediction.ndr_regions:
            html += """    <table>
        <tr>
            <th>Region</th>
            <th>Size</th>
            <th>Avg pLDDT</th>
            <th>Type</th>
        </tr>
"""
            for ndr in prediction.ndr_regions:
                html += f"""        <tr>
            <td>{ndr.start}-{ndr.end}</td>
            <td>{ndr.size}</td>
            <td>{ndr.avg_plddt:.1f}</td>
            <td>{ndr.reason}</td>
        </tr>
"""
            html += "    </table>\n"
        else:
            html += "    <p>No NDR regions identified.</p>\n"

        html += f"""
    <h2>Parameters Used</h2>
    <div class="params">
        <pre>{self._format_params(prediction.parameters)}</pre>
    </div>

    <h2>ChimeraX Commands</h2>
    <pre style="background-color: #f5f5f5; padding: 10px; border-radius: 4px;">
{self.generate_chimerax_commands(prediction)}
    </pre>

</body>
</html>
"""
        with open(output_path, "w") as f:
            f.write(html)

    def _segments_to_chimerax_selection(self, segments: List[Tuple[int, int]]) -> str:
        """Convert segments to ChimeraX selection string."""
        parts = []
        for start, end in segments:
            parts.append(f":{start}-{end}")
        return ",".join(parts)

    def _segments_to_pymol_selection(self, segments: List[Tuple[int, int]]) -> str:
        """Convert segments to PyMOL selection string."""
        parts = []
        for start, end in segments:
            parts.append(f"resi {start}-{end}")
        return " or ".join(parts)

    def _hex_to_rgb(self, hex_color: str) -> Tuple[float, float, float]:
        """Convert hex color to RGB (0-1 range)."""
        hex_color = hex_color.lstrip("#")
        r = int(hex_color[0:2], 16) / 255
        g = int(hex_color[2:4], 16) / 255
        b = int(hex_color[4:6], 16) / 255
        return r, g, b

    def _format_params(self, params: Dict) -> str:
        """Format parameters for display."""
        lines = []
        for key, value in params.items():
            lines.append(f"{key}: {value}")
        return "\n".join(lines)


class ComparisonVisualizer:
    """
    Visualize comparison between ABC predictions and other methods.
    """

    def __init__(self, colors: List[str] = None):
        self.colors = colors or DOMAIN_COLORS

    def compare_methods(
        self,
        predictions: Dict[str, "ABCPrediction"],
        output_path: str,
    ):
        """
        Create comparison visualization showing multiple method predictions.

        Parameters:
        -----------
        predictions : Dict[str, ABCPrediction]
            Mapping of method name to prediction
        output_path : str
            Path for output PNG
        """
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as mpatches
        except ImportError:
            return

        n_methods = len(predictions)
        fig, axes = plt.subplots(n_methods, 1, figsize=(14, 2 * n_methods), sharex=True)

        if n_methods == 1:
            axes = [axes]

        # Get sequence length from first prediction
        seq_len = list(predictions.values())[0].sequence_length

        for ax, (method, pred) in zip(axes, predictions.items()):
            # Draw backbone
            ax.axhline(y=0.5, color="gray", linewidth=2, zorder=1)

            # Draw domains
            for i, domain in enumerate(pred.domains):
                color = self.colors[i % len(self.colors)]
                for start, end in domain.segments:
                    rect = mpatches.FancyBboxPatch(
                        (start, 0.3), end - start + 1, 0.4,
                        boxstyle="round,pad=0.02,rounding_size=0.1",
                        facecolor=color, edgecolor="black", linewidth=1, zorder=2
                    )
                    ax.add_patch(rect)

            # Draw NDRs
            for ndr in pred.ndr_regions:
                rect = mpatches.Rectangle(
                    (ndr.start, 0.35), ndr.end - ndr.start + 1, 0.3,
                    facecolor="#CCCCCC", edgecolor="gray", linewidth=0.5,
                    alpha=0.7, zorder=1
                )
                ax.add_patch(rect)

            ax.set_xlim(0, seq_len + 1)
            ax.set_ylim(0, 1)
            ax.set_ylabel(method)
            ax.set_yticks([])

        axes[-1].set_xlabel("Residue Position")
        plt.suptitle("Method Comparison")
        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close()
