"""
Module for visualizing CQ filtered results across the genome.
Generates both static (Matplotlib) and interactive (Plotly) distribution plots
to visualize specific genomic regions identified through CQ analysis.
"""
import pandas as pd
import matplotlib.pyplot as plt
import logging
import os

# Try to import plotly for interactive plotting
try:
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# Initialize module-level logger
logger = logging.getLogger(__name__)

class CQPlotter:
    """
    Handles the generation of genomic distribution plots for CQ results.
    """
    def __init__(self, cq_paths, length_path, output_path, pdf_maker=False, html_maker=True):
        """
        Initialize the plotter with data paths and configuration.
        
        Args:
            cq_paths (str): Path to the filtered CQ result TSV file.
            length_path (str): Path to the chromosome length TSV file.
            output_path (str): Base path for output (extension will be appended/changed).
            pdf_maker (bool): If True, generate a PDF version of the static plot.
            html_maker (bool): If True, generate an interactive HTML plot using Plotly.
        """
        self.cq_paths = cq_paths
        self.length_path = length_path
        self.output_path = output_path
        self.pdf_maker: bool = pdf_maker
        self.html_maker: bool = html_maker

    def plot_interactive_cq(self, cq_df, chromosomes, chrom_len_dict):
        """
        Generate an interactive HTML plot using Plotly.
        
        Allows users to hover over data points to see specific window coordinates,
        read counts, CPM values, and calculated CQ results.
        
        Args:
            cq_df (pd.DataFrame): Dataframe containing filtered CQ results.
            chromosomes (list): Ordered list of chromosome names.
            chrom_len_dict (dict): Dictionary mapping chromosome names to their total lengths.
        """
        if not PLOTLY_AVAILABLE:
            logger.warning("Plotly is not installed. Skipping interactive plot.")
            return

        html_output = os.path.splitext(self.output_path)[0] + ".html"
        logger.info(f"Generating interactive HTML plot: {html_output}")

        # Create interactive scatter plot
        fig = px.scatter(
            cq_df, 
            x="start", 
            y="chr", 
            color="cq_result",
            color_continuous_scale="Magma", # Consistent with static plot
            hover_data={
                "chr": True,
                "start": True,
                "end": True,
                "cq_result": ":.4f", 
                "f_num_reads": True,
                "m_num_reads": True,
                "F_CPM": ":.2f",
                "M_CPM": ":.2f"
            },
            title="Genomic Distribution of Filtered CQ Regions (Interactive)",
            labels={"cq_result": "CQ Value", "start": "Position (bp)", "chr": "Chromosome"},
            height=max(500, len(chromosomes) * 35) # Dynamic height based on track count
        )

        # Add horizontal backbone lines for each chromosome
        for chrom in chromosomes:
            if chrom in chrom_len_dict:
                fig.add_shape(
                    type="line",
                    x0=0, x1=chrom_len_dict[chrom],
                    y0=chrom, y1=chrom,
                    line=dict(color="lightgrey", width=1),
                    opacity=0.8,
                    layer="below"
                )

        # Refine layout
        fig.update_layout(
            plot_bgcolor='white',
            xaxis=dict(showgrid=True, gridcolor='lightgrey', title="Genomic Position (bp)"),
            yaxis=dict(
                showgrid=False, 
                categoryorder='array', 
                categoryarray=chromosomes 
            ),
            coloraxis_colorbar=dict(title="CQ Value")
        )

        # Write to file
        fig.write_html(html_output)
        logger.info(f"Interactive plot successfully saved.")

    def plot_cq(self):
        """
        Primary plotting method. Generates both static and interactive plots.
        Filters data to only show chromosomes present in the filtered result set.
        """
        # 1. Load CQ results
        logger.info(f"Reading filtered CQ results from: {self.cq_paths}")
        if not os.path.exists(self.cq_paths):
            logger.warning(f"File not found: {self.cq_paths}. Skipping plot.")
            return

        try:
            cq_df = pd.read_csv(self.cq_paths, sep='\t')
        except pd.errors.EmptyDataError:
            logger.warning("CQ result file is empty.")
            return
        except Exception as e:
            logger.error(f"Error reading CQ file: {e}")
            return

        if cq_df.empty:
            logger.warning("No data points found. Skipping plot.")
            return

        # Identify chromosomes with significant regions
        active_chroms = cq_df['chr'].unique()

        # 2. Load and filter chromosome lengths
        logger.info(f"Reading chromosome lengths from: {self.length_path}")
        try:
            chrom_len_df = pd.read_csv(self.length_path, sep='\t', names=['chr', 'length'])
            # Keep only active tracks
            chrom_len_df = chrom_len_df[chrom_len_df['chr'].isin(active_chroms)]
            
            # Map chromosomes to Y-axis indices (inverted order for visual consistency)
            chromosomes = chrom_len_df['chr'].tolist()[::-1]
            y_map = {chrom: i for i, chrom in enumerate(chromosomes)}
            chrom_len_dict = dict(zip(chrom_len_df['chr'], chrom_len_df['length']))
        except Exception as e:
            logger.error(f"Error processing chromosome lengths: {e}")
            return

        # 3. Static Plotting (Matplotlib)
        # Compact track spacing
        plt.figure(figsize=(15, max(5, len(y_map) * 0.3))) 
        ax = plt.gca()

        # Draw chromosome backbones
        logger.info("Drawing track backbones...")
        for chrom, y in y_map.items():
            if chrom in chrom_len_dict:
                ax.hlines(y, 0, chrom_len_dict[chrom], colors='lightgrey', linewidth=1, alpha=0.8, zorder=1)

        # Map data points to Y-coordinates
        cq_df['y_pos'] = cq_df['chr'].map(y_map)
        cq_df = cq_df.dropna(subset=['y_pos'])

        if not cq_df.empty:
            # Draw scatter plot with color gradient
            scatter = ax.scatter(
                cq_df['start'], 
                cq_df['y_pos'], 
                c=cq_df['cq_result'], 
                cmap='magma', # High contrast scientific colormap
                s=35, 
                edgecolor='none', 
                zorder=2,
                alpha=0.8
            )

            # Add Colorbar and Labels
            cbar = plt.colorbar(scatter)
            cbar.set_label('CQ Value (Lower = More Specific)', rotation=270, labelpad=15)

            ax.set_yticks(list(y_map.values()))
            ax.set_yticklabels(list(y_map.keys()))
            ax.set_xlabel('Genomic Position (bp)')
            ax.set_ylabel('Chromosomes')
            ax.set_title('Genomic Distribution of Filtered CQ Regions')

            plt.grid(axis='x', linestyle='--', alpha=0.3)
            plt.tight_layout()

            # Save static image
            static_output = self.output_path
            if self.pdf_maker:
                static_output = os.path.splitext(self.output_path)[0] + ".pdf"
            
            logger.info(f"Saving static plot to: {static_output}")
            plt.savefig(static_output, dpi=300)
            plt.close()

        # 4. Interactive Plotting (Plotly)
        if self.html_maker:
            self.plot_interactive_cq(cq_df, chromosomes, chrom_len_dict)


if __name__ == "__main__":
    # Test block for module verification
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    test_cq = "F_M_CQ.filter.tsv"
    test_length = "chromosome_length.tsv"
    test_output = "CQ_test_plot.png"

    plotter = CQPlotter(test_cq, test_length, test_output, pdf_maker=False, html_maker=True)
    plotter.plot_cq()
