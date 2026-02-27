"""
Main entry point for CQ_tools.
Provides CLI commands for sequence alignment and CQ value calculation.
"""
import typer
import logging
import sys
import os
from CQAligent.align import run_bwa_samtools

from CQAnalysis.get_paths import get_paths
from CQAnalysis.windows_make import MakeWindows
from CQAnalysis.CQ_runner import CQRunner
from CQAnalysis.coverage_calculate import CoverageCalculate
from CQPlot.CQ_plot import CQPlotter

app = typer.Typer(
    help="CQ_tools: A tool for calculating CQ values from sequencing data.",
    rich_markup_mode="rich",
)

def setup_logging(log_file: str = None):
    """
    Configure unified logging for the entire project.
    Outputs to both console and a log file if provided.
    """
    # Create root logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # Format: Time - Name (Module) - Level - Message
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler (optional)
    if log_file:
        # Ensure log directory exists
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir)
        
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

# Get logger for the main module
logger = logging.getLogger("CQ_tools.main")

@app.command()
def align(
    fasta: str = typer.Option(..., help="Path to the reference genome FASTA file."),
    pair_1: str = typer.Option(..., help="Path to the first pair of FASTQ reads."),
    pair_2: str = typer.Option(..., help="Path to the second pair of FASTQ reads."),
    output_dir: str = typer.Option(..., help="Directory where alignment results will be saved."),
    cpu: int = typer.Option(..., help="Number of CPU cores to use for alignment."),
):
    """
    Perform sequence alignment using BWA and process with Samtools.
    """
    # Setup simple console logging for align command
    setup_logging()
    logger.info("Starting alignment workflow...")
    run_bwa_samtools(fasta=fasta, sample_1=pair_1, sample_2=pair_2, output_dir=output_dir, cpu=cpu)
    logger.info("Alignment workflow completed successfully.")

@app.command()
def cq(
        fasta: str = typer.Option(..., help="Path to the reference genome FASTA file."),
        f_bam: str = typer.Option(..., help="Path to the sorted BAM file for the female sample."),
        m_bam: str = typer.Option(..., help="Path to the sorted BAM file for the male sample."),
        output_dir: str = typer.Option(..., help="Directory where CQ results will be saved."),
        cq_value: float = typer.Option(0.3, help="CQ threshold for filtering (default: 0.3)."),
        threshold: int = typer.Option(0, help="Minimum M_CPM reads threshold for filtering (default: 0)."),
        parallel: int = typer.Option(2, help="Parallel processes: 1 (serial) or 2 (parallel). Default: 2."),
):
    """
    Calculate CQ values and filter windows based on coverage and CPM.
    """
    # 获取所有必要的文件路径
    paths = get_paths(f_bam, m_bam, fasta, output_dir)
    
    # 初始化统一日志系统
    setup_logging(paths['log_path'])
    logger.info("="*50)
    logger.info("CQ Analysis Pipeline Started")
    logger.info(f"Female BAM: {f_bam}")
    logger.info(f"Male BAM: {m_bam}")
    logger.info("="*50)

    # 1. 建立基因组索引并切分窗口 (Window creation)
    logger.info("Step 1: Genomic window creation...")
    windows = MakeWindows(paths)
    windows.executor()

    # 2. 计算雌雄样本在各窗口的覆盖度 (Coverage calculation)
    logger.info("Step 2: Coverage calculation...")
    coverage_calculator = CoverageCalculate(paths, parallel)
    coverage_calculator.execute()

    # 3. 计算 CPM 标准化值并计算 CQ 结果 (CQ calculation)
    logger.info("Step 3: CQ calculation and filtering...")
    cq_runner = CQRunner(paths, cq_value, threshold)
    cq_runner.cq_runner()

    logger.info("="*50)
    logger.info(f"CQ Analysis Pipeline Finished. Results: {paths['filtered_cq']}")
    logger.info("="*50)

@app.command()
def plot(
    cq_result: str = typer.Option(..., help="Path to the filtered CQ result TSV file."),
    chrom_length: str = typer.Option(..., help="Path to the chromosome length TSV file."),
    output: str = typer.Option("CQ_distribution.png", help="Path to save the output plot."),
    pdf: bool = typer.Option(False, help="Whether to generate a PDF version of the static plot."),
    html: bool = typer.Option(True, help="Whether to generate an interactive HTML plot."),
):
    """
    Generate visualization plots (Static and/or Interactive) for CQ results.
    """
    setup_logging()
    logger.info("Starting Plotting workflow...")
    
    plotter = CQPlotter(
        cq_paths=cq_result,
        length_path=chrom_length,
        output_path=output,
        pdf_maker=pdf,
        html_maker=html
    )
    
    plotter.plot_cq()
    logger.info("Plotting workflow completed.")


if __name__ == "__main__":
    app()