"""
Module for genome alignment and processing.
Handles BWA indexing, alignment, and Samtools processing (BAM conversion, sorting, indexing).
"""
import os
import subprocess
import logging

# Initialize module-level logger
logger = logging.getLogger(__name__)

def index_bwa(reference_genome, prefix="A"):
    """
    Create an index for the reference genome using BWA.
    
    Args:
        reference_genome (str): Path to the reference FASTA file.
        prefix (str): Prefix for the index files.
    """
    idx_bwa = ["bwa", "index", "-p", prefix, reference_genome]
    try:
        subprocess.run(idx_bwa, check=True, text=True)
        logger.info("BWA index finished successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"BWA index error: {e.stderr}")
        raise

def align_bwa(reference_prefix, sample_1, sample_2, bwa_output, cpu):
    """
    Align paired-end reads to the reference genome using BWA MEM.
    
    Args:
        reference_prefix (str): Prefix of the BWA index.
        sample_1 (str): Path to FASTQ R1.
        sample_2 (str): Path to FASTQ R2.
        bwa_output (str): Path to the output SAM file.
        cpu (int): Number of threads.
    """
    run_BWA = ["bwa", "mem", "-t", str(cpu), "-o", bwa_output, reference_prefix, sample_1, sample_2]
    try:
        subprocess.run(run_BWA, check=True, text=True)
        logger.info(f"BWA alignment finished. Output: {bwa_output}")
    except subprocess.CalledProcessError as e:
        logger.error(f"BWA alignment error: {e.stderr}")
        raise

def run_bwa(reference_genome, sample_1, sample_2, bwa_output, cpu, prefix="A"):
    """
    Full BWA workflow: indexing followed by alignment.
    """
    index_bwa(reference_genome, prefix)
    align_bwa(prefix, sample_1, sample_2, bwa_output, cpu)

def samtools(bwa_output, output_bam, output_sort, output_tsv):
    """
    Process SAM file into sorted and indexed BAM using Samtools.
    
    Args:
        bwa_output (str): Input SAM file path.
        output_bam (str): Unsorted BAM output path.
        output_sort (str): Sorted BAM output path.
        output_tsv (str): Flagstat report output path.
    """
    try:
        # Convert SAM to BAM
        logger.info(f"Converting SAM to BAM: {output_bam}")
        subprocess.run(["samtools", "view", "-bS", "-o", output_bam, bwa_output], check=True)

        # Sort BAM
        logger.info(f"Sorting BAM: {output_sort}")
        subprocess.run(["samtools", "sort", "-o", output_sort, output_bam], check=True)

        # Index BAM
        logger.info("Indexing BAM...")
        subprocess.run(["samtools", "index", output_sort], check=True)

        # Generate flagstat statistics
        logger.info(f"Generating flagstat report: {output_tsv}")
        with open(output_tsv, 'w') as f:
            subprocess.run(["samtools", "flagstat", output_sort], stdout=f, check=True)
        logger.info("Samtools processing finished successfully.")

    except subprocess.CalledProcessError as e:
        logger.error(f"Samtools processing error: {e}")
        raise

def run_bwa_samtools(fasta, sample_1, sample_2, output_dir, cpu):
    """
    Unified entry point for alignment and Samtools processing.
    """
    bwa_output = os.path.join(output_dir, "bwa_output.sam")
    output_bam = os.path.join(output_dir, "output.bam")
    output_sort = os.path.join(output_dir, "output.sort.bam")
    output_tsv = os.path.join(output_dir, "output.tsv")

    # 1. Run BWA alignment
    run_bwa(fasta, sample_1, sample_2, bwa_output, cpu)

    # 2. Convert, Sort, Index using Samtools
    samtools(bwa_output, output_bam, output_sort, output_tsv)
