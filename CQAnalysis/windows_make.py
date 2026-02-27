"""
Module for genomic window generation.
Handles FASTA indexing and slicing chromosomes into fixed-size windows.

Author: Yang Xiang
Email: georicl@outlook.com
"""
import logging
import pysam

# Initialize module-level logger
logger = logging.getLogger(__name__)

class MakeWindows:
    """
    Generates window files from a reference genome.
    """
    def __init__(self, paths):
        self.paths = paths

    def build_fai(self):
        """
        Create FASTA index file (FAI) using pysam.
        Ensures the reference genome is indexed before processing.
        """
        try:
            logger.info("Indexing FASTA file...")
            pysam.faidx(self.paths['fa_path'])
            logger.info("FAI file created successfully.")
        except Exception as e:
            logger.error(f"Error creating FAI file: {str(e)}")
            raise

    def calculate_chromosome_length(self):
        """
        Extract chromosome names and lengths from the FAI file.
        Outputs a TSV file used as a genome file for bedtools.
        """
        try:
            logger.info("Calculating chromosome lengths...")
            with open(self.paths['fai_path'], "r") as f:
                # FAI format: chr_name \t length \t ...
                lines = [(line.strip().split("\t")[0], line.strip().split("\t")[1]) for line in f]
            
            with open(self.paths['chromosome_length'], "w") as f:
                for chrom, length in lines:
                    f.write(f"{chrom}\t{length}\n")
            logger.info("Chromosome lengths written successfully.")
        except Exception as e:
            logger.error(f"Error calculating chromosome lengths: {str(e)}")
            raise

    def create_windows_py(self, windows_size=1000, step_size=500):
        """
        Create a BED-style table of windows across the genome.
        
        Args:
            windows_size (int): Size of each window in bp.
            step_size (int): Step size for sliding windows.
        """
        try:
            logger.info(f'Creating windows (size={windows_size}, step={step_size})...')
            with open(self.paths['chromosome_length'], "r") as f:
                chrom_data = {lines.strip().split()[0]: lines.strip().split()[1] for lines in f}
                
                windows = []
                for chrom, length in chrom_data.items():
                    start = 0
                    chrom_len = int(length)
                    while start < chrom_len:
                        # Ensure the end of the window doesn't exceed chromosome length
                        end = min(start + windows_size, chrom_len)
                        windows.append((chrom, start, end))
                        start += step_size
                
                with open(self.paths['windows_tsv'], 'w') as f:
                    for window in windows:
                        f.write(f'{window[0]}\t{window[1]}\t{window[2]}\n')
            logger.info('Windows created successfully.')
        except Exception as e:
            logger.error(f'Error creating windows: {str(e)}')
            raise

    def executor(self):
        """
        Runs the full window creation workflow.
        """
        self.build_fai()
        self.calculate_chromosome_length()
        self.create_windows_py()