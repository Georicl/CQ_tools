"""
Module for CPM (Counts Per Million) normalization.
Calculates normalization factors based on total mapped reads from BAM files.
"""
import subprocess
import pandas as pd
import logging

# Initialize module-level logger
logger = logging.getLogger(__name__)

class CPMCalculator:
    """
    Calculates CPM values for genomic windows to allow comparison between samples.
    """
    def __init__(self, m_bam_file, f_bam_file, reads_table_path):
        self.m_bam_path = m_bam_file
        self.f_bam_path = f_bam_file
        self.reads_table_path = reads_table_path

    def _calculate_sum_reads(self, bam_file: str) -> int:
        """
        Get the total number of mapped reads in a BAM file using Samtools.
        
        Uses 'samtools view -c -F 4' to count only mapped reads.
        """
        logger.info(f"Counting total mapped reads for: {bam_file}")
        cmd = f"samtools view -c -F 4 {bam_file}"
        results = subprocess.run(cmd, shell=True, text=True, capture_output=True)
        if results.returncode != 0:
            logger.error(f"Samtools count error: {results.stderr}")
            raise RuntimeError(f"Error running samtools: {results.stderr}")
        
        total_reads = int(results.stdout.strip())
        logger.info(f"Total mapped reads: {total_reads}")
        return total_reads

    def calculate_cpm_value(self) -> pd.DataFrame:
        """
        Perform CPM normalization on the input read counts table.
        
        Formula: CPM = (reads_in_window / total_mapped_reads) * 1,000,000
        """
        # 1. Fetch total mapped reads for both samples
        total_m = self._calculate_sum_reads(self.m_bam_path)
        total_f = self._calculate_sum_reads(self.f_bam_path)

        # 2. Read the merged coverage table
        logger.info(f"Reading merged coverage table from: {self.reads_table_path}")
        reads_table = pd.read_csv(self.reads_table_path, sep='\t')

        # 3. Calculate CPM values
        logger.info("Performing CPM normalization...")
        reads_table['F_CPM'] = (reads_table['f_num_reads'] / total_f) * 1e6
        reads_table['M_CPM'] = (reads_table['m_num_reads'] / total_m) * 1e6

        # 4. Generate a unique window identifier string
        reads_table['window'] = reads_table.apply(lambda x: f"{x['chr']}-{x['start']}-{x['end']}", axis=1)

        return reads_table
