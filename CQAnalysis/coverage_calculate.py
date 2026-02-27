"""
Module for calculating sequence coverage across genomic windows.
Uses bedtools coverage to count reads and calculate the fraction of coverage for each window.

Author: Yang Xiang
Email: georicl@outlook.com
"""
import logging
import os
import subprocess
import time
from concurrent.futures.process import ProcessPoolExecutor

# Initialize module-level logger
logger = logging.getLogger(__name__)

class CoverageCalculate:
    """
    Handles parallel coverage calculation for multiple BAM files.
    """
    def __init__(self, paths, num_parallel=2):
        self.paths = paths
        self.parallel = int(num_parallel)

    def _run_full_coverage(self, bam_input, output_path):
        """
        Execute bedtools coverage on the full windows file.
        
        Input windows_tsv is a 3-column BED file (chr, start, end).
        Bedtools coverage output format (standard 3-col BED input):
            1. Chromosome
            2. Start
            3. End
            4. Number of reads overlapping the window (Used for CPM)
            5. Number of bases with coverage
            6. Window length
            7. Fraction of window covered (Used as coverage value)
        """
        cmd = (
            f"bedtools coverage -a {self.paths['windows_tsv']} -b {bam_input} "
            f"-g {self.paths['chromosome_length']} -sorted > {output_path}"
        )
        try:
            subprocess.run(cmd, shell=True, check=True)
            logger.info(f"Coverage file generated: {output_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Coverage calculation failed: {e}")
            raise

    def m_coverage_calculate(self) -> None:
        """Calculate coverage for the male sample."""
        start_time = time.time()
        self._run_full_coverage(self.paths['m_bam'], self.paths['m_cov'])
        logger.info(f"Male coverage done in {time.time() - start_time:.2f}s")

    def f_coverage_calculate(self) -> None:
        """Calculate coverage for the female sample."""
        start_time = time.time()
        self._run_full_coverage(self.paths['f_bam'], self.paths['f_cov'])
        logger.info(f"Female coverage done in {time.time() - start_time:.2f}s")

    def _merge_final_result(self) -> None:
        """
        Merge female and male coverage results into a single TSV.
        Performs strict coordinate validation to ensure data integrity.
        """
        header = "chr\tstart\tend\tf_num_reads\tf_coverage\tm_num_reads\tm_coverage\n"
        mismatch_count = 0
        line_count = 0

        with open(self.paths['f_m_merge'], "w") as f_out:
            f_out.write(header)

            with open(self.paths['f_cov'], "r") as f_f, \
                 open(self.paths['m_cov'], "r") as f_m:
                
                for line_num, (line_f, line_m) in enumerate(zip(f_f, f_m), start=1):
                    parts_f = line_f.strip().split("\t")
                    parts_m = line_m.strip().split("\t")

                    # Coordinates: chr(0), start(1), end(2)
                    chr_f, start_f, end_f = parts_f[0], parts_f[1], parts_f[2]
                    chr_m, start_m, end_m = parts_m[0], parts_m[1], parts_m[2]

                    # Coordinate validation
                    if chr_f != chr_m or start_f != start_m or end_f != end_m:
                        logger.error(f"Row {line_num} coordinate mismatch: {chr_f}:{start_f} vs {chr_m}:{start_m}")
                        mismatch_count += 1
                        continue

                    # Extract count (col 4) and fraction (col 7)
                    f_num_reads, f_coverage = parts_f[3], parts_f[6]
                    m_num_reads, m_coverage = parts_m[3], parts_m[6]

                    f_out.write(f"{chr_f}\t{start_f}\t{end_f}\t{f_num_reads}\t{f_coverage}\t{m_num_reads}\t{m_coverage}\n")
                    line_count += 1

            if mismatch_count > 0:
                logger.warning(f"Merge finished with {mismatch_count} mismatches.")
            logger.info(f"Successfully merged {line_count} rows to {self.paths['f_m_merge']}")

    def execute(self) -> None:
        """
        Run coverage calculation and merge the results.
        Supports serial (1) or parallel (2) execution.
        """
        start_time = time.time()

        if self.parallel > 1:
            logger.info("Starting parallel coverage calculation (2 processes)...")
            with ProcessPoolExecutor(max_workers=2) as executor:
                future_m = executor.submit(self.m_coverage_calculate)
                future_f = executor.submit(self.f_coverage_calculate)
                future_m.result()
                future_f.result()
        else:
            logger.info("Starting serial coverage calculation...")
            self.m_coverage_calculate()
            self.f_coverage_calculate()

        self._merge_final_result()
        logger.info(f"Total coverage workflow completed in {time.time() - start_time:.2f}s")
