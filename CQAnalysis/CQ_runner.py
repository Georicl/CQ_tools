"""
Module for calculating and filtering CQ results.
Computes the CQ ratio (Female CPM / Male CPM) and filters windows based on user thresholds.
"""
import pandas as pd
import logging
from .standardization import CPMCalculator

# Initialize module-level logger
logger = logging.getLogger(__name__)

class CQRunner:
    """
    Orchestrates the CPM normalization and CQ calculation process.
    """
    def __init__(self, paths: dict, cq_value=0.3, m_reads_threshold=0):
        self.paths: dict = paths
        self.cq_value = float(cq_value)
        self.m_reads_threshold = int(m_reads_threshold)

    def _cq_calculate(self, combined_df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate CQ values and filter the DataFrame.
        
        CQ = F_CPM / M_CPM.
        Filters rows where CQ < cq_value and Male CPM > threshold.
        """
        # Vectorized CQ calculation
        combined_df['cq_result'] = combined_df['F_CPM'] / combined_df['M_CPM']
        
        # Handle division by zero (M_CPM = 0)
        combined_df.loc[combined_df['M_CPM'] == 0, 'cq_result'] = float('inf')

        # Apply filtering criteria
        filter_df = combined_df[
            (combined_df['cq_result'] < self.cq_value) & 
            (combined_df['M_CPM'] > self.m_reads_threshold)
        ].copy()

        return filter_df

    def cq_runner(self) -> pd.DataFrame:
        """
        Main execution logic for CQ analysis.
        Returns a filtered DataFrame and saves it to a TSV file.
        """
        # 1. Perform CPM normalization
        cpm_calc = CPMCalculator(
            m_bam_file=self.paths['m_bam'],
            f_bam_file=self.paths['f_bam'],
            reads_table_path=self.paths['f_m_merge']
        )
        result_df = cpm_calc.calculate_cpm_value()

        # 2. Calculate CQ and filter
        logger.info(f"Filtering results (CQ < {self.cq_value}, Male_CPM > {self.m_reads_threshold})...")
        final_df = self._cq_calculate(result_df)

        # 3. Save final results
        logger.info(f"Saving filtered CQ results to: {self.paths['filtered_cq']}")
        final_df.to_csv(self.paths['filtered_cq'], sep='\t', index=False)
        return final_df
